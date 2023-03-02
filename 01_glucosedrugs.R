# libraries ----
library(tidyr);
library(dplyr);
library(DBI);
library(rio);     # format-agnostic convenient file import
library(dint);    # date-conversion
library(zoo);     # more date conversion?!
library(digest);
library(forcats);
library(stringr);

# init  ----
options(datatable.integer64='numeric');

runReplace <- function(xx,fillwith,fillinto,length=6){
  if(length(setdiff(xx,c(fillinto,fillwith)))>0){#browser()
    stop('The xx argument must have two unique values, that correspond to the fillwith and fillinto arguments')};
    xxrle <- rle(xx);
    xxrle$values[with(xxrle,values==fillinto & lengths<=length)] <- fillwith;
    inverse.rle(xxrle);
}

# runtime duration: 11:22

source('default_config.R');
#' The local path names for the data files should be stored in a vector
#' named `inputdata` that get set in a script named `local_config.R`
if(file.exists('local_config.R')) source('local_config.R');

# load sql queries ----
source('99_sqlqueries.R');

# verify data sources ----
# Verify that the three original raw files are still the same and try to exit
# if they are not.
if(!file.exists('.ignorehashes') && digest::digest(inputdata['fullsqlfile'],file=T)!="6a16692a07c116e1559b09cbe0b92268"){
  rstudioapi::restartSession(command="warning('The original raw data files seem to be corrupted! You should not proceed further. Ask Alex to restore them from backup!')");
  stop('The original raw data files seem to be corrupted! You should not proceed further. Ask Alex to restore them from backup!');
  NOFUNCTION();
};
# Make sure the fullsql database exists and is true to the original, otherwise
# get a fresh copy
if(!file.exists(fullsql) || digest::digest(fullsql,file=T) != "6a16692a07c116e1559b09cbe0b92268"){
  file.copy(inputdata['fullsqlfile'],fullsql<-file.path(tempdir,'fullsql.db'),overwrite = T)};

# prepare db connections ----
fullsql <- file.path(tempdir,'fullsql.db');
# scriptmp is where we write temporary tables without cluttering up the fullsql
# source data. Maybe now the hashes will stop drifting.
scriptmp <- file.path(tempdir,'scriptmp.db');
if(file.exists(scriptmp)) unlink(scriptmp);

# Create connections for the two databases
fullsqlcon <- dbConnect(RSQLite::SQLite(), fullsql,read_only=TRUE);
scriptmpcon <- dbConnect(RSQLite::SQLite(), scriptmp);

# patient mapping ----
# Import the i2b2<->Clarity mapping table
id_patmap <- import(inputdata['patmap']) %>%
  mutate(across(any_of(c('PATIENT_NUM','PAT_MRN_ID')),as.character)) %>%
  select(PATIENT_NUM,PAT_MRN_ID,DATE_SHIFT,PATIENT_IDE_UPDATED) %>% unique;
# Load this table into scriptmp
dbWriteTable(scriptmpcon,'id_patmap',id_patmap,overwrite=T);
dbDisconnect(scriptmpcon);

# prepare temp SQLite tables ----
# Now attach scriptmp from within the read-only fullsql db
dbSendQuery(fullsqlcon,sprintf("ATTACH '%s' AS scriptmp",scriptmp));
# Create the temporary obs table used by other queries
dbSendQuery(fullsqlcon,.drugobstable);
# Create the combos table
dbSendQuery(fullsqlcon,.drugcombostable);
# release the temp db
dbSendQuery(fullsqlcon,"DETACH scriptmp");
# Done with fullsqlcon
dbDisconnect(fullsqlcon);
# Now we reconnect to the temp db and get whatever we need from there
scriptmpcon <- dbConnect(RSQLite::SQLite(), scriptmp);

# pull in the data ----
# All glucose-lowering drug ocurrences
dat0 <- dbGetQuery(scriptmpcon,.drugresultquery) %>%
  mutate(across(any_of(c('patient_num','PAT_MRN_ID')),as.character)
         ,across(ends_with('_date'),as.Date)) %>%
  group_by(patient_num) %>% arrange(PAT_MRN_ID,start_date);

# patient history, can be joined to main data
dat1 <- mutate(
  dat0
  ,across(any_of(names(.drugnames)),~runReplace(.x,cur_column(),''))
  ,Secretagogues=ifelse(Sulfonylureas!='' | Glinides!='','Secretagogues','')
  ,AnyOther=ifelse(SGLT2I!=''|DDP4I!=''|GLP1A!=''|TZD!=''|AGI!=''|OtherDrugs!='','AnyOther','')
  ,SecretagoguesMono=ifelse(Secretagogues!=''&AnyOther==''&Metformin==''&Insulin=='','SecretagoguesMono','')
  ,MetforminMono=ifelse(Metformin!=''&AnyOther==''&Secretagogues==''&Insulin=='','MetforminMono','')
  ,InsulinMono=ifelse(Metformin==''&AnyOther==''&Secretagogues==''&Insulin!='','InsulinMono','')
  ,None = ifelse(Secretagogues == '' & AnyOther == '' & Metformin == '' & Insulin == '','None','')
  ,CohortDetail=paste(Glinides,SGLT2I,DDP4I,GLP1A,TZD,AGI,OtherDrugs,Sulfonylureas,Metformin,Insulin,None,sep='+') %>%
    gsub('[+]+','+',.) %>% gsub('^[+]|[+]$','',.)
  #,CohortFactor=paste(Secretagogues,Metformin,Insulin,AnyOther,None,sep='+')
  #,CohortFactor=paste0(None,AnyOther,SecretagoguesMono,MetforminMono) #%>% gsub('[+]+','+',.)
  ,MonthFactor = paste(Secretagogues,Metformin,Insulin,AnyOther,None,sep='+') %>%
    gsub('[+]+','+',.) %>% gsub('^[+]|[+]$','',.)
  # ,CohortDetail=if(all(MonthFactor=='None')) 'None' else {
  #   paste0(sort(setdiff(MonthFactor,'None')),collapse=';')}
  # ,CohortFactor=if(any(grepl(';',CohortDetail))) 'Multi' else CohortDetail
  ) %>% arrange(PAT_MRN_ID,start_date) %>% ungroup %>% group_by(patient_num) %>%
  mutate(
    # MFLag = lag(MonthFactor,default='None'), MFLead = lead(MonthFactor,default='None')
    # Months where the preceding and next month are the same as each other and
    # involve multiple drugs, and this current (intervening) month has a
    # strict subset of those drugs
    # ,MF_MFLead= MFLag == MFLead & grepl('\\+',MFLead) & mapply(function(aa,bb){
    #   all(str_split_1(aa,'\\+') %in% str_split_1(bb,'\\+'))
    #   },MonthFactor,MFLead)
    # ...or, the intervening month is 'None' which is a subset of any treatment
    #,sandwiched = (MonthFactor=='None' & MonthFactor != MFLead) | MF_MFLead
    # If either of these conditions for sandwichedness are met, copy over the
    # preceding month's status
    #,MonthFactor = ifelse(sandwiched,MFLag,MonthFactor)
     CohortDetail=lapply(CohortDetail,str_split_1,'\\+') %>% unlist %>%
           setdiff('None') %>% sort %>% paste0(collapse='+') %>%
           ifelse(.=='','None',.)
         ,CohortFactor = lapply(MonthFactor,str_split_1,'\\+') %>% unlist %>%
           setdiff('None') %>% sort %>% paste0(collapse='+') %>%
           ifelse(.=='','None',.)
         ) %>% ungroup;
  # ungroup(start_date) %>%
  # mutate(CohortDetail=if(all(MonthFactor=='None')) 'None' else {
  #   paste0(sort(setdiff(MonthFactor,'None')),collapse=';')}
  #   ,CohortFactor=if(any(grepl(';',CohortDetail))) 'Multi' else CohortDetail
  #   # identify months of 'None' sandwiched between months of identical treatments
  #   # and make them part of the same interval
  # ) %>%
  # arrange(PAT_MRN_ID,start_date);

uniform_periods <- with(dat1,rle(paste0(PAT_MRN_ID,':',MonthFactor)));
uniform_periods$values <- seq_along(uniform_periods$values);
dat1$rownumber <- inverse.rle(uniform_periods);

dat1lds <- select(dat1,!any_of(c(idfields,'start_date'))) %>%
  rename(start_month=shifted_date);

dat2<-group_by(dat1,rownumber,PAT_MRN_ID,MonthFactor) %>%
  summarise( randomgroup =substring(patient_num,nchar(patient_num)-2)
            ,from_date=min(start_date),to_date=max(start_date)
            ,duration= (as.yearmon(to_date)-as.yearmon(from_date))*12
            ,min_since_pc=min(months_since_pcvisit)
            ,max_since_pc=max(months_since_pcvisit)
            # All the individual drug T/F columns are created by the next line
            ,across(any_of(c(names(.drugnames),'Secretagogues','None')),~any(.x!=''))
            ,N_UHS_Records=max(N_UHS_Records)
            ,N_UTMed_Records=max(N_UTMed_Records)
            ,name_chars=paste0(unique(coalesce(name_chars,'')),collapse=';')
            ,concept_cds=paste0(unique(coalesce(concept_cds,'')),collapse=';')
            ,CohortFactor=max(CohortFactor)
            ,CohortDetail=max(CohortDetail)
            ,patient_num=paste0(unique(patient_num),collapse=';')
            ,PATIENT_IDE_UPDATED=max(PATIENT_IDE_UPDATED)
            ,DATE_SHIFT=paste0(unique(DATE_SHIFT),collapse=';')
            ) %>% unique %>% ungroup %>% group_by(patient_num) %>%
  mutate(
   # duration= difftime(to_date,from_date,units='weeks') %>% as.numeric
    MFLag = lag(MonthFactor,default='None'), MFLead = lead(MonthFactor,default='None')
    # Intervals where the preceding and next month are the same as each other and
    # involve multiple drugs, and this current (intervening) month has a
    # strict subset of those drugs
    ,MF_MFLead= MFLag == MFLead & grepl('\\+',MFLead) & mapply(function(aa,bb){
      all(str_split_1(aa,'\\+') %in% str_split_1(bb,'\\+'))
    },MonthFactor,MFLead)
    ,sandwiched = duration <= 6 & ((MonthFactor=='None' & MonthFactor != MFLead) | MF_MFLead )
    ,MonthFactor2 = ifelse(sandwiched,MFLag,MonthFactor)
  );

# uniform_periods2 <- with(dat2,rle(paste0(PAT_MRN_ID,':',MonthFactor2)));
# uniform_periods2$values <- seq_along(uniform_periods2$values);
# dat2$rownumber2 <- inverse.rle(uniform_periods2);

# DONE: (not needed after all) make the long-form version of dat2 (i.e. no aggregation, just
#       rearranging the columns in the above order)
# DONE: make the deid version of dat1, for direct incorporation into the
#       analysis
# # check for uniqueness of patient_num,CohortDetail,CohortFactor
# # make sure None, Metformin, Secretagogue
#   summarise(
#     across(any_of(c(names(.drugnames),'None')),~{xx<- .x; c(setdiff(xx,''),'')[1]})
#     ,Cohort_SecretagoguesMono = if(all(paste0(SecretagoguesMono,None)!='') && any(SecretagoguesMono!='')){
#       'Cohort_SecretagoguesMono'} else ''
#     ,Cohort_MetforminMono = if(all(paste0(MetforminMono,None)!='') && any(MetforminMono!='')){
#       'Cohort_MetforminMono'} else ''
#     ,Cohort_None = if(all(None!='')) 'Cohort_None' else ''
#     #,dbg = if(any(patient_num=='1000090534')) {browser(); ''} else {''}
#     ,Cohort_Other = if(all(Cohort_MetforminMono==''&Cohort_SecretagoguesMono==''&Cohort_None=='')){
#       'Cohort_Other'} else ''
#       );

# dat1 <- summarise(dat0a
#                   ,Glinides=eval(hasString('Glinides'))
#                   ,SGLT2I=eval(hasString('SGLT2I'))
#                   ,DDP4I=eval(hasString('DDP4I'))
#                   ,GLP1A=eval(hasString('GLP1A'))
#                   ,TZD=eval(hasString('TZD'))
#                   ,Metformin=eval(hasString('Metformin'))
#                   ,Sulfonylureas=eval(hasString('Sulfonylureas'))
#                   ,Secretagogues=Sulfonylureas|Glinides
#                   ,AnyOther=SGLT2I|DDP4I|GLP1A|TZD
#                   ,Secretagogues_Mono=Secretagogues & !AnyOther & !Metformin
#                   ,Metformin_Mono=Metformin & !AnyOther & !Secretagogues
#                   ,None=!Metformin & !Secretagogues & !AnyOther
# ) %>% ungroup(start_date) %>%
#   mutate(Cohort_None=all(None)
#          ,Cohort_MetforminMono=any(Metformin_Mono)&all(None|Metformin_Mono)
#          ,Cohort_SecretagoguesMono=any(Secretagogues_Mono)&all(None|Secretagogues_Mono)
#          ,Cohort_Other=!(Cohort_MetforminMono|Cohort_SecretagoguesMono|Cohort_None)
#          ,across(any_of(.drugcols),~ifelse(.x,cur_column(),''),.names = 'temp_{.col}')
# #         ,across(starts_with('Cohort_'),~ifelse(.x,cur_column(),''),.names = 'temp_{.col}')
#   );
# dat1$CohortFactor <- ungroup(dat1) %>% select(starts_with('temp_Cohort_')) %>%
#   interaction(drop=T,sep='+') %>%
#   fct_relabel(~gsub('\\++','+',.x) %>% gsub('^\\+|\\+$','',.));
# dat1 <- select(dat1,!starts_with('temp_Cohort_'));
# CohortDetail <- ungroup(dat1) %>%
#   select(starts_with('temp_') & !matches('Glinides|Sulfonylureas')) %>%
#   interaction(drop=T,sep='+') %>%
#   fct_relabel(~gsub('\\++','+',.x) %>% gsub('^\\+|\\+$','',.));
# dat1$CohortDetail <- ifelse(dat1$CohortFactor=='Cohort_Other'
#                             ,as.character(CohortDetail)
#                             ,as.character(dat1$CohortFactor));
# dat1 <- select(dat1,!starts_with('temp_'));
#
# # start
# dat2 <- select(dat1,any_of(c('patient_num','CohortFactor'
#                                           ,'start_date',cDrugGroups,'None'))) %>%
#   # these steps are to organize the data by patient, drug-group, and contiguous
#   # sequence of encounters during which they were prescribed drugs from that
#   # group
#   pivot_longer(cols=any_of(c(cDrugGroups,'None'))
#                ,names_to='Drug',values_to='Active') %>%
#   group_by(patient_num,CohortFactor,Drug) %>% arrange(start_date) %>%
#   mutate(interval=with(rle(Active),rep(seq_along(lengths),lengths))) %>% # 24030610 rows
#   subset(Active) %>%                                                     #  3596758 rows
#   group_by(interval,.add=T) %>%
#   summarize(FromDate=min(start_date),ToDate=max(start_date)) %>%
#   # merging in MRNs and date-shifts to enable chart review
#   left_join(id_patmap) %>%
#   # mutate(FromDate=FromDate+DATE_SHIFT,ToDate=ToDate+DATE_SHIFT) %>%
#   # now, cleaning up for readability
#   select(-c('interval','DATE_SHIFT')) %>%
#   arrange(patient_num,FromDate) %>%
#   mutate(CohortFactor=gsub('Cohort_','',CohortFactor)) %>%
#   rename(Cohort=CohortFactor) %>%
#   relocate(PAT_MRN_ID,Drug,FromDate,ToDate,Cohort);

# Possibly useful summary:
# cf <- group_by(dat1,CohortFactor) %>% summarise(N=length(unique(patient_num)))
# cd <- group_by(dat1,CohortDetail) %>% summarise(N=length(unique(patient_num)))
# cfd <- full_join(cf,cd,by=c(CohortFactor="CohortDetail")) %>% mutate(NN=coalesce(N.y,N.x));
# View(cfd);


export(dat1,file='PHI_GLUDRUGS.tsv.zip');
export(dat1lds,file='DEID_GLUDRUGS.tsv.zip');
export(dat2,file='PHI_GLUDRUG_DATES.xlsx',overwrite=T,keepNA=F,firstRow=T);
