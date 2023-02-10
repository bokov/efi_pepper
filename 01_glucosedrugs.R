library(tidyr);
library(dplyr);
library(DBI);
library(rio);     # format-agnostic convenient file import
library(dint);    # date-conversion
library(zoo);     # more date conversion?!
library(digest);
library(forcats);
library(stringr);

options(datatable.integer64='numeric');

runReplace <- function(xx,fillwith,fillinto,length=6){
  if(length(setdiff(xx,c(fillinto,fillwith)))>0){browser()
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

# Verify that the three original raw files are still the same and try to exit
# if they are not.

if(!file.exists('.ignorehashes') && digest::digest(inputdata['fullsqlfile'],file=T)!="6a16692a07c116e1559b09cbe0b92268"){
  rstudioapi::restartSession(command="warning('The original raw data files seem to be corrupted! You should not proceed further. Ask Alex to restore them from backup!')");
  stop('The original raw data files seem to be corrupted! You should not proceed further. Ask Alex to restore them from backup!');
  NOFUNCTION();
};

fullsql <- file.path(tempdir,'fullsql.db');
if(!file.exists(fullsql) || digest::digest(fullsql,file=T) != "6a16692a07c116e1559b09cbe0b92268"){
  file.copy(inputdata['fullsqlfile'],fullsql<-file.path(tempdir,'fullsql.db'),overwrite = T)};
fullsqlcon <- dbConnect(RSQLite::SQLite(), fullsql);

.drugnames <- list(
    Glinides = c("Repaglinide", "Nateglinide"),
    SGLT2I = c("Empagliflozin","Canagliflozin","Dapagliflozin","Ertugliflozin"),
    DDP4I = c("Sitagliptin","Alogliptin", "Saxagliptin", "Linagliptin", "Vildagliptin"),
    GLP1A = c("Albiglutide","Dulaglutide","Liraglutide","Semaglutide","Exenatide","Lixisenatide"),
    TZD = c("Rosiglitazone", "Pioglitazone", "Troglitazone"),
    Sulfonylureas = c("Glipizide","Glimepiride", "Glyburide", "Tolazamide","Tolbutamide","Chlorpropamide","Glibornuride","Gliclazide","Gliquidone"),
    Metformin='Metformin',
    Insulin='Insulin'
  );

.drugmatches <- lapply(.drugnames,function(xx){
  sprintf(" name_char LIKE '%%%s%%' ",xx) %>% paste0(collapse = ' OR ')});
.druggroups <- unlist(.drugmatches) %>%
  paste0('\n, CASE WHEN ',.," THEN '",names(.),"' ELSE '' END ",names(.)
         ,collapse='') %>%
  gsub('name_char','name_chars',.) %>%
  paste0("\n, CASE WHEN  name_chars IS NULL THEN 'None' ELSE '' END None");
.drugrows <- unlist(.drugmatches) %>% paste(collapse=' OR ');
.pcviscodes <- paste0('concept_cd LIKE "%SPEC:',c(7,9,17,81,84,183),'%"'
                      ,collapse=' OR ');

.drugsqltemplate <- "WITH RECURSIVE dates(date) AS (
  VALUES('2005-01-01')
  UNION ALL
  SELECT date(date, '+1 month')
  FROM dates
  WHERE date <= '2022-12-01'
)
-- the first and last month for each patient, so we don't have a lot of empty
-- months at either end
,obs AS (
  SELECT DISTINCT -- 37,914,003 vs 42,618,747
    observation_fact.patient_num,concept_cd,sourcesystem_cd
    ,start_date start_date_orig ,COALESCE(end_date,start_date) end_date_orig
    ,DATE(DATE(start_date,-DATE_SHIFT||' days'),'start of month') start_date
    ,DATE(DATE(COALESCE(end_date,start_date),-DATE_SHIFT||' days'),'start of month') end_date
    ,DATE_SHIFT,PATIENT_IDE_UPDATED,PAT_MRN_ID
  FROM observation_fact
  INNER JOIN xwalk
  ON observation_fact.patient_num = xwalk.PATIENT_NUM
)
,q0 AS (
  SELECT patient_num
    ,MIN(start_date) dfrom
    ,MAX(end_date) dto
    ,DATE_SHIFT,PATIENT_IDE_UPDATED,PAT_MRN_ID
  FROM obs GROUP BY patient_num
)
-- all included months for all patients
,q1 AS (
  SELECT q0.patient_num,dates.date,DATE_SHIFT,PATIENT_IDE_UPDATED,PAT_MRN_ID
  FROM q0 INNER JOIN dates ON dates.date BETWEEN q0.dfrom AND q0.dto
)
-- unique combos of patients, start/end dates, drug codes, and meds
-- where the meds are the antihyperglycemic drugs of interest
,q2 AS (
  SELECT DISTINCT patient_num
    ,start_date
    ,end_date
    ,obs.concept_cd,cd.name_char,obs.sourcesystem_cd
  FROM obs
  INNER JOIN concept_dimension cd
      ON obs.concept_cd = cd.concept_cd
      AND (cd.concept_cd LIKE 'NDC:%%' OR cd.concept_cd LIKE 'RXCUI:%%')
      AND ( %s )
)
-- primary care visits
,qv0 AS (
  SELECT DISTINCT patient_num,start_date
  FROM obs WHERE  %s
)
-- collapse down to unique patient-months and for each patient-month concatenate
-- all unique drug codes and names
,q3 AS (
  SELECT q1.patient_num,q1.date start_date
  ,REPLACE(GROUP_CONCAT(DISTINCT concept_cd),',','; ') concept_cds
  -- because name_chars can have embedded commas
  ,REPLACE(REPLACE(GROUP_CONCAT(DISTINCT REPLACE(name_char,',','=')),',','; '),'=',',') name_chars
  -- time elapsed since last visit
  ,ROUND((JULIANDAY(q1.date) - JULIANDAY(MAX(qv0.start_date)))/(365/12)) months_since_pcvisit
  ,MAX(q1.DATE_SHIFT) DATE_SHIFT
  ,MAX(PATIENT_IDE_UPDATED) PATIENT_IDE_UPDATED
  ,MAX(PAT_MRN_ID) PAT_MRN_ID
  ,SUM(CASE sourcesystem_cd WHEN 'Sunrise@UHS' THEN 1 ELSE 0 END) N_UHS_Records
  ,SUM(CASE sourcesystem_cd WHEN 'Clarity@UTMed' THEN 1 ELSE 0 END) N_UTMed_Records
  FROM q1 LEFT JOIN q2
  ON q1.patient_num = q2.patient_num
  AND q1.date BETWEEN q2.start_date AND q2.end_date
  LEFT JOIN qv0
  ON q1.patient_num = qv0.patient_num AND qv0.start_date <= q1.date
  GROUP BY q1.patient_num,q1.date
)
-- finally, a pivot table of drug categories vs patient-months
-- with a column for None and columns for all unique codes/meds for that patient-month
-- (from above)
SELECT patient_num,PAT_MRN_ID,start_date
  ,N_UHS_Records,N_UTMed_Records
  ,DATE(DATE(start_date,DATE_SHIFT||' days'),'start of month') shifted_date
  %s
  , concept_cds
  , name_chars
  , months_since_pcvisit
  , DATE_SHIFT,PATIENT_IDE_UPDATED
FROM q3
";
.drugsql <- sprintf(.drugsqltemplate,.drugrows,.pcviscodes,.druggroups);

# All glucose-lowering drug ocurrences
id_patmap <- import(inputdata['patmap']) %>%
  mutate(across(any_of(c('PATIENT_NUM','PAT_MRN_ID')),as.character)) %>%
  select(PATIENT_NUM,PAT_MRN_ID,DATE_SHIFT,PATIENT_IDE_UPDATED) %>% unique;

dbWriteTable(fullsqlcon,'xwalk',id_patmap
             #,unique(import(inputdata['patmap'])[,c('PATIENT_NUM','DATE_SHIFT','PATIENT_IDE_UPDATED','PAT_MRN_ID')])
             ,overwrite=T);
dat0 <- dbGetQuery(fullsqlcon,.drugsql) %>%
  mutate(across(any_of(c('patient_num','PAT_MRN_ID')),as.character)
         ,across(ends_with('_date'),as.Date)) %>%
  group_by(patient_num) %>% arrange(PAT_MRN_ID,start_date);
#dbExecute(fullsqlcon,'DROP TABLE xwalk');
#dbDisconnect(fullsqlcon);

# patient history, can be joined to main data
dat1 <- mutate(
  dat0
  ,across(any_of(c('Glinides', 'SGLT2I', 'DDP4I', 'GLP1A', 'TZD', 'Sulfonylureas', 'Metformin', 'Insulin')),~runReplace(.x,cur_column(),''))
  ,Secretagogues=ifelse(Sulfonylureas!='' | Glinides!='','Secretagogues','')
  ,AnyOther=ifelse(SGLT2I!=''|DDP4I!=''|GLP1A!=''|TZD!='','AnyOther','')
  ,SecretagoguesMono=ifelse(Secretagogues!=''&AnyOther==''&Metformin==''&Insulin=='','SecretagoguesMono','')
  ,MetforminMono=ifelse(Metformin!=''&AnyOther==''&Secretagogues==''&Insulin=='','MetforminMono','')
  ,InsulinMono=ifelse(Metformin==''&AnyOther==''&Secretagogues==''&Insulin!='','InsulinMono','')
  ,None = ifelse(Secretagogues == '' & AnyOther == '' & Metformin == '' & Insulin == '','None','')
  ,CohortDetail=paste(Glinides,SGLT2I,DDP4I,GLP1A,TZD,Sulfonylureas,Metformin,Insulin,None,sep='+') %>%
    gsub('[+]+','+',.) %>% gsub('^[+]|[+]$','',.)
  #,CohortFactor=paste(Secretagogues,Metformin,Insulin,AnyOther,None,sep='+')
  #,CohortFactor=paste0(None,AnyOther,SecretagoguesMono,MetforminMono) #%>% gsub('[+]+','+',.)
  ,MonthFactor = paste(Secretagogues,Metformin,Insulin,AnyOther,None,sep='+') %>%
    gsub('[+]+','+',.) %>% gsub('^[+]|[+]$','',.)
  # ,CohortDetail=if(all(MonthFactor=='None')) 'None' else {
  #   paste0(sort(setdiff(MonthFactor,'None')),collapse=';')}
  # ,CohortFactor=if(any(grepl(';',CohortDetail))) 'Multi' else CohortDetail
  ) %>% arrange(PAT_MRN_ID,start_date) %>% ungroup %>% group_by(patient_num) %>%
  mutate(CohortDetail=lapply(CohortDetail,str_split_1,'\\+') %>% unlist %>%
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
  #   ,MFLag = lag(MonthFactor,default='None'), MFLead = lead(MonthFactor,default='None')
  #   # Months where the preceding and next month are the same as each other and
  #   # involve multiple drugs, and this current (intervening) month has a
  #   # strict subset of those drugs
  #   ,MF_MFLead= MFLag == MFLead & grepl('\\+',MFLead) & mapply(function(aa,bb){
  #     all(str_split_1(aa,'\\+') %in% str_split_1(bb,'\\+'))
  #     },MonthFactor,MFLead)
  #   # ...or, the intervening month is 'None' which is a subset of any treatment
  #   ,sandwiched = (MonthFactor=='None' & MonthFactor != MFLead) | MF_MFLead
  #   # If either of these conditions for sandwichedness are met, copy over the
  #   # preceding month's status
  #   ,MonthFactor = ifelse(sandwiched,MFLag,MonthFactor)
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
            ,Metformin=any(Metformin!='')
            ,Insulin=any(Insulin!='')
            ,Secretagogues=any(Secretagogues!='')
            ,SGLT2I=any(SGLT2I!=''),DDP4I=any(DDP4I!=''),GLP1A=any(GLP1A!='')
            ,TZD=any(TZD!=''),None=all(None!='')
            ,N_UHS_Records=max(N_UHS_Records)
            ,N_UTMed_Records=max(N_UTMed_Records)
            ,name_chars=paste0(unique(coalesce(name_chars,'')),collapse=';')
            ,concept_cds=paste0(unique(coalesce(concept_cds,'')),collapse=';')
            ,CohortFactor=max(CohortFactor)
            ,CohortDetail=max(CohortDetail)
            ,patient_num=paste0(unique(patient_num),collapse=';')
            ,PATIENT_IDE_UPDATED=max(PATIENT_IDE_UPDATED)
            ,DATE_SHIFT=paste0(unique(DATE_SHIFT),collapse=';')
            ) %>% unique %>% ungroup #%>%
  # mutate(
  #   duration= difftime(to_date,from_date,units='weeks') %>% as.numeric
  #   ,MFLag = lag(MonthFactor,default='None'), MFLead = lead(MonthFactor,default='None')
  #   ,MF_MFLead= MFLag == MFLead & grepl('\\+',MFLead) & mapply(function(aa,bb){
  #     all(str_split_1(aa,'\\+') %in% str_split_1(bb,'\\+'))
  #   },MonthFactor,MFLead)
  #   ,sandwiched = duration <= 24 & ((MonthFactor=='None' & MonthFactor != MFLead) | MF_MFLead )
  #   ,MonthFactor2 = ifelse(sandwiched,MFLag,MonthFactor)
  # );

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
