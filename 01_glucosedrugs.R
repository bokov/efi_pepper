library(tidyr);
library(dplyr);
library(DBI);
library(rio);     # format-agnostic convenient file import
library(dint);    # date-conversion
library(digest);
library(forcats);

# runtime duration: 11:22

source('default_config.R');
#' The local path names for the data files should be stored in a vector
#' named `inputdata` that get set in a script named `local_config.R`
if(file.exists('local_config.R')) source('local_config.R');

# Verify that the three original raw files are still the same and try to exit
# if they are not.

if(digest::digest(inputdata['fullsqlfile'],file=T)!="6a16692a07c116e1559b09cbe0b92268"){
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
    DDP4I = c("Sitagliptin","Alogliptin", "Saxagliptin", "Linagliptin"),
    GLP1A = c("Albiglutide","Dulaglutide","Liraglutide","Semaglutide","Exenatide","Lixisenatide"),
    TZD = c("Rosiglitazone", "Pioglitazone"),
    Sulfonylureas = c("Glipizide","Glimepiride", "Glyburide"),
    Metformin='Metformin'
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
  WHERE date <= '2021-12-01'
)
-- the first and last month for each patient, so we don't have a lot of empty
-- months at either end
,q0 AS (
  SELECT DISTINCT patient_num
    ,DATE(MIN(start_date),'start of month') dfrom
    ,DATE(MAX(COALESCE(end_date,start_date)),'start of month') dto
  FROM observation_fact GROUP BY patient_num
)
-- all included months for all patients
,q1 AS (
  SELECT q0.patient_num,dates.date
  FROM q0 INNER JOIN dates ON dates.date BETWEEN q0.dfrom AND q0.dto
)
-- unique combos of patients, start/end dates, drug codes, and meds
-- where the meds are the antihyperglycemic drugs of interest
,q2 AS (
  SELECT DISTINCT patient_num
    ,DATE(start_date,'start of month') start_date
    ,DATE(COALESCE(end_date,start_date),'start of month') end_date
    ,ob.concept_cd,cd.name_char
  FROM observation_fact ob
  INNER JOIN concept_dimension cd
      ON ob.concept_cd = cd.concept_cd
      AND (cd.concept_cd LIKE 'NDC:%%' OR cd.concept_cd LIKE 'RXCUI:%%')
      AND (%s)
)
,q3 AS (
  SELECT DISTINCT patient_num,DATE(start_date,'start of month') start_date
  FROM observation_fact WHERE  %s
)
-- collapse down to unique patient-months and for each patient-month concatenate
-- all unique drug codes and names
,q4 AS (
  SELECT q1.patient_num,q1.date
  ,REPLACE(GROUP_CONCAT(DISTINCT concept_cd),',','; ') concept_cds
  -- because name_chars can have embedded commas
  ,REPLACE(REPLACE(GROUP_CONCAT(DISTINCT REPLACE(name_char,',','=')),',','; '),'=',',') name_chars
  -- time elapsed since last visit
  ,ROUND((JULIANDAY(date) - JULIANDAY(MAX(q3.start_date)))/(365/12)) months_since_pcvisit
  FROM q1
  LEFT JOIN q2
  ON q1.patient_num = q2.patient_num
    AND q1.date BETWEEN q2.start_date AND q2.end_date
  LEFT JOIN q3
  ON q1.patient_num = q3.patient_num
    AND q3.start_date <= q1.date
  GROUP BY q1.patient_num,q1.date
)
-- finally, a pivot table of drug categories vs patient-months
-- with a column for None and columns for all unique codes/meds for that patient-month
-- (from above)
SELECT patient_num,date
%s
, concept_cds
, name_chars
, months_since_pcvisit
FROM q4
;";
.drugsql <- sprintf(.drugsqltemplate,.drugrows,.pcviscodes,.druggroups);

# #.drugwhere <- .drugnames %>% sapply(function(xx) sprintf(' name_char LIKE "%%%s%%" ',xx) %>% paste0(collapse='OR'));
# .drugsql <- .drugnames %>% unlist %>% sprintf(' name_char LIKE "%%%s%%" ',.) %>% paste0(collapse = 'OR') %>%
#   paste0('
# WITH q0 AS (
#   SELECT DISTINCT patient_num,cd.concept_cd
#     ,DATE(start_date,"start of month") start_date
#     ,DATE(COALESCE(end_date,start_date),"start of month") end_date,cd.name_char
#   FROM observation_fact ob
#   LEFT JOIN concept_dimension cd
#     ON ob.concept_cd = cd.concept_cd
#     AND (cd.concept_cd LIKE "NDC:%" OR cd.concept_cd LIKE "RXCUI:%")
#     AND (',.,'))
# SELECT q0.patient_num,q1.concept_cd,q0.start_date,q1.name_char
# FROM q0 LEFT JOIN q0 q1
#   ON q0.patient_num = q1.patient_num
#   AND q0.start_date BETWEEN q1.start_date AND q1.end_date
#   AND q1.name_char IS NOT NULL
#          ');

.drugcols <- c(paste0('Cohort_',c(paste0(cDrugGroupsMain,'Mono'),'None','Other')),cDrugGroups);
.drugregexps <- sapply(.drugnames,paste0,collapse='|');

# Local function to avoid a lot of repetitive code when searching for strings in a column
hasString <- function(xx) substitute(any(grepl(.drugregexps[xx],name_char,ignore.case=T)));

# All glucose-lowering drug ocurrences
id_patmap <- import(inputdata['patmap']) %>%
  mutate(patient_num=as.character(PATIENT_NUM)) %>%
  select(patient_num,PAT_MRN_ID,DATE_SHIFT,PATIENT_IDE_UPDATED) %>% unique;

dat0 <- dbGetQuery(fullsqlcon,.drugsql) %>%
  mutate(patient_num=as.character(patient_num),date=as.Date(date)) %>%
  group_by(patient_num,date);
dbDisconnect(fullsqlcon);
# dat0a <- mutate(dat0,patient_num=as.character(patient_num),start_date=as.Date(start_date)) %>%
#   group_by(patient_num,start_date);


# patient history, can be joined to main data
dat1 <- mutate(dat0
               ,Secretagogues=ifelse(Sulfonylureas!='' | Glinides!='','Secretagogues','')
               ,AnyOther=ifelse(SGLT2I!=''|DDP4I!=''|GLP1A!=''|TZD!='','AnyOther','')
               ,SecretagoguesMono=ifelse(Secretagogues!=''&AnyOther==''&Metformin=='','SecretagoguesMono','')
               ,MetforminMono=ifelse(Metformin!=''&AnyOther==''&Secretagogues=='','MetforminMono','')
               #,CohortFactor=paste0(None,AnyOther,SecretagoguesMono,MetforminMono) #%>% gsub('[+]+','+',.)
               ,MonthFactor = paste(Secretagogues,Metformin,AnyOther,None,sep='+') %>%
                 gsub('[+]+','+',.) %>% gsub('^[+]|[+]$','',.)
) %>%
  ungroup(date) %>%
  mutate(CohortDetail=if(all(MonthFactor=='None')) 'None' else {
    paste0(sort(setdiff(MonthFactor,'None')),collapse=';')}
    ,CohortFactor=if(any(grepl(';',CohortDetail))) 'Multi' else CohortDetail);

# #nrow(subset(dat1,CohortFactor=='Metformin' & paste0()))
#
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

# TODO: paused here.
dat1 <- summarise(dat0a
                  ,Glinides=eval(hasString('Glinides'))
                  ,SGLT2I=eval(hasString('SGLT2I'))
                  ,DDP4I=eval(hasString('DDP4I'))
                  ,GLP1A=eval(hasString('GLP1A'))
                  ,TZD=eval(hasString('TZD'))
                  ,Metformin=eval(hasString('Metformin'))
                  ,Sulfonylureas=eval(hasString('Sulfonylureas'))
                  ,Secretagogues=Sulfonylureas|Glinides
                  ,AnyOther=SGLT2I|DDP4I|GLP1A|TZD
                  ,Secretagogues_Mono=Secretagogues & !AnyOther & !Metformin
                  ,Metformin_Mono=Metformin & !AnyOther & !Secretagogues
                  ,None=!Metformin & !Secretagogues & !AnyOther
) %>% ungroup(start_date) %>%
  mutate(Cohort_None=all(None)
         ,Cohort_MetforminMono=any(Metformin_Mono)&all(None|Metformin_Mono)
         ,Cohort_SecretagoguesMono=any(Secretagogues_Mono)&all(None|Secretagogues_Mono)
         ,Cohort_Other=!(Cohort_MetforminMono|Cohort_SecretagoguesMono|Cohort_None)
         ,across(any_of(.drugcols),~ifelse(.x,cur_column(),''),.names = 'temp_{.col}')
#         ,across(starts_with('Cohort_'),~ifelse(.x,cur_column(),''),.names = 'temp_{.col}')
  );
dat1$CohortFactor <- ungroup(dat1) %>% select(starts_with('temp_Cohort_')) %>%
  interaction(drop=T,sep='+') %>%
  fct_relabel(~gsub('\\++','+',.x) %>% gsub('^\\+|\\+$','',.));
dat1 <- select(dat1,!starts_with('temp_Cohort_'));
CohortDetail <- ungroup(dat1) %>%
  select(starts_with('temp_') & !matches('Glinides|Sulfonylureas')) %>%
  interaction(drop=T,sep='+') %>%
  fct_relabel(~gsub('\\++','+',.x) %>% gsub('^\\+|\\+$','',.));
dat1$CohortDetail <- ifelse(dat1$CohortFactor=='Cohort_Other'
                            ,as.character(CohortDetail)
                            ,as.character(dat1$CohortFactor));
dat1 <- select(dat1,!starts_with('temp_'));

# start
dat2 <- select(dat1,any_of(c('patient_num','CohortFactor'
                                          ,'start_date',cDrugGroups,'None'))) %>%
  # these steps are to organize the data by patient, drug-group, and contiguous
  # sequence of encounters during which they were prescribed drugs from that
  # group
  pivot_longer(cols=any_of(c(cDrugGroups,'None'))
               ,names_to='Drug',values_to='Active') %>%
  group_by(patient_num,CohortFactor,Drug) %>% arrange(start_date) %>%
  mutate(interval=with(rle(Active),rep(seq_along(lengths),lengths))) %>% # 24030610 rows
  subset(Active) %>%                                                     #  3596758 rows
  group_by(interval,.add=T) %>%
  summarize(FromDate=min(start_date),ToDate=max(start_date)) %>%
  # merging in MRNs and date-shifts to enable chart review
  left_join(id_patmap) %>%
  # mutate(FromDate=FromDate+DATE_SHIFT,ToDate=ToDate+DATE_SHIFT) %>%
  # now, cleaning up for readability
  select(-c('interval','DATE_SHIFT')) %>%
  arrange(patient_num,FromDate) %>%
  mutate(CohortFactor=gsub('Cohort_','',CohortFactor)) %>%
  rename(Cohort=CohortFactor) %>%
  relocate(PAT_MRN_ID,Drug,FromDate,ToDate,Cohort);


export(dat1,file='DEID_GLUDRUGS.tsv');
export(dat2,file='PHI_GLUDRUG_DATES.xlsx');
