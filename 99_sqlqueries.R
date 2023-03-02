# drugnames ####
# This is where drugs AND their categories are assigned for the entire project
.drugnames <- list(
  Glinides = c("Repaglinide", "Nateglinide"),
  SGLT2I = c("Empagliflozin","Canagliflozin","Dapagliflozin","Ertugliflozin"),
  DDP4I = c("Sitagliptin","Alogliptin", "Saxagliptin", "Linagliptin", "Vildagliptin"),
  GLP1A = c("Albiglutide","Dulaglutide","Liraglutide","Semaglutide","Exenatide","Lixisenatide"),
  TZD = c("Rosiglitazone", "Pioglitazone", "Troglitazone"),
  Sulfonylureas = c("Glipizide","Glimepiride", "Glyburide", "Tolazamide","Tolbutamide","Chlorpropamide","Glibornuride","Gliclazide","Gliquidone"),
  Metformin='Metformin',
  Insulin='Insulin',
  AGI=c('Acarbose','Miglitol'),
  OtherDrugs=c('Pramlintide','Ropinirole','Cholecalciferol')
);

# This is an intermediate obect use to create both snippets of drug-related SQL.
.drugmatches <- lapply(.drugnames,function(xx){
  sprintf(" name_char LIKE '%%%s%%' ",xx) %>% paste0(collapse = ' OR ')});



# Drug-category columns ####
.druggroups <- unlist(.drugmatches) %>%
  paste0('\n, CASE WHEN ',.," THEN '",names(.),"' ELSE '' END ",names(.)
         ,collapse='') %>%
  gsub('name_char','name_chars',.) %>%
  paste0("\n, CASE WHEN  name_chars IS NULL THEN 'None' ELSE '' END None");

# Drug selection criteria ####
.drugrows <- unlist(.drugmatches) %>% paste(collapse=' OR ');
# These are codes for PCP visits
.pcviscodes <- paste0("concept_cd LIKE '%SPEC:",c(7,9,17,81,84,183),"%'"
                      ,collapse=' OR ');

# Drug obs table ####
.drugobstable <- "CREATE TABLE scriptmp.obs AS
SELECT DISTINCT -- 37,914,003 vs 42,618,747
  observation_fact.patient_num,concept_cd,sourcesystem_cd
  ,start_date start_date_orig ,COALESCE(end_date,start_date) end_date_orig
  ,DATE(DATE(start_date,-DATE_SHIFT||' days'),'start of month') start_date
  ,DATE(DATE(COALESCE(end_date,start_date),-DATE_SHIFT||' days'),'start of month') end_date
  ,DATE_SHIFT,PATIENT_IDE_UPDATED,PAT_MRN_ID
FROM observation_fact
INNER JOIN scriptmp.id_patmap
ON observation_fact.patient_num = id_patmap.PATIENT_NUM
";

# Drug combos table ####
.drugcombostable <- sprintf("CREATE TABLE scriptmp.combos AS
SELECT DISTINCT patient_num
    ,start_date
    ,end_date
    ,obs.concept_cd,cd.name_char,obs.sourcesystem_cd
  FROM scriptmp.obs
  INNER JOIN concept_dimension cd
      ON obs.concept_cd = cd.concept_cd
      AND (cd.concept_cd LIKE 'NDC:%%' OR cd.concept_cd LIKE 'RXCUI:%%')
      AND ( %s )",.drugrows);

.drugresultquery <- sprintf("WITH RECURSIVE dates(date) AS (
  VALUES('2005-01-01')
  UNION ALL
  SELECT date(date, '+1 month')
  FROM dates
  WHERE date <= '2022-12-01'
)
, q0 AS (
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
-- primary care visits
,qv0 AS (
  SELECT DISTINCT patient_num,start_date
  FROM obs
  -- .pcviscodes go here
  WHERE %s
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
  FROM q1 LEFT JOIN combos q2
  ON q1.patient_num = q2.patient_num
  AND q1.date BETWEEN q2.start_date AND q2.end_date
  LEFT JOIN qv0
  ON q1.patient_num = qv0.patient_num AND qv0.start_date <= q1.date
  GROUP BY q1.patient_num,q1.date
)
-- finally, a pivot table of drug categories vs patient-months with a column for
-- None and columns for all unique codes/meds for that patient-month (from above)
SELECT patient_num,PAT_MRN_ID,start_date
  ,N_UHS_Records,N_UTMed_Records
  ,DATE(DATE(start_date,DATE_SHIFT||' days'),'start of month') shifted_date
  -- .druggroups go here
  %s
  , concept_cds
  , name_chars
  , months_since_pcvisit
  , DATE_SHIFT,PATIENT_IDE_UPDATED
FROM q3",.pcviscodes,.druggroups)

# Old stuff ####


# WITH RECURSIVE dates(date) AS (
#   VALUES('2005-01-01')
#   UNION ALL
#   SELECT date(date, '+1 month')
#   FROM dates
#   WHERE date <= '2022-12-01'
# )
# -- the first and last month for each patient, so we don't have a lot of empty
# -- months at either end





#### The old version of the main query that brough SQLite engine to its knees
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
