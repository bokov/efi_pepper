library(dplyr);
library(DBI);
library(rio);     # format-agnostic convenient file import
library(dint);    # date-conversion
library(digest);

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
  file.copy(inputdata['fullsqlfile'],fullsql<-file.path(tempdir,'fullsql.db'))};
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

#.drugwhere <- .drugnames %>% sapply(function(xx) sprintf(' name_char LIKE "%%%s%%" ',xx) %>% paste0(collapse='OR'));
.drugsql <- .drugnames %>% unlist %>% sprintf(' name_char LIKE "%%%s%%" ',.) %>% paste0(collapse = 'OR') %>%
  paste0('SELECT DISTINCT patient_num,cd.concept_cd,start_date,cd.name_char
         FROM observation_fact ob
         LEFT JOIN concept_dimension cd
          ON ob.concept_cd = cd.concept_cd
          AND (cd.concept_cd LIKE "NDC:%" OR cd.concept_cd LIKE "RXCUI:%")
         AND (',.,')');
.drugregexps <- sapply(.drugnames,paste0,collapse='|');

# Local function to avoid a lot of repetitive code when searching for strings in a column
hasString <- function(xx) substitute(any(grepl(.drugregexps[xx],name_char,ignore.case=T)));

# All glucose-lowering drug ocurrences
dat0 <- dbGetQuery(fullsqlcon,.drugsql) %>%
  mutate(patient_num=as.character(patient_num),start_date=as.Date(start_date)) %>%
  group_by(patient_num,start_date) %>% arrange(patient_num,start_date);
dbDisconnect(fullsqlcon);

# patient history, can be joined to main data
dat1 <- summarise(dat0
                  ,Glinides=eval(hasString('Glinides'))
                  ,SGLT2I=eval(hasString('SGLT2I'))
                  ,DDP4I=eval(hasString('DDP4I'))
                  ,GLP1A=eval(hasString('GLP1A'))
                  ,TZD=eval(hasString('TZD'))
                  ,AnyNonSorM = Glinides|SGLT2I|DDP4I|GLP1A|TZD
                  ,Sulfonylureas=eval(hasString('Sulfonylureas'))
                  ,Metformin=eval(hasString('Metformin'))
);


export(dat1,file='DEID_GLUDRUGS.tsv');