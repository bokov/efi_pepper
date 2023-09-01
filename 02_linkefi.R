library(dplyr);
library(tidyr);
#library(DBI);
library(rio);     # format-agnostic convenient file import
library(dint);    # date-conversion
library(stringr);

source('default_config.R');
#' The local path names for the data files should be stored in a vector
#' named `inputdata` that get set in a script named `local_config.R`
if(file.exists('local_config.R')) source('local_config.R');

# If developing the datasources themselves, give preference to locally generated
# versions
if(.devmode){
  inputdata <- ifelse(file.exists(basename(inputdata))
                      ,basename(inputdata),inputdata) %>%
    setNames(names(inputdata))};


savename <- gsub('[A-Za-z]+_[A-Za-z]+','DEID_EFI',basename(inputdata['samplecsv'])); #%>% gsub('.csv.zip','.tsv',.);
#' for automated rebuilding of some or all data files, create a file named
#' .usecachedfiles in the project folder
if(file.exists('.usecachedfiles')){
  inputdata <- sapply(inputdata,function(xx) if(file.exists(basename(xx))) basename(xx) else xx);
}

#' Import the main data and the EFI mappings
dat0 <- import(inputdata['samplecsv'],colClasses=cClasses);

#' Columns containing nested data that takes up a lot of memory and isn't really
#' analysis-ready
json_cols <- names(dat0) %>% grep('_cd$|_mn$|_tf$',.,inv=T,val=T) %>% grep('^v[0-9]',.,val=T);

efi <- import(inputdata['efixwalk'],colClasses=cClasses);
#' Import various de-identified data elements to link
hba1c <- import(inputdata['hba1c'],colClasses=cClasses);

gludrugs <- import(inputdata['gludrugs'],colClasses=cClasses) %>%
  select(patient_num,start_month,MonthFactor,CohortFactor,CohortDetail,months_since_pcvisit) %>%
  unique() %>%
  group_by(patient_num,start_month,CohortFactor,CohortDetail) %>%
  summarize(MonthFactor = MonthFactor[which.max(str_length(MonthFactor))]
            ,months_since_pcvisit = max(months_since_pcvisit,na.rm=T)) %>%
  mutate(gludrugs_begin=min(start_month),gludrugs_end=max(start_month)) #%>%
  #subset(is.finite(months_since_pcvisit));

#' Filter out the patients who don't have valid EFIs within the range of the data
#' TODO: merge records for the two Epic IDs that each have two patient_num
#' DONE: how many patients lost due to dropping infinite months_since_pcvisit -- 88
#' DONE: 1140622 total rows in gludrugs and 1140622 unique patient-date combos
#' DONE: patients in gludrugs vs in dat1 -- 8966 vs 2523 and all dat1 patients are in gludrugs
#' TODO: Document the deduplication effort and counts of CohortFactors with and
#'       without exclusion due to too few PC visits
#' DONE: Does it get better if we only apply the PC visit exclusion for
#'       encounters in 2016-2019? Yes. Do that.
dat1 <- subset(dat0,patient_num %in% efi$patient_num) %>%
  mutate(start_month=first_of_month(start_date));
#' But for the rest, keep all available dates because some of them contain special data elements
dat2 <- left_join(dat1,efi) %>%
  left_join(gludrugs) %>%
  left_join(hba1c[,c('patient_num','start_date','medhba1c','vfhba1c')]) %>%
  fill(medhba1c,vfhba1c,FRAIL6MO,FRAIL12MO,FRAIL24MO) %>%
  select(!any_of(idfields)) %>%
  # this is the fix for missing CohortFactor
  subset(between(start_month,gludrugs_begin,gludrugs_end));
dat3 <- select(dat2,!any_of(json_cols));

message('Saving full data as ',savename);
export(dat2,file=savename);
message('Saving analytic-only data as ',nojsonsavename <- gsub('^DEID_EFI_','DEID_EFI_NOJSON_',savename));
export(dat3,file=nojsonsavename);

