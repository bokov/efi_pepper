#' ## Script for deploying newly created tables
#'
#' The scripts in this project _write_ data to the project directory but
#' _read_ data from the locations specified in the `inputdata` variable by
#' `default_config.R` and `local_config.R` if it exists. This is on purpose as
#' a safety speedbump. The assumption is that after reviewing the updated tables
#' and confirming that they are okay, the researcher can run this script to
#' create a backup copy of old data files that will be overwritten, copy new
#' files to the PHI and LDS folder locations overwriting old versions if any,
#' and copy those same files to a staging directory for deployment to the server
#'
#' NOTE: This script assumes that all PHI files are in one folder and all LDS
#'       files are in another. Furthermore, it assumes that the PHI and LDS
#'       folders are both in the same parent folder.
#' NOTE:
library(dplyr);

source('default_config.R');
#' The local path names for the data files should be stored in a vector
#' named `inputdata` that get set in a script named `local_config.R`
if(file.exists('local_config.R')) source('local_config.R');
basedir <- unique(dirname(dirname(inputdata)))[1];
stagedir <- file.path('/tmp/efi_stage');
backupdir <- file.path(paste0('/tmp/efi_old_',as.numeric(Sys.time())));
ldsdir <- '01 LDS Derived Data 220624';
phidir <- '01 PHI Derived Data 220624';
dir.create(file.path(backupdir,ldsdir),recursive=T,showWarnings = F);
dir.create(file.path(backupdir,phidir),recursive=T,showWarnings = F);
dir.create(file.path(stagedir,ldsdir),recursive=T,showWarnings = F);
dir.create(file.path(stagedir,phidir),recursive=T,showWarnings = F);
new_files <- list.files(pattern='^DEID_|^ID_|^PHI_');
# zip the large files
old_files <- list.files(basedir,recursive = T,full.names=T);
changed_files <- intersect(list.files(),basename(old_files));
new_files <- setdiff(new_files,changed_files);
stage_files <- c()
for(xx in old_files){
  if(basename(xx) %in% changed_files){
    if(!file.exists(xx)||tools::md5sum(xx)!=tools::md5sum(basename(xx))){
      stage_files <- rbind(stage_files,c(basename(xx),xx));
    }
  }
}
for(xx in new_files){
  .xxdest <- if(grepl('^DEID_',xx)) ldsdir else phidir;
  stage_files <- rbind(stage_files,c(xx,file.path(basedir,.xxdest,xx)));
}
# back up the old files
file.copy(stage_files[,2],sub(basedir,backupdir,stage_files[,2]),overwrite = T,copy.date = T);
# copy the new files to the staging dir
file.copy(stage_files[,1],sub(basedir,stagedir,stage_files[,2]),overwrite = T,copy.date = T)
# also overwrite the old files with the new files
file.copy(stage_files[,1],stage_files[,2],overwrite = T,copy.date = T);
