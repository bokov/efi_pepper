#' ## Script for deploying newly created tables
#'
#' The scripts in this project _write_ data to the project directory but
#' _read_ data from the locations specified in the `inputdata` variable by
#' `default_config.R` and `local_config.R` if it exists. This is on purpose as
#' a safety speedbump. The assumption is that after reviewing the updated tables
#' and confirming that they are okay, the researcher can run this script to
#' create a backup copy of old data files that will be overwritten, copy new
#' files to the PHI and LDS folder locations OVERWRITING OLD VERSIONS IF ANY,
#' and copy those same files to a staging directory for deployment to the server
#'
#' NOTE: This script assumes that all PHI files are in one folder and all LDS
#'       files are in another. Furthermore, it assumes that the PHI and LDS
#'       folders are both in the same parent folder.
#' NOTE:
library(dplyr);

# init ----
source('default_config.R');
#' The local path names for the data files should be stored in a vector
#' named `inputdata` that get set in a script named `local_config.R`
if(file.exists('local_config.R')) source('local_config.R');
#.filecopy_dryrun <- FALSE;
basedir <- unique(dirname(dirname(inputdata)))[1];
# stagedir <- file.path('/tmp/efi_stage');
# backupdir <- file.path(paste0('/tmp/efi_old_',as.numeric(Sys.time())));
# ldsdir <- '01 LDS Derived Data 220624';
# phidir <- '01 PHI Derived Data 220624';
dir.create(file.path(backupdir,ldsdir),recursive=T,showWarnings = F);
dir.create(file.path(backupdir,phidir),recursive=T,showWarnings = F);
dir.create(file.path(stagedir,ldsdir),recursive=T,showWarnings = F);
dir.create(file.path(stagedir,phidir),recursive=T,showWarnings = F);

# Every created file has to have one of these prefixes.
new_files <- list.files(pattern='^DEID_|^ID_|^PHI_|^consort_');
message('New files:\n',paste0(new_files,collapse='\n'));
message('They will be placed in: ',stagedir);
message('They will also be merged into: ',basedir);
message('Any old versions will be backed up to: ',backupdir);

# TODO: zip the large files

# File compare ----
# See which files are updates of existing ones (changed_files) and which
# are completely new (new_files)
old_files <- list.files(basedir,recursive = T,full.names=T);
changed_files <- intersect(list.files(),basename(old_files));
new_files <- setdiff(new_files,changed_files);
stage_files <- c();

# old files from basedir to backupdir  ----
# For changed_files, we copy them from the basedir to the backupdir in
# anticipation of basedir later being (manually) overwritten with the updated
# files we are about to copy to stagedir. So here we collect the names and paths
# of files to back up.
for(xx in old_files){
  if(basename(xx) %in% changed_files){
    if(!file.exists(xx)||tools::md5sum(xx)!=tools::md5sum(basename(xx))){
      stage_files <- rbind(stage_files,c(basename(xx),xx));
    }
  }
}
for(xx in new_files){
  .xxdest <- if(grepl('^(DEID_|consort_)',xx)) ldsdir else phidir;
  stage_files <- rbind(stage_files,c(xx,file.path(basedir,.xxdest,xx)));
}
# back up the old files
.junk<-sapply(sprintf('Copying %s to %s',stage_files[,2],
                      backup_targets<-sub(basedir,backupdir,stage_files[,2]))
              ,message);
if(!.filecopy_dryrun) file.copy(stage_files[,2],backup_targets,overwrite = T,copy.date = T);

# copy the new files to the staging dir ----
.junk<-sapply(sprintf('Copying %s to %s',stage_files[,1],
                      stage_targets<-sub(basedir,stagedir,stage_files[,2]))
              ,message);
if(!.filecopy_dryrun) file.copy(stage_files[,1],stage_targets,overwrite = T,copy.date = T);

# overwrite old files in basedir ----
# Also overwrite the old files (in basedir) with the new files. This is why we
# backed them up. Why did we bother copying them to stagedir if we're updating
# basedir directly from this script? So that we can also copy stagedir to the
# server, for example, or see exactly what got updated.
.junk<-sapply(sprintf('Copying %s to %s',stage_files[,1],stage_files[,2])
              ,message);
if(!.filecopy_dryrun) file.copy(stage_files[,1],stage_files[,2],overwrite = T,copy.date = T);

# git hash ----
.thisgithash <- system('git rev-parse --verify HEAD',intern=T);
if(system('git',intern=F)==1 &&
   system('if [ -z "$(git status --porcelain=v1 2>/dev/null)" ] ; then echo "No changes"; else echo "Changes"; fi',intern=T) == 'No changes' &&
   system('git ls-remote --head --exit-code origin main | cut -f 1 ',intern=T)== .thisgithash){
  .githashfile <- file.path(stagedir,paste0('git_hash_',as.numeric(Sys.time()),'.txt'));
  message('Writing value ',.thisgithash,' to ',.githashfile);
  if(!.filecopy_dryrun) write(.thisgithash,file=.githashfile);
}
