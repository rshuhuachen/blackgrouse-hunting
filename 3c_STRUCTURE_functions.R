#make directories
make_directories <- function(path){
  system(paste0("mkdir ", paste0(path, "Run_files")))
  system(paste0("mv ", paste0(path, paste0("*_f ", paste0(path, "Run_files")))))
  system(paste0("mkdir ", paste0(path, "Run_files/evanno")))
  system(paste0("mkdir ", paste0(path, "Run_files/outputclumpp")))
  system(paste0("mkdir ", paste0(path, "Run_files/struct-plots")))
  system(paste0("mkdir ", paste0(path, "Run_files/outputclumpp/collectedruns")))
  system(paste0("chmod ugo+rwx  ", paste0(path, "Run_files/outputclumpp")))
}


## load with clumpp
load_and_clumpp <- function(path_to_structure_out){
  
  all_files <- list.files(path_to_structure_out, pattern = "^results")
  
  # creates character vector containing path to relevant files
  struc_out_paths <- paste0(path_to_structure_out, all_files)
  
  slist <- readQ(files=struc_out_paths, filetype = "structure")
  
  # export for CLUMPP later on
  clumppExport(qlist=slist, exportpath= paste0(path_to_structure_out, "outputclumpp/"))# if CLUMPP takes to long set parammode to 3, otherwise leave out
  # collect CLUMPP output
  collectClumppOutput(filetype="merged", runsdir = paste0(path_to_structure_out, "outputclumpp/"),
                      newdir = paste0(path_to_structure_out, "outputclumpp/collectedruns/")) # aligned
  system("rm -r pop_K*")
  
  # move clump files to correct directory
  system(paste0("mv ", path_to_structure_out, "outputclumpp/pop_* ", path_to_structure_out))
  
}


#### Function to get K summary stats from run files ####

load_and_K <- function(path_to_structure_out){
  
  all_files <- list.files(path_to_structure_out, pattern = "^results")
  
  # creates character vector containing path to relevant files
  struc_out_paths <- paste0(path_to_structure_out, all_files)
  
  slist <- readQ(files=struc_out_paths, filetype = "structure")
  
  # Get summary stats
  
  em <- evannoMethodStructure(summariseQ(tabulateQ(slist)))
  plot(em$deltaK)
  
  evannoMethodStructure(data=em, exportplot=T, writetable = T, exportpath = paste0(path_to_structure_out, "outputclumpp"))
  system(paste0("mv ", path_to_structure_out, "outputclumpp/evannoMethodStructure* ", path_to_structure_out, "evanno/"))
  
  em <- em
  
}

### optimal k ###

optimal_k <- function(x) {
  elpd <- which.max(x$elpdmean)
  if (elpd == 1) {
    return(1)
  } else {
    delta_k <- which.max(x$deltaK)
    delta_k
  }
}


#### output best k's ####

both_k <- function(x) {
  elpd <- which.max(x$elpdmean)
  deltaK <- which.max(x$deltaK)
  ks <- c(lnk = x$k[elpd], deltak = x$k[deltaK])
  # write outfile
  #write.table(ks, paste0(path_to_structure_out, "Ks.txt"))
  ks <- ks
  ks
}
