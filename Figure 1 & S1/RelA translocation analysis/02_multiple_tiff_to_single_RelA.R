library(ijtiff)
library(dplyr)
library(stringr)
library(parallel)

# Folder containing all images
PATH = "~"
# File ending of unmerged images
ENDING = "-zFinal.tif"

# Search in the given folder for all images
Files <- list.files(PATH, full.names = T, recursive = T, include.dirs = F)

# Make table
FileList <- NULL
FileList$Path <- Files
FileList <- as_tibble(FileList)

# Make images list
FileList <-
  FileList %>%
  filter(
    grepl(ENDING, Path)
  ) %>% 
  mutate(
    Image = str_sub(Path, end = -nchar(ENDING)-5),
    Image = paste0(Image, ".tif"),
    Channel = str_sub(Path, nchar(Image)),
    Channel = str_sub(Channel, end=1)
  ) 

# Pull and merge images
CombineFx <- function(ImageX){
  # Pull images list
  TempFileList <-
    FileList %>% 
    filter(
      Image == Images[ImageX]
    )
  
  # Pull images
  Imgs <- lapply(TempFileList$Path, read_tif)
  
  # Make it a 4D array
  Imgs <- array(unlist(Imgs), dim = c(NROW(Imgs[[1]]), NCOL(Imgs[[1]]), NROW(TempFileList), count_frames(TempFileList$Path[1])[1]))
  
  # Save image
  write_tif(Imgs, TempFileList$Image[1], overwrite = T)
  
}
# Loop
Images <- unique(FileList$Image)
mclapply(1:NROW(Images), CombineFx)

