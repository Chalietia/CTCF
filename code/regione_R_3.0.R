## ---------------------------
##
## Script name:regione_R_3.0.R
##
## Loop through all track file (bed format, 20kb size bin) and target file (bed format)
##
## Author: RY
##
## Date Created: 2021-04-14

library(tidyr)
library(regioneR)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("USAGE: regionR PathToBEDDirectory, PathToTargetDirectory, permutation ntime, output_dir", call.=FALSE)
}


AllGene <- read.delim ("/Users/ruofany/Downloads/seq/gr2020/Allgene.bed.100k.bed", header=FALSE)

pwd_bed <- as.string (args[1])
pwd_target <- as.string (args[2])
n = args[3]
pwd_outdir <- as.string (args[4])


#pwd_bed = "/Users/ruofany/Downloads/seq/PermTest/SON_TSA_track"
#pwd_target = "/Users/ruofany/Downloads/seq/PermTest/Target_profile_100k"
#n=1000
#pwd_outdir <- "/Users/ruofany/Downloads/seq/PermTest/Results/100kb_04142021"

dir.create(pwd_outdir)
file_bed <- list.files(path= pwd_bed)
setwd(pwd_bed)

#looping through all bed and target files.
for (i in 1:length(file_bed))
  {
    setwd(pwd_bed)
    bedname<- toString (file_bed[i])
    bedtf <- read.table (file_bed[i])
    print(bedname)
    setwd(pwd_target)
    file_target <- list.files(path=pwd_target)
    
    for (j in 1:length(file_target))
      {
      setwd(pwd_target)
      print(file_target[j])
      targetname<- toString(file_target[j])
      Target <- read.delim(file_target[j], header=FALSE)
      pt<- permTest(A= Target, randomize.function=resampleRegions, universe= AllGene,ntime= n, evaluate.function=meanInRegions, x= bedtf)
      setwd(pwd_outdir)
      cat(capture.output(print(pt), file=paste (targetname,bedname,n,".txt")))
      png (filename = paste (targetname,bedname,n,".png"))
      plot(pt)
      dev.off()
      }
    
    }




