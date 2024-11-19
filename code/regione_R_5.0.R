# Load necessary libraries
library(tidyr)         # For data tidying and manipulation
library(regioneR)      # For region permutation testing and analysis

# Capture command-line arguments
args = commandArgs(trailingOnly=TRUE)

# Ensure that the required arguments are provided
if (length(args) == 0) {
  stop("USAGE: regionR PathToBEDDirectory, PathToTargetDirectory, permutation ntime, output_dir", call.=FALSE)
}

# Uncomment these lines to accept paths and parameters as arguments
# pwd_bed <- as.string(args[1])           # Path to BED directory
# pwd_target <- as.string(args[2])        # Path to target region directory
# n = args[3]                             # Number of permutations
# pwd_outdir <- as.string(args[4])        # Path to output directory
# pwd_univ <- as.string(args[5])          # Path to universal background BED file

# Define the universal background BED file path (edit as needed)
pwd_univ = "/Users/ruofany/Downloads/seq/PermTest/background_100k/hg19genic100k.bed"

# Load the universal background BED file
AllGene <- read.delim(pwd_univ, header=FALSE)

# Define the paths to BED and target directories, number of permutations, and output directory
pwd_bed = "/Users/ruofany/Downloads/seq/PermTest/TSA_track_test"
pwd_target = "/Users/ruofany/Downloads/seq/PermTest/profile"
n = 3000    # Number of permutations
pwd_outdir <- "/Users/ruofany/Downloads/seq/PermTest/Results/05242023REDO"

# Create the output directory if it doesn't already exist
dir.create(pwd_outdir, showWarnings=FALSE)

# List all BED files in the specified BED directory
file_bed <- list.files(path=pwd_bed)

# Set the working directory to the BED directory
setwd(pwd_bed)

# Loop through each BED file
for (i in 1:length(file_bed)) {
    setwd(pwd_bed)
    bedname <- toString(file_bed[i])      # Get the name of the current BED file
    bedtf <- read.table(file_bed[i])     # Read the BED file into a data frame
    bedtf <- na.omit(bedtf)              # Remove any rows with NA values
    print(bedname)                       # Print the BED file name for tracking

    # Switch to the target directory
    setwd(pwd_target)
    file_target <- list.files(path=pwd_target, pattern="\\.bed$")  # List all target BED files

    # Loop through each target BED file
    for (j in 1:length(file_target)) {
        setwd(pwd_target)
        print(file_target[j])            # Print the target file name
        targetname <- toString(file_target[j])  # Get the target file name
        Target <- read.delim(file_target[j], header=FALSE)  # Read the target BED file

        # Perform permutation testing
        pt <- permTest(
            A=Target, 
            randomize.function=resampleRegions, 
            universe=AllGene, 
            ntime=n, 
            evaluate.function=meanInRegions, 
            x=bedtf
        )

        # Save the results
        setwd(pwd_outdir)
        # Save permutation test results to a text file
        cat(capture.output(print(pt), file=paste(targetname, bedname, n, ".txt")))

        # Save a plot of permutation test results as a PNG
        png(filename=paste(targetname, bedname, n, ".png"), width=4.2, height=5, units="in", res=1200)
        lapply(pt, plot.permTestResults)  # Plot results
        dev.off()
    }
}

# Switch to the output directory
setwd(pwd_outdir)

# Read all .txt files in the output directory
file_txt <- list.files(pattern=".*.txt")
df <- data.frame("sample", "pvalue", "zscore", "ntimes")  # Initialize a data frame for summary statistics

# Loop through each result file and extract statistics
for (i in 1:length(file_txt)) {
    x <- scan(file_txt[i], what="", sep="\n")   # Read the text file line-by-line
    nametxt <- toString(file_txt[i])           # Extract the file name
    x[1] <- nametxt                            # Replace the first line with the file name
    subx <- substring(x[1:4], regexpr(":", x[1:4]) + 1)  # Extract relevant statistics (e.g., p-value, z-score)
    df <- rbind(df, subx)                      # Append extracted statistics to the data frame
}

# Write the summary statistics to a file
write.table(df, file="stats.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
