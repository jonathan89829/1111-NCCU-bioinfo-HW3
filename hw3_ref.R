######################################
# the reference code of program3 
######################################

######################################
# initial
######################################
library("Biostrings",verbose=F,quietly=T)

# read parameters
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("USAGE: Rscript hw3_<your student ID>.R --input test.fasta --score pam250.txt --aln global --gap -10 --output test_output.fasta", call.=FALSE)
}
