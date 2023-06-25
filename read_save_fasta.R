library("Biostrings",verbose=F,quietly=T)

# read fasta file
ff <- readAAStringSet("test.fasta")
seq_name = names(ff)
sequence = paste(ff)
print(sequence[1])
print(sequence[2])


# save fasta file
ff[[1]] <- AAString(x=seq1, start=1, nchar=NA)
ff[[2]] <- AAString(x=seq2, start=1, nchar=NA)
writeXStringSet(ff, "output_path")
