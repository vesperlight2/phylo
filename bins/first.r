#genome rename and gene prediction for all methods
library(stringr)
library(seqinr)

arg1 <- commandArgs(T)
sample_list <- arg1[1]
out_dir <- str_c(arg1[2],"genomes/")

sample_tab <- read.csv(sample_list,row.names = 1)
spe <- sample_tab[,2]
ids <- sample_tab[,1]
sample_num <- length(ids)
file1 <- as.character(sample_tab[,4])
dir.create(str_c(out_dir,"genome_files_named_with_ids"))
dir.create(str_c(out_dir,"gene_prediction"))
dir.create(str_c(out_dir,"gene_prediction/protein_seqs"))
dir.create(str_c(out_dir,"gene_prediction/nucleotide_seqs"))  
for(i in seq(length(file1))){
  if(file1[i] %in% dir(out_dir)==FALSE){
    stop(str_c(i,"......",file1[i],"  does not exsit, please check sample list and sequence files!"))
  }
  file.copy(str_c(out_dir,file1[i]),str_c(out_dir,"genome_files_named_with_ids/"))
  file.rename(str_c(out_dir,"genome_files_named_with_ids/",file1[i]),str_c(out_dir,"genome_files_named_with_ids/",ids[i],".fasta"))
  seqs1 <- read.fasta(str_c(out_dir,"genome_files_named_with_ids/",ids[i],".fasta"))
  name1 <- as.list(paste(ids[i],seq(length(names(seqs1))),sep = "_"))
  write.fasta(sequences = seqs1,names = name1,file.out = str_c(out_dir,"genome_files_named_with_ids/",ids[i],".fasta"))
  print(str_c("Gene prediction[ ",file1[i], " ]...................",i,"/",length(file1)))
  system(str_c("prodigal -a ",out_dir,"gene_prediction/protein_seqs/aa_orf_",ids[i],".fasta -d ",out_dir,"gene_prediction/nucleotide_seqs/",ids[i],"_nucleotide.fasta -i ",out_dir,"genome_files_named_with_ids/", ids[i],".fasta -q >> ",out_dir,"gene_prediction/prodigal.md 2>&1"))
}
