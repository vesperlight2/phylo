library(stringr)
arg1 <- commandArgs(T)
sample_list <- arg1[1]
out_dir <- str_c(arg1[2],"single_copy_gene/")
phylophlan_dir <- arg1[3]
path <- arg1[4]

project <- rnorm(1)

sample_tab <- read.csv(sample_list,row.names = 1)
sample_list <- sample_tab[,1]
seq_list <- str_c("aa_orf_",sample_list,".fasta")
seq_name <- str_c("aa_orf_",sample_list)

system(str_c("rm -rf ",phylophlan_dir,"input/phylo_res",project))
system(str_c("mkdir ",phylophlan_dir,"input/phylo_res",project))

for(i in seq(length(sample_list))){
  system(str_c("cp ",arg1[2],"genomes/gene_prediction/protein_seqs/",seq_list[i]," ",phylophlan_dir,"input/phylo_res",project,"/",sample_list[i],".faa "))
}

system(str_c("sh ",path,"bins/phylophlan_run.sh phylo_res", project," ", phylophlan_dir))

#system(str_c("cd ",phylophlan_dir," && python phylophlan.py -c phylo_res",project))
#system(str_c("python phylophlan.py --nproc 10 -u phylo_res",project," >> ",out_dir,"phylo.log 2>&1"))
system(str_c("cp ",phylophlan_dir,"output/phylo_res",project,"/* ",out_dir))

