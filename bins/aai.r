library(stringr)
arg1 <- commandArgs(T)
sample_list <- arg1[1]
out_dir <- str_c(arg1[2],"aai/")

sample_tab <- read.csv(sample_list,row.names = 1)
sample_list <- sample_tab[,1]

system(str_c("comparem aai_wf -x fasta -c 10 ",arg1[2],"genomes/genome_files_named_with_ids ",out_dir))
aai <- read.csv(str_c(out_dir,"aai/aai_summary.tsv"),sep = "\t")
aai1 <- as.data.frame(aai[aai[,1]%in%sample_list & aai[,3]%in%sample_list,c(1,3,6)])
aai2 <- cbind(c(as.character(aai1[,1]),as.character(aai1[,2])),c(as.character(aai1[,2]),as.character(aai1[,1])),c(aai1[,3],aai1[,3]))
a2 <- matrix(ncol=length(sample_list),nrow=length(sample_list))
rownames(a2) <- colnames(a2) <- sample_list

for(i in seq(length(sample_list))){
  for(j in seq(length(sample_list))){
    if(i==j){
      a2[i,j] <- "100" 
    } else{
      a2[i,j] <- a2[j,i] <- aai2[aai2[,1] %in% rownames(a2)[i] & aai2[,2] %in% colnames(a2)[j],][3]
    }
  }
}
a3 <- a2[lower.tri(a2)]
query <- NULL
ref <- NULL
for(i in seq(from=1,to=length(sample_list))){
  if(i<length(sample_list)){
    ref <- c(ref,rownames(a2)[(i+1):length(sample_list)])
    query <- c(query,rep(colnames(a2)[i],(length(sample_list)-i))) 
  }
}
a4 <- cbind(query,ref,a3)
colnames(a4) <- c("query","ref","AAI")

a5 <- NULL
query1 <- NULL
ref1 <- NULL

for(i in seq(from=1,to=length(sample_list))){
	a5 <- c(a5,a2[i,])
	query1 <- c(query1,rep(rownames(a2)[i],length(sample_list)))
	ref1 <- c(ref1,colnames(a2))
}
a6 <- cbind(query1,ref1,as.numeric(a5))
colnames(a6) <- c("query","ref","AAI")

write.csv(a4,str_c(out_dir,"aai_list.csv"))
write.csv(a6,str_c(out_dir,"aai_list_full.csv"))
write.csv(a2,str_c(out_dir,"aai_matrix.csv"))
