library(stringr)
library(seqinr)

arg1 <- commandArgs(T)
sample_list <- arg1[1]
out_dir <- str_c(arg1[2],"pocp/")
seq_dir <- arg1[3]


if(str_sub(seq_dir,-1,-1)!="/"){
  seq_dir <- str_c(seq_dir,"/")
}

sample_tab <- read.csv(sample_list,row.names = 1)
list1 <- sample_tab[,1]
dir.create(str_c(out_dir,"blast_res"))
for(i in seq(length(list1))){
system(str_c("makeblastdb -in ",seq_dir,"aa_orf_",list1[i],".fasta -dbtype prot -out ",out_dir,list1[i]))
for(j in seq(length(list1))){
 system(str_c("blastp -query ",seq_dir,"aa_orf_", list1[j],".fasta -db ",out_dir,list1[i]," -outfmt 6 -evalue 1e-5 -num_threads 10 -out ",out_dir,"blast_res/", list1[i],".",list1[j],".txt"))
}
print(str_c(i,Sys.time()))
}
dis1 <- matrix(nrow=length(list1),ncol=length(list1))
lenall <- NULL
for(i in seq(list1)){
  seq1 <- read.fasta(str_c(seq_dir,"aa_orf_",list1[i],".fasta"))
  len1 <- lapply(seq1,length)
  lenall <- c(lenall,length(seq1))
  for(j in seq(list1)){
    a <- read.table(str_c(out_dir,"blast_res/",list1[j],".",list1[i],".txt"))
    a1 <- as.data.frame(cbind(len1,names(len1)))
    a2 <- merge(x=a,y=a1,by.x = "V1",by.y="V2",all.y = T)
    a2 <- cbind(a2[,1:12],as.numeric(a2[,13]))
    a3 <- a2[a2[,3]>40 & a2[,4] >= (0.5*a2[,13]),]
    a4 <- a3[!duplicated(a3[,1]),]
    dis1[j,i] <- nrow(a4)
  }
  print(i)
}

lenall <- as.numeric(lenall)
pocp <- matrix(nrow=length(list1),ncol=length(list1))
for(i in seq(list1)){
  for(j in seq(list1)){
    pocp[i,j] <- (dis1[i,j]+dis1[j,i])/(lenall[i]+lenall[j])
  }
}
rownames(pocp) <- colnames(pocp) <- list1
a3 <- pocp[lower.tri(pocp)]
query <- NULL
ref <- NULL
for(i in seq(from=1,to=length(list1))){
  if(i<length(list1)){
    ref <- c(ref,rownames(pocp)[(i+1):length(list1)])
    query <- c(query,rep(colnames(pocp)[i],(length(list1)-i))) 
  }
}
a4 <- cbind(query,ref,round(as.numeric(a3)*100,3))
colnames(a4) <- c("query","ref","POCP")

a5 <- NULL
query1 <- NULL
ref1 <- NULL

for(i in seq(from=1,to=length(list1))){
	a5 <- c(a5,pocp[i,])
	query1 <- c(query1,rep(rownames(pocp)[i],length(list1)))
	ref1 <- c(ref1,colnames(pocp))
}
a6 <- cbind(query1,ref1,as.numeric(a5)*100)
colnames(a6) <- c("query","ref","POCP")

write.csv(pocp,str_c(out_dir,"pocp_matrix.csv"))
write.csv(a6,str_c(out_dir,"pocp_list_full.csv"))
write.csv(a4,str_c(out_dir,"pocp_list.csv"))