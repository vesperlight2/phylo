#16s identity

#import library
library(stringr)
library(seqinr)

#arguements
arg1 <- commandArgs(T)
sample_list <- arg1[1]
out_dir <- str_c(arg1[2],"16s_identity/")


#read sample table
sample_tab <- read.csv(sample_list,row.names = 1)
list1 <- sample_tab[,1]
dir.create(str_c(out_dir,"renamed_16s"))
dir.create(str_c(out_dir,"blast_res"))
file1 <- as.character(sample_tab[,5])

#rename
for(i in seq(length(file1))){
  file.rename(str_c(out_dir,"seqs/",file1[i]),str_c(out_dir,"renamed_16s/",list1[i],".fasta"))
}

a1 <- matrix(nrow=length(list1),ncol=length(list1))
#pairwised 16s blastn
for(i in seq(length(list1))){
  system(str_c("makeblastdb -in ",out_dir,"renamed_16s/",list1[i],".fasta -dbtype nucl -out ",out_dir,list1[i]))
  for(j in seq(length(list1))){
    system(str_c("blastn -query ",out_dir,"renamed_16s/", list1[j],".fasta -db ",out_dir,list1[i]," -outfmt 6 -num_threads 10 -out ",out_dir,"blast_res/", list1[i],".",list1[j],".txt"))
    if(file.size(str_c(out_dir,"blast_res/", list1[i],".",list1[j],".txt"))>0){
      a <- read.table(str_c(out_dir,"blast_res/",list1[i],".",list1[j],".txt"))
      a1[i,j] <- a[1,3]
    }else{
      a1[i,j] <- 0
    }
    }
  print(str_c(i,"/",length(list1),"..........",Sys.time()))
}


rownames(a1) <- colnames(a1) <- list1
a3 <- a1[lower.tri(a1)]
query <- NULL
ref <- NULL
for(i in seq(from=1,to=length(list1))){
  if(i<length(list1)){
    ref <- c(ref,rownames(a1)[(i+1):length(list1)])
    query <- c(query,rep(colnames(a1)[i],(length(list1)-i))) 
  }
}
a4 <- cbind(query,ref,a3)
colnames(a4) <- c("query","ref","16S_identity")

write.csv(a1,str_c(out_dir,"16s_matrix.csv"))
write.csv(a4,str_c(out_dir,"16s_list.csv"))
