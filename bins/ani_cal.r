library(stringr)
arg1 <- commandArgs(T)
sample_list <- arg1[1]
out_dir <- str_c(arg1[2],"ani/")
seq_dir <- arg1[3]

if(str_sub(seq_dir,-1,-1)!="/"){
  seq_dir <- str_c(seq_dir,"/")
}

sample_tab <- read.csv(sample_list,row.names = 1)
f1 <- sample_tab[,1]
a1 <- NULL
a2 <- matrix(ncol=length(f1),nrow=length(f1))
dir.create(str_c(out_dir,"res"))
for(i in seq(length(f1))){
  for(j in seq(length(f1))){
    system(str_c("java -jar /home/sam/software/OAT_cmd.jar -fasta1 ",seq_dir,f1[i],".fasta -fasta2 ",seq_dir,f1[j],
                 ".fasta -method ani -blastplus_dir /usr/bin | cat > ",out_dir,"res/",f1[i],"_",f1[j],".txt"))
    a <- readLines(str_c(out_dir,"res/",f1[i],"_",f1[j],".txt"))
    a <- str_sub(a[length(a)-1],12,-5)
    a1 <- rbind(a1,c(as.character(f1[i]),as.character(f1[j]),a))
    a2[i,j] <- a
  }
  print(str_c(i,Sys.time()))
}  
colnames(a2) <- rownames(a2) <- f1
a3 <- a2[lower.tri(a2)]
query <- NULL
ref <- NULL
for(i in seq(from=1,to=length(f1))){
  if(i<length(f1)){
    ref <- c(ref,rownames(a2)[(i+1):length(f1)])
    query <- c(query,rep(colnames(a2)[i],(length(f1)-i))) 
  }
}
a4 <- cbind(query,ref,a3)
colnames(a4) <- c("query","ref","ANI")

a5 <- NULL
query1 <- NULL
ref1 <- NULL

for(i in seq(from=1,to=length(f1))){
	a5 <- c(a5,a2[i,])
	query1 <- c(query1,rep(rownames(a2)[i],length(f1)))
	ref1 <- c(ref1,colnames(a2))
}
a6 <- cbind(query1,ref1,as.numeric(a5))
colnames(a6) <- c("query","ref","ANI")


write.csv(a4,str_c(out_dir,"ani_list.csv"))
write.csv(a6,str_c(out_dir,"ani_list_full.csv"))
write.csv(a2,str_c(out_dir,"ani_matrix.csv"))
