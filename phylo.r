#genome phylogenetic analysis tool kit 

#Version: 1.1

#Author: Bo Tu

#E-mail: tubo@caas.cn




#requirements
if(require(optparse)==F){
print("install require package: optparse")
install.packages("optparse")
}
if(require(stringr)==F){
print("install require package: stringr")
install.packages("stringr")
}
if(require(seqinr)==F){
print("install require package: seqinr")
install.packages("seqinr")
}

if(require(ape)==F){
print("install require package: ape")
install.packages("ape")
}

path1 <- "/home/sam/software/script/phylo/"
phylophlan_path <- readLines(str_c(path1,"config"))[4]

#options
option_list=list(
    make_option(c("-m","--method"),type="character",default=NULL,help="
	可使用的方法:
	cgf: core protein family method
	ribo: risomal method
	mlsa: MLSA method
	scg：single copy gene
	pocp: percetage of conserved protein
	ani:  anverage nucletide identity
	aai: anverage amino acid identity
	16s：16S identity
	std: for all genome-based analysis above except MLSA
	all: all analysis pipeline except MLSA"
	,metavar="character"),
    make_option(c("-s","--sample_list"),type="character",default=NULL,help="
	指定sample_list.csv文件，sample_list表格格式需要遵循示例文件规则
	No.	IDs	Species	Genome_IDs	Genome_files	16S_files
	1	A_o01_f01_g01_s01	Tepidanaerobacter acetatoxydans	GCF_000328765.2	GCA_000328765.2.fasta	
	2	A_o01_f02_g01_s01	Caldanaerovirga acetigignens	GCF_900142995.1	GCA_900142995.1.fasta	
",metavar="character"),
    make_option(c("-i","--input"),type="character",default=NULL,help="
	输入基因组文件所在位置",metavar="character"),
    make_option(c("-r","--ribo"),type="character",default=NULL,help="
	输入的16S序列文件位置，如果选择std方法或16s方法，该参数必须指定",metavar="character"),
    make_option(c("-o","--output"),type="character",default=NULL,help="
	输出文件位置",metavar="character"),
	make_option(c("-c","--cpus"),type="numeric",default=10,help="
	使用的cpu数目，默认10",metavar="INT")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if(is.null(opt$m)){
stop("Error!A method is needed.")
}

if(is.null(opt$s)){
stop("Error!A csv format sample list file is needed.")
}

if(is.null(opt$i)&opt$m!="16s"){
stop("Error!You need to specify the genome files or 16S sequence path")
}

if(is.null(opt$o)){
stop("Error!You need to specify the output direction")
}

if((opt$m=="16s")&is.null(opt$r)){
stop("Error!You need to specify the input 16S sequence file direction")
}


#dir format
if(opt$m!="16s"){
	if(str_sub(opt$i,-1,-1)!="/"){
		input_dir <- str_c(opt$i,"/")
			}else{
			input_dir <- opt$i
			}
}else if((opt$m=="16s")&(str_sub(opt$r,-1,-1)!="/")){
			input16s_dir <- str_c(opt$r,"/")
				}else if(opt$m=="16s"){
				input16s_dir <- opt$r
				}

if(str_sub(opt$o,-1,-1)!="/"){
output_dir <- str_c(opt$o,"/")
}else{output_dir <- opt$o}



#rename and gene prediction
if(opt$m!="16s"){
system(str_c("mkdir ",output_dir))
system(str_c("mkdir ",output_dir,"genomes"))
system(str_c("cp ",input_dir,"*.f* ",output_dir,"genomes/"))
system(str_c("Rscript ",path1,"bins/first.r ", opt$s, " ", output_dir))
}else{
system(str_c("mkdir ",output_dir))
}

#main
# core protein family
if(opt$m=="cgf"){
system(str_c("mkdir ",output_dir,"core_protein_fam/"))
system(str_c("Rscript ",path1,"bins/phylo_genome.r ", opt$s, " ", output_dir))
}

#scg
#single copy gene
if(opt$m=="scg"){
system(str_c("mkdir ",output_dir,"single_copy_gene"))
system(str_c("Rscript ",path1,"bins/scg.r ", opt$s, " ", output_dir," ",phylophlan_path," ",path1))
}

#ribo
if(opt$m=="ribo"){
system(str_c("mkdir ",output_dir))
system(str_c("mkdir ",output_dir,"ribo"))
system(str_c("mkdir ",output_dir,"ribo/ribo_hmm_files/"))
system(str_c("cp ",path1,"database/ribo_hmm_files/* ",output_dir,"ribo/ribo_hmm_files/"))
system(str_c("Rscript ",path1,"bins/phylo_ribo.r ", opt$s, " ", output_dir," ",output_dir,"genomes/gene_prediction/protein_seqs"))
}

# MLSA
#if(opt$m=="mlsa"&&is.null(opt$g)==F){
#system(str_c("mkdir ",output_dir))
#system(str_c("cp ",input_dir,"*fa* ",output_dir))
#system(str_c("Rscript ",path1,"/bins/phylo_genome.r ", opt$s, " ", output_dir))
#}

#pocp
if(opt$m=="pocp"){
system(str_c("mkdir ",output_dir))
system(str_c("mkdir ",output_dir,"pocp"))
system(str_c("Rscript ",path1,"bins/pocp_cal.r ", opt$s, " ", output_dir," ",output_dir,"genomes/gene_prediction/protein_seqs"))
}

#ani
if(opt$m=="ani"){
system(str_c("mkdir ",output_dir))
system(str_c("mkdir ",output_dir,"ani"))
system(str_c("Rscript ",path1,"/bins/ani_cal.r ", opt$s, " ", output_dir," ",output_dir,"genomes/genome_files_named_with_ids"))
}

#aai
if(opt$m=="aai"){
system(str_c("mkdir ",output_dir))
system(str_c("Rscript ",path1,"/bins/aai.r ", opt$s, " ", output_dir))
}

#16s identity
if(opt$m=="16s"){
system(str_c("mkdir ",output_dir,"16s_identity/"))
system(str_c("mkdir ",output_dir,"16s_identity/seqs"))
system(str_c("cp ",input16s_dir,"*.f* ",output_dir,"16s_identity/seqs"))
system(str_c("Rscript ",path1,"bins/16s.r ", opt$s, " ", output_dir))
}


#std
if(opt$m=="std"){
#cpf
system(str_c("mkdir ",output_dir,"core_protein_fam/"))
system(str_c("Rscript ",path1,"bins/phylo_genome.r ", opt$s, " ", output_dir))
#scg
system(str_c("mkdir ",output_dir,"single_copy_gene"))
system(str_c("Rscript ",path1,"bins/scg.r ", opt$s, " ", output_dir," ",phylophlan_path))
#ribo
system(str_c("mkdir ",output_dir,"ribo"))
system(str_c("mkdir ",output_dir,"ribo/ribo_hmm_files/"))
system(str_c("cp ",path1,"database/ribo_hmm_files/* ",output_dir,"ribo/ribo_hmm_files/"))
system(str_c("Rscript ",path1,"bins/phylo_ribo.r ", opt$s, " ", output_dir," ",output_dir,"genomes/gene_prediction/protein_seqs"))
#aai
system(str_c("mkdir ",output_dir))
system(str_c("Rscript ",path1,"/bins/aai.r ", opt$s, " ", output_dir))
#ani
system(str_c("mkdir ",output_dir,"ani"))
system(str_c("Rscript ",path1,"/bins/ani_cal.r ", opt$s, " ", output_dir," ",output_dir,"genomes/genome_files_named_with_ids"))
#pocp
system(str_c("mkdir ",output_dir,"pocp"))
system(str_c("Rscript ",path1,"bins/pocp_cal.r ", opt$s, " ", output_dir," ",output_dir,"genomes/gene_prediction/protein_seqs"))

aai <- read.csv(str_c(output_dir,"aai/aai_list.csv"),row.names=1)
ani <- read.csv(str_c(output_dir,"ani/ani_list.csv"),row.names=1)
pocp <- read.csv(str_c(output_dir,"pocp/pocp_list.csv"),row.names=1)
res <- cbind(aai,ani,pocp)
res1 <- res[,c(1,2,3,6,9)]
write.csv(res1,str_c(output_dir,"table_res.csv"))

aai1 <- read.csv(str_c(output_dir,"aai/aai_list_full.csv"),row.names=1)
ani1 <- read.csv(str_c(output_dir,"ani/ani_list_full.csv"),row.names=1)
pocp1 <- read.csv(str_c(output_dir,"pocp/pocp_list_full.csv"),row.names=1)
res2 <- cbind(aai1,ani1,pocp1)
res3 <- res2[,c(1,2,3,6,9)]
write.csv(res3,str_c(output_dir,"table_res_full.csv"))


print("===================ALL DONE!!!!!!!!!!!!===================================")
}


if(opt$m=="all"){
#cpf
system(str_c("mkdir ",output_dir,"core_protein_fam/"))
system(str_c("Rscript ",path1,"bins/phylo_genome.r ", opt$s, " ", output_dir))
#scg
system(str_c("mkdir ",output_dir,"single_copy_gene"))
system(str_c("Rscript ",path1,"bins/scg.r ", opt$s, " ", output_dir," ",phylophlan_path))
#ribo
system(str_c("mkdir ",output_dir,"ribo"))
system(str_c("mkdir ",output_dir,"ribo/ribo_hmm_files/"))
system(str_c("cp ",path1,"database/ribo_hmm_files/* ",output_dir,"ribo/ribo_hmm_files/"))
system(str_c("Rscript ",path1,"bins/phylo_ribo.r ", opt$s, " ", output_dir," ",output_dir,"genomes/gene_prediction/protein_seqs"))
#aai
system(str_c("mkdir ",output_dir))
system(str_c("Rscript ",path1,"/bins/aai.r ", opt$s, " ", output_dir))
#16s
system(str_c("mkdir ",output_dir,"16s_identity/"))
system(str_c("mkdir ",output_dir,"16s_identity/seqs"))
system(str_c("cp ",input16s_dir,"*.f* ",output_dir,"16s_identity/seqs"))
system(str_c("Rscript ",path1,"bins/16s.r ", opt$s, " ", output_dir))
#ani
system(str_c("mkdir ",output_dir,"ani"))
system(str_c("Rscript ",path1,"/bins/ani_cal.r ", opt$s, " ", output_dir," ",output_dir,"genomes/genome_files_named_with_ids"))
#pocp
system(str_c("mkdir ",output_dir,"pocp"))
system(str_c("Rscript ",path1,"bins/pocp_cal.r ", opt$s, " ", output_dir," ",output_dir,"genomes/gene_prediction/protein_seqs"))

#table combine

aai <- read.csv(str_c(output_dir,"aai/aai_list.csv"),row.names=1)
ani <- read.csv(str_c(output_dir,"ani/ani_list.csv"),row.names=1)
pocp <- read.csv(str_c(output_dir,"pocp/pocp_list.csv"),row.names=1)
s16 <- read.csv(str_c(output_dir,"16s_identity/pocp_list.csv"),row.names=1)
res <- cbind(aai,ani,pocp,s16)
res1 <- res[,c(1,2,3,6,9,12)]
write.csv(res1,str_c(output_dir,"table_res.csv"))

print("===================ALL DONE!!!!!!!!!!!!===================================")
}
