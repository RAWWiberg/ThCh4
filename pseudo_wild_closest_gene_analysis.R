##
# D.pseudoobscura project wild linse
# Closest Genes
# Last Modified: April 2017
##
#clean environment
rm(list = ls(all = TRUE))
#
# Load libraries and define functions
# #####
library(plyr)
library(ggplot2)
library(scales)
#install.packages("grid")
library(gridExtra)
library(grid)
#install.packages("reshape")
library(reshape)
library(dplyr)
library(qvalue)
se <- function(x) {sqrt(var(x, na.rm = TRUE))/sqrt(length(x))}
source("~/Desktop/Data/RData/RScripts/ggplot_theme.R")
options(scipen=6)
# ####
setwd("~/Desktop/Data/RData/pseudo_wild/")
#


# Load Dpse annotation
dpse_ann<-read.table("~/Desktop/Data/pseudo_project/evol_lines/genomics/SNP_analysis/Dpse_genes_dmelnames.gtf",header = FALSE,sep="\t")
colnames(dpse_ann)<-c("chr","source","type","start","end",
                      "-","-","-","gene")
head(dpse_ann)
# Load the lists of SNPs
fixed_SNPs<-read.table(
  "/home/axel/Desktop/Data/pseudo_project/evol_lines/genomics/SNP_analysis/pseudo_evol_glm_top_bonf_snps.list",
  sep = "\t", 
  header=FALSE)

fixed_SNPs_reg<-read.table(
  "/home/axel/Desktop/Data/pseudo_project/evol_lines/genomics/SNP_analysis/pseudo_evol_glm_top_q-value_snps.list",
  sep = "\t", 
  header=FALSE)

#
#
#------------------------------------------#
# Look at distribution of SNPs in relation #
# to annotated genes.                      #
#------------------------------------------#
# Load Dpse annotation vs SNP gtf overlaps files
# ####
fixed_SNPs_genes<-read.table("/home/axel/Desktop/Data/pseudo_project/wild_lines/genomics/SNP_analysis/pseudo_wild_fixed_snps_closest_genes_reduced.tab",
                           header = FALSE, sep = "\t")
fixed_SNPs_genes<-fixed_SNPs_genes[,c(1,2,3,11,12)]
colnames(fixed_SNPs_genes) <- c("Chr","start","end","gene_id","dist")
fixed_SNPs_genes <- fixed_SNPs_genes[!(fixed_SNPs_genes$gene_id == "."),]
nrow(fixed_SNPs_genes)

fixed_SNPs_reg_genes<-read.table("/home/axel/Desktop/Data/pseudo_project/wild_lines/genomics/SNP_analysis/pseudo_wild_fixed_snps_reg_closest_genes_reduced.tab",
                        header = FALSE, sep = "\t")
fixed_SNPs_reg_genes<-fixed_SNPs_reg_genes[,c(1,2,3,11,12)]
colnames(fixed_SNPs_reg_genes) <- c("Chr","start","end","gene_id","dist")
fixed_SNPs_reg_genes <- fixed_SNPs_reg_genes[!(fixed_SNPs_reg_genes$gene_id == "."),]
nrow(fixed_SNPs_reg_genes)

all_genes <-read.table("/home/axel/Desktop/Data/pseudo_project/wild_lines/genomics/SNP_analysis/pseudo_wild_all_snps.tab",
                       header = FALSE, sep = "\t")
colnames(all_genes) <- c("SNP","gene_id","dist")
all_genes <- all_genes[
  !(all_genes$gene_id == "." | all_genes$gene_id == "-1"),]
nrow(all_genes)

# distribution of "absolute" distances from gene
fixed_SNPs_genes_dist_hist <- ggplot()+
  geom_histogram(data = fixed_SNPs_genes, aes(abs(dist)/1000))+
  xlab("|Distance from Gene| (kb)")+
  ylab("Count")
fixed_SNPs_genes_dist_hist + my.theme

fixed_SNPs_reg_genes_dist_hist <- ggplot()+
  geom_histogram(data = fixed_SNPs_reg_genes, aes(abs(dist)/1000))+
  xlab("|Distance from Gene| (kb)")+
  ylab("Count")
fixed_SNPs_reg_genes_dist_hist + my.theme

all_dist_hist <- ggplot()+
  geom_histogram(data = all_genes, aes(abs(dist)/1000))+
  xlab("|Distance from Gene| (kb)")+
  ylab("Count")
all_dist_hist + my.theme
# ####


# How many SNPs outside a gene
nrow(fixed_SNPs_genes[abs(fixed_SNPs_genes$dist) > 0,])
nrow(fixed_SNPs_genes[abs(fixed_SNPs_genes$dist) > 0,])/nrow(fixed_SNPs_genes)

nrow(fixed_SNPs_reg_genes[abs(fixed_SNPs_reg_genes$dist) > 0,])
nrow(fixed_SNPs_reg_genes[abs(fixed_SNPs_reg_genes$dist) > 0,])/nrow(fixed_SNPs_reg_genes)

nrow(all_genes[abs(all_genes$dist) > 0,])
nrow(all_genes[abs(all_genes$dist) > 0,])/nrow(all_genes)

# How many SNPs WITHIN coding region of a gene
nrow(fixed_SNPs_genes[fixed_SNPs_genes$dist == 0,])
nrow(fixed_SNPs_genes[fixed_SNPs_genes$dist == 0,])/nrow(fixed_SNPs_genes)

nrow(fixed_SNPs_reg_genes[fixed_SNPs_reg_genes$dist == 0,])
nrow(fixed_SNPs_reg_genes[fixed_SNPs_reg_genes$dist == 0,])/nrow(fixed_SNPs_reg_genes)

nrow(all_genes[all_genes$dist == 0,])
nrow(all_genes[all_genes$dist == 0,])/nrow(all_genes)

# Subset to include only genes with SNP within a coding region
fixed_SNPs_genes_withingene <- fixed_SNPs_genes[
  abs(fixed_SNPs_genes$dist) == 0,]

fixed_SNPs_reg_genes_withingene <- fixed_SNPs_reg_genes[
  abs(fixed_SNPs_reg_genes$dist) == 0,]

all_genes_withingene <- all_genes[
  abs(all_genes$dist) == 0,]
# get unique genes
length(unique(fixed_SNPs_genes_withingene$gene_id))

length(unique(fixed_SNPs_reg_genes_withingene$gene_id))

length(unique(all_genes_withingene$gene_id))


#------------------------------------------------------------#
# How many SNPs within XX (1kb or 1Mb) of a gene
#------------------------------------------------------------#
XX<-1000000
#------------------------------------------------------------#
# UP or DOWNSTREAM, 
# *includes* SNPs that occur WITHIN a gene
#------------------------------------------------------------#
nrow(fixed_SNPs_genes[
  abs(fixed_SNPs_genes$dist) <= XX,])
nrow(fixed_SNPs_genes[
  abs(fixed_SNPs_genes$dist) <= XX,])/nrow(fixed_SNPs_genes)

nrow(fixed_SNPs_reg_genes[
  abs(fixed_SNPs_reg_genes$dist) <= XX,])
nrow(fixed_SNPs_reg_genes[
  abs(fixed_SNPs_reg_genes$dist) <= XX,])/nrow(fixed_SNPs_reg_genes)

nrow(all_genes[abs(all_genes$dist) <= XX,])
nrow(all_genes[abs(all_genes$dist) <= XX,])/nrow(all_genes)

# Subset to include only genes with SNP within XX
# of a SNP *including* SNPs within gene region
fixed_SNPs_genes_XX <- fixed_SNPs_genes[
  abs(fixed_SNPs_genes$dist) <= XX,]

fixed_SNPs_reg_genes_XX <- fixed_SNPs_reg_genes[
  abs(fixed_SNPs_reg_genes$dist) <= XX,]

all_genes_XX <- all_genes[
  abs(all_genes$dist) <= XX,]

length(unique(fixed_SNPs_genes_XX$gene_id))

length(unique(fixed_SNPs_reg_genes_XX$gene_id))
gsub("_.","",gsub("FBGN","FBgn",unique(fixed_SNPs_reg_genes_XX$gene_id)))

length(unique(all_genes_XX$gene_id))

# write Gene IDs as lists
write.table(gsub("_.","",
                 gsub("FBGN","FBgn",
                      unique(fixed_SNPs_genes_XX$gene_id))),
            paste("/home/axel/Desktop/Data/pseudo_project/wild_lines/genomics/SNP_analysis/pseudo_wild_fixed_snps_closest_genes_",
                  XX,".list",
                  sep = ""),
            quote = FALSE,col.names=FALSE,row.names=FALSE)

write.table(gsub("_.","",
                 gsub("FBGN","FBgn",
                      unique(fixed_SNPs_reg_genes_XX$gene_id))),
            paste("/home/axel/Desktop/Data/pseudo_project/wild_lines/genomics/SNP_analysis/pseudo_wild_fixed_snps_reg_closest_genes_",
                  XX,".list",
                  sep = ""),
            quote = FALSE,col.names=FALSE,row.names=FALSE)


write.table(unique(all_genes_XX$gene_id),
            paste("/home/axel/Desktop/Data/pseudo_project/wild_lines/genomics/SNP_analysis/pseudo_wild_all_snps_closest_genes_",
                  XX,".list",
                  sep = ""),
            quote = FALSE,col.names=FALSE,row.names=FALSE)


#------------------------------------------------------------#
# UP or DOWNSTREAM, 
# *excludes* SNPs that occur WITHIN a gene
#------------------------------------------------------------#
nrow(fixed_SNPs_genes[
  abs(fixed_SNPs_genes$dist) <= XX & 
    fixed_SNPs_genes$dist != 0,])
nrow(fixed_SNPs_genes[
  abs(fixed_SNPs_genes$dist) <= XX &
    fixed_SNPs_genes$dist != 0,])/nrow(fixed_SNPs_genes)

nrow(fixed_SNPs_reg_genes[
  abs(fixed_SNPs_reg_genes$dist) <= XX & 
    fixed_SNPs_reg_genes$dist != 0,])
nrow(fixed_SNPs_reg_genes[
  abs(fixed_SNPs_reg_genes$dist) <= XX &
    fixed_SNPs_reg_genes$dist != 0,])/nrow(fixed_SNPs_reg_genes)

nrow(all_genes[
  abs(all_genes$dist) <= XX &
    all_genes$dist != 0,])
nrow(all_genes[
  abs(all_genes$dist) <= XX &
    all_genes$dist != 0,])/nrow(all_genes)

# Subset to include only genes with SNP within XX
# of a SNP *and* not within gene region.
fixed_SNPs_genes_XX_out <- fixed_SNPs_genes[
  abs(fixed_SNPs_genes$dist) <= XX & 
    fixed_SNPs_genes$dist != 0,]

fixed_SNPs_reg_genes_XX_out <- fixed_SNPs_reg_genes[
  abs(fixed_SNPs_reg_genes$dist) <= XX & 
    fixed_SNPs_reg_genes$dist != 0,]

all_genes_XX_out <- all_genes[
  abs(all_genes$dist) <= XX & 
    all_genes$dist != 0,]

length(unique(fixed_SNPs_genes_XX_out$gene_id))

length(unique(fixed_SNPs_reg_genes_XX_out$gene_id))

length(unique(all_genes_XX_out$gene_id))

# Write Gene IDs as lists
write.table(unique(fixed_SNPs_genes_XX_out$gene_id),
            paste("/home/axel/Desktop/Data/pseudo_project/evol_lines/genomics/SNP_analysis/pseudo_wild_fixe_SNPs_closest_genes_",
                  XX,"_out.list",
                  sep=""),
            quote = FALSE,col.names=FALSE,row.names=FALSE)

write.table(unique(fixed_SNPs_reg_genes_XX_out$gene_id),
            paste("/home/axel/Desktop/Data/pseudo_project/evol_lines/genomics/SNP_analysis/pseudo_wild_fixed_SNPs_reg_closest_genes_",
                  XX,"_out.list",
                  sep=""),
            quote = FALSE,col.names=FALSE,row.names=FALSE)

write.table(unique(all_genes_XX_out$gene_id),
            paste("/home/axel/Desktop/Data/pseudo_project/evol_lines/genomics/SNP_analysis/pseudo_wild_all_closest_genes_",
                  XX,"_out.list",
                  sep=""),
            quote = FALSE,col.names=FALSE,row.names=FALSE)
# Write SNPs as .gtf file

#FIXED SNPS
fixed_SNPs_genes_XX_out
fixed_SNPs_genes_XX_out_gtf <- cbind(
  as.character(gsub("chrom","",fixed_SNPs_genes_XX_out$Chr)),
  rep("PoPoolation2",nrow(fixed_SNPs_genes_XX_out)),
  rep("CDS",nrow(fixed_SNPs_genes_XX_out)),
  (fixed_SNPs_genes_XX_out$start-30),
  (fixed_SNPs_genes_XX_out$end+30),
  rep(".",nrow(fixed_SNPs_genes_XX_out)),
  rep(".",nrow(fixed_SNPs_genes_XX_out)),
  rep(".",nrow(fixed_SNPs_genes_XX_out)),
  paste(
    paste('gene_id "SNP_',
          seq(1,nrow(fixed_SNPs_genes_XX_out),1),
          '"',sep=""),
    paste('transcript_id "SNP_',
          seq(1,nrow(fixed_SNPs_genes_XX_out),1),
          '"',sep=""),
    paste('nearest_gene ',fixed_SNPs_genes_XX_out$gene_id),
    sep = "; "))
fixed_SNPs_genes_XX_out_gtf<-as.data.frame(fixed_SNPs_genes_XX_out_gtf)
fixed_SNPs_genes_XX_out_gtf$V4 <- as.integer(
  as.character(fixed_SNPs_genes_XX_out_gtf$V4)) 
fixed_SNPs_genes_XX_out_gtf$V5 <- as.integer(
  as.character(fixed_SNPs_genes_XX_out_gtf$V5))
head(fixed_SNPs_genes_XX_out_gtf)
write.table(fixed_SNPs_genes_XX_out_gtf, 
            file = paste(
              "/home/axel/Desktop/Data/pseudo_project/wild_lines/",
              "genomics/SNP_analysis/pseudo_wild_fixed_SNPs_out",
              XX,".gtf",sep=""),
            sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
#FIXED SNPS in REGIONS
fixed_SNPs_reg_genes_XX_out
fixed_SNPs_reg_genes_XX_out_gtf <- cbind(
  as.character(gsub("chrom","",fixed_SNPs_reg_genes_XX_out$Chr)),
  rep("PoPoolation2",nrow(fixed_SNPs_reg_genes_XX_out)),
  rep("CDS",nrow(fixed_SNPs_reg_genes_XX_out)),
  (fixed_SNPs_reg_genes_XX_out$start-30),
  (fixed_SNPs_reg_genes_XX_out$end+30),
  rep(".",nrow(fixed_SNPs_reg_genes_XX_out)),
  rep(".",nrow(fixed_SNPs_reg_genes_XX_out)),
  rep(".",nrow(fixed_SNPs_reg_genes_XX_out)),
  paste(
    paste('gene_id "SNP_',
          seq(1,nrow(fixed_SNPs_reg_genes_XX_out),1),
          '"',sep=""),
    paste('transcript_id "SNP_',
          seq(1,nrow(fixed_SNPs_reg_genes_XX_out),1),
          '"',sep=""),
    paste('nearest_gene ',fixed_SNPs_reg_genes_XX_out$gene_id),
    sep = "; "))
fixed_SNPs_reg_genes_XX_out_gtf<-as.data.frame(fixed_SNPs_reg_genes_XX_out_gtf)
fixed_SNPs_reg_genes_XX_out_gtf$V4 <- as.integer(
  as.character(fixed_SNPs_reg_genes_XX_out_gtf$V4)) 
fixed_SNPs_reg_genes_XX_out_gtf$V5 <- as.integer(
  as.character(fixed_SNPs_reg_genes_XX_out_gtf$V5))
head(fixed_SNPs_reg_genes_XX_out_gtf)
write.table(fixed_SNPs_reg_genes_XX_out_gtf, 
            file = paste(
              "/home/axel/Desktop/Data/pseudo_project/wild_lines/",
              "genomics/SNP_analysis/pseudo_wild_fixed_SNPs_reg_out",
              XX,".gtf",sep=""),
            sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

#------------------------------------------------------------#
# How many SNPs within XX UPSTREAM of a gene
# *Excludes* SNPs that occur WITHIN a gene
#------------------------------------------------------------#
nrow(top_bonf_genes[top_bonf_genes$dist > -XX & 
                      top_bonf_genes$dist < 0,])
nrow(top_bonf_genes[top_bonf_genes$dist > -XX & 
                      top_bonf_genes$dist < 0,])/nrow(top_bonf_genes)

nrow(top_q_genes[top_q_genes$dist > -XX & 
                   top_q_genes$dist < 0,])
nrow(top_q_genes[top_q_genes$dist > -XX & 
                   top_q_genes$dist < 0,])/nrow(top_q_genes)

nrow(all_genes[all_genes$dist > -XX & 
                 all_genes$dist < 0,])
nrow(all_genes[all_genes$dist > -XX & 
                 all_genes$dist < 0,])/nrow(all_genes)

nrow(top_bonf_genes[top_bonf_genes$dist > -1000,])
nrow(top_bonf_genes[top_bonf_genes$dist > -1000,])/nrow(top_bonf_genes)
# Subset to include only genes with SNP within XX
# of a SNP
top_bonf_genes_XX_ups <- top_bonf_genes[
  top_bonf_genes$dist > -XX & 
    top_bonf_genes$dist < 0,]

top_q_genes_XX_ups <- top_q_genes[
  top_q_genes$dist > -XX & 
    top_q_genes$dist < 0,]


