##
# D.pseudoobscura project wild linse
# PoPoolation2 allele frequencies data analysis
# Last Modified: March 2017
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
se <- function(x) {sqrt(var(x, na.rm = TRUE))/sqrt(length(x))}
setwd("~/Desktop/Data/RData/pseudo_wild/")
source("../RScripts/ggplot_theme.R")
options(scipen = 6)

# Read in data
alfreq<-read.table("pseudo_wild_SHA-SLOB7subs_trcomp_pwc_red_chr.tab",
                  sep = "\t")
#colnames(alfreq)<-c(names(alfreq)[1:8],"R1","R2","R3","R4")
colnames(alfreq)<-c("chrom","pos","ref","alleles","states",
                    "V6","type","maj","SLo","Lew","Sha")
alfreq$SLo<-as.numeric(as.character(alfreq$SLo))
alfreq$Lew<-as.numeric(as.character(alfreq$Lew))
alfreq$Sha<-as.numeric(as.character(alfreq$Sha))
str(alfreq)
names(alfreq)
head(alfreq)
tail(alfreq)

# How many multiallelic
nrow(alfreq[alfreq$alleles > 2,])
nrow(alfreq[alfreq$alleles == 2,])
# Subset to biallelic
alfreq<-alfreq[alfreq$alleles == 2,]
# How many sites
nrow(alfreq)
# How many sites show fixed differences in all comparisons? (diff = 1)
nrow(alfreq[which(alfreq$SLo == 1 & 
              alfreq$Lew == 1 & 
              alfreq$Sha == 1),])
# As a percentage of all SNPs
(nrow(alfreq[which(alfreq$SLo == 1 & 
                    alfreq$Lew == 1 & 
                    alfreq$Sha == 1),])/nrow(alfreq))*100

head(alfreq[which(alfreq$SLo == 1 & 
              alfreq$Lew == 1 & 
              alfreq$Sha == 1),])

# Subset to those SNPs that are fixed in all populations
alfreq_fixed<-alfreq[which(alfreq$SLo == 1 & 
               alfreq$Lew == 1 & 
               alfreq$Sha == 1),]
alfreq_fixed$Chr2<-gsub("chrom","",alfreq_fixed$chrom)
alfreq_fixed$Chr2<-gsub("_group",".",alfreq_fixed$Chr2)
alfreq_fixed$Chr3<-gsub("\\..*","",alfreq_fixed$Chr2)
head(alfreq_fixed)
nrow(alfreq_fixed)
nrow(alfreq)
tapply(alfreq_fixed$chrom,INDEX = list(alfreq_fixed$Chr3),length)/nrow(alfreq_fixed)

# Plot the fixed SNPs
# Regions from tom price
reg<-read.table("price_BC_remating_regions.tab", header = TRUE, sep = "\t")
reg$Chr2<-gsub("chrom","",reg$chr)
reg$Chr2<-gsub("_group",".",reg$Chr)
# Chromosome lengths
main_chroms <- read.table(file = "/home/axel/Desktop/Data/pseudo_project/reference/dpse_r31_FB2013_02/pseudo_main_chrom_data.tab",
                          header = TRUE, sep = ",")
main_chroms$Chr2<-gsub("_group",".",main_chroms$Chr)
head(main_chroms)
# Make an ideogrammatic chromosome plot with fixed SNPs
ggplot()+
  geom_rect(data=main_chroms,aes(xmin=0,xmax=length,ymin=0,ymax=1),
            fill="white",colour="black")+
  geom_vline(data=alfreq_fixed,aes(xintercept = pos))+
  geom_rect(data=reg,aes(xmin = start, xmax=end,ymin=0,ymax=1),
            fill = "red",alpha = 1/2)+
  geom_vline(data=reg,aes(xintercept = start),
            colour = "red",linetype="dashed")+
  geom_vline(data=reg,aes(xintercept = end),
             colour = "red",linetype="dashed")+
  ylim(0,1)+
  facet_grid(Chr2~.)+
  my.theme+
  theme(panel.background = element_rect(fill="grey"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        strip.text.y = element_text(angle = 0,size = 10))

tapply(alfreq_fixed$chrom,INDEX = list(alfreq_fixed$Chr3),length)

# How many sites show at least one fixation
nrow(alfreq[alfreq$SLo == 1 | 
              alfreq$Lew == 1 | 
              alfreq$Sha == 1,])

# Count the number of fixed differences in windows
# Load chromosome info
chrom_dat<-read.table("~/Desktop/Data/pseudo_project/reference/dpse_r31_FB2013_02/pseudo_chrom_data.tab",
                      header = TRUE, sep = ",")

head(chrom_dat[chrom_dat$Chr == "2",])
head(alfreq)

# Melt the data
alfreq_m<-melt(data = alfreq,id.vars = c("chrom","pos"),
               measure.vars = c("SLo","Lew","Sha"))
colnames(alfreq_m)<-c("chr","pos","rep","diff")
head(alfreq_m)

nrow(alfreq[alfreq$pos >= 1 & alfreq$pos <= 100000 & alfreq$SLo == 1,])
nrow(alfreq[alfreq$pos >= 1 & alfreq$pos <= 100000 & alfreq$Lew == 1,])
nrow(alfreq[alfreq$pos >= 1 & alfreq$pos <= 100000 & alfreq$Sha == 1,])

# Make some vectors
dxy_v<-vector()
nsubs_v<-vector()
npoly_v<-vector()
avg_diff_v<-vector()
wind_v<-vector()
chr_v<-vector()
start_v<-vector()
end_v<-vector()
nsnps_v<-vector()
rep_v<-vector()
win_len<-50000
win_ovlap<-0
for(rep in unique(alfreq_m$rep)){
  for(chr in unique(alfreq_m$chr)){
    max_l<-chrom_dat$length[as.character(chrom_dat$Chr) == chr]
    chr_dat<-alfreq_m[as.character(alfreq_m$chr) == chr &
                        as.character(alfreq_m$rep) == rep,]
    start <- 1
    end <- start + win_len
    wind <- 1
    while(end < max_l){
      chr_v <- c(chr_v,chr)
      start_v <- c(start_v,start)
      end_v <- c(end_v,end)
      wind_v <- c(wind_v, wind)
      nsnps<-nrow(chr_dat[chr_dat$pos >= start & 
                            chr_dat$pos <= end,])
      nsubs <- nrow(chr_dat[chr_dat$pos >= start & 
                            chr_dat$pos <= end & 
                            chr_dat$diff == 1,])
      npoly <- nrow(chr_dat[chr_dat$pos >= start & 
                            chr_dat$pos <= end &
                            chr_dat$diff != 1,])
      avg_diff<-mean(chr_dat$diff[chr_dat$pos >= start & 
                               chr_dat$pos <= end ],na.rm = TRUE)
      dxy <- nsubs/win_len
      dxy_v <- c(dxy_v,dxy)
      nsubs_v <- c(nsubs_v,nsubs)
      npoly_v <- c(npoly_v,npoly)
      avg_diff_v<-c(avg_diff_v,avg_diff)
      nsnps_v<-c(nsnps_v,nsnps)
      rep_v<-c(rep_v,rep)
      start<-end
      end <- start+win_len
      if(end > max_l){
        end <- max_l
      }
      wind<-wind+1
    }
  }
}
d_theta_dat<-data.frame("chr"=chr_v,
                     "start"=start_v,
                     "end"=end_v,
                     "wind"=wind_v,
                     "nsubs"=nsubs_v,
                     "dxy"=dxy_v,
                     "npoly"=npoly_v,
                     "avg_diff"=avg_diff_v,
                     "nsnps"=nsnps_v,
                     "rep"=rep_v)
d_theta_dat$WinPos <- (d_theta_dat$start+
                     ((d_theta_dat$end-d_theta_dat$start)/2))
d_theta_dat$chr<-gsub("_group",".",d_theta_dat$chr)
head(d_theta_dat)
summary(d_theta_dat)

pal <- c('#d7191c','#fdae61','#2c7bb6')
# Plot
# Plot dXY
adiff_plot_dxy<-ggplot(data=d_theta_dat)+
  geom_line(aes(x = (start+((end-start)/2)),y=dxy,colour=rep))+
  #scale_y_continuous(limits=c(0,1),breaks = c(0,0.5,1))+
  xlab("")+
  ylab(expression(paste(italic(d)[XY])))+
  scale_y_continuous(breaks = c(0,0.015,0.03))+
  scale_colour_manual("Replicate",values=pal)+
  facet_grid(chr~.)

adiff_plot_dxy + my.theme + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10,face = "bold"),
        strip.text.y = element_text(size = 10,angle = 0),
        axis.ticks.x = element_blank(),
        legend.position = "top",
        legend.text = element_text(size=10,face="bold"),
        legend.title = element_text(size=10,face="bold"))



# Plot allele frequency differences, for "top SNPs"
# Get mean alfreq difference
meandata<-tapply(d_theta_dat$avg_diff,
                 INDEX = list(d_theta_dat$WinPos,d_theta_dat$chr),mean)

str(meandata)
winpos<-dimnames(meandata)[[1]]
chrs<-dimnames(meandata)[[2]]

combmeandata<-data.frame(winpos=vector(),adiff=vector,chr=vector())
for(chr in chrs){
  adiff<-meandata[,chr]
  chr<-rep(chr,length(adiff))
  winpos<-as.numeric(names(adiff))
  data<-data.frame(winpos=winpos,adiff=adiff,chr=chr)
  combmeandata<-rbind(combmeandata,data)
}
combmeandata<-na.omit(combmeandata)
rownames(combmeandata)<-seq(1,nrow(combmeandata))
str(combmeandata)
head(combmeandata)

adiff_plot<-ggplot(data=d_theta_dat)+
  geom_line(aes(x=(start+((end-start)/2)),
                 y=avg_diff,colour=rep),alpha=1/2)+
  geom_line(data=combmeandata,aes(x=winpos,
                y=adiff),colour="black",size=1)+
  #scale_y_continuous(limits=c(0,1),breaks = c(0,0.5,1))+
  xlab("")+
  ylab("Allele Frequency Differences")+
  scale_colour_manual("Replicate",values=pal)+
  scale_y_continuous(limits=c(0,1),breaks=c(0,0.5,1))+
  facet_grid(chr~.)

adiff_plot + my.theme + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12,face = "bold"),
        strip.text.y = element_text(size = 12,angle = 0),
        axis.ticks.x = element_blank(),
        legend.position = "top",
        legend.text = element_text(size=10,face="bold"),
        legend.title = element_text(size=10,face="bold"))


# Plot raw difference by SNP: histogram
head(alfreq_m)
rawdiff_plot<-ggplot(data=alfreq_m)+
  geom_histogram(aes(diff,fill=rep),position="dodge",alpha=1/2)+
  xlab("")+
  ylab("Allele Frequency Differences")+
  scale_fill_manual("Replicate",values=pal)

rawdiff_plot + my.theme + 
  theme(#axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12,face = "bold"),
        strip.text.y = element_text(size = 12,angle = 0),
        #axis.ticks.x = element_blank(),
        legend.position = "top",
        legend.text = element_text(size=10,face="bold"),
        legend.title = element_text(size=10,face="bold"))
#
#
#
#
#
##------------------------------------------------------------##
# SIMULATIONS
# What proportion of sites do we expect to be fixed by chance
##------------------------------------------------------------##

# Describe the population
### H. melpomene:    ##
# Keightley et al., 2015 Molecular Biology and Evolution 32:239-243; Ne = 2*10^6
### D. melanogaster: ##
# Keightley et al., 2014 Genetics 196:313-320; Ne = 1.4*10^6
### D. pseudoobscura:##
# Schaeffer 1995 Genetics of Natural Populations: The Continuing Importance of Theodosius Dobzhansky; Ne = ??? (NOT AVAILABLE TO ME)
# Noor et al., 2000 Genetic Research 75:25-35; Ne = 141,000-512,000 (microsats)
# Jensen & Bachtrog 2011 Genome Biology and Evolution 3:687-701; Ne = 4.5*10^6 (nucleotide variation at 123 coding sequence fragments)
### Ne IS DEFINED IN THE SIMULATION: RUN FOR Ne = 1,2,3,and 4*10^6

### D. melanogaster: ##
# Haag-Liautard et al., 2007 Nature 445:82-85; u = 8.4*10^-9
# Keightley et al., 2014 Genetics 196:313-320; u = 2.8*10^-9 (95% CI = 1.0*10^-9 6.1*10^-9)
### H. melpomene: ##
# Keightley et al., 2015 Molecular Biology and Evolution 32:239-243; u = 2.9*10^-9 (95% CI = 1.3*10^-9 5.5*10^-9)

# Simulation parameters
# Simulate XX draws of n individuals from a population (do this XXXX times [sims])
sims<-100
draws<-10000
nes<-c(1000000,2000000,3000000,4000000) # Range of Ne
ns<-c(2,3,4,5,7,8,9,10) # Number of HIGH-LOW comparisons
MU <-c((1.0*10^-9),(2.8*10^-9),(2.9*10^-9),(5.5*10^-9),(8.4*10^-9)) # Range of mutation rates

diffs<-vector(length=sims*length(ns)*length(nes)*length(MU))
Ns<-vector(length=sims*length(ns)*length(nes)*length(MU))
pA<-vector(length=sims*length(ns)*length(nes)*length(MU))
Nes<-vector(length=sims*length(ns)*length(nes)*length(MU))
mus<-vector(length=sims*length(ns)*length(nes)*length(MU))
i<-1
for(Ne in nes){
  for(mu in MU){
    # Charlesworth and Charlesworth 2008 Elements of Population Genetics
    alph<-4*Ne*mu # the alpha parameter for the beta distribution
    # This distribution describes a population at mutation-drift balance
    bet<-alph # alpha and beta are the same
    dist<-rbeta(100000,alph,bet) # Define an allele frequency distribution
    dist<-dist[which(dist != 1 | dist != 0)] # Remove fixed sites
    for(n in ns){
      sim<-1
      while(sim<=sims){
        cat("Ne: ",Ne," - mu: ",mu," - N: ",n," - sim:",sim,"\n")
        draw<-1
        draw_diffs<-vector(length=draws)
        while(draw <= draws){
          # Get the original allele frequency (pA)
          p<-sample(x = dist,size = 1,replace=TRUE)
          # Get haplotype probabilities
          # Assume HWE (p^2 + 2pq + q^2)
          prob_hom_p<-p^2
          prob_hom_q<-(1-p)^2
          prob_het <- 2*p*(1-p)
          probs<-c(prob_hom_p,prob_het,prob_hom_q)
          # Draw n individuals from the "HIGH" lines population
          high_haps<-rmultinom(n=n,prob=probs,size=1)
          # If individuals drawn are hets then p = 0.5; 
          # else p = 1 (hom_p) or 0 (hom_q)
          high<-vector(length=n)
          high[which(high_haps[1,]==1)] <- 1 # hom_p
          high[which(high_haps[3,]==1)] <- 0 # hom_q
          high[which(high_haps[2,]==1)] <- 0.5 # hets
          # Draw n individuals from the "LOW" lines population
          low_haps<-rmultinom(n=n,prob=probs,size=1)
          # If individuals drawn are hets then p = 0.5; 
          # else p = 1 (hom_p) or 0 (hom_q)
          low<-vector(length=n)
          low[which(low_haps[1,]==1)] <- 1 # hom_p
          low[which(low_haps[3,]==1)] <- 0 # hom_q
          low[which(low_haps[2,]==1)] <- 0.5 # hets
          # Calculate the differences between the lines
          diff<-high-low
          # Are all line comparisons fixed for the opposite allele
          if(all(diff==1)){
            draw_diffs[draw]<-"Fixed"
          }else{
            draw_diffs[draw]<-"Variable"
          }
          draw<-draw+1
        }
        sim<-sim+1
        if(is.na((table(draw_diffs)["Fixed"]/length(draw_diffs))*100)){
          diffs[i]<-0
        } else{
        diffs[i]<-(table(draw_diffs)["Fixed"]/length(draw_diffs))*100 
        }
        Ns[i]<-paste("n = ",n,sep="")
        Nes[i]<-paste("Ne = ",Ne,sep="")
        mus[i]<-mu
        i<-i+1
      }
    }
  }
}
p_dat<-data.frame(diffs=diffs,n=Ns,Ne=Nes,mu=mus)
head(p_dat)
p_dat$n<-factor(p_dat$n,levels=paste("n = ",ns,sep=""),
                labels=paste("n = ",ns,sep=""))
p_dat$Ne<-factor(p_dat$Ne,labels=paste("Ne = ",nes,sep=""))
p_dat$mu<-factor(p_dat$mu,
                 levels=MU,
                 labels=c(expression(paste("1.0x",10^-9,sep="")),
                          expression(paste("2.8x",10^-9,sep="")),
                          expression(paste("2.9x",10^-9,sep="")),
                          expression(paste("5.5x",10^-9,sep="")),
                          expression(paste("8.4x",10^-9,sep=""))))
tail(p_dat)
# What percentage of the time do we expect all line comparisons
# to be fixed for the alternative?
tapply(p_dat$diffs,INDEX = list(p_dat$n,p_dat$Ne,p_dat$mu),mean)
# What is the 95th percentile of the distribution of number of expected fixations
tapply(p_dat$diffs,INDEX = list(p_dat$n,p_dat$Ne,p_dat$mu),quantile, probs=0.95)
tapply(p_dat$diffs,INDEX = list(p_dat$n,p_dat$Ne,p_dat$mu),quantile, probs=0.99)
tapply(p_dat$diffs,INDEX = list(p_dat$n,p_dat$Ne,p_dat$mu),max)

# What is the relationship between this case and pA
head(p_dat)
str(p_dat)
ggplot()+
  geom_boxplot(data=p_dat,aes(n,diffs,fill=mu),
               position=position_dodge(0.9))+
#  geom_point(data=p_dat,aes(n,diffs,colour=Ne),alpha = 1/2)+
  geom_hline(yintercept = 0.022,linetype="dashed",size=0.5)+
  xlab("")+
  ylab("Proportion of All Fixed Differences (%)")+
  scale_fill_manual("Mutation Rate",
                    labels=c(expression(paste("1.0x",10^-9,sep="")),
                             expression(paste("2.8x",10^-9,sep="")),
                             expression(paste("2.9x",10^-9,sep="")),
                             expression(paste("5.5x",10^-9,sep="")),
                             expression(paste("8.4x",10^-9,sep=""))),
                    values=c("darkred","red","darkorange","yellow","white"))+
  scale_colour_discrete("")+
  my.theme+
  facet_grid(Ne~.)+
  theme(axis.text.x = element_text(size=10,face="bold"),
    axis.text.y = element_text(size = 10,face = "bold"),
    strip.text.y = element_text(size = 10,angle = 0),
    axis.title.y = element_text(size=10,face="bold"), 
    #axis.ticks.x = element_blank(),
    legend.position = "right",
    legend.key = element_rect(fill="white"),
    legend.text = element_text(size=10,face="bold"),
    legend.title = element_text(size=10,face="bold"))




##------------------------------------------------------------##
# PLOT CARTOON
##------------------------------------------------------------##
x<-rnorm(1000000,0)
range(x)
h<-quantile(x,probs = 0.95)
l<-quantile(x,probs = 0.05)
ggplot()+geom_histogram(aes(x),fill="white",colour="black")+
  xlab("Female Re-Mating Rate")+
  xlim(-5,5)+
  geom_rect(aes(xmin = h,xmax=h+0.1,ymax=100000,ymin=0),fill="red")+
  geom_rect(aes(xmin = l,xmax=l-0.1,ymax=100000,ymin=0),fill="blue")+
  theme(axis.title.x = element_text(size=10,angle=0),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_line(NULL),
        axis.title.y = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

##------------------------------------------------------------##
# ALLELE FREQUENCY DISTRIBUTIONS IN EACH LINE
##------------------------------------------------------------##
# Load data
al_freq_dat <- read.table("pseudo_wild_SHA-SLOB7subs_filt_alfreq.tab",
                       sep = "\t", header = FALSE)
colnames(al_freq_dat) <- c("Chr","Pos",
                           "ref","maj",
                           "SLOC9","SLOB7",
                           "LEW17","LEW23",
                           "SHAC1","SHAA10")
head(al_freq_dat)
al_freq_dat_m<-melt(al_freq_dat[,c(5,6,7,8,9)])
ggplot()+
  geom_histogram(data=al_freq_dat,aes(SHAC1))+
  my.theme+
  theme(axis.text.x = element_text(size=12,face="bold"),
        axis.text.y = element_text(size = 12,face = "bold"),
        strip.text.y = element_text(size = 12,angle = 0),
        #axis.ticks.x = element_blank(),
        legend.position = "top",
        legend.text = element_text(size=10,face="bold"),
        legend.title = element_text(size=10,face="bold"))






