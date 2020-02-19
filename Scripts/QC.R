### Title : shRNA screening (DGG Library) - TMZ drug study.
### Author : Perales-Paton, Javier
### Date/Version : July17 - v1.0

# SOLVED -> OLD: WARNING THERE ARE DUPLICATED TAGS (due to multiple-used hairpins across pool)


#---- Intern Variables : Defined by the user #####
  lib.tag <- "DDR_miRE"
  #

Proj.dir <- "/drives/slave/jperales/Projects/SquatritoM_U251_sgMSH6_May17/HighTMZt1w_HighDMSOt1w";

#---- The output files ####
Out_plots.dir <- paste0(Proj.dir,"/Results/Plots/")
Out_tables.dir <- paste0(Proj.dir,"/Results/Tables/")


setwd(paste0(Proj.dir,"/Ranalysis/"))

#---- Load libs and wd
source(paste0(Proj.dir,"/Build/shRNA_func.R"))
library(edgeR)
# library(survcomp) # the combine.test function
library(gplots)

#---- The input files ####
# Proc log files
Proc.dir <- paste0(Proj.dir,"/Processed_FastQ/")
proc.log_fls <- list.files(Proc.dir,full.names = FALSE,pattern = ".log2")
samplesNames.tmp <- gsub(pattern = ".log2","",proc.log_fls)
names(proc.log_fls) <- samplesNames.tmp
rm(samplesNames.tmp)

# Lib
shRNA_lib.fl <- paste0(Proj.dir,"/../REFERENCES/shRNA_libs/",lib.tag,"/",lib.tag,".rev_compSPECIFIC_STEMs.tsv")
shRNA_lib <- read.table(shRNA_lib.fl, col.names = c("ID","seq","batch"),
                        sep="\t",stringsAsFactors = FALSE)
shRNA_lib.tmp4sortingtags <- data.frame(V1=shRNA_lib$batch,
                                        V2=sapply(shRNA_lib$ID,function(z) strsplit(z,split="_")[[1]][1]),
                                        V3=as.numeric(sapply(shRNA_lib$ID,function(z) strsplit(z,split="_")[[1]][2])),
                                        stringsAsFactors = FALSE)
shRNA_lib <- shRNA_lib[order(shRNA_lib.tmp4sortingtags$V1,shRNA_lib.tmp4sortingtags$V2,shRNA_lib.tmp4sortingtags$V3),]


#hgnc lib
hgnc.obj <- read.table(file="/local/jperales/REFERENCES/HGNC/hgnc_complete_set__v170516_1009.txt",sep="\t",
                       comment.char = "",header = TRUE,quote = "",fill = TRUE,
                       stringsAsFactors = FALSE)

# Read count files
Counts.dir <- paste0(Proj.dir,"/Count/")
fls <- list.files(Counts.dir,full.names = FALSE)

# Reorder files (put together the biological group)
fls <- fls[c(grep("DMSO",fls),
             grep("TMZ",fls))]

#---- Experimental Design (order was established in previous line)
  # Samples
  samples.names <- gsub(".counts.tsv","",fls)
  samples.gr <- as.factor(sapply(samples.names, function(x) strsplit(x,split="_")[[1]][5]))
  samples.replicates <- as.factor(sapply(samples.names, function(x) strsplit(x,split="_")[[1]][7]))

  # Hairpins per pool
  unq <- unique(shRNA_lib$batch)
  pool.idx <- list()
  for (i in unq) {
    pool.idx[[i]] = shRNA_lib[which(shRNA_lib$batch==i),"ID"]
  }
  rm(i,unq)

  # Create a tmp lib
  shRNA_lib2 <- unique(shRNA_lib)
  rownames(shRNA_lib2) <- shRNA_lib2$ID

  # Code of colors
  col.LIBpools <- topo.colors(length(unique(shRNA_lib$batch))); names(col.LIBpools) <- unique(shRNA_lib$batch)


#---- Load data #######
# Processing Log files
proc.mt <- sapply(proc.log_fls,function(z) as.numeric(strsplit(readLines(paste0(Proc.dir,"/",z))[4],split=" ")[[1]][c(3,5,8)]))

proc.mt <- as.data.frame(proc.mt)
rownames(proc.mt) <- c("Total","Surviours","Dropped")

# Counts files
count.fl_list <- as.list(fls)
names(count.fl_list)  <- samples.names

# This function comes from the source file, which need list of tags. We know that
#   there are all hairpin IDs from the lib, but also other tag called "NOTaligned"
#   which contains those reads without any hit against a hairpin
cnts.df <- readCounts(listOfFiles = count.fl_list,sampleNames = samples.names,
           listOfTags = c(unique(shRNA_lib$ID),"NOTaligned"),pathDIR = Counts.dir)

# Sanity: eparate cnts of hairpins from notalinged 
cnts.NOTaligned <- cnts.df[which(rownames(cnts.df)=="NOTaligned"),]
cnts.df <- cnts.df[which(rownames(cnts.df)!="NOTaligned"),]

# Normalized counts -> Counts Per Million (we use cpm function from edgeR package)
# cpm.df <- cpm(cnts.df,normalized.lib.sizes=FALSE,lib.size=colSums(cnts.df),log=FALSE)
cpm.df <- cpm.ByPool(x = cnts.df,shRNA_lib = shRNA_lib)

#---- Exploratory data analysis ######
# Processing raw data:
proc.mt2 <- rbind(proc.mt,Hit=proc.mt[2,]-cnts.NOTaligned,NOhit=cnts.NOTaligned)
proc.mt2 <- proc.mt2[,colnames(cnts.df)]
stopifnot(all(proc.mt2["Surviours",samples.names]==(colSums(cnts.df)+cnts.NOTaligned)))

png(filename = paste0(Out_plots.dir,lib.tag,"_processed.png"),width = 1000,height = 1200,res=180)
par(mar=c(19,4,2,7),xpd = TRUE)
barplot(prop.table(as.matrix(proc.mt2[3:5,]),margin=2)*100,main="Lost reads during Processing raw data",
        col=c("grey4","lightblue","grey"),ylab="(%) Percentage total reads",las=2)
legend("topright",legend = c("Bad reads (discarded)","Hits on shRNA","no hit"),
       fill=c("grey4","lightblue","grey"), cex=0.6,inset = c(-0.40, 0))
dev.off()

# Plot number of hairpins that could be matched per sample
# F(x): visualize any samples under represented by the sequencing run
png(filename = paste0(Out_plots.dir,lib.tag,"_cntsPerSample.png"),width = 1000,height = 1200,res=180)
par(mar=c(19,4,1,7),xpd = TRUE)
barplot(as.matrix(rbind(hits=colSums(cnts.df)*1e-6,
                        no_hits=cnts.NOTaligned*1e-6)),
        beside = F, las = 2, 
        main = "Effective Read counts per Sample", ylab="Library size (millions)",
        cex.names = 1, cex.axis = 1.2)
legend("topright",legend = c("Hits on shRNA","no hit"),fill=c("grey4","grey"), cex=0.6,inset = c(-0.40, 0))
dev.off()

# Plot number of percentage of hairpins in each pool
# F(x): bad pooling of hairpins
png(filename = paste0(Out_plots.dir,lib.tag,"_propPool.png"),width = 1500,height = 1200,res=180)
cntBypool <- matrix(NA,nrow=length(unique(shRNA_lib$batch))+1,ncol=length(samples.names),
                    dimnames=list(c(unique(shRNA_lib$batch),"NOhit"),samples.names))
cntBypool <- as.data.frame(cntBypool)

for(id in unique(shRNA_lib$batch)) {
  cntBypool[id,samples.names] <- colSums(cnts.df[pool.idx[[id]],samples.names])
}
cntBypool["NOhit",samples.names] <- cnts.NOTaligned[,samples.names]
cntBypool.perc <- prop.table(as.matrix(cntBypool),margin=2)*100
par(mar=c(19,4,4,12),xpd = TRUE)
barplot(cntBypool.perc,las=2,main=paste0(lib.tag," : Proportion of total hairpins sequenced per pool"),
        col = c(topo.colors(nrow(cntBypool.perc)-1),"grey"))
legend("topright", inset = c(-0.30, 0), fill = c(topo.colors(nrow(cntBypool.perc)-1),"grey"), 
       legend = c(paste0(rownames(cntBypool.perc)[1:length(pool.idx)]," (n=",
                         unlist(lapply(pool.idx,function(z) length(z))),")"),"No hit"))
dev.off()

# Plot per hairpin totals across groups
col.features <- sapply(shRNA_lib2[rownames(cnts.df),"batch"], function(y) col.LIBpools[y])

png(filename = paste0(Out_plots.dir,lib.tag,"_totalCntGr.png"),width = 2400,height = 1800,res=200)
par(mar=c(8,5,2,2),mfrow=c(2,1))
for(gr in levels(samples.gr)) {


if(sum(samples.gr==gr)==1) {
  top15 <- rownames(cnts.df)[order(cnts.df[,samples.gr==gr],decreasing = TRUE)][1:15]
  topNames <- rownames(cnts.df)
  topNames[!topNames%in%top15] <- ""
  barplot(cnts.df[,samples.gr==gr],border = col.features,
          las = 2, main=paste0(gr,": Total counts per hairpin within group",
                          "(n=",sum(samples.gr==gr),")"), cex.names = 1,
          #         sub=paste(sapply(c(1,5,10),function(x) {
          #           paste(top15[seq(from = x,to=x+5,by=1)],collapse = "; ")}
          #                         ),collapse = "\n"),cex.sub=0.7,
          axisnames = TRUE,names.arg = topNames
  )
} else if (sum(samples.gr==gr)>1) {
  top15 <- names(sort(rowSums(cnts.df[,samples.gr==gr]),decreasing=TRUE)[1:15])
  topNames <- rownames(cnts.df[,samples.gr==gr])
  topNames[!topNames%in%top15] <- ""
  barplot(rowSums(cnts.df[,samples.gr==gr]),border = col.features,
          las = 2, main = paste0(gr,": Total counts per hairpin within group",
                                 "(n=",sum(samples.gr==gr),")"), cex.names = 1,
          axisnames = TRUE,names.arg = topNames
  ) 
} else {
  stop("ERROR")
}

}
rm(gr,top15,topNames)
dev.off()

#---- Analyses of ratios: FOLD-change

# Normalization TRIMMED-MEAN Of M-values
# y <- readDGE(files = count.fl_list,path = "../Count/",group = samples.gr,labels=names(count.fl_list))
# y <- y[rownames(y)!="NOTaligned", , keep.lib.sizes=FALSE]
# cnts.df <- y$counts
# y <- calcNormFactors(y)
# cpm.df <- cpm(y)
cnts.df <- readCounts(listOfFiles = count.fl_list,sampleNames = samples.names,
                      listOfTags = c(unique(shRNA_lib$ID),"NOTaligned"),pathDIR = Counts.dir)
cnts.df <- cnts.df[which(rownames(cnts.df)!="NOTaligned"),]
cpm.df <- cpm.ByPool(x = cnts.df,shRNA_lib = shRNA_lib)


#---- Build raw and normalized tables:
# if(all(rownames(cnts.df)==names(sel))) {
  cnts.out <- cbind(hairpinID=rownames(cnts.df),
                    pool=shRNA_lib2[rownames(cnts.df),"batch"],
                    as.data.frame(cnts.df)
                    )
  write.table(cnts.out,file="./counts.tsv",sep="\t",row.names = FALSE,col.names=TRUE)
# } else {
#   stop("ERROR: #1")
# }

# if(all(rownames(cpm.df)==names(sel))) {
  cpm.out <- cbind(hairpinID=rownames(cpm.df),
                   pool=shRNA_lib2[rownames(cpm.df),"batch"],
                   as.data.frame(cpm.df)
  )
  write.table(cnts.out,file="./cpm.tsv",sep="\t",row.names = FALSE,col.names=TRUE)
  
# # } else {
#   stop("ERROR: #2")
# }

#3) correlation between groups
colors.batch <- sapply(shRNA_lib2[rownames(cnts.df),"batch"], function(y) col.LIBpools[y])


# Intra-group
# DMSO
png(paste0(Out_plots.dir,"/consistency_intragroupDMSO.png"),width = 1200*3,height = 415*3,res=230)
par(mfrow=c(1,3))

plot(log2(cpm.df[,samples.names[samples.gr=="DMSO" & samples.replicates=="E1"]]+1),
     log2(rowMeans(cpm.df[,samples.names[samples.gr=="DMSO" & samples.replicates%in%c("E2","E3")]])+1),
     col="black",pch=21,cex=0.9,
     main="DMSO E1",
     xlab="log2(E1_CPM+1)",ylab="log2(avg E23_CPM+1)",
     bg=colors.batch)
plot(log2(cpm.df[,samples.names[samples.gr=="DMSO" & samples.replicates=="E2"]]+1),
     log2(rowMeans(cpm.df[,samples.names[samples.gr=="DMSO" & samples.replicates%in%c("E1","E3")]])+1),
     col="black",pch=21,cex=0.9,
     main="DMSO E2",
     xlab="log2(E2_CPM+1)",ylab="log2(avg E13_CPM+1)",
     bg=colors.batch)
plot(log2(cpm.df[,samples.names[samples.gr=="DMSO" & samples.replicates=="E3"]]+1),
     log2(rowMeans(cpm.df[,samples.names[samples.gr=="DMSO" & samples.replicates%in%c("E1","E2")]])+1),
     col="black",pch=21,cex=0.9,
     main="DMSO E3",
     xlab="log2(E3_CPM+1)",ylab="log2(avg E12_CPM+1)",
     bg=colors.batch)
dev.off()
# TMZ100
png(paste0(Out_plots.dir,"/consistency_intragroupTMZ.png"),width = 1200*3,height = 415*3,res=230)
par(mfrow=c(1,3))
plot(log2(cpm.df[,samples.names[samples.gr=="TMZ" & samples.replicates=="E1"]]+1),
     log2(rowMeans(cpm.df[,samples.names[samples.gr=="TMZ" & samples.replicates%in%c("E2","E3")]])+1),
     col="black",pch=21,cex=0.9,
     main="TMZ E1",
     xlab="log2(E1_CPM+1)",ylab="log2(avg E23_CPM+1)",
     bg=colors.batch)
plot(log2(cpm.df[,samples.names[samples.gr=="TMZ" & samples.replicates=="E2"]]+1),
     log2(rowMeans(cpm.df[,samples.names[samples.gr=="TMZ" & samples.replicates%in%c("E1","E3")]])+1),
     col="black",pch=21,cex=0.9,
     main="TMZ E2",
     xlab="log2(E2_CPM+1)",ylab="log2(avg E13_CPM+1)",
     bg=colors.batch)
plot(log2(cpm.df[,samples.names[samples.gr=="TMZ" & samples.replicates=="E3"]]+1),
     log2(rowMeans(cpm.df[,samples.names[samples.gr=="TMZ" & samples.replicates%in%c("E1","E2")]])+1),
     col="black",pch=21,cex=0.9,
     main="TMZ E3",
     xlab="log2(E3_CPM+1)",ylab="log2(avg E12_CPM+1)",
     bg=colors.batch)
dev.off()
# dev.off()
# # Intra-exp
# png(paste0(Out_plots.dir,"/consistency_repl.png"),width = 1200*3,
#     height = (830/2)*3,res=230)
# par(mfrow=c(1,3))
# plot(log2(cpm.df[,samples.names[samples.gr=="t1w" & samples.replicates=="E1"]]+1),
#      log2(cpm.df[,samples.names[samples.gr=="t1wTMZ100" & samples.replicates=="E1"]]+1),
#      col="black",pch=21,cex=0.9,
#      main="E1",
#      xlab="DMSO - log2(E1_CPM+1)",ylab="TMZ100uM - log2(E1_CPM+1)",
#      bg=colors.batch)
# 
# plot(log2(cpm.df[,samples.names[samples.gr=="t1w" & samples.replicates=="E2"]]+1),
#      log2(cpm.df[,samples.names[samples.gr=="t1wTMZ100" & samples.replicates=="E2"]]+1),
#      col="black",pch=21,cex=0.9,
#      main="E2",
#      xlab="DMSO - log2(E2_CPM+1)",ylab="TMZ100uM - log2(E2_CPM+1)",
#      bg=colors.batch)
# 
# plot(log2(cpm.df[,samples.names[samples.gr=="t1w" & samples.replicates=="E3"]]+1),
#      log2(cpm.df[,samples.names[samples.gr=="t1wTMZ100" & samples.replicates=="E3"]]+1),
#      col="black",pch=21,cex=0.9,
#      main="E3",
#      xlab="DMSO - log2(E3_CPM+1)",ylab="TMZ100uM - log2(E3_CPM+1)",
#      bg=colors.batch)
# dev.off()

# png(paste0(Out_plots.dir,"/consistency_replTMZ200.png"),width = (1200/4)*3,
#     height = (830/3)*3,res=230/2)
# par(mfrow=c(1,1))
# plot(log2(cpm.df[,samples.names[samples.gr=="t1w" & samples.replicates=="E3"]]+1),
#      log2(cpm.df[,samples.names[samples.gr=="t1wTMZ200" & samples.replicates=="E3"]]+1),
#      col="black",pch=21,cex=0.9,
#      main="E3",
#      xlab="DMSO - log2(E3_CPM+1)",ylab=" TMZ200 log2(E3_CPM+1)",
#      bg=colors.batch)
# dev.off()

# PCA
# library("scatterplot3d")
png(paste0(Out_plots.dir,"/PCA_PC12.png"),width = 600*3,height = 425*3,res=230)
samples4PCA <- colnames(cnts.df)
base.pca2 <- prcomp(t(cpm(cnts.df[,samples4PCA],log=TRUE)))
eig.base.pca <- data.frame(eig = (base.pca2$sdev)^2,
                           variance_expl = ((base.pca2$sdev)^2)*100/sum((base.pca2$sdev)^2),
                           cumvariance = cumsum(((base.pca2$sdev)^2)*100/sum((base.pca2$sdev)^2))
)
stopifnot(all(names(samples.gr)==colnames(cnts.df)))
par(mar=c(4,4,4,10),xpd = TRUE)
plot(PC2~PC1, data=as.data.frame(base.pca2$x)[samples4PCA,],
     cex=1.3,
     col=c("coral4","olivedrab3","skyblue2","orange2")[as.integer(samples.gr[samples4PCA])],
     # bg=c("coral4","olivedrab3","skyblue2","orange2")[as.integer(samples.gr[samples4PCA])], 
     # pch=c(7,21,22,23,24,25)[as.integer(samples.replicates[samples4PCA])],
     pch=19,
     main="PCA",
     xlab=paste0("PC1 (",round(eig.base.pca$variance_expl[1],digits=1)," variance explained)"),
     ylab=paste0("PC2 (",round(eig.base.pca$variance_expl[2],digits=1)," variance explained)"))
text(as.data.frame(base.pca2$x)[samples4PCA,"PC1"],
     as.data.frame(base.pca2$x)[samples4PCA,"PC2"],
     labels=as.character(samples.replicates[samples4PCA]),cex=1.3)

# legend("topright",
#        legend = unique(samples.replicates),
#        pch=c(7,21,22,23,24,25),
#        col=c("black"), cex=1,inset = c(-0.20, 0))
legend("right",title="Condition",
       legend = unique(levels(samples.gr))[unique(as.integer(samples.gr[samples4PCA]))],
       pch=c(21),
       pt.bg=c("coral4","olivedrab3","skyblue2","orange2")[as.integer(unique(samples.gr[samples4PCA]))],
       cex=1.3,inset = c(-0.35, 0))
dev.off()

png(paste0(Out_plots.dir,"/PCA_PC13.png"),width = 600*3,height = 425*3,res=230)
par(mar=c(4,4,4,10),xpd = TRUE)
plot(PC3~PC1, data=as.data.frame(base.pca2$x)[samples4PCA,],
     cex=1.3,
     col=c("coral4","olivedrab3","skyblue2","orange2")[as.integer(samples.gr[samples4PCA])],
     # bg=c("coral4","olivedrab3","skyblue2","orange2")[as.integer(samples.gr[samples4PCA])], 
     # pch=c(7,21,22,23,24,25)[as.integer(samples.replicates[samples4PCA])],
     pch=19,
     main="PCA",
     xlab=paste0("PC1 (",round(eig.base.pca$variance_expl[1],digits=1)," variance explained)"),
     ylab=paste0("PC3 (",round(eig.base.pca$variance_expl[3],digits=1)," variance explained)"))
text(as.data.frame(base.pca2$x)[samples4PCA,"PC1"],
     as.data.frame(base.pca2$x)[samples4PCA,"PC3"],
     labels=as.character(samples.replicates[samples4PCA]),cex=1.3)
legend("right",title="Condition",
       legend = unique(levels(samples.gr))[unique(as.integer(samples.gr[samples4PCA]))],
       pch=c(21),
       pt.bg=c("coral4","olivedrab3","skyblue2","orange2")[as.integer(unique(samples.gr[samples4PCA]))],
       cex=1.3,inset = c(-0.35, 0))
dev.off()


png(paste0(Out_plots.dir,"/MGMT_boxplots.png"),width = 800*3,height = 600*3,res=230)
MGMT.hairpins <- sort(grep("^MGMT_",rownames(cpm.df),value=TRUE))
par(mfrow=c(2,3),mar=c(8,4,4,4))
for(MGMT.hairpin in MGMT.hairpins) {
# MGMT.case <- sapply(levels(y$samples$group),function(z) as.matrix(log2(cpm.df[MGMT.hairpin,y$samples$group==z]+1)))
  MGMT.case <- sapply(levels(samples.gr),function(z) as.matrix(cpm.df[MGMT.hairpin,samples.gr==z]))
  
boxplot(MGMT.case,main=MGMT.hairpin,ylab="CPM",las=2)
stripchart(MGMT.case, vertical = TRUE,cex=1.5, 
           method = "jitter", add = TRUE, pch = 20, col=c("coral4","olivedrab3","skyblue2","orange2"))
}
dev.off()

# library(rgl)
# plot3d(as.data.frame(base.pca2$x)[,1:3],
#        cex=1.5,
#        col=c("coral4","olivedrab3","skyblue2")[as.integer(samples.gr[samples4PCA])],
#        pch=19,main="PCA",
#        xlab=paste0("PC1 (",round(eig.base.pca$variance_expl[1],digits=1)," variance explained)"),
#        ylab=paste0("PC2 (",round(eig.base.pca$variance_expl[2],digits=1)," variance explained)"),
#        zlab =paste0("PC3 (",round(eig.base.pca$variance_expl[3],digits=1)," variance explained)"))


library("scatterplot3d")
png(paste0(Out_plots.dir,"/PCA.png"),width = 600*3,height = 600*3,res=230)
# par(mar=c(4,4,4,10),xpd = TRUE)
DF <- as.data.frame(base.pca2$x)[samples4PCA,]
s3d <- scatterplot3d(x = DF$PC1,y= DF$PC2,z = DF$PC3,pch=19,type="h",
              color=c("coral4","olivedrab3","skyblue2","orange2")[as.integer(samples.gr[samples4PCA])],
              xlab=paste0("PC1 (",round(eig.base.pca$variance_expl[1],digits=1)," variance explained)"),
              ylab=paste0("PC2 (",round(eig.base.pca$variance_expl[2],digits=1)," variance explained)"),
              zlab=paste0("PC3 (",round(eig.base.pca$variance_expl[3],digits=1)," variance explained)"),
              main="Principal Component Analysis")
text(s3d$xyz.convert(DF$PC1, DF$PC2, DF$PC3), labels=as.character(samples.replicates[samples4PCA]), pos=1)
legend("topright",title="Condition",
       legend = unique(levels(samples.gr))[unique(as.integer(samples.gr[samples4PCA]))],
       pch=c(21),
       pt.bg=c("coral4","olivedrab3","skyblue2","orange2")[as.integer(unique(samples.gr[samples4PCA]))],
       cex=1)
dev.off()

