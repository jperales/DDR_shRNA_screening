### Title : shRNA screening (DGG Library) - TMZ drug study.
### Author : Perales-Paton, Javier
### Date/Version : July2017 - v6.0


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

# Reorder samples
fls <- fls[c(grep("_DMSO",fls),grep("_TMZ",fls))]


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
# Counts files
count.fl_list <- as.list(fls)
names(count.fl_list)  <- samples.names


#---2) limma
# cnts.df <- readCounts(listOfFiles = count.fl_list,sampleNames = samples.names,
#                       listOfTags = c(unique(shRNA_lib$ID),"NOTaligned"),pathDIR = Counts.dir)
# cnts.df <- cnts.df[which(rownames(cnts.df)!="NOTaligned"),]

y<- readDGE(fls,path = paste0(Proj.dir,"/Count"),group = samples.gr,labels = samples.names)
# y <- DGEList(counts = cnts.df,lib.size = colSums(cnts.df),group = samples.gr)
y <- y[which(rownames(y)!="NOTaligned"), , keep.lib.sizes=FALSE]
des <-  model.matrix(~samples.replicates+samples.gr)

# y <- calcNormFactors(y)
# 
# y <- estimateCommonDisp(y)

v <- voom(y,design = des)
# Fit the Linear model for each gene
fit.v <- lmFit(v,des)


# Please check that the order of biological groups are as expected they must be: WT > KD > LOX
fit.v$design
eBayes.v <- eBayes(fit.v)

TMZvsDMSO <-  topTable(eBayes.v, coef="samples.grTMZ",number = Inf)
TMZvsDMSO <- TMZvsDMSO[order(TMZvsDMSO$logFC,decreasing=FALSE),]

# ROAST
# Hairpins per gene
unq <- sapply(rownames(v),function(z) strsplit(z,split="_")[[1]][1])
targetGSlike.idx <- list()
for (i in unq) {
  targetGSlike.idx[[i]] = grep(paste0("^",i,"_"),rownames(v))
}
rm(i,unq)

GS <- mroast(v,index=targetGSlike.idx,design=des,contrast=ncol(des),nrot = 9999)
GS$score <- apply(GS,1,function(z) ifelse(z["Direction"]=="Down",log10(as.numeric(z["PValue"])),-log10(as.numeric(z["PValue"]))))

score <- vector(mode = "list")
score$RiScore <- GS$score
names(score$RiScore) <- rownames(GS)

RM <- romer(v,index=targetGSlike.idx,design=des,contrast=ncol(des),nrot = 9999)
tmp <- topRomer(RM,alternative = "down",n = Inf)
score$sScore <- log10(tmp[,"Down"])
names(score$sScore) <- rownames(tmp)



#---3) TMZ vs DMSO ####
# This function comes from the source file, which need list of tags. We know that
#   there are all hairpin IDs from the lib, but also other tag called "NOTaligned"
#   which contains those reads without any hit against a hairpin
cnts.df <- readCounts(listOfFiles = count.fl_list[samples.names],sampleNames = samples.names,
                      listOfTags = c(unique(shRNA_lib$ID),"NOTaligned"),pathDIR = Counts.dir)
cnts.df <- cnts.df[which(rownames(cnts.df)!="NOTaligned"),]
cpm.df <- cpm.ByPool(x = cnts.df,shRNA_lib = shRNA_lib)

depletion.log2cutoff <- log2(1/5)
depletion.tag <- "1/5"

# Fold change = log2 ( replicate j TMZ / replicate j DMSO)
TMZvsDMSO.log2E1 <- log2((cpm.df[,samples.names[samples.gr=="TMZ" & samples.replicates=="E1"]]+1)/
                               (cpm.df[,samples.names[samples.gr=="DMSO" & samples.replicates=="E1"]]+1))

TMZvsDMSO.log2E2 <- log2((cpm.df[,samples.names[samples.gr=="TMZ" & samples.replicates=="E2"]]+1)/
                               (cpm.df[,samples.names[samples.gr=="DMSO" & samples.replicates=="E2"]]+1))

TMZvsDMSO.log2E3 <- log2((cpm.df[,samples.names[samples.gr=="TMZ" & samples.replicates=="E3"]]+1)/
                               (cpm.df[,samples.names[samples.gr=="DMSO" & samples.replicates=="E3"]]+1))


TMZ_avgCPM <- rowMeans(cpm.df[,samples.names[samples.gr=="TMZ"]])
DMSO_avg_CPM <- rowMeans(cpm.df[,samples.names[samples.gr=="DMSO"]])

TMZ.E123.avg_log2FC <- rowMeans(cbind(TMZvsDMSO.log2E1,TMZvsDMSO.log2E2,TMZvsDMSO.log2E3))




TMZvsDMSO <- data.frame(Symbol=sapply(rownames(cpm.df),function(z) strsplit(z,split="_")[[1]][1]),
                            sScore=sapply(rownames(cpm.df),function(z) score$sScore[strsplit(z,split="_")[[1]][1]]),
                            RiScore=sapply(rownames(cpm.df),function(z) score$RiScore[strsplit(z,split="_")[[1]][1]]),
                            hairpinID=rownames(cpm.df),
                            pool=shRNA_lib2[rownames(cpm.df),"batch"],
                            E123_avgLog2FC=TMZ.E123.avg_log2FC,
                            E1_log2FC=TMZvsDMSO.log2E1,
                            E2_log2FC=TMZvsDMSO.log2E2,
                            E3_log2FC=TMZvsDMSO.log2E3,
                            DMSO_avg_CPM=DMSO_avg_CPM,
                            TMZ_avgCPM=TMZ_avgCPM
)

TMZvsDMSO$Name <- sapply(TMZvsDMSO$Symbol,function(z) {
  if(length(which(hgnc.obj$symbol==z)) > 0) {
    hgnc.obj[which(hgnc.obj$symbol==z),"name"];
  } else {
    "";
  }
})


# Ranking
sorted.RnkByFC <- 1:length(TMZ.E123.avg_log2FC)
names(sorted.RnkByFC) <- names(sort(TMZ.E123.avg_log2FC))

TMZvsDMSO$Gene.RankingBy_avgLog2FC  <- sapply(TMZvsDMSO$Symbol,function(z) {
  paste(sort(sorted.RnkByFC[rownames(subset(TMZvsDMSO,Symbol==z))]),collapse = "; ")
})

TMZvsDMSO$Gene.NoHairpins_Below1Percentile <- sapply(TMZvsDMSO$Symbol,function(z) {
  sum(subset(TMZvsDMSO,Symbol==z)$E123_avgLog2FC< quantile(TMZvsDMSO$E123_avgLog2FC,0.10))
})

TMZvsDMSO$Gene.NoHairpins_Above99Percentile <- sapply(TMZvsDMSO$Symbol,function(z) {
  sum(subset(TMZvsDMSO,Symbol==z)$E123_avgLog2FC > quantile(TMZvsDMSO$E123_avgLog2FC,0.90))
})


# Clean global_env
rm(TMZ.E123.avg_log2FC,TMZvsDMSO.log2E1,TMZvsDMSO.log2E2,TMZvsDMSO.log2E3)
rm(TMZ_avgCPM,DMSO_avg_CPM)

# Stardard table for shRNA: avg FC
TMZvsDMSO.tmp <- TMZvsDMSO[order(TMZvsDMSO$E123_avgLog2FC,
                                 -TMZvsDMSO$Gene.NoHairpins_Below1Percentile,
                                 decreasing = FALSE),]
TMZvsDMSO.stard_rnk <- TMZvsDMSO.tmp[,c("Symbol",
                                                "Gene.RankingBy_avgLog2FC",
                                                "sScore",
                                                "RiScore",
                                                "Gene.NoHairpins_Below1Percentile",
                                                "Gene.NoHairpins_Above99Percentile",
                                                "hairpinID","pool",
                                                "E123_avgLog2FC","E1_log2FC","E2_log2FC","E3_log2FC",
                                                "DMSO_avg_CPM",
                                                "TMZ_avgCPM",
                                                "Name")]

tmz_vs_t1w.stard_genernk <- unique(TMZvsDMSO.stard_rnk[,c("Symbol",
                                                              "Gene.RankingBy_avgLog2FC",
                                                              "sScore",
                                                              "RiScore",
                                                              "Gene.NoHairpins_Below1Percentile",
                                                              "Gene.NoHairpins_Above99Percentile")])
cnts.out <- cbind(hairpinID=rownames(cnts.df),
                  pool=shRNA_lib2[rownames(cnts.df),"batch"],
                  as.data.frame(cnts.df)
)
cpm.out <- cbind(hairpinID=rownames(cpm.df),
                 pool=shRNA_lib2[rownames(cpm.df),"batch"],
                 as.data.frame(cpm.df)
)

#---- OUTPUT #### 
# CNTs
write.table(x=cnts.out,file=paste0(Out_tables.dir,lib.tag,"_cnts.txt"),sep="\t",row.names=FALSE,quote=FALSE)
write.table(x=cpm.out,file=paste0(Out_tables.dir,lib.tag,"_cpm.txt"),sep="\t",row.names=FALSE,quote=FALSE)

# TMZ100 vs DMSO
write.table(x = TMZvsDMSO.stard_rnk,file = paste0(Out_tables.dir,lib.tag,"_TMZ_standard_rnk.txt"),
            row.names = FALSE,sep = "\t",quote = FALSE)
write.table(x = tmz_vs_t1w.stard_genernk,file = paste0(Out_tables.dir,lib.tag,"_TMZ_standard_GeneRNK.txt"),
            row.names = FALSE,sep = "\t",quote = FALSE)



pdf(file = paste0(Out_tables.dir,lib.tag,"_TMZ_top25RNK.pdf"),onefile = T)
par(mfrow=c(4,1),cex.lab=0.001,cex=0.8,oma=c(3,1,5,0.5),mar=c(0,0,0,0))
positionInRanking <- 0

for(g in unique(as.character((TMZvsDMSO.stard_rnk$Symbol)))[1:25]) {
  positionInRanking <- positionInRanking+1
  g.ids <- grep(paste0("^",g,"_"),rownames(TMZvsDMSO))
  
  col.features <- as.list(topo.colors(length(g.ids)))
  names(col.features) <- rownames(TMZvsDMSO)[g.ids]
  
  textplot(format(TMZvsDMSO[g.ids,c("hairpinID","TMZ_avgCPM","DMSO_avg_CPM")],
                  digits=3),col.colnames="black",
           show.rownames = FALSE,show.colnames = TRUE,
           halign = "left",
           mar = c(0,1,1,1),cex = 0.8)
  par(new=TRUE,mar=c(1,1,1,3))
  #   plot(1000,1000, type="n", axes=F, xlab="", ylab="",)
  #   legend("right",legend = names(col.features),
  #          col=unlist(col.features),pch = 20,cex=0.7,y.intersp = ifelse(nchar(g)>5,1.5,0.5));
  legend("right",legend = names(col.features),
         col=unlist(col.features),pch = 20,cex=0.8);
  
  #   legend("bottomleft",legend = names(col.features),col=unlist(col.features),pch = 20,cex=0.5,horiz = T);
  
  par(mar=c(3,2,3,2))
  sapply(c("E1_log2FC","E2_log2FC","E3_log2FC"),function(s) {
    tmp <- rownames(TMZvsDMSO[g.ids,c("Symbol",s)])[order(TMZvsDMSO[g.ids,s],na.last = TRUE, decreasing = TRUE)]
    barcodeplot2(TMZvsDMSO[,s], index = g.ids,
                 col = unlist(col.features[tmp]),
                 main = gsub("_log2FC","",s),
                 labels = c("Positive log2FC", "Negative log2FC"))
  })
  title(g,outer=T,cex.main=1.5);
  title(paste0("\n\n\n\n\n","#",positionInRanking,
               " Ranking ::   ","Lost hairpins (TMZ Vs. DMSO) ",
               "\n\n","sScore: ",round(score$sScore[g],2),
               "; ","RiScore: ",round(score$RiScore[g],2)
  ),
  outer=T,cex.main=0.8);
}
dev.off()
rm(col.features,g)
# All
pdf(file = paste0(Out_tables.dir,lib.tag,"_TMZ_visual_ranking.pdf"),onefile = T)
par(mfrow=c(4,1),cex.lab=0.001,cex=0.8,oma=c(3,1,5,0.5),mar=c(0,0,0,0))
positionInRanking <- 0

for(g in unique(as.character((TMZvsDMSO.stard_rnk$Symbol)))) {
  positionInRanking <- positionInRanking+1
  g.ids <- grep(paste0("^",g,"_"),rownames(TMZvsDMSO))
  
  col.features <- as.list(topo.colors(length(g.ids)))
  names(col.features) <- rownames(TMZvsDMSO)[g.ids]
  
  textplot(format(TMZvsDMSO[g.ids,c("hairpinID","TMZ_avgCPM","DMSO_avg_CPM")],
                  digits=3),col.colnames="black",
           show.rownames = FALSE,show.colnames = TRUE,
           halign = "left",
           mar = c(0,1,1,1),cex = 0.8)
  par(new=TRUE,mar=c(1,1,1,3))
  #   plot(1000,1000, type="n", axes=F, xlab="", ylab="",)
  #   legend("right",legend = names(col.features),
  #          col=unlist(col.features),pch = 20,cex=0.7,y.intersp = ifelse(nchar(g)>5,1.5,0.5));
  legend("right",legend = names(col.features),
         col=unlist(col.features),pch = 20,cex=0.8);
  
  #   legend("bottomleft",legend = names(col.features),col=unlist(col.features),pch = 20,cex=0.5,horiz = T);
  
  par(mar=c(3,2,3,2))
  sapply(c("E1_log2FC","E2_log2FC","E3_log2FC"),function(s) {
    tmp <- rownames(TMZvsDMSO[g.ids,c("Symbol",s)])[order(TMZvsDMSO[g.ids,s],na.last = TRUE, decreasing = TRUE)]
    barcodeplot2(TMZvsDMSO[,s], index = g.ids,
                 col = unlist(col.features[tmp]),
                 main = gsub("_log2FC","",s),
                 labels = c("Positive log2FC", "Negative log2FC"))
  })
  title(g,outer=T,cex.main=1.5);
  title(paste0("\n\n\n\n\n","#",positionInRanking,
               " Ranking ::   ","Lost hairpins (TMZ Vs. DMSO) ",
               "\n\n","sScore: ",round(score$sScore[g],2),
               "; ","RiScore: ",round(score$RiScore[g],2)
  ),
  outer=T,cex.main=0.8);
}
dev.off()
rm(col.features,g)



# Selected target
png(filename = paste0(Out_plots.dir,lib.tag,"_FCdistrib2_selectedtargets.png"),width = 1800,height = 900,res=180)

selected.tar <- c("E2F1","PSMB7")
colors_selected.tar <- c("E2F1"="brown",
                         "PSMB7"="darkblue")

tmp <- TMZvsDMSO[order(TMZvsDMSO$E123_avgLog2FC),]
colors.bar <- rep("grey",nrow(tmp))

bin.bar <- 2
for(selected.gene in selected.tar) {
  row.id <- as.numeric(which(tmp$Symbol==selected.gene))
  for(row.i in row.id) {
    row.i <- ifelse(row.i<bin.bar,bin.bar,row.i)
    colors.bar[(row.i-bin.bar):(row.i+bin.bar)] <- colors_selected.tar[selected.gene]
  }
}
barplot(tmp$E123_avgLog2FC,border = colors.bar,ylab="log2(FC)",
        axisnames = FALSE,names.arg = FALSE,
        main="TMZ Vs. DMSO : lost hairpins");
legend("topleft",
       legend = sapply(selected.tar,function(z) paste0(z," (",
                                                       " sScore= ",round(score$sScore[z],digits=1)," ; ",
                                                       # "RiScore= ",round(score$RiScore[z],digits=1),
                                                       # ") "," (",
                                                       sum(grepl(paste0("^",z,"_"),shRNA_lib2$ID)),
                                                       " hairpins; ",
                                                       unique(shRNA_lib[grepl(paste0("^",z,"_"),shRNA_lib2$ID),"batch"]),
                                                       ")")),
       fill = colors_selected.tar)
abline(h=depletion.log2cutoff,col="blue");
text(x = length(tmp$E123_avgLog2FC)/3,y = depletion.log2cutoff+0.2,
     labels = paste0(sum(tmp$E123_avgLog2FC < depletion.log2cutoff),
                     " hairpins with <",depletion.tag,"-Fold-Change"),col="blue")
dev.off()

# Distribution of scores
png(filename = paste0(Out_plots.dir,lib.tag,"_scoresDistrib.png"),width = 700*1.5*3,height = 450*3,res=230)
par(mfrow=c(1,2),las=1)
hist(score$sScore,main="Distribution of Sensitizer Score (sScore)",ylab="Frequency",xlab="sScore",xlim=c(-4,0));
abline(v=log10(0.01),lwd=2,col="blue");
text(x = log10(0.01), y= 250,labels = "Significant features",cex=1.2,col="blue",pos = 2)
text(x = log10(0.01), y= 250,labels = paste0("\n\n\n(#",sum(score$sScore<log10(0.01))," Target genes)"),cex=0.8,col="blue",pos = 2)
text(x = log10(0.01), y= 250,labels = "Non-significant features",cex=1.2,col="black",pos = 4)

hist(score$RiScore,main="Distribution of RNA interference Score (RiScore)",
     ylab="Frequency",xlab="RiScore",xlim=c(-4,4));
abline(v=log10(0.05),lwd=2,col="blue");
abline(v=-log10(0.05),lwd=2,col="red");
text(x = log10(0.05), y= 80,labels = "Significant features",cex=1.2,col="blue",pos = 2)
text(x = log10(0.05), y= 80,labels = paste0("\n\n\n(#",sum(score$RiScore<log10(0.05))," Target genes)"),cex=0.8,col="blue",pos = 2)
text(x = -log10(0.05), y= 80,labels = "Significant features",cex=1,col="red",pos = 4)
text(x = -log10(0.05), y= 80,labels = paste0("\n\n\n(#",sum(score$RiScore> -log10(0.05))," Target genes)"),cex=0.8,col="red",pos = 4)
text(x = mean(c(log10(0.05),-log10(0.05))), y= 90,labels = "Non-significant features",cex=1.1,col="black")
dev.off()


# modeling
score.df <- data.frame(sScore=score$sScore[names(score$RiScore)],
                       RiScore=score$RiScore)
score.df$Significant <- apply(score.df,1,function(z) ifelse(z[1]<log10(0.01) & z[2]<log10(0.05),"YES","NO"))
score.df$Gene <- rownames(score.df)

# plot(score.df$sScore,score.df$RiScore,xlab="sScore",ylab="RiScore")
# rect(xleft=log10(0.01),ybottom=-10,xright=-10,ytop=log10(0.05),border = "red")
# text(subset(score.df,sScore<log10(0.01) & RiScore<log10(0.05))$sScore,
#      subset(score.df,sScore<log10(0.01) & RiScore<log10(0.05))$RiScore,
#      labels = rownames(subset(score.df,sScore<log10(0.01) & RiScore<log10(0.05))),cex=0.8,col="red2")
png(filename = paste0(Out_plots.dir,lib.tag,"_scatterplot_scores.png"),width = 600*3,height = 450*3,res=240)
library(ggrepel)
p = ggplot(score.df, aes(sScore, RiScore)) +
  geom_point(aes(col=Significant),size=1) +
  scale_color_manual(values=c("black", "blue")) + stat_smooth(level=0.99999)

p+geom_text_repel(data=subset(score.df, Significant=="YES"),size=4,
                  aes(label=Gene)) + theme_bw() +
  ggtitle("Scatter plot :: sScore ~ RiScore")
dev.off()


