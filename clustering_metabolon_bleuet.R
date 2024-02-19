######################################################################
#                     METABOLON CLUSTERING ANALYSIS                  # 
######################################################################

library(mixOmics)
library(tidyverse)
library(pvclust)
library(ggplot2)
library(pvclust)

set.seed(66)

#setwd("/Users/juan/Library/CloudStorage/OneDrive-UniversitéLaval/INAF/VALENTIN/scripts valentin/bleuet/")

meta.bleuet<-read.csv("data/meta_bleuets.csv", row.names = 1)
metabolites.bleuet<-read.csv("data/chem.csv", row.names = 1)
log.bleuet<-read.csv("data/log_bleuets.csv", row.names = 1)

colnames(log.bleuet) <- metabolites.bleuet$PLOT_NAME
rownames(log.bleuet)<-meta.bleuet$SUBJECT_OR_ANIMAL_ID

### IN CASE WE HAVE TO REMOVE OUTLIERS WE HAVE TO DO IT *BEFORE* UNIQUE AND 10% VARIANCE REMOVAL
outliers <- c("LARC-00295","LARC-00296","LARC-00323","LARC-00324")#,"LARC-00279", "LARC-00280")
#outliers <- c("NONE")
meta.bleuet<-subset(meta.bleuet,!rownames(meta.bleuet) %in% outliers)
log.bleuet<-subset(log.bleuet,rownames(log.bleuet) %in% meta.bleuet$SUBJECT_OR_ANIMAL_ID)
dim(log.bleuet)

### REMOVING METABOLITES WITH UNIQUE VALUES
varx <- function(dat) {
  out <- lapply(dat, function(x) length(unique(x)))
  want <- which(!out > 1)
  unlist(want)}

log.bleuet.dropsx<-as.vector(varx(log.bleuet))
log.bleuetii <- log.bleuet[-(log.bleuet.dropsx)]

### REMOVING METABOLITES WITHOUT NAME DESCRIPTION
log.bleuetii <-log.bleuetii[, !grepl("X-", names(log.bleuetii))]

dim(log.bleuet)
dim(log.bleuetii)

### REMOVING METABOLITES WITH LESS THAN 10% VARIANCE 
# calculate variance per column
variances <- apply(X=log.bleuetii, MARGIN=2, FUN=var)
vardf<-as.data.frame(variances)
vardf0.1 <- subset(vardf,variances>0.1)
dim(vardf0.1)
log.bleuet2 <- log.bleuetii[colnames(log.bleuetii) %in% rownames(vardf0.1)]
dim(log.bleuet2)
mean(sapply(log.bleuet2,var)) #MEAN VARIANCE OF THE RESULTING DB

# REMOVING METABOLITES DEFINED AS Xenobiotics
xeno <- subset(metabolites.bleuet,metabolites.bleuet$SUPER_PATHWAY == "Xenobiotics")
log.bleuet3 <-log.bleuet2[, !names(log.bleuet2) %in% xeno$PLOT_NAME]
dim(log.bleuet3)

# REMOVING METABOLITES DEFINED AS Partially Characterized Molecules
partial <- subset(metabolites.bleuet,metabolites.bleuet$SUPER_PATHWAY == "Partially Characterized Molecules")
log.bleuet4 <-log.bleuet3[, !names(log.bleuet3) %in% partial$PLOT_NAME]
dim(log.bleuet4)


####### PLSDA AND HIERARCHICAL CLUSTERING #####################

# multilevel PLS-DA
plsda.multilevel <- plsda(log.bleuet4, scale = F,
                          meta.bleuet$TIME_POINT, 
                          multilevel = meta.bleuet$SUBJECT_ID, 
                          ncomp = 10)

plotIndiv(plsda.multilevel, ind.names = T, legend=TRUE,
          ellipse = TRUE)#, style="3d", ncomp=3)

# undergo performance evaluation in order to tune the number of components to use
set.seed(66)  
perf.plsda.multilevel<- perf(plsda.multilevel, validation = "Mfold", 
                           folds = 5, nrepeat = 20, # use repeated cross-validation
                           progressBar = T, auc = TRUE) # include AUC values

# plot the outcome of performance evaluation across all ten components
#pdf("cluster/bleuet/only_0.2+_variance/plsda_perf_bleuet.pdf",width=6, height=5)
plot(perf.plsda.multilevel, sd = TRUE,
     legend.position = "vertical",
     dist = c("max.dist"),
     measure = c("overall"))
#dev.off()

#pdf("cluster/bleuet/only_0.2+_variance/plsda_scree_bleuet.pdf",width=6, height=5)
expvar<-plsda.multilevel$prop_expl_var$X
names(expvar) <- gsub("comp","",names(expvar))
bar<-barplot(expvar,ylim=c(0,0.15),ylab="Proportion of explained variance",xlab="Component")
lines(x = bar, y = plsda.multilevel$prop_expl_var$X)
points(x = bar, y = plsda.multilevel$prop_expl_var$X)
text(x = bar, y = round(expvar,2), label = round(expvar,2), pos = 3, cex = 0.8, col = "red")
box()
#dev.off()

#pdf("cluster/bleuet/only_0.2+_variance/plsda_bleuet.names.xeno.pdf",width=6, height=5)
plotIndiv(plsda.multilevel,style = "graphics",ellipse.level = 0.95,#pch=19,
          ind.names = meta.bleuet$SUBJECT_OR_ANIMAL_ID, ellipse=T,
          title = "",
          X.label = paste("Component 1 ", "(",round(plsda.multilevel$prop_expl_var$X[1:1]*100,0),"% )"), 
          Y.label = paste("Component 2 ", "(",round(plsda.multilevel$prop_expl_var$X[2:2]*100,0),"% )"),
          group = meta.bleuet$TIME_POINT, # col.per.group=c("black","darkorange1"),
          legend=T,point.lwd=1, cex=1, size.xlabel = rel(1.2),
          size.ylabel = rel(1.2),size.axis = rel(1.2))
#dev.off()

# LOADINGS
plotLoadings(plsda.multilevel, comp=1,ndisplay = 10)
plotLoadings(plsda.multilevel, comp=2,ndisplay = 10)


# ON UTILISE LES VARIANTES DES DEUX PREMIERS COMPOSANTS DU PLSDA POUR FAIRE LE HIERARCHICAL CLUSTERING
vars1<-plsda.multilevel$variates$X[,1:2] 

#Generate distance matrix -> find distance matrix # euclidean
d <- dist(as.matrix(vars1),method="euclidean")

##Hierarchical clustering -> find clustering method
hc <- hclust(d,method="ward.D")     
plot(hc)
rhc <- rect.hclust(hc, k = 4, border = c("red","blue"))
NONRESPONDERS <- as.numeric(gsub("v.*","",names(rhc[[1]])))

#### ombrage ---> https://arxiv.org/pdf/1411.5259.pdf
str(vars1)
comp11<-rep((0), times = length(meta.bleuet$SUBJECT_OR_ANIMAL_ID)) #dummy variable
comp12<-rep((0), times = length(meta.bleuet$SUBJECT_OR_ANIMAL_ID)) #dummy variable
vars2<-cbind(vars1,comp11,comp12)
vars3<-t(vars2[c(1,2,3,4)])
vars3<-t(vars2)
fit <- pvclust(vars3, method.hclust="ward.D",method.dist="euclidean")

# dendogram with p values

#pdf("cluster/bleuet/only_0.2+_variance/pvclust_bleuet.pdf",width=6, height=5)
plot(fit, pch=19, print.pv=c("au"),print.num=F,
     xlab = "", ylab = "Euclidean distance", sub = "", main="",
     cex.main=1, cex.lab=1, cex.axis=1,cex.pv=1,lwd=1)
#dev.off()

plot(fit, pch=19, print.pv=c("au"),print.num=F,
     xlab = "", ylab = "Euclidean distance", sub = "", main="",
     cex.main=1.2, cex.lab=1, cex.axis=1,cex.pv=0.5,lwd=1, cex = 0.6)
rect.hclust(hc, k = 4, border = c('#E41A1C','#377EB8'))#,"red","orange","purple","orange"))


NONRESPONDERS #keepx: 160 10 200 when removing outliers, 40 140 10 when keeping them
nonresponders.log.bleuet6.1.VALE <- c("58","71","28","68","107","43","106","96","111")
nonresponders.log.bleuet6.1.JUAN <- c("58","71","28","68","107","46","102","4","88")
nonresponders.log.bleuet6.1.JUAN2 <- c("71","28","68","107","102")
group.noscale.outliers <- c(59,79,43,78,57,23,110) # 60 60 40 / 20 10 100 ? nope erreur
group.noscale.no.outliers <- c(61,28,71,68,106) # 200 20 120

nonresponders8 <- c(102,107,28,68,71,46,58,4,88)
nonresponders9 <- c(102,107,28,68,71,46,58,4,88,23,59)

# group.noscale.no.outliers <- NONRESPONDERS


# multilevel SPLS-DA
meta.bleuet3 <- meta.bleuet
meta.bleuet3$SUBJECT_ID <- gsub("B","",meta.bleuet3$SUBJECT_ID)

# nonresponders<-nonresponders.log.bleuet6.1.VALE
meta.bleuet3$CLUSTER <- NA
meta.bleuet3[meta.bleuet3$SUBJECT_ID %in% group.noscale.no.outliers, "CLUSTER"] <- 1
meta.bleuet3[!meta.bleuet3$SUBJECT_ID %in% group.noscale.no.outliers, "CLUSTER"] <- 2
meta.bleuet3$CLUSTER_POINT <- paste0(meta.bleuet3$CLUSTER,meta.bleuet3$TIME_POINT)


set.seed(66) # for reproducbility in this vignette, otherwise increase nrepeat
list.keepX <- c(10,  seq(20, 200, 20))

tune.splsda.clust <- tune.splsda(log.bleuet4, 
                                 as.factor(meta.bleuet3$CLUSTER_POINT), #visite X cluster
                                 # as.factor(meta.bleuet3$CLUSTER), #cluster only
                                 ncomp = 3,
                                 validation = 'Mfold',
                                 folds = 10, dist = 'max.dist', progressBar = T,
                                 measure = "BER", test.keepX = list.keepX,
                                 nrepeat = 50, multilevel=meta.bleuet3$SUBJECT_ID)   # we suggest nrepeat = 50

error <- tune.splsda.clust$error.rate
ncomp <- tune.splsda.clust$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
select.keepX <- tune.splsda.clust$choice.keepX  # optimal number of variables to select
select.keepX

plot(tune.splsda.clust, col = color.jet(3))

MyResult.splsda.final <- splsda(log.bleuet4,
                                as.factor(meta.bleuet3$CLUSTER_POINT), #visite X cluster
                                # as.factor(meta.bleuet3$CLUSTER), #cluster only
                                ncomp = 3, keepX = c(200,20),#,120),#select.keepX,
                                multilevel=meta.bleuet3$SUBJECT_ID)

plotIndiv(MyResult.splsda.final, ind.names = T, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA MODEL - CLUSTERING")#, style="3d")

perf.splsda.multilevel<- perf(MyResult.splsda.final, validation = "Mfold", 
                              folds = 5, nrepeat = 50, # use repeated cross-validation
                              progressBar = T, auc = TRUE) # include AUC values


plot(perf.splsda.multilevel, sd = TRUE,
     legend.position = "vertical",
     dist = c("max.dist"),
     measure = c("overall"))


# LOADINGS
plotLoadings(MyResult.splsda.final, comp=1,ndisplay = 20)
plotLoadings(MyResult.splsda.final, comp=2,ndisplay = 20)


#Loadings comp 1
plots.splsda1<-plotLoadings(MyResult.splsda.final,comp=1)#,ndisplay = 10)
loads.splsda <-as.data.frame(MyResult.splsda.final$loadings$X)
loadsDA1 <- subset(loads.splsda,loads.splsda$comp1 %in% plots.splsda1$importance,comp1, drop=F)
loadsDA1 <-loadsDA1[order(abs(loadsDA1$comp1)), ,drop=F]

loads.splsda.fam <- subset(metabolites.bleuet,metabolites.bleuet$PLOT_NAME %in% rownames(loadsDA1),
                           c("PLOT_NAME","SUPER_PATHWAY","SUB_PATHWAY","HMDB", "KEGG","PUBCHEM"))

loads.splsda.fam <-loads.splsda.fam[order(loads.splsda.fam$PLOT_NAME),]
loadsDA1 <-loadsDA1[order(rownames(loadsDA1)), ,drop=F]
loadsDA1.fam <- cbind(loadsDA1,loads.splsda.fam[,-1])
loadsDA1.fam <-loadsDA1.fam[order(loadsDA1.fam$SUPER_PATHWAY),]
loadsDA1.fam


# Classer loads en liste par ordre de valeur absolue ####################
df <- loadsDA1
original_row_names <- rownames(df)
ordered_df <- df[order(abs(df$comp1)), , drop = FALSE]
rownames(ordered_df) <- original_row_names[order(abs(df$comp1))]
ordered_df
# write.csv(ordered_df,"ordered_df.csv")
# END ######################

#loadings comp 2
plots.splsda2<-plotLoadings(MyResult.splsda.final,comp=2)#,ndisplay = 10)
loads.splsda <-as.data.frame(MyResult.splsda.final$loadings$X)
loadsDA2 <- subset(loads.splsda,loads.splsda$comp2 %in% plots.splsda2$importance,comp2, drop=F)
loadsDA2 <-loadsDA2[order(abs(loadsDA2$comp2)), ,drop=F]

loads.splsda.fam2 <- subset(metabolites.bleuet,metabolites.bleuet$PLOT_NAME %in% rownames(loadsDA2),
                            c("PLOT_NAME","SUPER_PATHWAY","SUB_PATHWAY","HMDB", "KEGG","PUBCHEM"))

loads.splsda.fam2 <-loads.splsda.fam2[order(loads.splsda.fam2$PLOT_NAME),]
loadsDA2 <-loadsDA2[order(rownames(loadsDA2)), ,drop=F]
loadsDA2.fam <- cbind(loadsDA2,loads.splsda.fam2[,-1])
loadsDA2.fam <-loadsDA2.fam[order(loadsDA2.fam$SUPER_PATHWAY),]
loadsDA2.fam

# Classer loads en liste par ordre de valeur absolue ####################
df <- loadsDA2
original_row_names <- rownames(df)
ordered_df <- df[order(abs(df$comp2)), , drop = FALSE]
rownames(ordered_df) <- original_row_names[order(abs(df$comp2))]
ordered_df

# write.csv(ordered_df,"ordered_df.csv")
# END ######################


names.loadsDA <- cbind(as.data.frame(rownames(loadsDA1.fam)),as.data.frame(rownames(loadsDA2.fam)))
names(names.loadsDA) <- c('loadsDA1','loadsDA2')
names.loadsDA

### LOADINGS FOR CORRELATION ANALYSIS
loadsDA1.cor <-loadsDA1[order(abs(loadsDA1$comp1), decreasing=T), ,drop=F]
loadsDA2.cor <-loadsDA2[order(abs(loadsDA2$comp2), decreasing=T), ,drop=F]
loads <- cbind(rownames(loadsDA1.cor), rownames(loadsDA2.cor))
colnames(loads) <- c("comp1","comp2")
rownames(loads) <- c(1:length(loads[,1]))

### PERFORMANCE ERROR RATE
perf.splsda.multilevel$error.rate$overall[2,1]


### HEATMAPS
# set the colours used for the subject assocaited 
# with each sample (left-most column)

col.ID <- c("orange", "red","purple","#377EB8")[meta.bleuet3$CLUSTER_POINT]
cim.colors <- as.data.frame(meta.bleuet3$CLUSTER_POINT)
colnames(cim.colors) <- "Cluster"
cim.colors$color <- ifelse(cim.colors$Cluster == "12","orange",ifelse(
  cim.colors$Cluster == "11","red",ifelse(
    cim.colors$Cluster == "22","purple","#377EB8"
  )
))


#MyResult.splsda.final$loadings
cim(MyResult.splsda.final, comp=1:2,
    row.sideColors = cim.colors$color, 
    row.names = F,#row.names = meta.rami3$CLUSTER_POINT2,
    col.names = FALSE, legend=list(legend = c("R1 Pre","R1 Post","R2 Pre","R2 Post"),
                                   col = c("red", "orange","#377EB8","purple"),
                                   title = "Cluster timepoint", cex = 1.0))
