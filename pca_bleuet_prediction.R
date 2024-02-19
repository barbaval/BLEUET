library(caret) #important de loader caret en premier
library(mixOmics)

#setwd("your_work_directory_here")

meta.rami<-read.csv("data//meta_bleuets.csv", row.names = 1)

metabolites.rami<-read.csv("data//chem.csv", row.names = 1)
metabolites.rami1 <- metabolites.rami[1:5,1:12]

log.rami<-read.csv("data//log_bleuets.csv", row.names = 1)
log.rami1 <- log.rami[1:5,1:100]

colnames(log.rami) <- metabolites.rami$PLOT_NAME

var <- function(dat) {
  out <- lapply(dat, function(x) length(unique(x)))
  want <- which(!out > 1)
  unlist(want)}

log.rami.drops<-as.vector(var(log.rami))
log.ramiii <- log.rami[-(log.rami.drops)]
log.ramiii <-log.ramiii[, !grepl("X-", names(log.ramiii))]
#log.ramiii <- log.rami[ -c(177,427,489,501,502,508,551,552,779,804,849,860) ]

dim(log.rami)
dim(log.ramiii)



# undergo pca and plot samples without any multilevel decomposition
pca.result1 <- pca(log.ramiii, scale = T, center = T, ncomp=10)
plotIndiv(pca.result1, 
          ind.names = rownames(log.ramiii), 
          group = meta.rami$TIME_POINT, 
          title = 'Figure 1: log rami')
plot(pca.result1, ncomp=10)
plotLoadings(pca.result1, ndisplay = 20)

# undergo pca and plot samples with multilevel decomposition
pca.result2 <- pca(log.ramiii, multilevel = meta.rami$SUBJECT_ID, scale = T, center = T, ncomp=10)
plotIndiv(pca.result2, 
          ind.names = rownames(log.ramiii), 
          group = meta.rami$TIME_POINT, 
          legend = TRUE, 
          legend.title = "Visit", 
          title = 'Figure 2: Multilevel PCA on log rami')
plot(pca.result2, ncomp=10)
plots<-plotLoadings(pca.result2, ndisplay = 20)

loads<-as.data.frame(pca.result2$loadings$X)
loadsPC1 <- subset(loads,loads$PC1 %in% plots$importance,PC1, drop=F)
loadsPC1 <-loadsPC1[order(abs(loadsPC1$PC1)), ,drop=F]

loads.fam <- subset(metabolites.rami,metabolites.rami$PLOT_NAME %in% rownames(loadsPC1),
                  c("PLOT_NAME","SUPER_PATHWAY","SUB_PATHWAY","HMDB", "KEGG","PUBCHEM"))

loads.fam <-loads.fam[order(loads.fam$PLOT_NAME),]
loadsPC1 <-loadsPC1[order(rownames(loadsPC1)), ,drop=F]
loadsPC1.fam <- cbind(loadsPC1,loads.fam[,-1])
loadsPC1.fam <-loadsPC1.fam[order(loadsPC1.fam$SUPER_PATHWAY),]


########## PLS-DA

MyResult.plsda <- plsda(log.ramiii, meta.rami$TIME_POINT) # 1 Run the method
plotIndiv(MyResult.plsda, title = "PLSDA on log rami")     
plotVar(MyResult.plsda)                 
plotVar(MyResult.plsda, cutoff = 0.75)
selectVar(MyResult.splsda, comp=1)$name 


########## SPLS-DA
#remove the outliers
outliers <- c("LARC-00295","LARC-00296","LARC-00323","LARC-00324")
log.ramiii2 <- subset(log.ramiii,!rownames(log.ramiii) %in% outliers)
dim(log.ramiii2)
meta.rami2 <- subset(meta.rami,!rownames(meta.rami) %in% outliers)
dim(meta.rami2)

# log.ramiii2.drops<-as.vector(var(log.ramiii2))
# log.ramiii2 <- log.ramiii2 [-(log.ramiii2.drops)]
# dim(log.ramiii2)


set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
list.keepX <- c(5:10,  seq(20, 100, 10))

tune.splsda <- tune.splsda(log.ramiii, as.factor(meta.rami$TIME_POINT), ncomp = 3,
                                 validation = 'Mfold',
                                 folds = 3, dist = 'max.dist', progressBar = T,
                                 measure = "BER", test.keepX = list.keepX,
                                 nrepeat = 50, multilevel=meta.rami$SUBJECT_ID)   # we suggest nrepeat = 50

error <- tune.splsda$error.rate
ncomp <- tune.splsda$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
select.keepX <- tune.splsda$choice.keepX  # optimal number of variables to select
select.keepX
#ncomp <- 3
plot(tune.splsda, col = color.jet(3))


MyResult.splsda.final <- splsda(log.ramiii2, meta.rami2$TIME_POINT, ncomp = 2, keepX = c(5,100),
                                multilevel=meta.rami2$SUBJECT_ID)

plotIndiv(MyResult.splsda.final, ind.names = T, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA - final result")

plots.splsda<-plotLoadings(MyResult.splsda.final,comp=1)

loads.splsda<-as.data.frame(MyResult.splsda.final$loadings$X)
loadsDA1 <- subset(loads.splsda,loads.splsda$comp1 %in% plots.splsda$importance,comp1, drop=F)
loadsDA1 <-loadsDA1[order(abs(loadsDA1$comp1)), ,drop=F]

loads.splsda.fam <- subset(metabolites.rami,metabolites.rami$PLOT_NAME %in% rownames(loadsDA1),
                    c("PLOT_NAME","SUPER_PATHWAY","SUB_PATHWAY","HMDB", "KEGG","PUBCHEM"))

loads.splsda.fam <-loads.splsda.fam[order(loads.splsda.fam$PLOT_NAME),]
loadsDA1 <-loadsDA1[order(rownames(loadsDA1)), ,drop=F]
loadsDA1.fam <- cbind(loadsDA1,loads.splsda.fam[,-1])
loadsDA1.fam <-loadsDA1.fam[order(loadsDA1.fam$SUPER_PATHWAY),]

###methyl glucopyranoside
#https://content.iospress.com/articles/journal-of-berry-research/jbr135
#For example, a metabolic evaluation of liver tissue shows accumulation of raspberry-related metabolites including methyl glucopyranoside 
#and S-methylmethionine only when whole raspberry products were provided to mice fed a HF diet 

###pentose acid ---> valeric acid
#https://pubs.rsc.org/en/content/articlelanding/2018/fo/c7fo02061a/unauth
#Our data indicate that RA consumption significantly increases the concentration of fecal
#acetic acid, propionic acid, and valeric acid and remarkably elevates the butyric acid level.



####################################
#meta.rami<-read.csv("meta_rami.csv", row.names = 1)
outliers <- c("LARC-00295","LARC-00296","LARC-00323","LARC-00324")
meta.rami <- subset(meta.rami,!rownames(meta.rami) %in% outliers)

meta.rami.before <- subset(meta.rami,meta.rami$TIME_POINT == 1)

set.seed(30)
inTrain <- createDataPartition(
  y = meta.rami.before$GENDER,
  p = .5,
  list = FALSE,
  times = 1)

rami.train <- meta.rami.before[ inTrain,]
rami.test  <- meta.rami.before[-inTrain,]

prop.table(table(rami.train$GENDER))
prop.table(table(rami.test$GENDER))

prop<-prop.table(table(meta.rami.before$GENDER))
prop

meta.rami.train <- subset(meta.rami,meta.rami$SUBJECT_ID %in% rami.train$SUBJECT_ID)
meta.rami.test <- subset(meta.rami,meta.rami$SUBJECT_ID %in% rami.test$SUBJECT_ID)


log.ramiii.train <- subset(log.ramiii,rownames(log.ramiii) %in% rownames(meta.rami.train))
log.ramiii.test <- subset(log.ramiii,rownames(log.ramiii) %in% rownames(meta.rami.test))

log.rami.train.drops<-as.vector(var(log.ramiii.train))
log.rami.test.drops<-as.vector(var(log.ramiii.test))
log.rami.all.drops <- c(log.rami.train.drops,log.rami.test.drops)

log.ramiii.train <- log.ramiii.train[-(log.rami.all.drops)]
log.ramiii.test <- log.ramiii.test[-(log.rami.all.drops)]

dim(log.ramiii.train)


set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
list.keepX <- c(5:10,  seq(20, 100, 10))

tune.splsda.train <- tune.splsda(log.ramiii.train, as.factor(meta.rami.train$TIME_POINT), ncomp = 3,
                                validation = 'Mfold',
                                folds = 3, dist = 'max.dist', progressBar = T,
                                measure = "BER", test.keepX = list.keepX,
                                nrepeat = 50, multilevel=meta.rami.train$SUBJECT_ID)   # we suggest nrepeat = 50

error <- tune.splsda.train$error.rate
ncomp.train <- tune.splsda.train$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
select.keepX.train <- tune.splsda.train$choice.keepX  # optimal number of variables to select


# train the model

train.splsda <- splsda(log.ramiii.train, meta.rami.train$TIME_POINT, ncomp = 2, keepX = c(20,60),
                       multilevel=meta.rami.train$SUBJECT_ID)

plotIndiv(train.splsda, ind.names = T, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA - TRAINING")


plotLoadings(train.splsda,comp=1)
plotLoadings(train.splsda,comp=2)


# use the model on the Xtest set
predict.splsda <- predict(train.splsda, log.ramiii.test,dist = "mahalanobis.dist")
                         # multilevel=meta.rami.test$SUBJECT_ID)


# evaluate the prediction accuracy for the first two components
predict.comp1 <- predict.splsda$class$mahalanobis.dist[,1]
conf.matrix<-table(factor(predict.comp1, levels=levels(time)), meta.rami.test$TIME_POINT)         
conf.matrix
BER<-get.BER(conf.matrix)
(Accuracy<-round(1-BER,2))

plot(log.ramiii$`guaiacol sulfate` ~ meta.rami$TIME_POINT)


