library(mixOmics)

#setwd("/Users/juan/Library/CloudStorage/OneDrive-UniversitéLaval/INAF/VALENTIN/pca")

meta.bleuet<-read.csv("data/meta_bleuets.csv", row.names = 1)

metabolites.bleuet<-read.csv("data/chem.csv", row.names = 1)
metabolites.bleuet1 <- metabolites.bleuet[1:5,1:12]

# bleuets<-read.csv("data/log_bleuets.csv", row.names = 1)
# bleuets1 <- bleuets[1:5,1:5]

log.bleuets<-read.csv("data/log_bleuets.csv", row.names = 1)
# log.bleuets1 <- log.bleuets[1:5,1:100]

colnames(log.bleuets) <- metabolites.bleuet$PLOT_NAME

var <- function(dat) {
  out <- lapply(dat, function(x) length(unique(x)))
  want <- which(!out > 1)
  unlist(want)}

log.rami.drops<-as.vector(var(log.bleuets))
log.ramiii <- log.bleuets[-(log.rami.drops)]
log.ramiii <-log.ramiii[, !grepl("X-", names(log.ramiii))]
#log.ramiii <- log.rami[ -c(177,427,489,501,502,508,551,552,779,804,849,860) ]

dim(log.bleuets)
dim(log.ramiii)

log.bleuetsx <-log.ramiii[, names(log.ramiii) %in% xeno$PLOT_NAME]

# undergo pca and plot samples without any multilevel decomposition
pca.result1 <- pca(log.ramiii, scale = T, center = T, ncomp=10)
plotIndiv(pca.result1, 
          ind.names = rownames(log.ramiii), 
          group = meta.bleuet$TIME_POINT, 
          title = 'Figure 1: log rami')
plot(pca.result1, ncomp=10)
plotLoadings(pca.result1, ndisplay = 20)

# undergo pca and plot samples with multilevel decomposition
pca.result2 <- pca(log.ramiii, multilevel = meta.bleuet$SUBJECT_ID, scale = T, center = T, ncomp=10)
plotIndiv(pca.result2, 
          ind.names = rownames(log.ramiii), 
          group = meta.bleuet$TIME_POINT, 
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

MyResult.plsda <- plsda(log.ramix, meta.rami$TIME_POINT) # 1 Run the method
plotIndiv(MyResult.plsda, ellipse=T,
          legend=F)     
plotVar(MyResult.plsda)                 
plotVar(MyResult.plsda, cutoff = 0.75)
selectVar(MyResult.plsda, comp=1)$name 

plotLoadings(MyResult.plsda,comp=1,ndisplay=5)
# ,contrib = 'max', method = 'median'
########## SPLS-DA

outliers <- c("LARC-00307","LARC-00308","LARC-00347","LARC-00348")
log.ramiii2 <- subset(log.ramiii,!rownames(log.ramiii) %in% outliers)
dim(log.ramiii2)
meta.rami2 <- subset(meta.rami,!rownames(meta.rami) %in% outliers)
dim(meta.rami2)

log.ramiii2.drops<-as.vector(var(log.ramiii2))
log.ramiii2 <- log.ramiii2 [-(log.ramiii2.drops)]
dim(log.ramiii2)


set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
list.keepX <- c(5:10,  seq(20, 100, 10))

tune.splsda <- tune.splsda(log.ramix, as.factor(meta.rami$TIME_POINT), ncomp = 3,
                                 validation = 'Mfold',
                                 folds = 3, dist = 'max.dist', progressBar = T,
                                 measure = "BER", test.keepX = list.keepX,
                                 nrepeat = 50, multilevel=meta.rami$SUBJECT_ID)   # we suggest nrepeat = 50

error <- tune.splsda$error.rate
ncomp <- tune.splsda$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
select.keepX <- tune.splsda$choice.keepX  # optimal number of variables to select
select.keepX
#ncomp <- 3
plot(tune.splsda, col = c("red","blue","green"))

# SWITCHED TO INCLUDE OUTLIERS meta.rami vs meta.rami2
MyResult.splsda.final <- splsda(log.ramix, meta.rami$TIME_POINT, ncomp = 2, keepX = c(5,5),
                                multilevel=meta.rami$SUBJECT_ID)

plotIndiv(MyResult.splsda.final,style = "graphics",ellipse.level = 0.99, 
          X.label =  paste('Component 1 ', '(',round(MyResult.splsda.final$prop_expl_var$X['comp1']*100,0),'% )'),
          Y.label = paste('Component 2 ', '(',round(MyResult.splsda.final$prop_expl_var$X['comp2']*100,0),'% )'),
          ind.names = FALSE, ellipse=T, title = "", cex = 2, point.lwd=0.5, pch=16,
          col.per.group=c('#377EB8','#E41A1C'), legend = T, legend.title = "Timepoint")

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




# OUTLIER DETECTION
# https://privefl.github.io/blog/detecting-outlier-samples-in-pca/
library(bigutilsr)

X <- log.bleuet4
pca <- prcomp(X, scale. = TRUE, rank. = 10)
U <- pca$x

library(ggplot2)
theme_set(bigstatsr::theme_bigstatsr(0.8))
qplot(U[, 1], U[, 2]) + coord_equal()

apply(U, 2, function(x) which( abs(x - mean(x)) > (6 * sd(x)) ))


ind.out <- apply(U, 2, function(x) which( (abs(x - median(x)) / mad(x)) > 6 )) %>%
  Reduce(union, .) %>%
  print()

col <- rep("black", nrow(U)); col[ind.out] <- "red"
qplot(U[, 1], U[, 3], color = I(col), size = I(2)) + coord_equal()

dist <- apply(U, 2, function(x) abs(x - median(x)) / mad(x)) %>%
  apply(1, max)
qplot(U[, 1], U[, 3], color = dist, size = I(3)) + coord_equal() + 
  scale_color_viridis_c(trans = "log", breaks = c(1, 3, 6))

qplot(y = sort(dist, decreasing = TRUE)) +
  geom_hline(yintercept = 6, color = "red")

dist2 <- covRob(U, estim = "pairwiseGK")$dist
qplot(dist, sqrt(dist2))

cowplot::plot_grid(
  qplot(U[, 1], U[, 2], color = dist2, size = I(2)) + coord_equal() + 
    scale_color_viridis_c(trans = "log", breaks = NULL),
  qplot(U[, 3], U[, 7], color = dist2, size = I(2)) + coord_equal() + 
    scale_color_viridis_c(trans = "log", breaks = NULL),
  rel_widths = c(0.7, 0.4), scale = 0.95
)

pval <- pchisq(dist2, df = 10, lower.tail = FALSE)
hist(pval)      

is.out <- (pval < (0.05 / length(dist2)))  # Bonferroni correction
sum(is.out)

qplot(U[, 3], U[, 7], color = is.out, size = I(3)) + coord_equal()


llof <- LOF(U)  # log(LOF) by default
qplot(dist2, llof)

cowplot::plot_grid(
  qplot(U[, 1], U[, 2], color = llof, size = I(3)) + coord_equal() + 
    scale_color_viridis_c(breaks = NULL),
  qplot(U[, 3], U[, 7], color = llof, size = I(3)) + coord_equal() + 
    scale_color_viridis_c(breaks = NULL),
  rel_widths = c(0.7, 0.4), scale = 0.95
)


hist(dist, breaks = nclass.scottRob)
str(hist_out(dist))
abline(v = hist_out(dist)$lim[2], col = "red")

hist(dist2, breaks = nclass.scottRob, labels = T)
abline(v = hist_out(dist2)$lim[2], col = "red")

hist(llof, breaks = nclass.scottRob)
abline(v = hist_out(llof)$lim[2], col = "red")

eigval <- pca$sdev^2
hist(eigval, breaks = "FD")  # "FD" gives a bit more bins than scottRob
abline(v = hist_out(eigval, breaks = "FD")$lim[2], col = "red")

sum(eigval > hist_out(eigval, breaks = "FD")$lim[2])
pca_nspike(eigval)  # directly implemented in {bigutilsr}
