library("mixOmics")
library("zoo")
library("R.utils")
library(nlme)
library(emmeans)
library(ggpubr)
library(GGally)
library(gridExtra)
# library(formattable)


#setwd("/Users/juan/Library/CloudStorage/OneDrive-UniversitéLaval/INAF/VALENTIN/pca")
data.r <- read.csv("data//bleuet_data.csv")
#data.r$nomat <- sapply(data.r$nomat, as.character)
#data.r$sample_id <- paste(data.r$nomat,data.r$Visites)
outliers <- c(23,59)#,33)#,33)
data.r <- subset(data.r, !nomat %in% outliers)

# data.r$isMissingInsGluc <- is.na(data.r$Glucosefast) | is.na(data.r$Insulinefast)
# 
# data.r <- data.r %>%
#   mutate(FIRI = ifelse(isMissingInsGluc == "TRUE", NA, (Glucosefast*Insulinefast)/25))
# 
# data.r <- data.r %>%
#   mutate(QUICKI = ifelse(isMissingInsGluc == "TRUE", NA, (1/(log10(Insulinefast/6)+log10(Glucosefast*18.02)))))


data.r$sex[data.r$sex == ""] <- NA
data.r$sex <- na.locf(data.r$sex)
data.r$age[data.r$age == ""] <- NA
data.r$age <- na.locf(data.r$age)
bleuets.lme <- subset(data.r, !Visites == "V3" & !tx == "B")
# nonresponders<-nonresponders5
nonresponders<-group.noscale.no.outliers
# nonresponders <- group.noscale.outliers
# nonresponders <- nonresponders.log.bleuet6.1.VALE

bleuets.lme$cluster <- NA
bleuets.lme[bleuets.lme$nomat %in% nonresponders, "cluster"] <- "NR" #Groupe détecté (le + petit des 2)
bleuets.lme[!bleuets.lme$nomat %in% nonresponders, "cluster"] <- "R" #L'autre groupe
bleuets.lme$cluster
bleuets.lme[,c("nomat","cluster")]


data <-as.data.frame(bleuets.lme[,c("weight","bmi","waistc","hipc","whr","sbp","dbp","apob_neph","chol","tg",
                               "hdlc", "ldlc","chol_hdlc","creatinine","hba1c", "crp_neph","glucosefast","insulinefast","hr","matsuda5t","homa",
                               "glucosem","glucose0","glucose15","glucose30","glucose60","glucose90","glucose120","insulinem","insuline0",
                               "insuline15","insuline30","insuline60","insuline90","insuline120","glucosemoy5t","insulinemoy5t","plaq","vgm","na","k")]) # ,"FIRI","QUICKI"

#data <-as.data.frame(bleuets.lme[,c("gb","gr","hb","ht","plaq","vgm","tgmh","cgmh","dve","vpm","neutro_rel","neutro_abs","lymph_rel","lymph_abs","mono_rel","mono_abs","eosi_rel","eosi_abs","baso_rel","baso_abs","apob_neph","chol","tg","hdlc","ldlc","chol_hdlc","creatinine","na","k","chlore2018","ast","alt2018","calcium_ser","mg_serum","phosphore_ser","gluc_endocrin","insfast2018","insfast_hemo","hba1c","folate_serum","vitb12_serum2018","crp_neph","date_visite","weight","height","bmi","waistc","hipc","whr","sbp1","sbp2","sbp3","sbp","dbp1","dbp2","dbp3","dbp","hr","Glucosem","Glucose0","Glucose15","Glucose30","Glucose60","Glucose90","Glucose120","Insulinem","Insuline0","Insuline15","Insuline30","Insuline60","Insuline90","Insuline120","Glucosefast","Insulinefast","Glucosemoy5t","Insulinemoy5t","Matsuda5t","Homa"
#)]) # 

######linear mixed model ###############


list<-list()
list.gr<-list()
list.vi<-list()
list.grvi<-list()

list.emm.v2.NR<-list()
list.emm.v4.NR<-list()
list.se.v2.NR<-list()
list.se.v4.NR<-list()

list.emm.v2.R<-list()
list.emm.v4.R<-list()
list.se.v2.R<-list()
list.se.v4.R<-list()

for (i in 1:length(data)) { 
  variable <- data[,i]
  lme_cer <- lme(variable ~ sex+age+cluster*Visites, random = ~1 | nomat, method = "REML",data = bleuets.lme, na.action=na.omit)
  ano_cer <- anova.lme(lme_cer, type = "sequential", adjustSigma = F)
  gr <- ano_cer[4,4]
  vi <- ano_cer[5,4]
  grvi <- ano_cer[6,4]
  # print(summary(lme_cer)$tTable)
  list[[i]] <-ano_cer
  list.gr[[i]] <-gr
  list.vi[[i]] <-vi
  list.grvi[[i]] <-grvi
  
  emm  <- emmeans(lme_cer, ~ Visites|cluster)
  emm.df <- as.data.frame(emm[1:4])
  
  # emt <- emtrends(lmm1, ~ Visites|cluster, var = bleuets.lme$bmi)
  
  list.emm.v2.NR[[i]]  <- emm.df[1,3]
  list.emm.v4.NR[[i]]  <- emm.df[2,3]
  list.se.v2.NR[[i]]  <- emm.df[1,4]
  list.se.v4.NR[[i]]  <- emm.df[2,4]
  
  list.emm.v2.R[[i]]  <- emm.df[3,3]
  list.emm.v4.R[[i]]  <- emm.df[4,3]
  list.se.v2.R[[i]]  <- emm.df[3,4]
  list.se.v4.R[[i]]  <- emm.df[4,4]
  
}

df.gr<-data.frame(matrix(unlist(list.gr)))
df.vi<-data.frame(as.numeric(matrix(unlist(list.vi))))
df.grvi<-data.frame(matrix(unlist(list.grvi)))
df.emm.v2.NR<-data.frame(matrix(unlist(list.emm.v2.NR)))
df.se.v2.NR <- data.frame(matrix(unlist(list.se.v2.NR)))
df.emm.v4.NR<-data.frame(matrix(unlist(list.emm.v4.NR)))
df.se.v4.NR <- data.frame(matrix(unlist(list.se.v4.NR)))
df.emm.v2.R<-data.frame(matrix(unlist(list.emm.v2.R)))
df.se.v2.R <- data.frame(matrix(unlist(list.se.v2.R)))
df.emm.v4.R<-data.frame(matrix(unlist(list.emm.v4.R)))
df.se.v4.R <- data.frame(matrix(unlist(list.se.v4.R)))
df.names <- names(data)

results <- cbind(df.names,df.emm.v2.NR,df.se.v2.NR,df.emm.v4.NR,df.se.v4.NR,df.emm.v2.R,df.se.v2.R,df.emm.v4.R,df.se.v4.R,df.gr,df.vi,df.grvi)
names(results) <- c("variables","NRv2","se","NRv4","se","Rv2","se","Rv4","se","group","visite","interac") #c("variables","NRv2","NRv4","Rv2","Rv4","group","visite","interac")
results[2:12] <- round(results[2:12],4)
results$sig <- ifelse(results$interac < 0.05,"*", "")
results


# Average and tsnadard deviation
cols <- c("nomat","Visites","hdlc","ldlc","chol","chol_hdlc")
temp <- bleuets.lme

temp <-as.data.frame(bleuets.lme[,cols])
cols2 <- c("hdlc","ldlc","chol","chol_hdlc")
temp.v2 <- subset(temp, temp$Visites == "V2")
row.names(temp.v2) <- temp.v2$nomat
temp.v2 <- temp.v2[,cols2]
temp.v4 <- subset(temp, temp$Visites == "V4")
row.names(temp.v4) <- temp.v4$nomat
temp.v4 <- temp.v4[,cols2]

delta <- temp.v4 - temp.v2

group <- c(61,28,71,68,106)

delta.NR <- subset(delta, row.names(delta) %in% group)
delta.R <- subset(delta, !row.names(delta) %in% group)


averages <- colMeans(delta.NR)
standard_deviations <- apply(delta.NR, 2, sd)
results.NR <- data.frame(Average = averages, Standard_Deviation = standard_deviations)
results.NR

averages <- colMeans(delta.R)
standard_deviations <- apply(delta.R, 2, sd)
results.R <- data.frame(Average = averages, Standard_Deviation = standard_deviations)
results.R

#library(ggboxplot)
#library(ggpubr)
bxp <- ggboxplot(
  bleuets.lme, x = "Visites", y = "chol",
  color = "cluster", palette = "jco",
  facet.by = "cluster", short.panel.labs = FALSE,add = "jitter"
  #bxp.errorbar = T#,
  #ggtheme = theme_pubr()
)
bxp



#### BARPLOTS ----
# http://www.sthda.com/english/wiki/ggplot2-barplots-quick-start-guide-r-software-and-data-visualization
nonresponders<-nonresponders2
# nonresponders.c1
all.nonresponders <- nonresponders2

#all.nonresponders <- c(102,107,28,68,71,46,58,04,88)
lme.plots <- subset(data.r, !Visites == "V3" & !tx == "B")

lme.plots.pre <- subset(lme.plots, lme.plots$Visites == "V2")
lme.plots.post <- subset(lme.plots, lme.plots$Visites == "V4")

rownames(lme.plots.pre) <- lme.plots.pre$nomat
rownames(lme.plots.post) <- lme.plots.post$nomat

clinic.variables <- subset(results$variables, results$interac < 0.12, drop = T)

lme.plots.pre <- lme.plots.pre[,c(clinic.variables)]
lme.plots.post <- lme.plots.post[,c(clinic.variables)]

lme.plots.delta <- lme.plots.post - lme.plots.pre

lme.plots.delta$cluster <- NA
lme.plots.delta[row.names(lme.plots.delta) %in% nonresponders, "cluster"] <- "R" #répondeurs sur composante 2
lme.plots.delta[!row.names(lme.plots.delta) %in% nonresponders, "cluster"] <- "NR" #non répondeurs
# lme.plots.delta[row.names(lme.plots.delta) %in% nonresponders.c1, "cluster"] <- "R1" #répondeurs sur composante 1

lme.plots.delta$cluster
lme.plots.delta$nomat <- row.names(lme.plots.delta)

lme.plots.delta <- as.data.frame(lme.plots.delta)


hdlc.p <- ggplot(data=lme.plots.delta, aes(x = reorder(nomat, -hdlc), y = hdlc, fill=cluster))+
  geom_bar(stat="identity", color="black") +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_minimal()
hdlc.p

chol.p <- ggplot(data=lme.plots.delta, aes(x = reorder(nomat, -chol), y = chol, fill=cluster))+
  geom_bar(stat="identity", color="black") +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_minimal()
chol.p

chol_hdlc.p <- ggplot(data=lme.plots.delta, aes(x = reorder(nomat, -chol_hdlc), y = chol_hdlc, fill=cluster))+
  geom_bar(stat="identity", color="black") +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_minimal()
chol_hdlc.p

ldlc.p <- ggplot(data=lme.plots.delta, aes(x = reorder(nomat, -ldlc), y = ldlc, fill=cluster))+
  geom_bar(stat="identity", color="black") +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_minimal()
ldlc.p

tg.p <- ggplot(data=lme.plots.delta, aes(x = reorder(nomat, -tg), y = tg, fill=cluster))+
  geom_bar(stat="identity", color="black") +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_minimal()
tg.p

### CONTRASTS
# Insuline
lmm1 <- lme(Insulinefast ~ sex+age+cluster*Visites, random = ~1 | nomat, method = "REML",data = bleuets.lme, na.action=na.omit)
emm1 <- emmeans(lmm1, ~ Visites|cluster)
contrast(emm1, interaction = "pairwise", by = NULL)

# Matsuda
lmm2 <- lme(Matsuda5t ~ sex+age+cluster*Visites, random = ~1 | nomat, method = "REML",data = bleuets.lme, na.action=na.omit)
emm2 <- emmeans(lmm2, ~ Visites|cluster)
contrast(emm2, interaction = "pairwise", by = NULL)

# whr
lmm3 <- lme(whr ~ sex+age+cluster*Visites, random = ~1 | nomat, method = "REML",data = bleuets.lme, na.action=na.omit)
emm3 <- emmeans(lmm3, ~ Visites|cluster)
contrast(emm3, interaction = "pairwise", by = NULL)

# apob
lmm4 <- lme(apob_neph ~ sex+age+cluster*Visites, random = ~1 | nomat, method = "REML",data = bleuets.lme, na.action=na.omit)
emm4 <- emmeans(lmm4, ~ Visites|cluster)
contrast(emm4, interaction = "pairwise", by = NULL)

# crp
lmm5 <- lme(crp_neph ~ sex+age+cluster*Visites, random = ~1 | nomat, method = "REML",data = bleuets.lme, na.action=na.omit)
emm5 <- emmeans(lmm5, ~ Visites|cluster)
contrast(emm5, interaction = "pairwise", by = NULL)

# waistc
lmm5 <- lme(waistc ~ sex+age+cluster*Visites, random = ~1 | nomat, method = "REML",data = bleuets.lme, na.action=na.omit)
emm5 <- emmeans(lmm5, ~ Visites|cluster)
contrast(emm5, interaction = "pairwise", by = NULL)

# homa
lmm5 <- lme(Homa ~ sex+age+cluster*Visites, random = ~1 | nomat, method = "REML",data = bleuets.lme, na.action=na.omit)
emm5 <- emmeans(lmm5, ~ Visites|cluster)
contrast(emm5, interaction = "pairwise", by = NULL)

ggpairs(bleuets.lme[c("waistc","whr", "apob_neph","crp_neph","Insulinefast","Matsuda5t","Homa","cluster")], aes(color = cluster), lower = list(continuous = "smooth"), columns = 1:8)


contrast(emm, interaction = "pairwise", by = NULL)


### STRIPCHART
# http://www.sthda.com/english/wiki/ggplot2-stripchart-jitter-quick-start-guide-r-software-and-data-visualization

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

# p <- ggplot(lme.plots.delta, aes(x=cluster, y=Insulinefast, color=cluster)) + 
#   geom_jitter()
# p + stat_summary(fun.data=data_summary, color="blue") + theme(legend.position="right")

p.waistc <- ggplot(bleuets.lme, aes(x=Visites, y=waistc, color=cluster, shape=cluster)) + 
  geom_jitter()
p.waistc + stat_summary(fun.data=data_summary, color="blue") + theme(legend.position="right") + scale_color_brewer(palette="Dark2") + theme_classic()

p.whr <- ggplot(bleuets.lme, aes(x=Visites, y=whr, color=cluster, shape=cluster)) + 
  geom_jitter()
p.whr + stat_summary(fun.data=data_summary, color="blue") + theme(legend.position="right") + scale_color_brewer(palette="Dark2") + theme_classic()

p.apob_neph <- ggplot(bleuets.lme, aes(x=Visites, y=apob_neph, color=cluster, shape=cluster)) + 
  geom_jitter()
p.apob_neph + stat_summary(fun.data=data_summary, color="blue") + theme(legend.position="right") + scale_color_brewer(palette="Dark2") + theme_classic()

p.crp_neph <- ggplot(bleuets.lme, aes(x=Visites, y=crp_neph, color=cluster, shape=cluster)) + 
  geom_jitter()
p.crp_neph + stat_summary(fun.data=data_summary, color="blue") + theme(legend.position="right") + scale_color_brewer(palette="Dark2") + theme_classic()

p.Insulinefast <- ggplot(bleuets.lme, aes(x=Visites, y=Insulinefast, color=cluster, shape=cluster)) + 
  geom_jitter(position=)
p.Insulinefast + stat_summary(fun.data=data_summary, color="blue") + theme(legend.position="right") + scale_color_brewer(palette="Dark2") + theme_classic()

p.Matsuda5t <- ggplot(bleuets.lme, aes(x=Visites, y=Matsuda5t, color=cluster, shape=cluster)) + 
  geom_jitter()
p.Matsuda5t + stat_summary(fun.data=data_summary, color="blue") + theme(legend.position="right") + scale_color_brewer(palette="Dark2") + theme_classic()

p.Homa <- ggplot(bleuets.lme, aes(x=Visites, y=Homa, color=cluster, shape=cluster)) + 
  geom_jitter()
p.Homa + stat_summary(fun.data=data_summary, color="blue") + theme(legend.position="right") + scale_color_brewer(palette="Dark2") + theme_classic()

p.Quicki <- ggplot(bleuets.lme, aes(x=Visites, y=QUICKI, color=cluster, shape=cluster)) + 
  geom_jitter()
p.Quicki + stat_summary(fun.data=data_summary, color="blue") + theme(legend.position="right") + scale_color_brewer(palette="Dark2") + theme_classic()

p.ins <- p.Insulinefast + stat_summary(fun.data=data_summary, color="blue") + theme(legend.position="right") + scale_color_brewer(palette="Dark2") + theme_classic()
p.mat <- p.Matsuda5t + stat_summary(fun.data=data_summary, color="blue") + theme(legend.position="right") + scale_color_brewer(palette="Dark2") + theme_classic()
p.hom <- p.Homa + stat_summary(fun.data=data_summary, color="blue") + theme(legend.position="right") + scale_color_brewer(palette="Dark2") + theme_classic()
p.wai <- p.waistc + stat_summary(fun.data=data_summary, color="blue") + theme(legend.position="right") + scale_color_brewer(palette="Dark2") + theme_classic()
p.whr2 <- p.whr + stat_summary(fun.data=data_summary, color="blue") + theme(legend.position="right") + scale_color_brewer(palette="Dark2") + theme_classic()
p.crp <- p.crp_neph + stat_summary(fun.data=data_summary, color="blue") + theme(legend.position="right") + scale_color_brewer(palette="Dark2") + theme_classic()
p.apo <- p.apob_neph + stat_summary(fun.data=data_summary, color="blue") + theme(legend.position="right") + scale_color_brewer(palette="Dark2") + theme_classic()
p.qui <- p.Quicki + stat_summary(fun.data=data_summary, color="blue") + theme(legend.position="right") + scale_color_brewer(palette="Dark2") + theme_classic()


# grid.arrange(p.ins,p.mat,p.hom,p.wai,p.whr2,p.crp,p.apo,p.qui, nrow = 3)

ggarrange(p.ins,p.qui,p.mat,p.hom,p.wai,p.whr2,p.crp,p.apo, nrow=3, ncol=3, levels=c("A","B","C","D","E","F","G"), common.legend = T, legend="bottom")

# dev.off()

# test jitter v2
p.Insulinefast2 <- ggplot(bleuets.lme, aes(x=Visites, y=Insulinefast, color=cluster, shape=cluster)) + 
  geom_jitter()
p.Insulinefast2 + stat_summary(fun.data=data_summary, color="blue") + theme(legend.position="right") + scale_color_brewer(palette="Dark2") + theme_classic()


p.ins3 <- ggplot(bleuets.lme, aes(x = Visites, y = Insulinefast, color = cluster)) +
  geom_point(position = position_jitterdodge(0.1) , size=1.5) +
  stat_summary(fun.data=data_summary, color="blue") +
  theme(legend.position = "right") +
  labs(title = "Fasting insulin changes per subgroup per visist", x = "Visits", y = "Fasting insulin", color = "Subgroup") +
  theme_fivethirtyeight() +
  theme(axis.title = element_text())
p.ins3

p.qui3 <- ggplot(bleuets.lme, aes(x = Visites, y = QUICKI, color = cluster, shape = cluster)) +
  geom_point(position = position_jitterdodge()) +
  stat_summary(fun.data=data_summary, color="blue") +
  theme(legend.position = "right") +
  labs(x = "x axis", y = "QUICKI", color = "subgroup") +
  scale_color_brewer(palette = "Dark2") +
  theme_classic()

p.wai3 <- ggplot(bleuets.lme, aes(x = Visites, y = waistc, color = cluster, shape = cluster)) +
  geom_point(position = position_jitterdodge()) +
  stat_summary(fun.data=data_summary, color="blue") +
  theme(legend.position = "right") +
  scale_color_brewer(palette = "Dark2") +
  theme_classic()

p.whr3 <- ggplot(bleuets.lme, aes(x = Visites, y = whr, color = cluster, shape = cluster)) +
  geom_point(position = position_jitterdodge()) +
  stat_summary(fun.data=data_summary, color="blue") +
  theme(legend.position = "right") +
  scale_color_brewer(palette = "Dark2") +
  theme_classic()

p.apob3 <- ggplot(bleuets.lme, aes(x = Visites, y = apob_neph, color = cluster, shape = cluster)) +
  geom_point(position = position_jitterdodge()) +
  stat_summary(fun.data=data_summary, color="blue") +
  theme(legend.position = "right") +
  scale_color_brewer(palette = "Dark2") +
  theme_classic()

p.crp3 <- ggplot(bleuets.lme, aes(x = Visites, y = crp_neph, color = cluster, shape = cluster)) +
  geom_point(position = position_jitterdodge()) +
  stat_summary(fun.data=data_summary, color="blue") +
  theme(legend.position = "right") +
  scale_color_brewer(palette = "Dark2") +
  theme_classic()

# ggarrange(p.ins3,p.qui3,p.wai3,p.whr3,p.crp3,p.apob3, nrow=2, ncol=3, levels=c("A","B","C","D","E"), common.legend = T, legend="bottom", font.label = list(size = 14, color = "black", face = "bold", family = NULL))
# ggarrange(p.ins3,p.qui3,p.wai3,p.whr3,p.crp3,p.apob3, nrow=3, ncol=2, levels=c("A","B","C","D","E"), common.legend = T, legend="bottom")
# 
# grid.arrange(p.ins3,p.qui3,p.wai3,p.whr3,p.crp3,p.apob3, nrow = 3, common.legend = T)
ggarrange(p.ins3,p.qui3,p.wai3,p.whr3,p.crp3,p.apob3, nrow=3, ncol=2, common.legend=T,legend="bottom", font.label = list(size = 14, color = "black", face = "bold", family = NULL))

# backup
# p.ins3 <- ggplot(bleuets.lme, aes(x = Visites, y = Insulinefast, color = cluster, shape = cluster)) +
#   geom_point(position = position_jitterdodge()) +
#   stat_summary(fun.data = mean_cl_boot, color = "blue") +
#   theme(legend.position = "right") +
#   scale_color_brewer(palette = "Dark2") +
#   theme_classic()



# SPLIT VIOLIN

GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}


p.violin.quicki <- ggplot(bleuets.lme, aes(x=Visites, y=QUICKI, fill = cluster, shape=cluster)) + geom_split_violin()
p.violin.insuline <- ggplot(bleuets.lme, aes(x=Visites, y=Insulinefast, fill = cluster, shape=cluster)) + geom_split_violin()
p.violin.matsuda <- ggplot(bleuets.lme, aes(x=Visites, y=Matsuda5t, fill = cluster, shape=cluster)) + geom_split_violin()
p.violin.Homa <- ggplot(bleuets.lme, aes(x=Visites, y=Homa, fill = cluster, shape=cluster)) + geom_split_violin()
p.violin.waistc <- ggplot(bleuets.lme, aes(x=Visites, y=waistc, fill = cluster, shape=cluster)) + geom_split_violin()
p.violin.whr <- ggplot(bleuets.lme, aes(x=Visites, y=whr, fill = cluster, shape=cluster)) + geom_split_violin()
p.violin.crp <- ggplot(bleuets.lme, aes(x=Visites, y=crp_neph, fill = cluster, shape=cluster)) + geom_split_violin()
p.violin.apob <- ggplot(bleuets.lme, aes(x=Visites, y=apob_neph, fill = cluster, shape=cluster)) + geom_split_violin()


p.ins2 <- p.violin.insuline + stat_summary(fun.data=data_summary, color="blue") + theme(legend.position="right") + scale_color_brewer(palette="Dark2") + theme_classic()
p.mat2 <- p.violin.matsuda + stat_summary(fun.data=data_summary, color="blue") + theme(legend.position="right") + scale_color_brewer(palette="Dark2") + theme_classic()
p.hom2 <- p.violin.Homa + stat_summary(fun.data=data_summary, color="blue") + theme(legend.position="right") + scale_color_brewer(palette="Dark2") + theme_classic()
p.wai2 <- p.violin.waistc + stat_summary(fun.data=data_summary, color="blue") + theme(legend.position="right") + scale_color_brewer(palette="Dark2") + theme_classic()
p.whr22 <- p.violin.whr + stat_summary(fun.data=data_summary, color="blue") + theme(legend.position="right") + scale_color_brewer(palette="Dark2") + theme_classic()
p.crp2 <- p.violin.crp + stat_summary(fun.data=data_summary, color="blue") + theme(legend.position="right") + scale_color_brewer(palette="Dark2") + theme_classic()
p.apo2 <- p.violin.apob + stat_summary(fun.data=data_summary, color="blue") + theme(legend.position="right") + scale_color_brewer(palette="Dark2") + theme_classic()
p.qui2 <- p.violin.quicki + stat_summary(fun.data=data_summary, color="blue") + theme(legend.position="right") + scale_color_brewer(palette="Dark2") + theme_classic()

ggarrange(p.ins2,p.qui2,p.mat2,p.hom2,p.wai2,p.whr22,p.crp2,p.apo2, nrow=3, ncol=3, levels=c("A","B","C","D","E","F","G"), common.legend = T, legend="bottom", font.label = list(size = 14, color = "black", face = "bold", family = NULL))


# Correlation LME
#Explore associations using anthropometrical variables
tab1vec_c <- c("waistc","whr","apob_neph","crp_neph","Insulinefast","Matsuda5t","Homa","QUICKI")
tab1vec_m <- c("chiro-inositol","xylose","3-ureidopropionate","N6,N6-dimethyllysine","creatine","adenosine","3-carboxy-4-methyl-5-pentyl-2-furanpropionate (3-CMPFP)**",
               "trigonelline (N'-methylnicotinate)","acisoga","octadecanedioylcarnitine (C18-DC)*",
               "arachidonoylcarnitine (C20:4)","carotene diol (2)","carotene diol (1)","3beta-hydroxy-5-cholestenoate",
               "octanoylcarnitine (C8)","tetradecadienoate (14:2)*","linoleate (18:2n6)","dodecadienoate (12:2)*",
               "nonanoylcarnitine (C9)","cis-4-decenoylcarnitine (C10:1)")

#Telomere length
#Generate list for store results
# bleuets_lme_t <- as.data.frame(t(bleuets.lme))
# colnames(bleuets_lme_t) <- bleuets_lme_t$subject_id
# 
# bleuets_lme2 <- cbind(bleuets.lme,log.bleuetsii)


# t1a_models <- list()
# t1a_models <- lapply(tab1vec_c,list)
# 
# t1a_summaries <- list()
# t1a_summaries <- lapply(tab1vec_c,list)
# 
# ##Model for differences
# for (j in seq_along(tab1vec_c)) {
#   for (i in seq_along(tab1vec_m)) {
#     t4a_models[[j]][[i]] <- lme(formula = paste0(tab1vec_c[[j]],
#                                                      " ~ sex + age + ",
#                                                      tab1vec_m[[i]],
#                                                      "*cluster*Visites"), random = ~1 | nomat, method = "REML",
#                                                      data = contrats_data, na.action = "na.exclude")
#                                     t1a_summaries[[j]][[i]] <- coef(summary(t1a_models[[j]][[i]]))[9,c(1,2,5)]
#   }
# }
# 
# 
# 
# ###Add labels and set colnames for each dependent variable
# names(t1a_summaries) <- tab1vec_c
# ###Adding names on nested lists
# for(i in seq_along(tab1vec_c)) {                                              # Run for-nested lists
#   names(t1a_summaries[[i]]) <- tab1vec_m
# }
# 
# ###Transform in dataframes
# for(i in seq_along(tab1vec_c)) {                                                # Run for-nested lists
#   t1a_summaries[[i]] <- round(do.call(rbind, t1a_summaries[[i]]), digits = 4)
# }
# 
# t1a_summaries 


# V2 #######################################################################

#Merge clinical variables and metabolites in the same DB
library(lme4)
library(nlme)
library(lmerTest)
library(emmeans)
library(margins)
library(effects)
library(janitor)

#Explore associations using clinical vars and metabolites variables
#NONRESPONDERS5
tab1vec_c <- c("tg","hdlc","chol_hdlc")
tab1vec_m <- c("pantoate","oleoyl-linoleoyl-glycerol (18:1/18:2) [2]",
               "1-palmitoyl-2-linoleoyl-GPI (16:0/18:2)",
               "1-(1-enyl-palmitoyl)-2-linoleoyl-GPE (P-16:0/18:2)*","linoleoyl-linoleoyl-glycerol (18:2/18:2) [1]*",
               "1-(1-enyl-stearoyl)-2-oleoyl-GPE (P-18:0/18:1)","andro steroid monosulfate C19H28O6S (1)*",
               "1-palmitoyl-2-palmitoleoyl-GPC (16:0/16:1)*",
               "N-succinyl-phenylalanine","oleoyl-linoleoyl-glycerol (18:1/18:2) [1]","1-linolenoylglycerol (18:3)",
               "S-adenosylhomocysteine (SAH)","inosine","gentisate","adenosine","linoleoylcarnitine (C18:2)*",
               "N-acetyl-3-methylhistidine*","eicosenedioate (C20:1-DC)*","indolepropionate","linoleoyl-linoleoyl-glycerol (18:2/18:2) [2]*",
               "ethylmalonate","1-palmitoyl-2-alpha-linolenoyl-GPC (16:0/18:3n3)*")          ##With the same name than your DB
tab1vec_m <- make_clean_names(tab1vec_m)

log.bleuetsii$subject_id <- row.names(log.bleuetsii)
contrats_data <- merge(bleuets.lme,log.bleuetsii, by="subject_id")

#Generate list for store results
#Merge clinical variables and metabolites in the same DB

colnames_data <- colnames(contrats_data)
colnames_data <- make_clean_names(colnames_data)
colnames(contrats_data) <- colnames_data

#Generate list for store results
t1a_models <- list()
t1a_models <- lapply(tab1vec_c,  list)

t1a_summaries <- list()
t1a_summaries <- lapply(tab1vec_c,  list)

##Model for differences
for (j in seq_along(tab1vec_c)) {
  for (i in seq_along(tab1vec_m)) {
    t1a_models[[j]][[i]]    <- lmer(formula = paste0(tab1vec_c[[j]], " ~ sex + age + ", 
                                                     tab1vec_m[[i]], "*cluster*visites + (1|nomat)"), 
                                    REML = TRUE, data = contrats_data, na.action = "na.exclude")
    t1a_summaries[[j]][[i]] <- coef(summary(t1a_models[[j]][[i]]))[10,c(1,2,5)]
    }
}

###Add labels and set colnames for each dependent variable
names(t1a_summaries) <- tab1vec_c
###Adding names on nested lists
for(i in seq_along(tab1vec_c)) {                                              # Run for-nested lists
  names(t1a_summaries[[i]]) <- tab1vec_m
}

###Transform in dataframes
for(i in seq_along(tab1vec_c)) {                                                # Run for-nested lists
  t1a_summaries[[i]] <- round(do.call(rbind, t1a_summaries[[i]]), digits = 4)
}

t1a_summaries   




# ## chol hdlc
# chol hdlc and homovanillate_hva - p = 0.0135
emtrends(t1a_models[[3]][[20]], consec ~ cluster*visites, var = "homovanillate_hva")

# chol hdlc and n_succinyl_phenylalanine - p = 0.0547
emtrends(t1a_models[[3]][[8]], consec ~ cluster*visites, var = "n_succinyl_phenylalanine")

# chol hdlc and inosine - p = 0.0099
emtrends(t1a_models[[3]][[12]], consec ~ cluster*visites, var = "inosine")


# ## HDLC
# HDLC and gentisate - p = 0.044
emtrends(t1a_models[[2]][[13]], consec ~ cluster*visites, var = "gentisate")

# HDLC and indolepropionate - p = 0.0233
emtrends(t1a_models[[2]][[17]], consec ~ cluster*visites, var = "indolepropionate")

# HDLC and x1_palmitoyl_2_alpha_linolenoyl_gpc_16_0_18_3n3 - p = 0.0233
emtrends(t1a_models[[2]][[22]], consec ~ cluster*visites, var = "x1_palmitoyl_2_alpha_linolenoyl_gpc_16_0_18_3n3")


# ## TG
# TG and gentisate - p = 0.0816
emtrends(t1a_models[[1]][[13]], consec ~ cluster*visites, var = "gentisate")

# TG and ethylmalonate - p = 0.0046
emtrends(t1a_models[[1]][[21]], consec ~ cluster*visites, var = "ethylmalonate")



# emtrends(t1a_models[[8]][[12]], c("consec", "consec") ~ cluster*visites, var = "carotene_diol_2")


# Boxplots figures for LME

bleuets.lme2 <- bleuets.lme
bleuets.lme2$visit <- ifelse(bleuets.lme2$Visites == "V2", "Pre","Post")
bleuets.lme2$cluster2 <- ifelse(bleuets.lme$cluster == "respc1", "R1", "R2")
bleuets.lme2$visit <- relevel(as.factor(bleuets.lme2$visit),"Pre","Post")

outlier <- c(72, 126) # CHECKER 3ème outlier et c("LARC-00307","LARC-00308","LARC-00347","LARC-00348")#,"LARC-00279", "LARC-00280")
bleuets.lme2 <- subset(bleuets.lme2, !bleuets.lme2$nomat %in% outlier)


# bleuets.lme2 %>% mutate(visit = case_when (Visites == "V2" ~ "Pre", Visites == "V4" ~ "Post"))
# bleuets.lme2 %>% mutate(cluster2 = case_when (cluster == "respc1" ~ "R1", cluster == "respc2" ~ "R2"))


p.ins4 <-  ggplot(bleuets.lme2, aes(x = reorder(visit, insulinefast), y = insulinefast, fill = cluster2)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_point(aes(fill = cluster2), size = 2, shape = 21, position = position_jitterdodge()) +
  facet_wrap(~cluster2) +
  theme(legend.position = "right") +
  labs(title = "A. Fasting insulin", x = "Visits", y = " Fasting insulin (pmol/L)", fill = "Group:") +
  theme_fivethirtyeight() +
  theme(axis.title = element_text(), plot.title = element_text(size=12, hjust = 0.5),
         plot.background = element_blank(),panel.background = element_blank(), legend.background = element_blank(), 
        legend.key = element_blank(), strip.background = element_blank(), strip.text.x = element_blank()) +
scale_fill_manual(values = c('#E41A1C', '#377EB8'))
p.ins4

p.qui4 <- ggplot(bleuets.lme2, aes(x = reorder(visit, QUICKI), y = QUICKI, fill = cluster2)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_point(aes(fill = cluster2), size = 2, shape = 21, position = position_jitterdodge()) +
  facet_wrap(~cluster2) +
  theme(legend.position = "right") +
  labs(title = "B. QUICKI", x = "Visits", y = "QUICKI", color = "Group:") +
  theme_fivethirtyeight() +
  theme(axis.title = element_text(), plot.title = element_text(size=12, hjust = 0.5),
        plot.background = element_blank(),panel.background = element_blank(), legend.background = element_blank(), 
        legend.key = element_blank(), strip.background = element_blank(), strip.text.x = element_blank()) +
  scale_fill_manual(values = c('#E41A1C', '#377EB8'))
p.qui4

p.wai4 <- ggplot(bleuets.lme2, aes(x = reorder(visit, waistc), y = waistc, fill = cluster2)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_point(aes(fill = cluster2), size = 2, shape = 21, position = position_jitterdodge()) +
  facet_wrap(~cluster2) +
  theme(legend.position = "right") +
  labs(title = "E. Waist circumference", x = "Visits", y = "Waist circumference (cm)", color = "Group:") +
  theme_fivethirtyeight() +
  theme(axis.title = element_text(), plot.title = element_text(size=12, hjust = 0.5),
        plot.background = element_blank(),panel.background = element_blank(), legend.background = element_blank(), 
        legend.key = element_blank(), strip.background = element_blank(), strip.text.x = element_blank()) +
  scale_fill_manual(values = c('#E41A1C', '#377EB8'))
p.wai4

p.whr4 <- ggplot(bleuets.lme2, aes(x = reorder(visit, whr), y = whr, fill = cluster2)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_point(aes(fill = cluster2), size = 2, shape = 21, position = position_jitterdodge()) +
  facet_wrap(~cluster2) +
  theme(legend.position = "right") +
  labs(title = "F. Waist-to-hip ratio", x = "Visits", y = "Waist-to-hip ratio", color = "Group:") +
  theme_fivethirtyeight() +
  theme(axis.title = element_text(), plot.title = element_text(size=12, hjust = 0.5),
        plot.background = element_blank(),panel.background = element_blank(), legend.background = element_blank(), 
        legend.key = element_blank(), strip.background = element_blank(), strip.text.x = element_blank()) +
  scale_fill_manual(values = c('#E41A1C', '#377EB8'))
p.whr4

p.apo4 <- ggplot(bleuets.lme2, aes(x = reorder(visit, apob_neph), y = apob_neph, fill = cluster2)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_point(aes(fill = cluster2), size = 2, shape = 21, position = position_jitterdodge()) +
  facet_wrap(~cluster2) +
  theme(legend.position = "right") +
  labs(title = "H. ApoB", x = "Visits", y = "Apo B (g/L)", color = "Group:") +
  theme_fivethirtyeight() +
  theme(axis.title = element_text(), plot.title = element_text(size=12, hjust = 0.5),
        plot.background = element_blank(),panel.background = element_blank(), legend.background = element_blank(), 
        legend.key = element_blank(), strip.background = element_blank(), strip.text.x = element_blank()) +
  scale_fill_manual(values = c('#E41A1C', '#377EB8'))
p.apo4

p.crp4 <- ggplot(bleuets.lme2, aes(x = reorder(visit, crp_neph), y = crp_neph, fill = cluster2)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_point(aes(fill = cluster2), size = 2, shape = 21, position = position_jitterdodge()) +
  facet_wrap(~cluster2) +
  theme(legend.position = "right") +
  labs(title = "G. CRP", x = "Visits", y = "CRP (mg/L)", color = "Group:") +
  theme_fivethirtyeight() +
  theme(axis.title = element_text(), plot.title = element_text(size=12, hjust = 0.5),
        plot.background = element_blank(),panel.background = element_blank(), legend.background = element_blank(), 
        legend.key = element_blank(), strip.background = element_blank(), strip.text.x = element_blank()) +
  scale_fill_manual(values = c('#E41A1C', '#377EB8'))
p.crp4

p.mat4 <- ggplot(bleuets.lme2, aes(x = reorder(visit, matsuda5t), y = matsuda5t, fill = cluster2)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_point(aes(fill = cluster2), size = 2, shape = 21, position = position_jitterdodge()) +
  facet_wrap(~cluster2) +
  theme(legend.position = "right") +
  labs(title = "C. Matsuda", x = "Visits", y = "Matsuda", color = "Group:") +
  theme_fivethirtyeight() +
  theme(axis.title = element_text(), plot.title = element_text(size=12, hjust = 0.5),
        plot.background = element_blank(),panel.background = element_blank(), legend.background = element_blank(), 
        legend.key = element_blank(), strip.background = element_blank(), strip.text.x = element_blank()) +
  scale_fill_manual(values = c('#E41A1C', '#377EB8'))
p.mat4

p.hom4 <- ggplot(bleuets.lme2, aes(x = reorder(visit, homa), y = homa, fill = cluster2)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_point(aes(fill = cluster2), size = 2, shape = 21, position = position_jitterdodge()) +
  facet_wrap(~cluster2) +
  theme(legend.position = "right") +
  labs(title = "D. HOMA-IR", x = "Visits", y = "HOMA-IR", color = "Group:") +
  theme_fivethirtyeight() +
  theme(axis.title = element_text(), plot.title = element_text(size=12, hjust = 0.5),
        plot.background = element_blank(),panel.background = element_blank(), legend.background = element_blank(), 
        legend.key = element_blank(), strip.background = element_blank(), strip.text.x = element_blank()) +
  scale_fill_manual(values = c('#E41A1C', '#377EB8')) #c("#C63F23", "#278F3D")
p.hom4

lme.figure <- ggarrange(p.ins4,p.qui4,p.mat4,p.hom4,p.wai4,p.whr4,p.crp4,p.apo4, nrow=2, ncol=4, common.legend=T,legend=NULL, 
          font.label = list(size = 14, color = "black", face = "bold", family = NULL))
# lme.figure

annotate_figure(lme.figure, top = text_grob("", color = "Black", face = "bold", size = 14))


check_apob <- subset(bleuets.lme, bleuets.lme$cluster == "respc2")
check_apob <- subset(check_apob, check_apob$Visites == "V2")
check_apob$apob_neph



# Boxplot with lines
# p.qui4 <- ggplot(bleuets.lme, aes(x = Visites, y = QUICKI, fill = cluster)) +
#   geom_boxplot() +
#   geom_point(aes(fill = cluster), size = 2, shape = 21) +
#   geom_line(aes(group=nomat)) +
#   facet_wrap(~cluster) +
#   theme(legend.position = "right") +
#   labs(title = "QUICKI changes per subgroup per visist", x = "Visits", y = "QUICKI", color = "Subgroup") +
#   theme_fivethirtyeight() +
#   theme(axis.title = element_text())
# p.qui4


# PLOT INTERACTIONS https://strengejacke.github.io/sjPlot/articles/plot_interactions.html
library(sjPlot)
library(sjmisc)
library(ggplot2)
library(ggalluvial)
data(efc)
theme_set(theme_sjplot())

# QUICKI and Carotene Diol 1
# fit <- lm(quicki ~  cluster * carotene_diol_1, data = contrats_data)
caro.in <- plot_model(t1a_models[[8]][[13]], type = "emm", terms = c("carotene_diol_1","cluster"),title="B.")

# Matsuda and x3_ureidopropionate
# fit2 <- lm(matsuda5t ~cluster * x3_ureidopropionate, data = contrats_data)
ure.mat <- plot_model(t1a_models[[6]][[3]], type = "emm", terms = c("x3_ureidopropionate","cluster"),title="D.")

# Insulinefast and carotene diol 1
# fit3 <- lm(insulinefast ~ visites * cluster * carotene_diol_1, data = contrats_data)
# fit3 <- lm(insulinefast ~ cluster * carotene_diol_1, data = contrats_data)
#plot_model(t1a_models[[5]][[13]], type = "pred", terms = c("carotene_diol_1","cluster","visites"))
caro.qui <- plot_model(t1a_models[[5]][[13]], type = "emm", terms = c("carotene_diol_1","cluster"),title="F.")

?plot_model
em.insu.caro

# data.aluv <- subset() VERIFIER SI ORDRE DES NOMAT OK
delta.clinic <- delta.bleuets
delta.clinic$nomat <- rownames(delta.clinic)

delta.met <- delta.bleuets2
delta.met$id <- rownames(delta.met)
delta.met$nomat <- delta.clinic$nomat
delta.met.short <- as.data.frame(c(delta.met$`carotene diol (1)`,delta.met$`3-ureidopropionate`,delta.met$nomat))
#extract 

id.match <- as.data.frame(meta.bleuets$SUBJECT_ID)
rownames(id.match) <- rownames(meta.bleuets)

delta.lme <- merge(delta.clinic,delta.met, by="nomat")

type2 <- c(39,69,90,96,98)

delta.lme$cluster <- NA
delta.lme[delta.lme$nomat %in% type2, "cluster"] <- "R2"
delta.lme[!delta.lme$nomat %in% type2, "cluster"] <- "R1"

delta.lme$QUICKI.delta <- ifelse(delta.lme$QUICK < 0, "neg","pos")

# PLOT `carotene diol (1)`
# ggplot(data = delta.lme,
#        aes(axis1 = cluster, axis2 = QUICKI.delta, y = `carotene diol (1)`)) +
#   geom_alluvium(aes(fill = cluster)) +
#   geom_stratum() +
#   geom_text(stat = "stratum",
#             aes(label = after_stat(stratum))) +
#   scale_x_discrete(limits = c("Survey", "Response"),
#                    expand = c(0.25, 0.05)) +
#   theme_void()


dim(bleuets.lme2)
bleuets.lme3 <- cbind(bleuets.lme2,ure=log.bleuets6$`3-ureidopropionate`,chiro=log.bleuets6$`chiro-inositol`, xylo=log.bleuets6$xylose, caro=log.bleuets6$`carotene diol (1)`)
p.ure <- ggplot(bleuets.lme3, aes(x = reorder(visit, ure), y = ure, fill = cluster2)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_point(aes(fill = cluster2), size = 2, shape = 21, position = position_jitterdodge()) +
  facet_wrap(~cluster2) +
  theme(legend.position = "right") +
  labs(title = "E.", x = "Visits", y = "3-ureidopropionate", color = "Subgroup") +
  theme_fivethirtyeight() +
  theme(axis.title = element_text(), plot.title = element_text(size=12, hjust = 0.5),
        plot.background = element_blank(),panel.background = element_blank(), legend.background = element_blank(), 
        legend.key = element_blank(), strip.background = element_blank(), strip.text.x = element_blank()) +
  scale_fill_manual(values = c('#E41A1C', '#377EB8')) #c("#C63F23", "#278F3D")
p.ure

p.caro1 <- ggplot(bleuets.lme3, aes(x = reorder(visit, caro), y = caro, fill = cluster2)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_point(aes(fill = cluster2), size = 2, shape = 21, position = position_jitterdodge()) +
  facet_wrap(~cluster2) +
  theme(legend.position = "right") +
  labs(title = "A.", x = "Visits", y = "Carotene diol 1", color = "Subgroup") +
  theme_fivethirtyeight() +
  theme(axis.title = element_text(), plot.title = element_text(size=12, hjust = 0.5),
        plot.background = element_blank(),panel.background = element_blank(), legend.background = element_blank(), 
        legend.key = element_blank(), strip.background = element_blank(), strip.text.x = element_blank()) +
  scale_fill_manual(values = c('#E41A1C', '#377EB8')) #c("#C63F23", "#278F3D")
p.caro1

p.caro12 <- ggplot(bleuets.lme3, aes(x = reorder(visit, caro), y = caro, fill = cluster2)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_point(aes(fill = cluster2), size = 2, shape = 21, position = position_jitterdodge()) +
  facet_wrap(~cluster2) +
  theme(legend.position = "right") +
  labs(title = "C.", x = "Visits", y = "Carotene diol 1", color = "Subgroup") +
  theme_fivethirtyeight() +
  theme(axis.title = element_text(), plot.title = element_text(size=12, hjust = 0.5),
        plot.background = element_blank(),panel.background = element_blank(), legend.background = element_blank(), 
        legend.key = element_blank(), strip.background = element_blank(), strip.text.x = element_blank()) +
  scale_fill_manual(values = c('#E41A1C', '#377EB8')) #c("#C63F23", "#278F3D")
p.caro12

lme.figure2 <- ggarrange(p.caro1,caro.in,p.caro12,caro.qui,p.ure,ure.mat, nrow=3, ncol=2, common.legend=T,legend=NULL, 
                        font.label = list(size = 14, color = "black", face = "bold", family = NULL))
# lme.figure

lme.figure3 <- ggarrange(p.caro1,caro.in,p.ure,ure.mat, nrow=2, ncol=2, common.legend=T,legend=NULL, 
                         font.label = list(size = 14, color = "black", face = "bold", family = NULL))

annotate_figure(lme.figure2, top = text_grob("", color = "Black", face = "bold", size = 14))
annotate_figure(lme.figure3, top = text_grob("", color = "Black", face = "bold", size = 14))

