#### LME2 for correlations ####
#Merge clinical variables and metabolites in the same DB
library(lme4)
library(nlme)
library(lmerTest)
library(emmeans)
library(margins)
library(effects)
library(janitor)
# On récupère les noms des agv et les variables cliniques significatives du 1er modèle
tab1vec_c <- c("chol","hdlc","chol_hdlc","ldlc")
# tab1vec_m <- c("chiro-inositol","xylose","3-ureidopropionate","N6,N6-dimethyllysine","creatine","adenosine","3-carboxy-4-methyl-5-pentyl-2-furanpropionate (3-CMPFP)**",
#                "trigonelline (N'-methylnicotinate)","acisoga","octadecanedioylcarnitine (C18-DC)*",
#                "arachidonoylcarnitine (C20:4)","carotene diol (2)","carotene diol (1)","3beta-hydroxy-5-cholestenoate",
#                "octanoylcarnitine (C8)","tetradecadienoate (14:2)*","linoleate (18:2n6)","dodecadienoate (12:2)*",
#                "nonanoylcarnitine (C9)","cis-4-decenoylcarnitine (C10:1)")          ##With the same name than your DB

df <- loadsDA1
original_row_names <- rownames(df)
ordered_df <- df[order(abs(df$comp1)), , drop = FALSE]
rownames(ordered_df) <- original_row_names[order(abs(df$comp1))]
ordered_df

comp1 <- ordered_df
comp1$metab <- row.names(comp1)
comp1 <- comp1[180:200,]
comp1.list <- c(comp1$metab) # TOP 50 mtabolites of comp1

df <- loadsDA2
original_row_names <- rownames(df)
ordered_df <- df[order(abs(df$comp2)), , drop = FALSE]
rownames(ordered_df) <- original_row_names[order(abs(df$comp2))]
ordered_df

comp2 <- ordered_df
comp2$metab <- row.names(comp2)
comp2.list <- c(comp2$metab) # all metabolites of comp2 (20)
# comp2.list <- comp2.list[-12] # remove indolepropionate (already in comp 1)

tab1vec_m <- c(comp1.list,comp2.list)

tab1vec_m <- make_clean_names(tab1vec_m)

log.bleuet4$subject_id <- row.names(log.bleuet4)
contrats_data <- merge(bleuets.lme,log.bleuet4, by="subject_id")

# tab1vec_m <- c("isovaleric_acid")
# tab1vec_c <- c("hdl")
colnames_data <- colnames(contrats_data)
colnames_data <- make_clean_names(colnames_data)
colnames(contrats_data) <- colnames_data

#Generate lists to store results
t1a_models <- list()
t1a_models <- lapply(tab1vec_c,  list)

t1a_summaries <- list()
t1a_summaries <- lapply(tab1vec_c,  list)

list.var.clin <-list
list.var.clin <- lapply(tab1vec_c,  list)
list.agv <- list()
list.agv <- lapply(tab1vec_c,  list)
list.contrast <- list()
list.contrast <- lapply(tab1vec_c,  list)
list.estimate <- list()
list.estimate <- lapply(tab1vec_c,  list)
list.se <- list()
list.se <- lapply(tab1vec_c,  list)
list.df <- list()
list.df <- lapply(tab1vec_c,  list)
list.t.ratio <- list()
list.t.ratio <- lapply(tab1vec_c,  list)
list.p.value <- list()
list.p.value <- lapply(tab1vec_c,  list)

##Model for differences
for (j in seq_along(tab1vec_c)) {
  for (i in seq_along(tab1vec_m)) {
    t1a_models[[j]][[i]]    <- lmer(formula = paste0(tab1vec_c[[j]], " ~ sex + age + ", 
                                                     tab1vec_m[[i]], "*cluster*visites + (1|nomat)"), 
                                    REML = TRUE, data = contrats_data, na.action = "na.exclude")
    t1a_summaries[[j]][[i]] <- coef(summary(t1a_models[[j]][[i]]))[7,c(1,2,5)]
    
    #check all contrasts
    emt <- emtrends(t1a_models[[j]][[i]], consec ~ cluster, var = names(t1a_models[[j]][[i]]@frame[4]))
    
    #build the results table
    list.var.clin[[j]][[i]] <- tab1vec_c[[j]]
    list.agv[[j]][[i]] <- tab1vec_m[[i]]
    list.contrast[[j]][[i]] <- unlist(emt$contrast@grid)
    list.estimate[[j]][[i]] <- as.data.frame(emt$contrasts)[1,2]
    list.se[[j]][[i]] <- as.data.frame(emt$contrasts)[1,3]
    list.df[[j]][[i]] <- as.data.frame(emt$contrasts)[1,4]
    list.t.ratio[[j]][[i]] <- as.data.frame(emt$contrasts)[1,5]
    list.p.value[[j]][[i]] <- as.data.frame(emt$contrasts)[1,6]
  }
}

df.var.clin <- data.frame(matrix(unlist(list.var.clin)))
df.agv <- data.frame(matrix(unlist(list.agv)))
df.contrast <- data.frame(matrix(unlist(list.contrast)))
df.estimate <- data.frame(matrix(unlist(list.estimate)))
df.se <- data.frame(matrix(unlist(list.se)))
df.df <- data.frame(matrix(unlist(list.df)))
df.t <- data.frame(matrix(unlist(list.t.ratio)))
df.p <- data.frame(matrix(unlist(list.p.value)))

results.2 <- cbind(df.var.clin,df.agv,df.contrast,df.estimate,df.se,df.df,df.t,df.p)
names(results.2) <- c("clinique","agv","contrast","estimate","se","df","t.ratio","p.value")
results.2[4:8] <- round(results.2[4:8],3)
results.2$sig <- ifelse(results.2$p.value < 0.05,"*", "")
results.2 # les résultats de tous les contrasts

results.2.sig <- subset(results.2, results.2$sig == "*")

###Add labels and set colnames for each dependent variable
names(t1a_summaries) <- tab1vec_c
###Adding names on nested lists
for(i in seq_along(tab1vec_c)) { # Run for-nested lists
  names(t1a_summaries[[i]]) <- tab1vec_m
}

###Transform in dataframes
for(i in seq_along(tab1vec_c)) { # Run for-nested lists
  t1a_summaries[[i]] <- round(do.call(rbind, t1a_summaries[[i]]), digits = 4)
}

t1a_summaries   

write.csv(t1a_summaries$ldlc, "ldlc.csv")

# contrats.sig <- 

#CHOL
emtrends(t1a_models[[1]][[2]], consec ~ cluster, var = "glutamate")
emtrends(t1a_models[[1]][[3]], consec ~ cluster, var = "glycine") #nope
emtrends(t1a_models[[1]][[5]], consec ~ cluster, var = "isovalerylglycine") #nope
emtrends(t1a_models[[1]][[7]], consec ~ cluster, var = "phenylacetylglutamine")
emtrends(t1a_models[[1]][[9]], consec ~ cluster, var = "x1_linolenoylglycerol_18_3") #nope
emtrends(t1a_models[[1]][[10]], consec ~ cluster, var = "glycosyl_ceramide_d18_1_23_1_d17_1_24_1") #nope
emtrends(t1a_models[[1]][[11]], consec ~ cluster, var = "x1_methyl_5_imidazoleacetate")
emtrends(t1a_models[[1]][[17]], consec ~ cluster, var = "x5_hydroxyindole_sulfate")
emtrends(t1a_models[[1]][[22]], consec ~ cluster, var = "x5_6_dihydrothymine") #nope
emtrends(t1a_models[[1]][[27]], consec ~ cluster, var = "maltose")
emtrends(t1a_models[[1]][[30]], consec ~ cluster, var = "x4_guanidinobutanoate") #NOPE
emtrends(t1a_models[[1]][[32]], consec ~ cluster, var = "x4_methoxyphenol_sulfate")
emtrends(t1a_models[[1]][[39]], consec ~ cluster, var = "heme")
# END CHOL
#HDLC
emtrends(t1a_models[[2]][[11]], consec ~ cluster, var = "x1_methyl_5_imidazoleacetate")
emtrends(t1a_models[[2]][[15]], consec ~ cluster, var = "x3_indoxyl_sulfate")
emtrends(t1a_models[[2]][[17]], consec ~ cluster, var = "x5_hydroxyindole_sulfate") #nope
emtrends(t1a_models[[2]][[32]], consec ~ cluster, var = "x4_methoxyphenol_sulfate") #nope
# END HDLC
# CHOL HDLC
emtrends(t1a_models[[3]][[2]], consec ~ cluster, var = "glutamate") #+++++
emtrends(t1a_models[[3]][[3]], consec ~ cluster, var = "glycine") #nope
emtrends(t1a_models[[3]][[4]], consec ~ cluster, var = "dihomo_linolenoyl_choline") #nope
emtrends(t1a_models[[3]][[7]], consec ~ cluster, var = "phenylacetylglutamine")
emtrends(t1a_models[[3]][[10]], consec ~ cluster, var = "glycosyl_ceramide_d18_1_23_1_d17_1_24_1") #nope
emtrends(t1a_models[[3]][[11]], consec ~ cluster, var = "x1_methyl_5_imidazoleacetate")
emtrends(t1a_models[[3]][[14]], consec ~ cluster, var = "x1_stearoyl_2_docosahexaenoyl_gpe_18_0_22_6") #nope
emtrends(t1a_models[[3]][[17]], consec ~ cluster, var = "x5_hydroxyindole_sulfate")
emtrends(t1a_models[[3]][[29]], consec ~ cluster, var = "x2_hydroxyoctanoate")
emtrends(t1a_models[[3]][[30]], consec ~ cluster, var = "x4_guanidinobutanoate")
emtrends(t1a_models[[3]][[32]], consec ~ cluster, var = "x4_methoxyphenol_sulfate")
emtrends(t1a_models[[3]][[34]], consec ~ cluster, var = "lithocholate_sulfate_1") #nope
emtrends(t1a_models[[3]][[37]], consec ~ cluster, var = "dopamine_3_o_sulfate") #nope
emtrends(t1a_models[[3]][[38]], consec ~ cluster, var = "ximenoylcarnitine_c26_1")
emtrends(t1a_models[[3]][[39]], consec ~ cluster, var = "heme") #nope
#END CHOL HDLC
#LDLC
emtrends(t1a_models[[4]][[7]], consec ~ cluster, var = "phenylacetylglutamine")
emtrends(t1a_models[[4]][[8]], consec ~ cluster, var = "phenylacetylglutamate") #nope
emtrends(t1a_models[[4]][[11]], consec ~ cluster, var = "x1_methyl_5_imidazoleacetate")
emtrends(t1a_models[[3]][[17]], consec ~ cluster, var = "x5_hydroxyindole_sulfate")
emtrends(t1a_models[[3]][[27]], consec ~ cluster, var = "maltose")
emtrends(t1a_models[[3]][[17]], consec ~ cluster, var = "x5_hydroxyindole_sulfate")
emtrends(t1a_models[[3]][[39]], consec ~ cluster, var = "heme") #nope
#END LDLC



# END NEW MODEL

library(sjPlot)
library(sjmisc)
library(ggplot2)
library(ggalluvial)
data(efc)
theme_set(theme_sjplot())

colors <- c("#E41A1C", "#377EB8")
# colors <- c("green", "yellow")
colors <- c("#FA407B","#46B553")
symbols <- c(19, 21)

chol.2 <- plot_model(t1a_models[[1]][[2]], type = "emm", terms = c("glutamate","cluster"), # BOF
                            title="Chol/Maltose, p=0.0468", axis.title =c("Glutamate","Cholesterol"), colors = colors)#, col = colors)
chol.2 

chol.7 <- plot_model(t1a_models[[1]][[7]], type = "emm", terms = c("phenylacetylglutamine","cluster"), # BOF
                         title="Chol/phenylacetylglutamine, p=0.0013", axis.title =c("phenylacetylglutamine","Cholesterol"), colors = colors)#, col = colors)
chol.7 #nope

chol.11 <- plot_model(t1a_models[[1]][[11]], type = "emm", terms = c("x1_methyl_5_imidazoleacetate","cluster"), # BOF
                         title="Chol/x1_methyl_5_imidazoleacetate, p=0.002", axis.title =c("x1_methyl_5_imidazoleacetate","Cholesterol"), colors = colors)#, col = colors)
chol.11

chol.17 <- plot_model(t1a_models[[1]][[17]], type = "emm", terms = c("x5_hydroxyindole_sulfate","cluster"), # BOF
                         title="Chol/x5_hydroxyindole_sulfate, p=0.0068", axis.title =c("x5_hydroxyindole_sulfate","Cholesterol"), colors = colors)#, col = colors)
chol.17

chol.27 <- plot_model(t1a_models[[1]][[27]], type = "emm", terms = c("maltose","cluster"), # BOF
                        title="Chol/maltose, p=0.005", axis.title =c("maltose","Cholesterol"), colors = colors)#, col = colors)
chol.27

chol.32  <- plot_model(t1a_models[[1]][[32]], type = "emm", terms = c("x4_methoxyphenol_sulfate ","cluster"), # BOF
                                    title="Chol/x4_methoxyphenol_sulfate , p=0.0317", axis.title =c("x4_methoxyphenol_sulfate ","Cholesterol"), colors = colors)#, col = colors)
chol.32


chol.39 <- plot_model(t1a_models[[1]][[39]], type = "emm", terms = c("heme","cluster"), # BOF
                        title="Chol/heme, p=0.0145", axis.title =c("heme","Cholesterol"), colors = colors)#, col = colors)
chol.39

hdlc.11 <- plot_model(t1a_models[[2]][[11]], type = "emm", terms = c("x1_methyl_5_imidazoleacetate","cluster"), # BOF
                      title="Hdlc/x1_methyl_5_imidazoleacetate, p=0.0426", axis.title =c("x1_methyl_5_imidazoleacetate","Hdlc"), colors = colors)#, col = colors)
hdlc.11

hdlc.15 <- plot_model(t1a_models[[2]][[15]], type = "emm", terms = c("x3_indoxyl_sulfate","cluster"), # BOF
                      title="Hdlc/x3_indoxyl_sulfate, p=", axis.title =c("x3_indoxyl_sulfate","Hdlc"), colors = colors)#, col = colors)
hdlc.15


chol_hdlc.2 <- plot_model(t1a_models[[3]][[2]], type = "emm", terms = c("glutamate","cluster"), # BOF
                      title="Chol_Hdlc/glutamate, p<0.0001", axis.title =c("glutamate","Chol_hdlc"), colors = colors)#, col = colors)
chol_hdlc.2

chol_hdlc.7 <- plot_model(t1a_models[[3]][[7]], type = "emm", terms = c("phenylacetylglutamine","cluster"), # BOF
                          title="Chol_Hdlc/phenylacetylglutamine, p=0.0013", axis.title =c("phenylacetylglutamine","Chol_hdlc"), colors = colors)#, col = colors)
chol_hdlc.7

chol_hdlc.11 <- plot_model(t1a_models[[3]][[11]], type = "emm", terms = c("x1_methyl_5_imidazoleacetate","cluster"), # BOF
                          title="Chol_Hdlc/x1_methyl_5_imidazoleacetate, p=0.0039", axis.title =c("x1_methyl_5_imidazoleacetate","Chol_hdlc"), colors = colors)#, col = colors)
chol_hdlc.11

chol_hdlc.17 <- plot_model(t1a_models[[3]][[17]], type = "emm", terms = c("x5_hydroxyindole_sulfate","cluster"), # BOF
                          title="Chol_Hdlc/x5_hydroxyindole_sulfate, p=0.0010", axis.title =c("x5_hydroxyindole_sulfate","Chol_hdlc"), colors = colors)#, col = colors)
chol_hdlc.17

chol_hdlc.29 <- plot_model(t1a_models[[3]][[29]], type = "emm", terms = c("x2_hydroxyoctanoate","cluster"), # BOF
                          title="Chol_Hdlc/x2_hydroxyoctanoate, p=0.0226", axis.title =c("x2_hydroxyoctanoate","Chol_hdlc"), colors = colors)#, col = colors)
chol_hdlc.29

chol_hdlc.30 <- plot_model(t1a_models[[3]][[30]], type = "emm", terms = c("x4_guanidinobutanoate","cluster"), # BOF
                          title="Chol_Hdlc/x4_guanidinobutanoate, p=0.0406", axis.title =c("x4_guanidinobutanoate","Chol_hdlc"), colors = colors)#, col = colors)
chol_hdlc.30

chol_hdlc.32 <- plot_model(t1a_models[[3]][[32]], type = "emm", terms = c("x4_methoxyphenol_sulfate","cluster"), # BOF
                           title="Chol_Hdlc/x4_methoxyphenol_sulfate, p=0.0317", axis.title =c("x4_methoxyphenol_sulfate","Chol_hdlc"), colors = colors)#, col = colors)
chol_hdlc.32

chol_hdlc.38 <- plot_model(t1a_models[[3]][[38]], type = "emm", terms = c("ximenoylcarnitine_c26_1","cluster"), # BOF
                           title="Chol_Hdlc/ximenoylcarnitine_c26_1, p=0.0048", axis.title =c("ximenoylcarnitine_c26_1","Chol_hdlc"), colors = colors)#, col = colors)
chol_hdlc.38

ldlc.7 <- plot_model(t1a_models[[4]][[7]], type = "emm", terms = c("phenylacetylglutamine","cluster"), # BOF
                           title="ldlc/phenylacetylglutamine, p=0.0201", axis.title =c("phenylacetylglutamine","ldlc"), colors = colors)#, col = colors)
ldlc.7

ldlc.11 <- plot_model(t1a_models[[4]][[11]], type = "emm", terms = c("x1_methyl_5_imidazoleacetate","cluster"), # BOF
                     title="ldlc/x1_methyl_5_imidazoleacetate, p=0.0228", axis.title =c("x1_methyl_5_imidazoleacetate","ldlc"), colors = colors)#, col = colors)
ldlc.11

ldlc.17 <- plot_model(t1a_models[[4]][[17]], type = "emm", terms = c("x5_hydroxyindole_sulfate","cluster"), # BOF
                     title="ldlc/x5_hydroxyindole_sulfate, p=0.0010", axis.title =c("x5_hydroxyindole_sulfate","ldlc"), colors = colors)#, col = colors)
ldlc.17

ldlc.27 <- plot_model(t1a_models[[4]][[27]], type = "emm", terms = c("maltose","cluster"), # BOF
                      title="ldlc/maltose, p=0.0146", axis.title =c("maltose","ldlc"), colors = colors)#, col = colors)
ldlc.27



# ADDED Rb here
rami.lme2 <- rami.lme
rami.lme2$visit <- ifelse(rami.lme2$Visites == "V2", "Pre Rb","Post Rb")
rami.lme2$cluster2 <- ifelse(rami.lme$cluster == "respc1", "R1", "R2")
rami.lme2$visit <- relevel(as.factor(rami.lme2$visit),"Pre Rb","Post Rb")

outlier <- c(72, 126) # CHECKER 3ème outlier et c("LARC-00307","LARC-00308","LARC-00347","LARC-00348")#,"LARC-00279", "LARC-00280")
rami.lme2 <- subset(rami.lme2, !rami.lme2$nomat %in% outlier)


rami.lme3 <- cbind(rami.lme2,ure=log.rami6$`3-ureidopropionate`,chiro=log.rami6$`chiro-inositol`, xylo=log.rami6$xylose, caro=log.rami6$`carotene diol (1)`, 
                   cmp=log.rami6$`3-carboxy-4-methyl-5-pentyl-2-furanpropionate (3-CMPFP)**`, creat=log.rami6$creatine,
                   nona=log.rami6$`nonanoylcarnitine (C9)`,octa=log.rami6$`octadecanedioylcarnitine (C18-DC)*`)
box.3CMPFP <- ggplot(rami.lme3, aes(x = reorder(visit, cmp), y = cmp, fill = cluster2)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_point(aes(fill = cluster2), size = 2, shape = 21, position = position_jitterdodge()) +
  facet_wrap(~cluster2) +
  theme(legend.position = "right") +
  labs(title = "A) 3-CMPFP, p=0.12", x = "Visits", y = "3-CMPFP", color = "Subgroup") +
  theme_fivethirtyeight() +
  theme(axis.title = element_text(), plot.title = element_text(size=15, hjust = 0.5),
        plot.background = element_blank(),panel.background = element_blank(), legend.background = element_blank(), 
        legend.key = element_blank(), strip.background = element_blank(), strip.text.x = element_blank()) +
  scale_fill_manual(values = c('#FA407B', '#46B553')) +#c("#C63F23", "#278F3D")+
  theme(axis.text = element_text(size = 14)) +
  theme(axis.title = element_text(size = 16)) 
box.3CMPFP

box.creat <- ggplot(rami.lme3, aes(x = visit, y = creat, fill = cluster2)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_point(aes(fill = cluster2), size = 2, shape = 21, position = position_jitterdodge()) +
  facet_wrap(~cluster2) +
  theme(legend.position = "right") +
  labs(title = "Creatine", x = "Visites", y = "Creatine", color = "Subgroup") +
  theme_fivethirtyeight() +
  theme(axis.title = element_text(), plot.title = element_text(size=12, hjust = 0.5),
        plot.background = element_blank(),panel.background = element_blank(), legend.background = element_blank(), 
        legend.key = element_blank(), strip.background = element_blank(), strip.text.x = element_blank()) +
  scale_fill_manual(values = c('#E41A1C', '#377EB8')) #c("#C63F23", "#278F3D")
box.creat

box.octa <- ggplot(rami.lme3, aes(x = visit, y = octa, fill = cluster2)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_point(aes(fill = cluster2), size = 2, shape = 21, position = position_jitterdodge()) +
  facet_wrap(~cluster2) +
  theme(legend.position = "right") +
  labs(title = "Octadecanedioylcarnitine (C18-DC)*", x = "Visits", y = "Octadecanedioylcarnitine", color = "Subgroup") +
  theme_fivethirtyeight() +
  theme(axis.title = element_text(), plot.title = element_text(size=12, hjust = 0.5),
        plot.background = element_blank(),panel.background = element_blank(), legend.background = element_blank(), 
        legend.key = element_blank(), strip.background = element_blank(), strip.text.x = element_blank()) +
  scale_fill_manual(values = c('#E41A1C', '#377EB8')) #c("#C63F23", "#278F3D")
box.octa

box.nona <- ggplot(rami.lme3, aes(x = visit, y = nona, fill = cluster2)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_point(aes(fill = cluster2), size = 2, shape = 21, position = position_jitterdodge()) +
  facet_wrap(~cluster2) +
  theme(legend.position = "right") +
  labs(title = "Nonanoylcarnitine (C9)", x = "Visits", y = "Nonanoylcarnitine", color = "Subgroup") +
  theme_fivethirtyeight() +
  theme(axis.title = element_text(), plot.title = element_text(size=12, hjust = 0.5),
        plot.background = element_blank(),panel.background = element_blank(), legend.background = element_blank(), 
        legend.key = element_blank(), strip.background = element_blank(), strip.text.x = element_blank()) +
  scale_fill_manual(values = c('#E41A1C', '#377EB8')) #c("#C63F23", "#278F3D")
box.nona

box.chiro <- ggplot(rami.lme3, aes(x = visit, y = chiro, fill = cluster2)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_point(aes(fill = cluster2), size = 2, shape = 21, position = position_jitterdodge()) +
  facet_wrap(~cluster2) +
  theme(legend.position = "right") +
  labs(title = "B) Chiro-Inositol, p=0.07", x = "Visits", y = "Chiro Inositol", color = "Subgroup") +
  theme_fivethirtyeight() +
  theme(axis.title = element_text(), plot.title = element_text(size=15, hjust = 0.5),
        plot.background = element_blank(),panel.background = element_blank(), legend.background = element_blank(), 
        legend.key = element_blank(), strip.background = element_blank(), strip.text.x = element_blank()) +
  scale_fill_manual(values = c('#FA407B', '#46B553')) +#c("#C63F23", "#278F3D")+
  theme(axis.text = element_text(size = 14)) +
  theme(axis.title = element_text(size = 16)) 
box.chiro

box.caro <- ggplot(rami.lme3, aes(x = reorder(visit, caro), y = caro, fill = cluster2)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_point(aes(fill = cluster2), size = 2, shape = 21, position = position_jitterdodge()) +
  facet_wrap(~cluster2) +
  theme(legend.position = "right") +
  labs(title = "C) Carotene diol 1, p=0.45", x = "Visits", y = "Carotene diol 1", color = "Subgroup") +
  theme_fivethirtyeight() +
  theme(axis.title = element_text(), plot.title = element_text(size=15, hjust = 0.5),
        plot.background = element_blank(),panel.background = element_blank(), legend.background = element_blank(), 
        legend.key = element_blank(), strip.background = element_blank(), strip.text.x = element_blank()) +
  scale_fill_manual(values = c('#FA407B', '#46B553')) +#c("#C63F23", "#278F3D")+
  theme(axis.text = element_text(size = 14)) +
  theme(axis.title = element_text(size = 16)) 
box.caro

box.ure <- ggplot(rami.lme3, aes(x = reorder(visit, ure), y = ure, fill = cluster2)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_point(aes(fill = cluster2), size = 2, shape = 21, position = position_jitterdodge()) +
  facet_wrap(~cluster2) +
  theme(legend.position = "right") +
  labs(title = "3-Ureidopropionate", x = "Visits", y = "3-Ureidopropionate", color = "Subgroup") +
  theme_fivethirtyeight() +
  theme(axis.title = element_text(), plot.title = element_text(size=12, hjust = 0.5),
        plot.background = element_blank(),panel.background = element_blank(), legend.background = element_blank(), 
        legend.key = element_blank(), strip.background = element_blank(), strip.text.x = element_blank()) +
  scale_fill_manual(values = c('#E41A1C', '#377EB8')) #c("#C63F23", "#278F3D")
box.ure

p.wai4.2 <- ggplot(rami.lme2, aes(x = reorder(visit, waistc), y = waistc, fill = cluster2)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_point(aes(fill = cluster2), size = 2, shape = 21, position = position_jitterdodge()) +
  facet_wrap(~cluster2) +
  theme(legend.position = "right") +
  labs(title = "Waist circumference", x = "Visits", y = "Waist circumference (cm)", color = "Group:") +
  theme_fivethirtyeight() +
  theme(axis.title = element_text(), plot.title = element_text(size=12, hjust = 0.5),
        plot.background = element_blank(),panel.background = element_blank(), legend.background = element_blank(), 
        legend.key = element_blank(), strip.background = element_blank(), strip.text.x = element_blank()) +
  scale_fill_manual(values = c('#FA407B', '#46B553'))
p.wai4.2

p.whr4.2 <- ggplot(rami.lme2, aes(x = reorder(visit, whr), y = whr, fill = cluster2)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_point(aes(fill = cluster2), size = 2, shape = 21, position = position_jitterdodge()) +
  facet_wrap(~cluster2) +
  theme(legend.position = "right") +
  labs(title = "Waist-to-hip ratio", x = "Visits", y = "Waist-to-hip ratio", color = "Group:") +
  theme_fivethirtyeight() +
  theme(axis.title = element_text(), plot.title = element_text(size=12, hjust = 0.5),
        plot.background = element_blank(),panel.background = element_blank(), legend.background = element_blank(), 
        legend.key = element_blank(), strip.background = element_blank(), strip.text.x = element_blank()) +
  scale_fill_manual(values = c('#FA407B', '#46B553'))
p.whr4.2

p.homa.2 <- ggplot(rami.lme2, aes(x = reorder(visit, Homa), y = Homa, fill = cluster2)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_point(aes(fill = cluster2), size = 2, shape = 21, position = position_jitterdodge()) +
  facet_wrap(~cluster2) +
  theme(legend.position = "right") +
  labs(title = "HOMA-IR", x = "Visits", y = "HOMA-IR", color = "Group:") +
  theme_fivethirtyeight() +
  theme(axis.title = element_text(), plot.title = element_text(size=12, hjust = 0.5),
        plot.background = element_blank(),panel.background = element_blank(), legend.background = element_blank(), 
        legend.key = element_blank(), strip.background = element_blank(), strip.text.x = element_blank()) +
  scale_fill_manual(values = c('#FA407B', '#46B553'))
p.homa.2

p.quicki.2 <- ggplot(rami.lme2, aes(x = reorder(visit, QUICKI), y = QUICKI, fill = cluster2)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_point(aes(fill = cluster2), size = 2, shape = 21, position = position_jitterdodge()) +
  facet_wrap(~cluster2) +
  theme(legend.position = "right") +
  labs(title = "QUICKI", x = "Visits", y = "QUICKI", color = "Group:") +
  theme_fivethirtyeight() +
  theme(axis.title = element_text(), plot.title = element_text(size=12, hjust = 0.5),
        plot.background = element_blank(),panel.background = element_blank(), legend.background = element_blank(), 
        legend.key = element_blank(), strip.background = element_blank(), strip.text.x = element_blank()) +
  scale_fill_manual(values = c('#FA407B', '#46B553'))
p.quicki.2

box.mat <- ggplot(rami.lme2, aes(x = reorder(visit, Matsuda5t), y = Matsuda5t, fill = cluster2)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_point(aes(fill = cluster2), size = 2, shape = 21, position = position_jitterdodge()) +
  facet_wrap(~cluster2) +
  theme(legend.position = "right") +
  labs(title = "Matsuda", x = "Visits", y = "Matsuda", color = "Group:") +
  theme_fivethirtyeight() +
  theme(axis.title = element_text(), plot.title = element_text(size=12, hjust = 0.5),
        plot.background = element_blank(),panel.background = element_blank(), legend.background = element_blank(), 
        legend.key = element_blank(), strip.background = element_blank(), strip.text.x = element_blank()) +
  scale_fill_manual(values = c('#FA407B', '#46B553'))
box.mat


contrasts.figure <- ggarrange(whr.x3carboxy, waistc.x3carboxy, quicki.creatine, quicki.chiro, 
                              homa.octa, homa.nona, caro.in, ure.mat, nrow=4, ncol=2, common.legend=T,legend=NULL, 
                              font.label = list(size = 14, color = "black", face = "bold", family = NULL))

# 2 rows 2 columns
contrasts.figure2 <- ggarrange(quicki.chiro, 
                               caro.in, nrow=1, ncol=2, common.legend=T,legend=NULL, 
                              font.label = list(size = 14, color = "black", face = "bold", family = NULL))

contrasts.figure2.1 <- ggarrange(waistc.x3carboxy,quicki.chiro,caro.in, nrow=1, ncol=3, common.legend=T,legend=NULL, 
                                 font.label = list(size = 14, color = "black", face = "bold", family = NULL))

# only 1 row
contrasts.figure3 <- ggarrange(waistc.x3carboxy, quicki.chiro, 
                               caro.in, nrow=1, ncol=3, common.legend=T,legend=NULL, 
                               font.label = list(size = 14, color = "black", face = "bold", family = NULL))


boxplots.figure.clin <- ggarrange(p.wai4.2, p.quicki.2, nrow=1, ncol=2, common.legend=T,legend=NULL, 
                             font.label = list(size = 14, color = "black", face = "bold", family = NULL))

boxplots.figure.meta <- ggarrange(box.3CMPFP, box.chiro, box.caro, 
                                  nrow=1, ncol=3, common.legend=T,legend=T, 
                                  font.label = list(size = 14, color = "black", face = "bold", family = NULL))

boxplots.figure <- ggarrange(boxplots.figure.clin, boxplots.figure.meta,
                             nrow=2, ncol=1, common.legend=T,legend=NULL, 
                             font.label = list(size = 14, color = "black", face = "bold", family = NULL))


# final.figure <- ggarrange(contrasts.figure,boxplots.figure, nrow=2, ncol=1, common.legend=T,legend=NULL, heights=c(2,3),
                          # font.label = list(size = 14, color = "black", face = "bold", family = NULL))


# png(filename = "contrasts2.png", width = 14, height = 18, units = 'in', res = 300)
# annotate_figure(final.figure, top = text_grob("", color = "Black", face = "bold", size = 14))
# dev.off()

# png(filename = "contrasts.png", width = 14, height = 18, units = 'in', res = 300)
annotate_figure(contrasts.figure, top = text_grob("", color = "Black", face = "bold", size = 14))
# dev.off()
# 
# png(filename = "clin.png", width = 14, height = 16, units = 'in', res = 300)
annotate_figure(boxplots.figure.clin, top = text_grob("", color = "Black", face = "bold", size = 14))
# dev.off()

png(filename = "metabo.box.png", width = 15, height = 10, units = 'in', res = 300)
annotate_figure(boxplots.figure.meta, top = text_grob("", color = "Black", face = "bold", size = 14))
dev.off()

png(filename = "contrasts2.1.png", width = 18, height = 5, units = 'in', res = 300)
annotate_figure(contrasts.figure2.1, top = text_grob("", color = "Black", face = "bold", size = 14))
dev.off()

# png(filename = "contrasts3.png", width = 18, height = 6, units = 'in', res = 300)
annotate_figure(contrasts.figure3, top = text_grob("", color = "Black", face = "bold", size = 14))
# dev.off()

# png(filename = "boxplots2.png", width = 14, height = 10, units = 'in', res = 300)
annotate_figure(boxplots.figure, top = text_grob("", color = "Black", face = "bold", size = 14))
# dev.off()
