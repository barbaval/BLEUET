library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
library("tableone")
library(fmsb)
library(tfplot)


flowchart <- grViz("digraph flowchart {
    #Graphics and text characteristics
    node [fontname = helvetica, fontsize = 24,
          shape = rectangle, style = rounded,
          penwidth = 2, width = 5]
    edge [penwidth = 2]

 

    #Definition of nodes and labels
    tab1 [label = 'Total analyzed metabolites \\n n=1132']
        blank1[label = '', width = 0.001, height = 0.001]
        excluded1[label = 'Metabolites showing no variance in samples \\n n=14']
    tab2 [label = 'Total metabolites showing variance in samples \\n n=1118']
        blank2[label = '', width = 0.001, height = 0.001]
        excluded2[label = 'Unknown metabolites \\n n=216']
    tab3 [label = 'Total known metabolites \\n n=902']
        blank3[label = '', width = 0.001, height = 0.001]
        excluded3[label = 'Total metabolites showing variance < 0.1 in samples \\n n=272']
    tab4 [label = 'Total remaining metabolites \\n n=630']
        blank4[label = '', width = 0.001, height = 0.001]
        excluded4[label = 'Xenobiotics \\n n=120']
    tab5 [label = 'Total remaining metabolites \\n n=510']
        blank5[label = '', width = 0.001, height = 0.001]
        excluded5[label = 'Partially described molecules \\n n=24']
    tab6 [label = 'Total remaining metabolites \\n n=486']
    

    #Edge definitions with the node IDs 
    tab1 -> blank1[dir = none];
    blank1 -> excluded1[minlen = 2];
    blank1 -> tab2;
    tab2 -> blank2[dir = none];
    blank2 -> excluded2[minlen = 2];
    blank2 -> tab3;
    tab3 -> blank3[dir = none];
    blank3 -> excluded3[minlen = 2];
    blank3 -> tab4;
    tab4 -> blank4[dir = none];
    blank4 -> excluded4[minlen = 2];
    blank4 -> tab5;
    tab5 -> blank5[dir = none];
    blank5 -> excluded5[minlen = 2];
    blank5 -> tab6;

    {rank = same; blank1 excluded1} 
    {rank = same; blank2 excluded2} 
    {rank = same; blank3 excluded3} 
    {rank = same; blank4 excluded4}
    {rank = same; blank5 excluded5}
  }
    ")



flowchart

flowchart <- export_svg(flowchart)
flowchart <- charToRaw(flowchart)
rsvg_svg(flowchart,
         file =  "Flowchart.svg",
         height = 10)




flowchart <- export_svg(flowchart)
flowchart <- charToRaw(flowchart)
rsvg_pdf(flowchart, "supp1.png")


# Table one
table.variables <- c("sex","age","weight","bmi", "waistc","hipc","whr","sbp","dbp","apob_neph","chol","tg",
                     "hdlc", "ldlc","chol_hdlc","hba1c", "crp_neph","Glucosefast","Insulinefast","Matsuda5t","Homa","FIRI","QUICKI","nomat","Visites", "tx")
summary.table <-as.data.frame(data.r[,c(table.variables)])
outliers <- c(72,126)

responders <- c(39,69,90,96,98)
respondersneg <- c(77,89,32,54)
summary.table <- subset(summary.table, !summary.table$nomat %in% outliers)
summary.table[summary.table$nomat %in% responders, "cluster"] <- "R" #répondeurs sur composante 2
summary.table[!summary.table$nomat %in% responders, "cluster"] <- "NR" #non répondeurs
# summary.table[summary.table$nomat %in% respondersneg, "cluster"] <- "R-" #répondeurs sur composante 2

summary.table.w0 <- subset(summary.table, summary.table$Visites == "V2" & summary.table$tx == "B")
summary.table.w1 <- subset(summary.table, summary.table$Visites == "V4" & summary.table$tx == "B")
  

summary.table.w0.R <- subset(summary.table.w0, summary.table.w0$cluster == "R")
# summary.table.w0.R2 <- subset(summary.table.w0, summary.table.w0$cluster == "R-") 
summary.table.w0.NR <- subset(summary.table.w0, summary.table.w0$cluster == "NR") 

summary.table.w1.R <- subset(summary.table.w1, summary.table.w1$cluster == "R")
# summary.table.w1.R2 <- subset(summary.table.w1, summary.table.w1$cluster == "R-")
summary.table.w1.NR <- subset(summary.table.w1, summary.table.w1$cluster == "NR") 

summary.table.w0.final <- as.data.frame(c(table.variables))



list.var <- c("sex","age","weight","bmi", "waistc","hipc","whr","sbp","dbp","apob_neph","chol","tg",
                     "hdlc", "ldlc","chol_hdlc","hba1c", "crp_neph","Glucosefast","Insulinefast","Matsuda5t","Homa", "FIRI", "QUICKI")
catVars <- c("cluster")

table1 <- CreateTableOne(list.var, summary.table.w0, catVars, strata = c("cluster"))
table1
?CreateTableOne
# Table one for the 6 clusters


rami_pre <- read.csv("data/rami_pre.csv", row.names = 1)
rami_post <- read.csv("data/rami_post.csv", row.names = 1)

outliers.larc <- c("LARC-00307","LARC-00308","LARC-00347","LARC-00348")
rami_pre <- subset(rami_pre, !row.names(rami_pre) %in% outliers.larc)

rc2 <- c("LARC-00337","LARC-00339","LARC-00341","LARC-00287","LARC-00305","LARC-00338","LARC-00340","LARC-00342","LARC-00288","LARC-00306")
rami_pre[row.names(rami_pre) %in% rc2, "cluster"] <- "R" #répondeurs sur composante 2
rami_pre[!row.names(rami_pre) %in% rc2, "cluster"] <- "NR" #non répondeurs

# summary.table.meta.w0.R <- subset(rami_pre, rami_pre$cluster == "R") 
# summary.table.meta.w0.NR <- subset(rami_pre, rami_pre$cluster == "NR") 

list.var.meta <- c(names(rami_pre))

table2 <- CreateTableOne(list.var.meta, rami_pre, catVars, strata = c("cluster"))
table2



### PATHWAY RADAR
# https://r-graph-gallery.com/web-circular-barplot-with-R-and-ggplot2.html
# https://r-charts.com/ranking/radar-chart/
# https://www.datanovia.com/en/blog/beautiful-radar-chart-in-r-using-fmsb-and-ggplot-packages/
# library(dplyr)
# library(ggplot2)
# library(stringr)
# 

# library(fmsb)
# pathway <- read.csv("data/pathway_comp_1.csv")
# # row.names(pathway) <- pathway[,1]
# pathway1 <- pathway[, 1:2]
# row.names(pathway1) <- pathway1[,1]
# 
# pathway2 <- as.data.frame(t(pathway1))
# pathway2 <- pathway2[2,]
# pathway2 <- as.data.frame(read.csv("data/pathway2.csv"))
# 
# pathway2$Fatty.acid.Metabolism  <- as.numeric(pathway2$Alpha.Linolenic.Acid.and.Linoleic.Acid.Metabolism)
# pathway2 <- pathway2[,-1]
# 
# radarchart(pathway2, cglty = 1, cglcol = "gray",
#            pcol = 4, plwd = 2,
#            pfcol = rgb(0, 0.4, 1, 0.25))
# 
# plot_df <- pathway %>%
#   group_by(row.names(pathway)) %>%
#   summarise(
#     sum_length = sum(length_num),
#     n = n()
#   ) %>%
#   mutate(mean_gain = round(mean_gain, digits = 0))
# 

library(fmsb)

# substrRight <- function(x, n){
#   substr(x, nchar(x)-n+1, nchar(x))
# }

delta.rami.radar <- rami.norm.pre - rami.norm.post
colnames(delta.rami.radar) <- metabolites.rami$PLOT_NAME
responders.radar <- c("LARC-00337","LARC-00339","LARC-00341","LARC-00287","LARC-00305","LARC-00338","LARC-00340","LARC-00342","LARC-00288","LARC-00306")

delta.rami.radar$cluster <- ifelse(row.names(delta.rami.radar) %in% responders.radar, "R", "NR")
delta.rami.radar.comp1 <- delta.rami.radar[, names(delta.rami.radar) %in% c(row.names(loadsDA1.cor)[1:10],"cluster")]
delta.rami.radar.comp2 <- delta.rami.radar[, names(delta.rami.radar) %in% c(row.names(loadsDA2.cor)[1:10],"cluster")]


delta.rami.radar.comp2.r <- subset(delta.rami.radar.comp2, delta.rami.radar.comp2$cluster == "R")
delta.rami.radar.comp2.r <- delta.rami.radar.comp2.r[,1:10]

delta.rami.radar.comp2.nr <- subset(delta.rami.radar.comp2, delta.rami.radar.comp2$cluster == "NR")
delta.rami.radar.comp2.nr <- delta.rami.radar.comp2.nr[,1:10]

# Means table for component 1, R
delta.rami.radar.comp1.r <- subset(delta.rami.radar.comp1, delta.rami.radar.comp1$cluster == "R")
delta.rami.radar.comp1.r <- delta.rami.radar.comp1.r[,1:10]
delta.rami.radar.comp1.r.mean <- as.data.frame(sapply(delta.rami.radar.comp1.r,mean))
delta.rami.radar.comp1.r.mean <- as.data.frame(t(delta.rami.radar.comp1.r.mean))
delta.rami.radar.comp1.r.mean <- rbind(delta.rami.radar.comp1.r.mean, c(-1))
delta.rami.radar.comp1.r.mean <- rbind(delta.rami.radar.comp1.r.mean, c(1))
delta.rami.radar.comp1.r.mean <- delta.rami.radar.comp1.r.mean[c(3,2,1),]
row.names(delta.rami.radar.comp1.r.mean) <- c("max","min","mean")

# Means table for component 1, NR
delta.rami.radar.comp1.nr <- subset(delta.rami.radar.comp1, delta.rami.radar.comp1$cluster == "NR")
delta.rami.radar.comp1.nr <- delta.rami.radar.comp1.nr[,1:10]
delta.rami.radar.comp1.nr.mean <- as.data.frame(sapply(delta.rami.radar.comp1.nr,mean))
delta.rami.radar.comp1.nr.mean <- as.data.frame(t(delta.rami.radar.comp1.nr.mean))
delta.rami.radar.comp1.nr.mean <- rbind(delta.rami.radar.comp1.nr.mean, c(-1))
delta.rami.radar.comp1.nr.mean <- rbind(delta.rami.radar.comp1.nr.mean, c(1))
delta.rami.radar.comp1.nr.mean <- delta.rami.radar.comp1.nr.mean[c(3,2,1),]
row.names(delta.rami.radar.comp1.nr.mean) <- c("max","min","mean")

# Means table for component 1, R + NR
delta.rami.radar.comp1.all <- rbind(delta.rami.radar.comp1.r.mean, delta.rami.radar.comp1.nr.mean[3,])
row.names(delta.rami.radar.comp1.all) <- c("max","min","R2","R1")


# Means table for component 2, R
delta.rami.radar.comp2.r <- subset(delta.rami.radar.comp2, delta.rami.radar.comp2$cluster == "R")
delta.rami.radar.comp2.r <- delta.rami.radar.comp2.r[,1:10]
delta.rami.radar.comp2.r.mean <- as.data.frame(sapply(delta.rami.radar.comp2.r,mean))
delta.rami.radar.comp2.r.mean <- as.data.frame(t(delta.rami.radar.comp2.r.mean))
delta.rami.radar.comp2.r.mean <- rbind(delta.rami.radar.comp2.r.mean, c(-1))
delta.rami.radar.comp2.r.mean <- rbind(delta.rami.radar.comp2.r.mean, c(1))
delta.rami.radar.comp2.r.mean <- delta.rami.radar.comp2.r.mean[c(3,2,1),]
row.names(delta.rami.radar.comp2.r.mean) <- c("max","min","mean")

# Means table for component 2, NR
delta.rami.radar.comp2.nr <- subset(delta.rami.radar.comp2, delta.rami.radar.comp2$cluster == "NR")
delta.rami.radar.comp2.nr <- delta.rami.radar.comp2.nr[,1:10]
delta.rami.radar.comp2.nr.mean <- as.data.frame(sapply(delta.rami.radar.comp2.nr,mean))
delta.rami.radar.comp2.nr.mean <- as.data.frame(t(delta.rami.radar.comp2.nr.mean))
delta.rami.radar.comp2.nr.mean <- rbind(delta.rami.radar.comp2.nr.mean, c(-1))
delta.rami.radar.comp2.nr.mean <- rbind(delta.rami.radar.comp2.nr.mean, c(1))
delta.rami.radar.comp2.nr.mean <- delta.rami.radar.comp2.nr.mean[c(3,2,1),]
row.names(delta.rami.radar.comp2.nr.mean) <- c("max","min","mean")

# Means table for component 1, R + NR
delta.rami.radar.comp2.all <- rbind(delta.rami.radar.comp2.r.mean, delta.rami.radar.comp2.nr.mean[3,])
row.names(delta.rami.radar.comp2.all) <- c("max","min","R2","R1")


# With percent change
rami.norm.pre.radar <- rami.norm.pre
rami.norm.post.radar <- rami.norm.post

colnames(rami.norm.pre.radar) <- metabolites.rami$PLOT_NAME
colnames(rami.norm.post.radar) <- metabolites.rami$PLOT_NAME
radar.rami.norm.pre.comp1 <- rami.norm.pre.radar[, names(rami.norm.pre.radar) %in% c(row.names(loadsDA1.cor)[1:10],"cluster")]
radar.rami.norm.post.comp1 <- rami.norm.post.radar[, names(rami.norm.post.radar) %in% c(row.names(loadsDA1.cor)[1:10],"cluster")]

radar.rami.norm.comp1 <- rbind(radar.rami.norm.pre.comp1, radar.rami.norm.post.comp1)

rami.norm.radar <- rami.norm
row.names(rami.norm.radar) <- rami.norm.radar[,1]
rami.norm.radar$visit <- ifelse(rami.norm.radar$bef == "RAMI BEFORE TREATMENT", "v2","v4")
rami.norm.radar <- rami.norm.radar[,-c(1:2)]
colnames(rami.norm.radar) <- c(metabolites.rami$PLOT_NAME, "visit")
responders.radar <- c("LARC-00337","LARC-00339","LARC-00341","LARC-00287","LARC-00305","LARC-00338","LARC-00340","LARC-00342","LARC-00288","LARC-00306")
rami.norm.radar$cluster <- ifelse(row.names(rami.norm.radar) %in% responders.radar, "R", "NR")
rami.norm.radar.comp1 <- rami.norm.radar[, names(rami.norm.radar) %in% c(row.names(loadsDA1.cor)[1:10],"visit","cluster")]

rami.norm.radar.comp1.pre <- subset(rami.norm.radar.comp1, rami.norm.radar.comp1$visit == "v2")
rami.norm.radar.comp1.pre <- rami.norm.radar.comp1.pre[, !names(rami.norm.radar.comp1.pre) %in% c("visit")]
rami.norm.radar.comp1.post <- subset(rami.norm.radar.comp1, rami.norm.radar.comp1$visit == "v2")
rami.norm.radar.comp1.post <- rami.norm.radar.comp1.post[, !names(rami.norm.radar.comp1.post) %in% c("visit")]

rami.norm.radar.comp1.pre.r <- subset(rami.norm.radar.comp1.pre, rami.norm.radar.comp1.pre$cluster == "R")
rami.norm.radar.comp1.pre.r <- rami.norm.radar.comp1.pre.r[,!names(rami.norm.radar.comp1.pre.r) %in% c("cluster")]
rami.norm.radar.comp1.pre.nr <- subset(rami.norm.radar.comp1.pre, rami.norm.radar.comp1.pre$cluster == "NR")
rami.norm.radar.comp1.pre.nr <- rami.norm.radar.comp1.pre.nr[,!names(rami.norm.radar.comp1.pre.nr) %in% c("cluster")]


# Calculation of % change
# Component 1 responders
delta.rami.radar.comp1.r.perc <- delta.rami.radar.comp1.r/rami.norm.radar.comp1.pre.r
delta.rami.radar.comp1.r.perc.mean <- as.data.frame(sapply(delta.rami.radar.comp1.r.perc,mean))
delta.rami.radar.comp1.r.perc.mean <- as.data.frame(t(delta.rami.radar.comp1.r.perc.mean))
delta.rami.radar.comp1.r.perc.mean <- rbind(delta.rami.radar.comp1.r.perc.mean, c(-0.6))
delta.rami.radar.comp1.r.perc.mean <- rbind(delta.rami.radar.comp1.r.perc.mean, c(0.6))
delta.rami.radar.comp1.r.perc.mean <- delta.rami.radar.comp1.r.perc.mean[c(3,2,1),]
row.names(delta.rami.radar.comp1.r.perc.mean) <- c("max","min","mean")

# Component 1 non-responders
delta.rami.radar.comp1.nr.perc <- delta.rami.radar.comp1.nr/rami.norm.radar.comp1.pre.nr
delta.rami.radar.comp1.nr.perc.mean <- as.data.frame(sapply(delta.rami.radar.comp1.nr.perc,mean))
delta.rami.radar.comp1.nr.perc.mean <- as.data.frame(t(delta.rami.radar.comp1.nr.perc.mean))
delta.rami.radar.comp1.nr.perc.mean <- rbind(delta.rami.radar.comp1.nr.perc.mean, c(-0.6))
delta.rami.radar.comp1.nr.perc.mean <- rbind(delta.rami.radar.comp1.nr.perc.mean, c(0.6))
delta.rami.radar.comp1.nr.perc.mean <- delta.rami.radar.comp1.nr.perc.mean[c(3,2,1),]
row.names(delta.rami.radar.comp1.nr.perc.mean) <- c("max","min","mean")

# Means table for component 1, R + NR
delta.rami.radar.comp1.perc.all <- rbind(delta.rami.radar.comp1.r.perc.mean, delta.rami.radar.comp1.nr.perc.mean[3,])
row.names(delta.rami.radar.comp1.perc.all) <- c("max","min","R2","R1")



### COMPONENT 2 ###
rami.norm.radar.comp2 <- rami.norm.radar[, names(rami.norm.radar) %in% c(row.names(loadsDA2.cor)[1:10],"visit","cluster")]

rami.norm.radar.comp2.pre <- subset(rami.norm.radar.comp2, rami.norm.radar.comp2$visit == "v2")
rami.norm.radar.comp2.pre <- rami.norm.radar.comp2.pre[, !names(rami.norm.radar.comp2.pre) %in% c("visit")]
rami.norm.radar.comp2.post <- subset(rami.norm.radar.comp2, rami.norm.radar.comp2$visit == "v2")
rami.norm.radar.comp2.post <- rami.norm.radar.comp2.post[, !names(rami.norm.radar.comp2.post) %in% c("visit")]

rami.norm.radar.comp2.pre.r <- subset(rami.norm.radar.comp2.pre, rami.norm.radar.comp2.pre$cluster == "R")
rami.norm.radar.comp2.pre.r <- rami.norm.radar.comp2.pre.r[,!names(rami.norm.radar.comp2.pre.r) %in% c("cluster")]
rami.norm.radar.comp2.pre.nr <- subset(rami.norm.radar.comp2.pre, rami.norm.radar.comp2.pre$cluster == "NR")
rami.norm.radar.comp2.pre.nr <- rami.norm.radar.comp2.pre.nr[,!names(rami.norm.radar.comp2.pre.nr) %in% c("cluster")]


# Calculation of % change
# Component 2 responders
delta.rami.radar.comp2.r.perc <- delta.rami.radar.comp2.r/rami.norm.radar.comp2.pre.r
delta.rami.radar.comp2.r.perc.mean <- as.data.frame(sapply(delta.rami.radar.comp2.r.perc,mean))
delta.rami.radar.comp2.r.perc.mean <- as.data.frame(t(delta.rami.radar.comp2.r.perc.mean))
delta.rami.radar.comp2.r.perc.mean <- rbind(delta.rami.radar.comp2.r.perc.mean, c(-8))
delta.rami.radar.comp2.r.perc.mean <- rbind(delta.rami.radar.comp2.r.perc.mean, c(1))
delta.rami.radar.comp2.r.perc.mean <- delta.rami.radar.comp2.r.perc.mean[c(3,2,1),]
row.names(delta.rami.radar.comp2.r.perc.mean) <- c("max","min","mean")

# Component 2 non-responders
delta.rami.radar.comp2.nr.perc <- delta.rami.radar.comp2.nr/rami.norm.radar.comp2.pre.nr
delta.rami.radar.comp2.nr.perc.mean <- as.data.frame(sapply(delta.rami.radar.comp2.nr.perc,mean))
delta.rami.radar.comp2.nr.perc.mean <- as.data.frame(t(delta.rami.radar.comp2.nr.perc.mean))
delta.rami.radar.comp2.nr.perc.mean <- rbind(delta.rami.radar.comp2.nr.perc.mean, c(-8))
delta.rami.radar.comp2.nr.perc.mean <- rbind(delta.rami.radar.comp2.nr.perc.mean, c(1))
delta.rami.radar.comp2.nr.perc.mean <- delta.rami.radar.comp2.nr.perc.mean[c(3,2,1),]
row.names(delta.rami.radar.comp2.nr.perc.mean) <- c("max","min","mean")

# Means table for component 2, R + NR
delta.rami.radar.comp2.perc.all <- rbind(delta.rami.radar.comp2.r.perc.mean, delta.rami.radar.comp2.nr.perc.mean[3,])
row.names(delta.rami.radar.comp2.perc.all) <- c("max","min","R2","R1")



library(tfplot)

# https://stackoverflow.com/questions/21117554/time-series-in-r-how-do-i-calculate-percent-change-from-a-fixed-year-for-multip

# charts LANCER CODE DE LA FONCTION + RADARCHART + PARAMETRES EN MEME TEMPS, ex ligne 338 à 370
### USING ONLY NORMALIZED DATA DELTA pre-post, NOT % CHANGE
# Radarchart for component 1
create_beautiful_radarchart <- function(data, color = "#00AFBB", 
                                        vlabels = colnames(data), vlcex = 0.7,
                                        caxislabels = NULL, title = NULL, ...){
  radarchart(
    data, axistype = 1,
    # Customize the polygon
    pcol = color, pfcol = scales::alpha(color, 0.5), plwd = 2, plty = 1,
    # Customize the grid
    cglcol = "grey", cglty = 1, cglwd = 0.8,
    # Customize the axis
    axislabcol = "black", 
    # Variable labels
    vlcex = vlcex, vlabels = vlabels,
    caxislabels = caxislabels, title = title, ...
  )
}


# Reduce plot margin using par()
op <- par(mar = c(1, 2, 2, 2))
# Create the radar charts
create_beautiful_radarchart(
  data = delta.rami.radar.comp1.all, caxislabels = c(-1, -0.5, 0, 0.5, 1),
  color = c('#377EB8','#E41A1C'), 
  title = ""
)
# Add an horizontal legend
legend(
  x = "topright", legend = rownames(delta.rami.radar.comp1.all[-c(1,2),]), horiz = F,
  bty = "n", pch = 20 , col = c('#377EB8','#E41A1C'),
  text.col = "black", cex = 1, pt.cex = 1.5
)
par(op)


# Radarchart for component 2
create_beautiful_radarchart2 <- function(data, color = "#00AFBB", 
                                        vlabels = colnames(data), vlcex = 0.7,
                                        caxislabels = NULL, title = NULL, ...){
  radarchart(
    data, axistype = 1,
    # Customize the polygon
    pcol = color, pfcol = scales::alpha(color, 0.5), plwd = 2, plty = 1,
    # Customize the grid
    cglcol = "grey", cglty = 1, cglwd = 0.8,
    # Customize the axis
    axislabcol = "black",
    # Variable labels
    vlcex = vlcex, vlabels = vlabels,
    caxislabels = caxislabels, title = title, ...
  )
}

# Reduce plot margin using par()
op <- par(mar = c(1, 2, 2, 2))
# Create the radar charts
create_beautiful_radarchart2(
  data = delta.rami.radar.comp2.all, caxislabels = c(-1, -0.5, 0, 0.5, 1),
  color = c('#377EB8','#E41A1C'),
  title = ""
)
# Add an horizontal legend
legend(
  x = "topright", legend = rownames(delta.rami.radar.comp2.all[-c(1,2),]), horiz = F,
  bty = "n", pch = 20 , col = c('#377EB8','#E41A1C'),
  title = "Group",
  text.col = "black", cex = 1, pt.cex = 1.5
)
par(op)


### CHARTS FOR percent change %%%%%%%

# Radarchart for component 1
create_beautiful_radarchart3 <- function(data, color = "#00AFBB", 
                                         vlabels = colnames(data), vlcex = 0.7,
                                         caxislabels = NULL, title = NULL, ...){
  radarchart(
    data, axistype = 1,
    # Customize the polygon
    pcol = color, pfcol = scales::alpha(color, 0.5), plwd = 2, plty = 1,
    # Customize the grid
    cglcol = "grey", cglty = 1, cglwd = 0.8,
    # Customize the axis
    axislabcol = "grey", 
    # Variable labels
    vlcex = vlcex, vlabels = vlabels,
    caxislabels = caxislabels, title = title, ...
  )
}

# Reduce plot margin using par()
op <- par(mar = c(1, 2, 2, 2))
# Create the radar charts
create_beautiful_radarchart3(
  data = delta.rami.radar.comp1.perc.all, caxislabels = c(-0.6, -0.3, 0, 0.3, 0.6),
  color = c("#00AFBB", "#FC4E07"),
  title = "Mean percent change in component 1 metabolites per group"
)
# Add an horizontal legend
legend(
  x = "topright", legend = rownames(delta.rami.radar.comp1.perc.all[-c(1,2),]), horiz = TRUE,
  bty = "n", pch = 20 , col = c("#00AFBB", "#FC4E07"),
  text.col = "black", cex = 1, pt.cex = 1.5
)
par(op)


# Radarchart for component 2
create_beautiful_radarchart4 <- function(data, color = "#00AFBB", 
                                         vlabels = colnames(data), vlcex = 0.7,
                                         caxislabels = NULL, title = NULL, ...){
  radarchart(
    data, axistype = 1,
    # Customize the polygon
    pcol = color, pfcol = scales::alpha(color, 0.5), plwd = 2, plty = 1,
    # Customize the grid
    cglcol = "grey", cglty = 1, cglwd = 0.8,
    # Customize the axis
    axislabcol = "grey", 
    # Variable labels
    vlcex = vlcex, vlabels = vlabels,
    caxislabels = caxislabels, title = title, ...
  )
}

# Reduce plot margin using par()
op <- par(mar = c(1, 2, 2, 2))
# Create the radar charts
create_beautiful_radarchart4(
  data = delta.rami.radar.comp2.perc.all, caxislabels = c(-2, 0.5, 3, 5.5, 8),
  color = c("#00AFBB", "#FC4E07"),
  title = "Mean percent change in component 2 metabolites per group"
)
# Add an horizontal legend
legend(
  x = "topright", legend = rownames(delta.rami.radar.comp2.perc.all[-c(1,2),]), horiz = TRUE,
  bty = "n", pch = 20 , col = c("#00AFBB", "#FC4E07"),
  text.col = "black", cex = 1, pt.cex = 1.5
)
par(op)


# Data for pathway/enrichment analysis
pathway.all <- as.data.frame(names(log.rami6))
names(pathway.all) <- c("PLOT_NAME")
metabo.pathway <- metabolites.rami[ ,c("PLOT_NAME","HMDB","KEGG", "CAS", "PUBCHEM")]
# pathway.all <- merge(pathway.all, metabolites.rami$HMDB, pathway.all$Metabolite,metabolites.rami$PLOT_NAME)
pathway.all <- subset(metabo.pathway, metabo.pathway$PLOT_NAME %in% pathway.all$PLOT_NAME)


# delta
# pathway.delta <- subset(delta.rami2, colnames(delta.rami2) %in% pathway.all$PLOT_NAME)
# delta.rami2

delta.rami3 <- (rami.norm.post - rami.norm.pre)/rami.norm.post

delta.mean <- as.data.frame(colMeans(delta.rami3)) #change to delta.rami2 for simple difference
row.names(delta.mean) <- metabolites.rami$PLOT_NAME
delta.mean <- subset(delta.mean, rownames(delta.mean) %in% pathway.all$PLOT_NAME)
pathway.all <- cbind(pathway.all,delta.mean)

loadingsc1 <- plotLoadings(plsda.multilevel, comp=1)
loadingsc2 <- plotLoadings(plsda.multilevel, comp=2)

loads.splsda.c1 <-as.data.frame(plsda.multilevel$loadings$X)
loadsDA1.c1 <- as.data.frame(loads.splsda.c1[,1])
row.names(loadsDA1.c1) <- row.names(loads.splsda.c1)
names(loadsDA1.c1) <- c("loadingWeight")
loadsDA1 <-loadsDA1[order(abs(loadsDA1$comp1)), ,drop=F]

loads.splsda.fam <- subset(metabolites.rami,metabolites.rami$PLOT_NAME %in% rownames(loadsDA1),
                            c("PLOT_NAME","SUPER_PATHWAY","SUB_PATHWAY","HMDB", "KEGG","PUBCHEM"))
row.names(pathway.all)

# pathway.all[1:5,1:5]
head(loadsDA1.c1)


 pathway.all <- cbind(pathway.all,loadsDA1.c1)

