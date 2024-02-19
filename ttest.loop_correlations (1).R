# setwd("/Users/juan/Library/CloudStorage/OneDrive-UniversitéLaval/INAF/JOSIANE")
# FAS138.REP.AGV.CLIN<-read.csv("FAS participants/FAS138.REP.AGV.CLIN.csv",row.names = 1)

data.v2 <- datat #subset(FAS138.REP.AGV.CLIN, visite ==2)
data <- as.data.frame(data.v2)

######################################################################
############ t.test PRE-INTERVENTION (CLINIQUE ET FFQ) ###############
######################################################################
p.list<-list()
for (i in 2:length(data)) {   
  variable <- data[,i]
  ttest <- t.test(variable ~ cluster, data = data, na.action=na.omit,paired = F)
  p <- ttest$p.value
  p.list[[i]] <-p
  p.df<-data.frame(matrix(unlist(p.list)))
}

mean<-aggregate(data[,-1], list(data$cluster), FUN=mean,na.rm=TRUE)
sd<-aggregate(data[,-1], list(data$cluster), FUN=sd,na.rm=TRUE)

tmean<-t(mean[,-1])
tsd<-t(sd[,-1])

results <- cbind("R"=tmean[,2],"sd"=tsd[,2],"NR"=tmean[,1],"sd"=tsd[,1],"p-value"=p.df[,1])
results <- as.data.frame(round(results,2))
results$sig <- ifelse(results$`p-value` < 0.05,"*", "")
results

  #########################################################
  ####################### correlation #####################
  ########################################################
  #http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software

  data.v2 <- subset(FAS138.REP.AGV.CLIN, visite ==2)
  data.v2c <- as.data.frame(data.v2[,c(4:9,15:52)])
  
  data.v3 <- subset(FAS138.REP.AGV.CLIN, visite ==3)
  data.v3c <- as.data.frame(data.v3[,c(4:9,15:52)])
  
  data.delta <- round(data.v3c-data.v2c,2)
  
  res <- cor(data.delta,use="complete.obs")
  res <- round(res, 2)
  res<-res[1:6,-c(1:6)]
  
  library("Hmisc")
  res2 <- rcorr(as.matrix(data.delta))
  res2
  res2.r<-round(t(res2$r[1:6,-c(1:6)]),2)
  res2.p<-round(t(res2$P[1:6,-c(1:6)]),2)
  res2.p[c("Pufap","pufag"),]

  library("PerformanceAnalytics") #pour les plots de correlation avec differents graphiques
  chart.Correlation(data.delta[,1:6], histogram=TRUE, pch=19)
  
  #isovaleric + pufap + pufag
  data.delta.isoval <- data.delta[,c("isovaleric_acid","Pufap","pufag")]
  chart.Correlation(data.delta.isoval, histogram=TRUE, pch=19)
  
  library("ggpubr") #pour les scatter plots
 
  #pufap all (more interesting)
  ggscatter(data.delta2, x = "isovaleric_acid", y = "Pufap", 
           add = "reg.line", conf.int = TRUE, 
           cor.coef = TRUE, cor.method = "pearson",
           xlab = "Isovaleric acid (µmol/L)", ylab = "Total PUFA (% calories)")+
   geom_point(aes(colour = resp),shape = 19) +
   theme(legend.title = element_blank(),legend.position = c(0.085, 0.85),
   panel.border = element_rect(colour = "black", fill=NA))
  
  #pufag all (more interesting)
  ggscatter(data.delta2, x = "isovaleric_acid", y = "pufag", 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "pearson",
            xlab = "Isovaleric acid (µmol/L)", ylab = "Total PUFA (g)")+
    geom_point(aes(colour = resp),shape = 19) +
    theme(legend.title = element_blank(),legend.position = c(0.085, 0.85),
    panel.border = element_rect(colour = "black", fill=NA))
        
  ###########################################################
  
  #resp (less interesting - just to compare with the "contrast" result)
  data.delta2 <- cbind(resp=data.v3$resp,  data.delta)
  data.delta.resp <-  subset(data.delta2, resp =="POS")
  res2.resp <- rcorr(as.matrix(data.delta.resp[,-1]))
  res2.resp.r<-round(t(res2.resp$r[1:6,-c(1:6)]),2)
  res2.resp.p<-round(t(res2.resp$P[1:6,-c(1:6)]),2)
  res2.resp.p[c("Pufap","pufag"),]
  
  #hdl resp vs no-resp 
  ggscatter(data.delta2, x = "isovaleric_acid", y = "hdl", 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "pearson",
            xlab = "Isovaleric acid (µmol/L)", ylab = "HDL",color="resp")+
    theme(legend.title = element_blank(),legend.position = c(0.085, 0.85),
          panel.border = element_rect(colour = "black", fill=NA))

  #no-resp (less interesting - just to compare with the "contrast" result)
  data.delta2 <- cbind(resp=data.v3$resp,  data.delta)
  data.delta.noresp <- subset(data.delta2, resp =="NEG")
  res2.noresp <- rcorr(as.matrix(data.delta.noresp[,-1]))
  res2.noresp.r<-round(t(res2.noresp$r[1:6,-c(1:6)]),2)
  res2.noresp.p<-round(t(res2.noresp$P[1:6,-c(1:6)]),2)
  
  #vitb resp vs no-resp 
  ggscatter(data.delta2, x = "isobutyric_acid", y = "hdl", 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "pearson",
            xlab = "Isovaleric acid (µmol/L)", ylab = "vitb",color="resp")+
    theme(legend.title = element_blank(),legend.position = c(0.085, 0.85),
          panel.border = element_rect(colour = "black", fill=NA))
  



  