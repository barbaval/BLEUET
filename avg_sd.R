data <- read.csv("data//bleuet_data.csv")
outliers <- c(23,59)#,33)#,33)
data <- subset(data, !nomat %in% outliers)
bleuets <- subset(data, !Visites == "V3" & !tx == "B")

# Average and tsnadard deviation
cols <- c("nomat","Visites","hdlc","ldlc","chol","chol_hdlc")
temp <- bleuets

temp <-as.data.frame(bleuets[,cols])
cols2 <- c("hdlc","ldlc","chol","chol_hdlc")
temp.v2 <- subset(temp, temp$Visites == "V2")
row.names(temp.v2) <- temp.v2$nomat
temp.v2 <- temp.v2[,cols2]
temp.v4 <- subset(temp, temp$Visites == "V4")
row.names(temp.v4) <- temp.v4$nomat
temp.v4 <- temp.v4[,cols2]

delta <- temp.v4 - temp.v2

group <- c(61,28,71,68,106) # Groupe "NR"

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