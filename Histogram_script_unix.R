#Usage: Rscript Histogram_script.R CPM_SLE1_3_MIC_trim.txt conditions.csv
#conditions.txt must have the groups ordered alphabetically



#library(edgeR)
#library(reshape)
library(reshape)
library(plyr)
library(ggplot2)
library(data.table)
library(plotrix)
library(RColorBrewer)




args <- commandArgs(TRUE)
file_name <- args[1]
conds_file <- read.csv(args[2], header=FALSE)
print(file_name)
print (conds_file)
file_prefix = strsplit(file_name,"[.]")


data <- fread(file_name, header=TRUE, drop=2)
data <- as.data.frame(data)
head(data)
names(data)

rownames(data) <- unlist(data[,1])
data[,1] <- NULL


data <- data[, order(names(data))]



df <- data + 1
df <- log2(df)


df <- df[rowSums(df) != 0, ]
df <- as.data.frame(df)
print(colnames(df))



conds_file <- conds_file[order(conds_file$V1), ]
conds_numeric <- as.factor(as.numeric(conds_file$V2))
conds_levels <- levels(conds_numeric)
conds_numeric <- as.data.frame(conds_numeric)


conds <- as.factor(conds_file$V2)
conds <- as.data.frame(conds)
print(conds)


conds["id"] <- NA
for (ii in 1:nrow(conds)){
    conds[ii,"id"] <- ii
}
print(conds)



conditions <- merge(conds,conds_numeric, by=0, sort=TRUE)
print(conditions)

conditions_ordered <- conditions[order(conditions$id), ]
print(conditions_ordered)

conditions_ordered$Row.names <- NULL
conditions_ordered$id <- NULL
print(conditions_ordered)



for(idx in 1:length(conds_levels)){

    assign(paste0("g",idx),colnames(df)[which(conditions_ordered$conds_numeric == conds_levels[idx])])

}



conditions_ordered["sample"] <- NA
for (idxx in 1:length(conds_levels)){

conditions_ordered$sample[which(conditions_ordered$conds_numeric == idxx)] <- unlist(get(paste0("g",idxx)))

}
conditions_ordered$conds_numeric <- NULL
print(conditions_ordered)







df_melt <- melt(df, variable.factor=TRUE, variable.name="variable", value.name="value")
df_melt <- as.data.frame(df_melt)
bins_pt5 <- table(cut(df_melt$value, breaks= seq(0, 16.5, by=0.5)), df_melt$variable)
bins_pt25 <- table(cut(df_melt$value, breaks= seq(0, 16.5, by=0.25)), df_melt$variable)
bins_pt1 <- table(cut(df_melt$value, breaks= seq(0, 16.5, by=0.1)), df_melt$variable)
write.csv(bins_pt5, file = paste0(file_prefix[[1]][1],"_bincounts_0.5.csv"))
write.csv(bins_pt25, file = paste0(file_prefix[[1]][1],"_bincounts_0.25.csv"))
write.csv(bins_pt1, file = paste0(file_prefix[[1]][1], "_bincounts_0.1.csv"))
# only want one histogram as output 12/1/2022





re_bins_pt5 <- data.frame(matrix(NA_integer_, nrow=nrow(bins_pt5), ncol=nrow(conditions_ordered)))


for (i in 1:nrow(conditions_ordered)){
    re_bins_pt5[,i] <- rev(bins_pt5[ ,i])
    re_bins_pt5[,i] <- cumsum(re_bins_pt5[ ,i])
    re_bins_pt5[,i] <- rev(re_bins_pt5[ ,i])

}

colnames(re_bins_pt5) <- colnames(bins_pt5)
rownames(re_bins_pt5) <- rownames(bins_pt5)


#rename rownames for incremental re_bins
for (j in 1:nrow(re_bins_pt5)){
    tmp = strsplit(rownames(bins_pt5)[j],',')
    rownames(re_bins_pt5)[j] <- paste(tmp[[1]][1],'<')
}

write.csv(re_bins_pt5, file= paste0(file_prefix[[1]][1],"_incremental_bincounts[0.5].csv"))

print(re_bins_pt5)

color_dev = c("#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46","#008941", "#006FA6", "#A30059", "#7A4900", "#0000A6", "#63FFAC", "#B79762","#004D43", "#8FB0FF", "#997D87","#5A0007", "#809693", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80","#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100","#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F","#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#001E09","#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66","#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C","#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81","#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00","#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700","#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329","#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C","#83AB58", "#001C1E", "#D1F7CE", "#004B28", "#C8D0F6", "#A3A489", "#806C66", "#222800","#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59", "#8ADBB4", "#1E0200", "#5B4E51","#C895C5", "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94", "#7ED379", "#012C58","#7A7BFF", "#D68E01", "#353339", "#78AFA1", "#FEB2C6", "#75797C", "#837393", "#943A4D","#B5F4FF", "#D2DCD5", "#9556BD", "#6A714A", "#001325", "#02525F", "#0AA3F7", "#E98176","#DBD5DD", "#5EBCD1", "#3D4F44", "#7E6405", "#02684E", "#962B75", "#8D8546", "#9695C5","#E773CE", "#D86A78", "#3E89BE", "#CA834E", "#518A87", "#5B113C", "#55813B", "#E704C4","#00005F", "#A97399", "#4B8160", "#59738A", "#FF5DA7", "#F7C9BF", "#643127", "#513A01","#6B94AA", "#51A058", "#A45B02", "#1D1702", "#E20027", "#E7AB63", "#4C6001", "#9C6966","#64547B", "#97979E", "#006A66", "#391406", "#F4D749", "#0045D2", "#006C31", "#DDB6D0","#7C6571", "#9FB2A4", "#00D891", "#15A08A", "#BC65E9", "#FFFFFE", "#C6DC99", "#203B3C","#671190", "#6B3A64", "#F5E1FF", "#FFA0F2", "#CCAA35", "#374527", "#8BB400", "#797868","#C6005A", "#3B000A", "#C86240", "#29607C", "#402334", "#7D5A44", "#CCB87C", "#B88183","#AA5199", "#B5D6C3", "#A38469", "#9F94F0", "#A74571", "#B894A6", "#71BB8C", "#00B433","#789EC9", "#6D80BA", "#953F00", "#5EFF03", "#E4FFFC", "#1BE177", "#BCB1E5", "#76912F","#003109", "#0060CD", "#D20096", "#895563", "#29201D", "#5B3213", "#A76F42", "#89412E","#1A3A2A", "#494B5A", "#A88C85", "#F4ABAA", "#A3F3AB", "#00C6C8", "#EA8B66", "#958A9F","#BDC9D2", "#9FA064", "#BE4700", "#658188", "#83A485", "#453C23", "#47675D", "#3A3F00","#061203", "#DFFB71", "#868E7E", "#98D058", "#6C8F7D", "#D7BFC2", "#3C3E6E", "#D83D66","#2F5D9B", "#6C5E46", "#D25B88", "#5B656C", "#00B57F", "#545C46", "#866097", "#365D25","#252F99", "#00CCFF", "#674E60", "#FC009C","#FFDBE5", "#92896B")

#Colored by samples
jpeg(paste0(file_prefix[[1]][1],"_incremental_freqplot[bin=0.5].jpg"), height=5, width=6, units="in", res=400)

mypalette = color_dev[1:ncol(re_bins_pt5)]
matplot(y=re_bins_pt5, xlab="Bins(binsize=0.5)", ylab="Incremental Gene Counts", type='l', lty=1, cex=0.3, col=mypalette, xaxt="n")
axis(1,at=1:nrow(re_bins_pt5), labels=NA)
staxlab(1, at=1:nrow(re_bins_pt5), labels=rownames(re_bins_pt5), nlines=1, top.line=0.1, line.spacing=0.1, srt=90, ticklen=0.05, cex=0.6, adj=1)
legend(x="topright", y="topright", legend=colnames(re_bins_pt5), lty=1, lwd=2, ncol=3, cex=0.5, col=mypalette)
dev.off()


#Colored by groups of samples

jpeg(paste0(file_prefix[[1]][1],"_grouped_incremental_freqplot[bin=0.5].jpg"), height=5, width=6, units="in", res=400)


matplot(y=re_bins_pt5, xlab="Bins(binsize=0.5)", ylab="Incremental Gene Counts", type='l', lty=1, cex=0.3, col=mypalette[conditions_ordered$conds], xaxt="n")
axis(1,at=1:nrow(re_bins_pt5), labels=NA)
staxlab(1, at=1:nrow(re_bins_pt5), labels=rownames(re_bins_pt5), nlines=1, top.line=0.1, line.spacing=0.1, srt=90, ticklen=0.05, cex=0.6, adj=1)
legend(x="topright", y="topright", legend=levels(conditions_ordered$conds), lty=1, lwd=2, ncol=3, cex=0.5, col=mypalette[1:length(levels(conditions_ordered$conds))])
dev.off()





re_bins_pt25 <- data.frame(matrix(NA_integer_, nrow=nrow(bins_pt25), ncol=nrow(conditions_ordered)))

for (i in 1:nrow(conditions_ordered)){
    re_bins_pt25[,i] <- rev(bins_pt25[,i])
    re_bins_pt25[,i] <- cumsum(re_bins_pt25[,i])
    re_bins_pt25[,i] <- rev(re_bins_pt25[,i])

}

colnames(re_bins_pt25) <- colnames(bins_pt25)
rownames(re_bins_pt25) <- rownames(bins_pt25)



#rename rownames for incremental re_bins
for (k in 1:nrow(re_bins_pt25)){
    tmp1 = strsplit(rownames(re_bins_pt25)[k],',')
    rownames(re_bins_pt25)[k] <- paste(tmp1[[1]][1],"<")
}

write.csv(re_bins_pt25,file=paste0(file_prefix[[1]][1],"_incremental_bincounts[0.25].csv"))


#Colored by samples
jpeg(paste0(file_prefix[[1]][1],"_incremental_freqplot[bin=0.25].jpg"), height=5, width=6, units="in", res=400)

mypalette = color_dev[1:ncol(re_bins_pt25)]
matplot(y=re_bins_pt25, xlab="Bins(binsize=0.25)", ylab="Incremental Gene Counts", type='l', lty=1, cex=0.3, col=mypalette, xaxt="n")
axis(1,at= c(1,5,10,15,20,25,30,35,40,45,50,55,60,63,nrow(re_bins_pt25)), labels=NA)
staxlab(1, at=c(1,5,10,15,20,25,30,35,40,45,50,55,60,63,nrow(re_bins_pt25)), labels=c(rownames(re_bins_pt25)[1],rownames(re_bins_pt25)[5],rownames(re_bins_pt25)[10],rownames(re_bins_pt25)[15],rownames(re_bins_pt25)[20],rownames(re_bins_pt25)[25],rownames(re_bins_pt25)[30],rownames(re_bins_pt25)[35],rownames(re_bins_pt25)[40],rownames(re_bins_pt25)[45],rownames(re_bins_pt25)[50],rownames(re_bins_pt25)[55],rownames(re_bins_pt25)[60],rownames(re_bins_pt25)[63],rownames(re_bins_pt25)[nrow(re_bins_pt25)]), nlines=1, top.line=0.1, line.spacing=0.1, srt=90, ticklen=0.05, cex=0.6, adj=1)
legend(x="topright", y="topright", legend=colnames(re_bins_pt25), lty=1, lwd=2, ncol=3, cex=0.5, col=mypalette)
dev.off()



#Colored by groups of samples

jpeg(paste0(file_prefix[[1]][1],"_grouped_incremental_freqplot[bin=0.25].jpg"), height=5, width=6, units="in", res=400)
#mypalette = rainbow(ncol(re_bins_pt25))
matplot(y=re_bins_pt25, xlab="Bins(binsize=0.25)", ylab="Incremental Gene Counts", type='l', lty=1, cex=0.3, col=mypalette[conditions_ordered$conds], xaxt="n")
axis(1,at= c(1,5,10,15,20,25,30,35,40,45,50,55,60,63,nrow(re_bins_pt25)), labels=NA)
staxlab(1, at=c(1,5,10,15,20,25,30,35,40,45,50,55,60,63,nrow(re_bins_pt25)), labels=c(rownames(re_bins_pt25)[1],rownames(re_bins_pt25)[5],rownames(re_bins_pt25)[10],rownames(re_bins_pt25)[15],rownames(re_bins_pt25)[20],rownames(re_bins_pt25)[25],rownames(re_bins_pt25)[30],rownames(re_bins_pt25)[35],rownames(re_bins_pt25)[40],rownames(re_bins_pt25)[45],rownames(re_bins_pt25)[50],rownames(re_bins_pt25)[55],rownames(re_bins_pt25)[60],rownames(re_bins_pt25)[63],rownames(re_bins_pt25)[nrow(re_bins_pt25)]), nlines=1, top.line=0.1, line.spacing=0.1, srt=90, ticklen=0.05, cex=0.6, adj=1)
legend(x="topright", y="topright", legend=levels(conditions_ordered$conds), lty=1, lwd=2, ncol=3, cex=0.5, col=mypalette[1:length(levels(conditions_ordered$conds))])
dev.off()











re_bins_pt1 <- data.frame(matrix(NA_integer_, nrow=nrow(bins_pt1), ncol=nrow(conditions_ordered)))

for (i in 1:nrow(conditions_ordered)){
    re_bins_pt1[,i] <- rev(bins_pt1[,i])
    re_bins_pt1[,i] <- cumsum(re_bins_pt1[,i])
    re_bins_pt1[,i] <- rev(re_bins_pt1[,i])

}

colnames(re_bins_pt1) <- colnames(bins_pt1)
rownames(re_bins_pt1) <- rownames(bins_pt1)





#rename rownames for incremental re_bins
for (l in 1:nrow(re_bins_pt1)){
    tmp2 = strsplit(rownames(re_bins_pt1)[l],',')
    rownames(re_bins_pt1)[l] <- paste(tmp2[[1]][1],"<")
}

write.csv(re_bins_pt1, file=paste0(file_prefix[[1]][1],"_incremental_bincounts[0.1].csv"))


#Colored by samples
jpeg(paste0(file_prefix[[1]][1],"_incremental_freqplot[bin=0.1].jpg"), height=5, width=6, units="in", res=400)

mypalette = color_dev[1:ncol(re_bins_pt1)]
matplot(y=re_bins_pt1, xlab="Bins(binsize=0.1)", ylab="Incremental Gene Counts", type='l', lty=1, cex=0.3, col=mypalette, xaxt="n")
axis(1,at= c(1,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,nrow(re_bins_pt1)), labels=NA)
staxlab(1, at=c(1,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,nrow(re_bins_pt1)), labels=c(rownames(re_bins_pt1)[1],rownames(re_bins_pt1)[10],rownames(re_bins_pt1)[20],rownames(re_bins_pt1)[30],rownames(re_bins_pt1)[40],rownames(re_bins_pt1)[50],rownames(re_bins_pt1)[60],rownames(re_bins_pt1)[70],rownames(re_bins_pt1)[80],rownames(re_bins_pt1)[90],rownames(re_bins_pt1)[100],rownames(re_bins_pt1)[110],rownames(re_bins_pt1)[120],rownames(re_bins_pt1)[130],rownames(re_bins_pt1)[140], rownames(re_bins_pt1)[150], rownames(re_bins_pt1)[160],rownames(re_bins_pt1)[nrow(re_bins_pt1)]), nlines=2, top.line=0.1, line.spacing=0.3, srt=90, ticklen=0.03, cex=0.6, adj=1)
legend(x="topright", y="topright", legend=colnames(re_bins_pt1), lty=1, lwd=2, ncol=3, cex=0.5, col=mypalette)
dev.off()


#Colored by groups of samples
jpeg(paste0(file_prefix[[1]][1],"_grouped_incremental_freqplot[bin=0.1].jpg"), height=5, width=6, units="in", res=400)
#mypalette = rainbow(ncol(re_bins_pt1))
matplot(y=re_bins_pt1, xlab="Bins(binsize=0.1)", ylab="Incremental Gene Counts", type='l', lty=1, cex=0.3, col=mypalette[conditions_ordered$conds], xaxt="n")
axis(1,at= c(1,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,nrow(re_bins_pt1)), labels=NA)
staxlab(1, at=c(1,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,nrow(re_bins_pt1)), labels=c(rownames(re_bins_pt1)[1],rownames(re_bins_pt1)[10],rownames(re_bins_pt1)[20],rownames(re_bins_pt1)[30],rownames(re_bins_pt1)[40],rownames(re_bins_pt1)[50],rownames(re_bins_pt1)[60],rownames(re_bins_pt1)[70],rownames(re_bins_pt1)[80],rownames(re_bins_pt1)[90],rownames(re_bins_pt1)[100],rownames(re_bins_pt1)[110],rownames(re_bins_pt1)[120],rownames(re_bins_pt1)[130],rownames(re_bins_pt1)[140], rownames(re_bins_pt1)[150], rownames(re_bins_pt1)[160],rownames(re_bins_pt1)[nrow(re_bins_pt1)]), nlines=2, top.line=0.1, line.spacing=0.3, srt=90, ticklen=0.03, cex=0.6, adj=1)
legend(x="topright", y="topright", legend=levels(conditions_ordered$conds), lty=1, lwd=2, ncol=3, cex=0.5, col=mypalette[1:length(levels(conditions_ordered$conds))])
dev.off()











df_melt["category"] <- NA
df_melt$category <- conditions_ordered$conds[match(df_melt$variable, conditions_ordered$sample)]
head(df_melt)







#Frequency Polygons with better labelling and staxlab

jpeg(paste0(file_prefix[[1]][1],"_frequencyPolygon[bin=0.5].jpg"), height=5, width=6, units="in", res=400)

palette = color_dev[1: ncol(bins_pt5)]
matplot(y=bins_pt5, xlab="Bins(binsize=0.5)", ylab="Gene Counts", type='l', lty=1, cex=0.3, col=palette, xaxt="n")
axis(1,at= c(1,5,10,15,20,25,30,nrow(bins_pt5)), labels=NA)
staxlab(1, at=1:nrow(bins_pt5), labels=rownames(bins_pt5), nlines=1, top.line=0.1, line.spacing=0.1, srt=90, ticklen=0.05, cex=0.6, adj=1)
legend(x="topright", y="topright", legend=colnames(bins_pt5), lty=1, lwd=2, ncol=3, cex=0.4, col=palette)
dev.off()


jpeg(paste0(file_prefix[[1]][1],"_grouped_frequencyPolygon[bin=0.5].jpg"), height=5, width=6, units="in", res=400)

matplot(y=bins_pt5, xlab="Bins(binsize=0.5)", ylab="Gene Counts", type='l', lty=1, cex=0.3, col=mypalette[conditions_ordered$conds], xaxt="n")
axis(1,at=1:nrow(bins_pt5), labels=NA)
staxlab(1, at=1:nrow(bins_pt5), labels=rownames(bins_pt5), nlines=1, top.line=0.1, line.spacing=0.1, srt=90, ticklen=0.05, cex=0.6, adj=1)
legend(x="topright", y="topright", legend=levels(conditions_ordered$conds), lty=1, lwd=2, ncol=3, cex=0.5, col=mypalette[1:length(levels(conditions_ordered$conds))])
dev.off()




















jpeg(paste0(file_prefix[[1]][1],"_frequencyPolygon[bin=0.25].jpg"), height=5, width=6, units="in", res=400)

palette = color_dev[1:ncol(bins_pt25)]
matplot(y=bins_pt25, xlab="Bins(binsize=0.25)", ylab="Gene Counts", type='l', lty=1, cex=0.3, col=palette, xaxt="n")
axis(1,at= c(1,5,10,15,20,25,30,35,40,45,50,55,60,63,nrow(bins_pt25)), labels=NA)
staxlab(1, at=c(1,5,10,15,20,25,30,35,40,45,50,55,60,63,nrow(bins_pt25)), labels=c(rownames(bins_pt25)[1],rownames(bins_pt25)[5],rownames(bins_pt25)[10],rownames(bins_pt25)[15],rownames(bins_pt25)[20],rownames(bins_pt25)[25],rownames(bins_pt25)[30],rownames(bins_pt25)[35],rownames(bins_pt25)[40],rownames(bins_pt25)[45],rownames(bins_pt25)[50],rownames(bins_pt25)[55],rownames(bins_pt25)[60],rownames(bins_pt25)[63],rownames(bins_pt25)[nrow(bins_pt25)]), nlines=1, top.line=0.1, line.spacing=0.1, srt=90, ticklen=0.05, cex=0.6, adj=1)
legend(x="topright", y="topright", legend=colnames(bins_pt25), lty=1, lwd=2, ncol=3, cex=0.5, col=palette)
dev.off()






jpeg(paste0(file_prefix[[1]][1],"_frequencyPolygon[bin=0.1].jpg"), height=5, width=6, units="in", res=400)

palette = color_dev[1:ncol(bins_pt1)]
matplot(y=bins_pt1, xlab="Bins(binsize=0.1)", ylab="Gene Counts", type='l', lty=1, cex=0.4, col=palette, xaxt="n")
axis(1,at= c(1,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,nrow(bins_pt1)), labels=NA)
staxlab(1, at=c(1,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,nrow(bins_pt1)), labels=c(rownames(bins_pt1)[1],rownames(bins_pt1)[10],rownames(bins_pt1)[20],rownames(bins_pt1)[30],rownames(bins_pt1)[40],rownames(bins_pt1)[50],rownames(bins_pt1)[60],rownames(bins_pt1)[70],rownames(bins_pt1)[80],rownames(bins_pt1)[90],rownames(bins_pt1)[100],rownames(bins_pt1)[110],rownames(bins_pt1)[120],rownames(bins_pt1)[130],rownames(bins_pt1)[140], rownames(bins_pt1)[150], rownames(bins_pt1)[160],rownames(bins_pt1)[nrow(bins_pt1)]), nlines=2, top.line=0.1, line.spacing=0.3, srt=90, ticklen=0.05, cex=0.6, adj=1)
legend(x="topright", y="topright", legend=colnames(bins_pt1), lty=1, lwd=2, ncol=3, cex=0.5, col=palette)
dev.off()
