# Oleg Moskvin, www.genepeak.com
HypExpo <- function(df, param.hyp=1, param.expo=1.8, param.diag=0, param.shift=0.17, param.xshift=0, baseName="test", preview.only=TRUE) {              
theCurve <- function(x) param.hyp / (log10(x)-param.xshift)^param.expo - param.diag * (log10(x)-param.xshift) + param.shift
paramChain <- paste(param.hyp, param.expo, param.diag, param.shift, param.xshift, sep="_")
imgName <- paste(baseName, paramChain, "png", sep=".")
df <- cbind(df, passed = as.numeric(log10(df[,2]) >= theCurve(df[,1])))
# vector of the selection results:
selVec <- df[,ncol(df)] # column number may differ depending on the number of selection runs with different parameters
selVec[is.na(selVec)] <- 0 
# dataframe of the selected gene's data:
df.sel <- df[selVec==1,] 
#
nGenes.initial <- nrow(df)
nGenes.selected <- nrow(df.sel)
exprRange <- log10(range(df[,1], na.rm = TRUE))
foldRange <- log10(range(df[,2], na.rm = TRUE))
#
if (preview.only) { 
plot(log10(df[,2]) ~ log10(df[,1]), pch=19, cex=0.5, col="grey", cex.lab=1.5, cex.axis=1.5, xlab="log10(expression)", ylab="log10(Rel. change)", main=paste(nGenes.selected, "out of", nGenes.initial, "genes"), xlim=exprRange, ylim=foldRange)
curve(param.hyp / (x-param.xshift)^param.expo - param.diag * (x-param.xshift) + param.shift, add=TRUE, col="green") }
#
if (!(preview.only)) {
png(imgName, width=960, height=480)
par(mfrow=c(1,2))
plot(log10(df[,2]) ~ log10(df[,1]), pch=19, cex=0.5, col="grey", cex.lab=1.5, cex.axis=1.5, xlab="log10(expression)", ylab="log10(Rel. change)", main=c(paste(nGenes.initial, "Genes"), baseName), xlim=exprRange, ylim=foldRange)
curve(param.hyp / (x-param.xshift)^param.expo - param.diag * (x-param.xshift) + param.shift, add=TRUE, col="green")
plot(log10(df.sel[,2]) ~ log10(df.sel[,1]), pch=19, cex=0.5, col="grey", cex.lab=1.5, cex.axis=1.5, xlab="log10(expression)", ylab="log10(Rel. change)", main=c(paste(nGenes.selected, "Genes"), paste("Hyp_Expo_Diag_Shift", paramChain)), xlim=exprRange, ylim=foldRange)
curve(param.hyp / (x-param.xshift)^param.expo - param.diag * (x-param.xshift) + param.shift, add=TRUE, col="green")
dev.off() 
out <- list(df, df.sel, rownames(df.sel))
out
} }
