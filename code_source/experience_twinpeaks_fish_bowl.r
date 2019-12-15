### twinpeaks####
n<-1200

x<-runif(n, -1, 1)
y<-runif(n, -1, 1)

twin_peaks<-data.frame(x = x,
                       y = y,
                       z = sin(pi * x) * tanh(3 * y))

scatterplot3d(twin_peaks$x, twin_peaks$y, twin_peaks$z, highlight.3d=TRUE,
              col.axis="blue", col.grid="lightblue",
              main="Twin_peaks n=1200", pch=20)
library(MASS)
library(rgl) # nice 3d graphs
library(vegan)

n<-900

x<-runif(n, -1, 1)
y<-runif(n, -1, 1)

twin_peaks<-data.frame(x = x,
                       y = y,
                       z = sin(pi * x) * tanh(3 * y))
#MDS + Isomap
d2 = dist(twin_peaks)
fit <- cmdscale(d2, eig = TRUE)
plot(fit$eig[1:20]/fit$eig[1], type = 'h', 
     main = "Normalized eigenvalues from MDS")
fit2 <- isomap(d2, ndim = 3,k=12)
plot(fit2$eig[1:20]/fit2$eig[1], type = 'h', 
     main = "Normalized eigenvalues from Isomap")
par(mfrow=c(2,2))
plot(fit$points[order(twin_peaks$x), ],ylab='',xlab='',col = rainbow(nrow(twin_peaks)),main="Plot of twinpeaks using method cmdscale", pch = 19)
plot(fit2$points[order(twin_peaks$z), ], ylab='',xlab='',col = rainbow(nrow(twin_peaks)),main="Plot of twinpeaks using method Isomap", pch = 19)
#KPCA try for rbfdot
library(kernlab)
skpca <- kpca(as.matrix(twin_peaks), kernel = "rbfdot", kpar = list(sigma = 0.001))
plot(eig(skpca),  xlim = c(0,10))
plot(pcv(skpca), col = rainbow(nrow(twin_peaks)), pch = 19)
#Sammon Non Linear Mapping
sam =sammon(d2)
#LLE 
library(lle)
calc_k(twin_peaks, m = 2,parallel=T, cpus=3)#13 optimal k
res = lle::lle(twin_peaks,m=2,k=13,reg=2,v=0.95)
par(mfrow=c(2,2))
plot(fit$points[order(twin_peaks$x), ],ylab='',xlab='',col = rainbow(nrow(twin_peaks)),main="Plot of twinpeaks using method cmdscale", pch = 19)
plot(fit2$points[order(twin_peaks$z), ], ylab='',xlab='',col = rainbow(nrow(twin_peaks)),main="Plot of twinpeaks using method Isomap", pch = 19)
plot( res$Y[order(twin_peaks$x), ], col = rainbow(nrow(twin_peaks)), pch = 19, main="Plot of twinpeaks using methode LLE", xlab=expression(y[1]), ylab=expression(y[2]) )
plot(sam$points[order(twin_peaks$z), ],ylab='',xlab='', col = rainbow(nrow(twin_peaks)), main="Plot of twinpeaks using methode Sammon",pch = 19)

d2 = dist(twin_peaks)
#RMSE CMDScale
mdsdist = cmdscale(d2, 2)
dist1 = dist(mdsdist)
RMSE_cmd = sqrt(mean((dist(twin_peaks) - dist1)^2))
#RMSE Sammon Non Linear Mapping
sam =sammon(d2)
dist2 = dist(sam$points)
RMSE_sam = sqrt(mean((dist(twin_peaks) - dist2)^2))
#RMSE MDS
fit <- cmdscale(d2, eig = TRUE)
dist3 = dist(fit$points)
RMSE_mds = sqrt(mean((dist(twin_peaks) - dist3)^2))
#RMSE Isomap
fit2 <- isomap(d2, ndim = 3,k=12)
dist4 = dist(fit2$points)
RMSE_isomap = sqrt(mean((dist(twin_peaks) - dist4)^2))

RMSE = data.frame(CMDScale=RMSE_cmd,Sammon_Non_Linear_Mapping=RMSE_sam,MDS=RMSE_mds,Isomap = RMSE_isomap)
View(RMSE)


#fish boll####
library(scatterplot3d)
library(MASS)
library(vegan)
library(lle)
library(FactoMineR)
library(kernlab)
n<- 1000
phi <- stats::runif(n, 0, 2 * pi)
psi <- acos(stats::runif(n, -1, 0.8))

fish_boll<-data.frame(x = cos(phi) * sin(psi),
                      y = sin(phi) * sin(psi),
                      z = cos(psi))




scatterplot3d(fish_boll[order(psi),3:1], 
              col.axis="blue", col.grid="lightblue",
              main="Fish_boll - 3", pch=20,color = rainbow(n))

scatterplot3d(fish_boll[order(psi),], 
              col.axis="blue", col.grid="lightblue",
              main="Fish_boll - 3", pch=20,color = rainbow(n))

scatterplot3d(fish_boll[order(psi),c(1,3,2)], 
              col.axis="blue", col.grid="lightblue",
              main="Fish_boll - 3", pch=20,color = rainbow(n))

plot(fish_boll[order(fish_boll),c(1,3)], main="Fish_boll - 2", pch=20,col = rainbow(n))
plot(fish_boll[order(fish_boll),c(1,2)], main="Fish_boll - 2", pch=20,col = rainbow(n))
# MDS
fit <- cmdscale(dist(fish_boll), k = 2, eig = TRUE)

plot(fit$eig[1:20]/fit$eig[1], type = 'h', 
     main = "Normalized eigenvalues from MDS for Fish_boll")
plot(fit$points[order(psi), ], col = rainbow(n), pch = 19,main="Plot of Fish_boll using methode MDS",xlab=expression(points[1]), ylab=expression(points[2]))

# Isomap
d <- dist(fish_boll)
fit2 <- isomap(d, ndim = 2, k = 13)
plot(fit2$eig[1:20]/fit2$eig[1], type = 'h', 
     main = "Normalized eigenvalues from Isomap for Fish_boll")
plot(fit2$points[order(psi), ], col = rainbow(n), pch = 19,main="Plot of Fish_boll using methode Isomap")
#sammon
fit3 <- sammon(d)
plot(fit3$points[order(psi), ], col = rainbow(n), pch = 19,main="Plot of Fish_boll using methode Sammon")
#### LLE
fit4<-lle( X=fish_boll, m=2, k=12, reg=2, ss=FALSE, id=TRUE, v=0.9 )
plot( fit4$Y[order(psi), ], col = rainbow(n), pch = 19, main="Plot of Fish_boll using methode LLE", xlab=expression(y[1]), ylab=expression(y[2]) )

### PCA lineair
fit5<-prcomp(fish_boll)
plot( fit5$Y[order(psi), ], col = rainbow(n), pch = 19, main="Plot of Fish_boll using methode LLE", xlab=expression(y[1]), ylab=expression(y[2]) )


fit5<-PCA(fish_boll,ncp = Inf)
plot(fit5[["svd"]][["vs"]],type = 'h')
plot( fit5[["ind"]][["coord"]][order(psi), ], col = rainbow(n), pch = 19, main="Plot of Fish_boll using methode Lineair PCA" )
d = dist(fish_boll)
#calcule de erreur
#RMSE LLE
fit4<-lle( X=fish_boll, m=2, k=12, reg=2, ss=FALSE, id=TRUE, v=0.9 )
dist1 = dist(fit4$Y)
RMSE_lle = sqrt(mean((dist(fish_boll) - dist1)^2))
#RMSE Sammon Non Linear Mapping
fit3 <- sammon(d)
dist2 = dist(fit3$points)
RMSE_sam = sqrt(mean((dist(fish_boll) - dist2)^2))
#RMSE MDS
fit <- cmdscale(dist(fish_boll), k = 2, eig = TRUE)
dist3 = dist(fit$points)
RMSE_mds = sqrt(mean((dist(fish_boll) - dist3)^2))
#RMSE Isomap
fit2 <- isomap(d, ndim = 2, k = 13)
dist4 = dist(fit2$points)
RMSE_isomap = sqrt(mean((dist(fish_boll) - dist4)^2))

RMSE = data.frame(LLE = RMSE_lle,Sammon_Non_Linear_Mapping=RMSE_sam,MDS=RMSE_mds,Isomap = RMSE_isomap)
View(RMSE)
#plot(MD[["ind"]][["coord"]],col=as.integer(cercles[,3]),xlab="1st Principal Component",ylab="2nd Principal Component")
###########""" Kernel ACP default
kpc <- kpca(~.,data=as.data.frame(fish_boll))
pcv(kpc)
plot(rotated(kpc)[order(psi),1:2],col=rainbow(n),
     xlab="1st Principal Component",ylab="2nd Principal Component",main="Plot of Fish_boll using default methode Kernel PCA")
###### Polynomial kernel function
kpc2 <- kpca(~.,data=as.data.frame(fish_boll), kernel = "polydot", kpar=list(degree = 2, scale = 1, offset = 1))
pcv(kpc2)
plot(rotated(kpc2)[order(psi),1:2],col=rainbow(n),
     xlab="1st Principal Component",ylab="2nd Principal Component",main="Plot of Fish_boll using methode Polinomiale 2Â° deg Kernel PCA - ")  

