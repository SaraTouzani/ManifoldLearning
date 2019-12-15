library(scatterplot3d)
library(MASS)
library(vegan)
library(lle)
library(FactoMineR)
library(kernlab)
################

n<-50

temp <- seq(-pi, 0, length = n)



hemisphere<-data.frame(x=c(rep(1, 50) %*% t(cos(temp))),
                       y=(1)*c(cos(temp) %*% t(sin(temp))),
                       z=(1)*c(sin(temp) %*% t(sin(temp))))
layout(matrix(1:6,ncol=3,byrow =TRUE))
scatterplot3d(hemisphere[order(rep(temp,length(temp))),], color = rainbow(n**2),
              col.axis="blue", col.grid="lightblue",
              main="hemisphaere - 3D", pch=20)

scatterplot3d(hemisphere[order(rep(temp,length(temp))),3:1], color = rainbow(n**2),
              col.axis="blue", col.grid="lightblue",
              main="hemisphaere - 3D", pch=20)

scatterplot3d(hemisphere[order(rep(temp,length(temp))),c(1,3,2)], color = rainbow(n**2),
              col.axis="blue", col.grid="lightblue",
              main="hemisphaere - 3D", pch=20)



plot(hemisphere[order(rep(temp,length(temp))),c(1,3)], main="hemisphere - 2D", pch=20,col = rainbow(n**2))
plot(hemisphere[order(rep(temp,length(temp))),c(2,3)], main="hemisphere - 2D", pch=20,col = rainbow(n**2))
plot(hemisphere[order(rep(temp,length(temp))),c(1,2)], main="hemisphere - 2D", pch=20,col = rainbow(n**2))

# MDS
layout(matrix(1:6, ncol=3, byrow=TRUE))
fit <- cmdscale(dist(hemisphere), k = 2, eig = TRUE)

plot(fit$eig[1:20]/fit$eig[1], type = 'h', 
     main = "Normalized eigenvalues from MDS for hemisphare")
plot(fit$points[order(rep(temp,length(temp))),], col = rainbow(n**2), pch = 19,main="Plot of hemisphere using methode MDS",xlab=expression(points[1]), ylab=expression(points[2]))

# Isomap
d <- dist(hemisphere)
fit2 <- isomap(d, ndim = 2,epsilon = 0.2)
plot(fit2$eig[1:20]/fit2$eig[1], type = 'h', 
     main = "Normalized eigenvalues from Isomap for hemisphere")
plot(fit2$points[order(rep(temp,length(temp))),], col = rainbow(n**2), pch = 19,main="Plot of hemisphere using methode Isomap")
#########sammon
#d <- dist(hemisphere)
#y'a des poins qui se répete à apartir de 2452 observation
d1 <- dist(hemisphere[-c(2452:2500),])
s<-length(temp)-1
fit3 <- sammon(d1,y = jitter(cmdscale(d, 2))) 
#plot(fit3$points[order(rep(temp,length(temp))),], col = rainbow(n**2), pch = 19,main="Plot of hemisphere using methode Sammon")
plot(fit3$points[order(rep(temp,s)),], col = rainbow(s**2), pch = 19,main="Plot of hemisphere using methode Sammon")

#### LLE
fit4<-lle( X=hemisphere[-c(2452:2500),], m=2, k=40, reg=1, ss=FALSE, id=TRUE, v=0.9 )
plot( fit4$Y, col = rainbow(n**2), pch = 19, main="Plot of hemisphere using methode LLE", xlab=expression(y[1]), ylab=expression(y[2]) )
#M Dim 
### PCA lineair
#fit5<-prcomp(hemisphere)
#plot( fit5$Y[order(rep(temp,length(temp))),], col = rainbow(n**2), pch = 19, main="Plot of hemisphere using methode LLE", xlab=expression(y[1]), ylab=expression(y[2]) )


fit5<-PCA(hemisphere,ncp = Inf)
plot(fit5[["svd"]][["vs"]],type = 'h', main="eigen values for PCA lineair")
plot( fit5[["ind"]][["coord"]][order(rep(temp,length(temp))),], col = rainbow(n**2), pch = 19, main="Plot of hemisphere using methode Lineair PCA" )

#plot(MD[["ind"]][["coord"]],col=as.integer(cercles[,3]),xlab="1st Principal Component",ylab="2nd Principal Component")
###########""" Kernel ACP default
kpc <- kpca(~.,data=hemisphere)
## quand hemispere est de grande dimention kernel function nous renvoie une matrice infinie
# n<-50
# temp <- seq(-pi, 0, length = n)
# 
# hemisphere<-data.frame(x=c(rep(1, 50) %*% t(cos(temp))),
#                        y=(1e20)*c(cos(temp) %*% t(sin(temp))),
#                        z=(1e20)*c(sin(temp) %*% t(sin(temp))))
##mais utiliser : 
# n<-50
# temp <- seq(-pi, 0, length = n)
# 
# hemisphere<-data.frame(x=c(rep(1, 50) %*% t(cos(temp))),
#                        y=c(cos(temp) %*% t(sin(temp))),
#                        z=c(sin(temp) %*% t(sin(temp))))
pcv(kpc)
plot(rotated(kpc)[order(rep(temp,length(temp))),], col = rainbow(n**2),
     xlab="1st Principal Component",ylab="2nd Principal Component",main="Plot of hemisphere using default methode Kernel PCA")
###### Polynomial kernel function
kpc2 <- kpca(~.,data=(hemisphere), kernel = "polydot", kpar=list(degree = 2, scale = 1, offset = 1))
## le PCA utilisant un kernel plonomiale est robuste pour les petite dimention et pour les grande
pcv(kpc2)
plot(rotated(kpc2)[order(rep(temp,length(temp))),], col = rainbow(n**2),
     xlab="1st Principal Component",ylab="2nd Principal Component",main="Plot of hemisphere using methode Polinomiale 2° deg Kernel PCA - ")  



kpc22 <- kpca(~.,data=(hemisphere), kernel = "polydot", kpar=list(degree = 5, scale = 1, offset = 1))
## le PCA utilisant un kernel plonomiale est robuste pour les petite dimention et pour les grande
pcv(kpc22)
plot(rotated(kpc22)[order(rep(temp,length(temp))),], col = rainbow(n**2),
     xlab="1st Principal Component",ylab="2nd Principal Component",main="Plot of hemisphere using methode Polinomiale 5° deg Kernel PCA - ")  

kpc26 <- kpca(~.,data=(hemisphere), kernel = "polydot", kpar=list(degree = 6, scale = 1, offset = 1))
## le PCA utilisant un kernel plonomiale est robuste pour les petite dimention et pour les grande
pcv(kpc26)
plot(rotated(kpc26)[order(rep(temp,length(temp))),], col = rainbow(n**2),
     xlab="1st Principal Component",ylab="2nd Principal Component",main="Plot of hemisphere using methode Polinomiale 6° deg Kernel PCA - ")  

kpc27 <- kpca(~.,data=(hemisphere), kernel = "polydot", kpar=list(degree = 7, scale = 1, offset = 1))
## le PCA utilisant un kernel plonomiale est robuste pour les petite dimention et pour les grande
pcv(kpc27)
plot(rotated(kpc27)[order(rep(temp,length(temp))),], col = rainbow(n**2),
     xlab="1st Principal Component",ylab="2nd Principal Component",main="Plot of hemisphere using methode Polinomiale 7° deg Kernel PCA - ")  

########################"Hyperbolic tangent kernel function

#kpc3 <- kpca(~.,data=as.data.frame(hemisphere), kernel = "tanhdot", kpar=list(scale = 1, offset = 1))
#pcv(kpc3)
#plot(rotated(kpc3)[order(v),1:2],col=rainbow(n),
#     xlab="1st Principal Component",ylab="2nd Principal Component") 

###################### laplacedot: Laplacian kernel function
kpc4 <- kpca(~.,data=(hemisphere), kernel = "laplacedot", kpar=list(sigma = 1))
# sensible au grande dimentions 
pcv(kpc4)
plot(rotated(kpc4)[order(rep(temp,length(temp))),], col = rainbow(n**2),
     xlab="1st Principal Component",ylab="2nd Principal Component",main="Plot of hemisphere using Laplacian Kernel PCA - ") 

################ besseldot Bessel kernel function
kpc5 <- kpca(~.,data=(hemisphere), kernel = "besseldot", kpar=list(sigma = 1, order = 1, degree = 1))
pcv(kpc5)
plot(rotated(kpc5)[order(rep(temp,length(temp))),], col = rainbow(n**2),
     xlab="1st Principal Component",ylab="2nd Principal Component",main="Plot of hemisphere using Bessel Kernel PCA  ") 

############ anovadot ANOVA RBF kernel function
kpc61 <- kpca(~.,data=(hemisphere), kernel = "anovadot", kpar=list(sigma = 1, degree = 1))
pcv(kpc61)
plot(rotated(kpc61)[order(rep(temp,length(temp))),], col = rainbow(n**2),
     xlab="1st Principal Component",ylab="2nd Principal Component",main="Plot of hemisphere using Anova RBF Kernel PCA sigma=1, d°1  ") 

kpc62 <- kpca(~.,data=(hemisphere), kernel = "anovadot", kpar=list(sigma = 0.5, degree = 3))
pcv(kpc62)
plot(rotated(kpc62)[order(rep(temp,length(temp))),], col = rainbow(n**2),
     xlab="1st Principal Component",ylab="2nd Principal Component",main="Plot of hemisphere using Anova RBF Kernel PCA sigma=0.5, d°3  ") 

kpc64 <- kpca(~.,data=(hemisphere), kernel = "anovadot", kpar=list(sigma = 0.5, degree = 10))
pcv(kpc64)
plot(rotated(kpc64)[order(rep(temp,length(temp))),], col = rainbow(n**2),
     xlab="1st Principal Component",ylab="2nd Principal Component",main="Plot of hemisphere using Anova RBF Kernel PCA sigma=0.5, d°10  ") 

kpc63 <- kpca(~.,data=(hemisphere), kernel = "anovadot", kpar=list(sigma = 0.25, degree = 5))
pcv(kpc63)
plot(rotated(kpc63)[order(rep(temp,length(temp))),], col = rainbow(n**2),
     xlab="1st Principal Component",ylab="2nd Principal Component",main="Plot of hemisphere using Anova RBF Kernel PCA sigma=0.25, d°5  ") 

kpc65 <- kpca(~.,data=(hemisphere), kernel = "anovadot", kpar=list(sigma = 0.25, degree = 10))
pcv(kpc65)
plot(rotated(kpc65)[order(rep(temp,length(temp))),], col = rainbow(n**2),
     xlab="1st Principal Component",ylab="2nd Principal Component",main="Plot of hemisphere using Anova RBF Kernel PCA sigma=0.25, d°10  ") 

kpc66 <- kpca(~.,data=(hemisphere), kernel = "anovadot", kpar=list(sigma = 0.25, degree = 15))
pcv(kpc66)
plot(rotated(kpc66)[order(rep(temp,length(temp))),], col = rainbow(n**2),
     xlab="1st Principal Component",ylab="2nd Principal Component",main="Plot of hemisphere using Anova RBF Kernel PCA sigma=0.25, d°15  ") 


########### splinedot pas de paramètre  Spline kernel
kpc7 <- kpca(~.,data=(hemisphere), kernel = "splinedot", kpar=list())
pcv(kpc7)
plot(rotated(kpc7)[order(rep(temp,length(temp))),], col = rainbow(n**2),
     xlab="1st Principal Component",ylab="2nd Principal Component",main="Plot of hemisphere using Spline Kernel PCA  ") 
##Erreur ####
Erreur<-data.frame(methode="MDS",
                   e=sqrt(mean((as.dist(d) - dist(fit$points))^2))
)

Erreur<-rbind(Erreur,
              data.frame(methode="ISOMAP",
                         e=sqrt(mean((as.dist(d) - dist(fit2$points))^2))))
Erreur<-rbind(Erreur,
              data.frame(methode="SAMMON",
                         e=sqrt(mean((as.dist(d1) - dist(fit3$points))^2))))
Erreur<-rbind(Erreur,
              data.frame(methode="LLE",
                         e=sqrt(mean((as.dist(d) - dist(fit4$Y))^2))))
Erreur<-rbind(Erreur,
              data.frame(methode="Lin-PCA",
                         e=sqrt(mean((as.dist(d) - dist( fit5[["ind"]][["coord"]]))^2))))
Erreur<-rbind(Erreur,
              data.frame(methode="Def-K-PCA",
                         e=sqrt(mean((as.dist(d) - dist(rotated(kpc)))^2))))
Erreur<-rbind(Erreur,
              data.frame(methode="Optimal_Pol_KPCA",
                         e=sqrt(mean((as.dist(d) - dist(rotated(kpc2)))^2))))
Erreur<-rbind(Erreur,
              data.frame(methode="Laplacian-KPCA",
                         e=sqrt(mean((as.dist(d) - dist(rotated(kpc4)))^2))))
Erreur<-rbind(Erreur,
              data.frame(methode="Bessel-KPCA",
                         e=sqrt(mean((as.dist(d) - dist(rotated(kpc5)))^2))))
Erreur<-rbind(Erreur,
              data.frame(methode="Spline-KPCA",
                         e=sqrt(mean((as.dist(d) - dist(rotated(kpc7)))^2))))

Erreur<-rbind(Erreur,
              data.frame(methode="Anova-RBF-sig1-d1-KPCA",
                         e=sqrt(mean((as.dist(d) - dist(rotated(kpc61)))^2))))
Erreur<-rbind(Erreur,
              data.frame(methode="Optimal-Anova-RBF-sig0.5-d3-KPCA",
                         e=sqrt(mean((as.dist(d) - dist(rotated(kpc62)))^2))))
Erreur<-rbind(Erreur,
              data.frame(methode="Optimal-Anova-RBF-sig0.5-d10-KPCA",
                         e=sqrt(mean((as.dist(d) - dist(rotated(kpc63)))^2))))
Erreur<-rbind(Erreur,
              data.frame(methode="Optimal-Anova-RBF-0.25-d5-KPCA",
                         e=sqrt(mean((as.dist(d) - dist(rotated(kpc64)))^2))))
Erreur<-rbind(Erreur,
              data.frame(methode="Optimal-Anova-RBF-0.25-d10-KPCA",
                         e=sqrt(mean((as.dist(d) - dist(rotated(kpc65)))^2))))
Erreur<-rbind(Erreur,
              data.frame(methode="Optimal-Anova-RBF-0.25-d15-KPCA",
                         e=sqrt(mean((as.dist(d) - dist(rotated(kpc66)))^2))))

####################

Erreur<-data.frame(methode="MDS",
                   e=sqrt(mean((as.dist(d) - dist(fit$points))^2))
)

Erreur<-rbind(Erreur,
              data.frame(methode="ISOMAP",
                         e=sqrt(mean((as.dist(d) - dist(fit2$points))^2))))
Erreur<-rbind(Erreur,
              data.frame(methode="SAMMON",
                         e=sqrt(mean((as.dist(d1) - dist(fit3$points))^2))))
Erreur<-rbind(Erreur,
              data.frame(methode="LLE",
                         e=sqrt(mean((as.dist(d) - dist(fit4$Y))^2))))
Erreur<-rbind(Erreur,
              data.frame(methode="Lin-PCA",
                         e=sqrt(mean((as.dist(d) - dist( fit5[["ind"]][["coord"]]))^2))))
Erreur<-rbind(Erreur,
              data.frame(methode="Def-K-PCA",
                         e=sqrt(mean((as.dist(d) - dist(rotated(kpc)))^2))))
Erreur<-rbind(Erreur,
              data.frame(methode="Optimal_Pol_KPCA",
                         e=sqrt(mean((as.dist(d) - dist(rotated(kpc2)))^2))))
Erreur<-rbind(Erreur,
              data.frame(methode="Laplacian-KPCA",
                         e=sqrt(mean((as.dist(d) - dist(rotated(kpc4)))^2))))
Erreur<-rbind(Erreur,
              data.frame(methode="Bessel-KPCA",
                         e=sqrt(mean((as.dist(d) - dist(rotated(kpc5)))^2))))
Erreur<-rbind(Erreur,
              data.frame(methode="Spline-KPCA",
                         e=sqrt(mean((as.dist(d) - dist(rotated(kpc7)))^2))))

Erreur<-rbind(Erreur,
              data.frame(methode="Anova-RBF-sig1-d1-KPCA",
                         e=sqrt(mean((as.dist(d) - dist(rotated(kpc61)))^2))))
Erreur<-rbind(Erreur,
              data.frame(methode="Optimal-Anova-RBF-sig0.5-d3-KPCA",
                         e=sqrt(mean((as.dist(d) - dist(rotated(kpc62)))^2))))
Erreur<-rbind(Erreur,
              data.frame(methode="Optimal-Anova-RBF-sig0.5-d10-KPCA",
                         e=sqrt(mean((as.dist(d) - dist(rotated(kpc63)))^2))))
Erreur<-rbind(Erreur,
              data.frame(methode="Optimal-Anova-RBF-0.25-d5-KPCA",
                         e=sqrt(mean((as.dist(d) - dist(rotated(kpc64)))^2))))
Erreur<-rbind(Erreur,
              data.frame(methode="Optimal-Anova-RBF-0.25-d10-KPCA",
                         e=sqrt(mean((as.dist(d) - dist(rotated(kpc65)))^2))))
Erreur<-rbind(Erreur,
              data.frame(methode="Optimal-Anova-RBF-0.25-d15-KPCA",
                         e=sqrt(mean((as.dist(d) - dist(rotated(kpc66)))^2))))

###################



