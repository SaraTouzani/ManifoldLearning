library(scatterplot3d)
library(MASS)
library(vegan)
library(lle)
library(FactoMineR)
library(kernlab)

n <- 2500 # Random position on the parameteric domain.
u <- matrix(runif(2 * n), ncol = 2)

v <- 3 * pi / 2 * (0.1 + 2 * u[, 1])

x <- -(v) 
y <- 20 * u[, 2]
z <- sin(v) 

sin_curv <- cbind(x, y , z)
#plot3d(sin_curv[order(v), ], col = rainbow(n), size = 10)

scatterplot3d(sin_curv[order(v),3:1 ], 
              col.axis="blue", col.grid="lightblue",
              main="S Curv - 3", pch=20,color = rainbow(n))

scatterplot3d(sin_curv[order(v), ], 
              col.axis="blue", col.grid="lightblue",
              main="S Curv - 3", pch=20,color = rainbow(n))

plot(sin_curv[order(v),c(1,3)], main="S Curv - 2", pch=20,col = rainbow(n))
plot(sin_curv[order(v),c(1,2)], main="S Curv - 2", pch=20,col = rainbow(n))
# MDS
fit <- cmdscale(dist(sin_curv), k = 2, eig = TRUE)
layout(matrix(1:6,ncol = 3,byrow = TRUE))
plot(fit$eig[1:20]/fit$eig[1], type = 'h', 
     main = "Normalized eigenvalues from MDS for S Curv")
plot(fit$points[order(v), ], col = rainbow(n), pch = 19,main="Plot of S_curv using methode MDS",xlab=expression(points[1]), ylab=expression(points[2]))

# Isomap
d <- dist(sin_curv)
fit2 <- isomap(d, ndim = 2, k = 13)
plot(fit2$eig[1:20]/fit2$eig[1], type = 'h', 
     main = "Normalized eigenvalues from Isomap for S Curv")
plot(fit2$points[order(v), ], col = rainbow(n), pch = 19,main="Plot of S_curv using methode Isomap")
#########sammon
fit3 <- sammon(d)
plot(fit3$points[order(v), ], col = rainbow(n), pch = 19,main="Plot of S_curv using methode Sammon")
#### LLE
fit4<-lle( X=sin_curv, m=2, k=12, reg=2, ss=FALSE, id=TRUE, v=0.9 )
plot( fit4$Y[order(v), ], col = rainbow(n), pch = 19, main="Plot of S_curv using methode LLE", xlab=expression(y[1]), ylab=expression(y[2]) )

### PCA lineair
fit5<-prcomp(sin_curv)
plot( fit5$Y[order(v), ], col = rainbow(n), pch = 19, main="Plot of S_curv using methode LLE", xlab=expression(y[1]), ylab=expression(y[2]) )


fit5<-PCA(sin_curv,ncp = Inf)
plot(fit5[["svd"]][["vs"]],type = 'h')
plot( fit5[["ind"]][["coord"]][order(v), ], col = rainbow(n), pch = 19, main="Plot of S_curv using methode Lineair PCA" )

#plot(MD[["ind"]][["coord"]],col=as.integer(cercles[,3]),xlab="1st Principal Component",ylab="2nd Principal Component")
###########""" Kernel ACP default
kpc <- kpca(~.,data=as.data.frame(sin_curv))
pcv(kpc)
plot(rotated(kpc)[order(v),1:2],col=rainbow(n),
     xlab="1st Principal Component",ylab="2nd Principal Component",main="Plot of S_curv using default methode Kernel PCA")
###### Polynomial kernel function
kpc2 <- kpca(~.,data=as.data.frame(sin_curv), kernel = "polydot", kpar=list(degree = 2, scale = 1, offset = 1))
pcv(kpc2)
plot(rotated(kpc2)[order(v),1:2],col=rainbow(n),
     xlab="1st Principal Component",ylab="2nd Principal Component",main="Plot of S_curv using methode Polinomiale 2° deg Kernel PCA - ")  

########################"Hyperbolic tangent kernel function

#kpc3 <- kpca(~.,data=as.data.frame(sin_curv), kernel = "tanhdot", kpar=list(scale = 1, offset = 1))
#pcv(kpc3)
#plot(rotated(kpc3)[order(v),1:2],col=rainbow(n),
#     xlab="1st Principal Component",ylab="2nd Principal Component") 

###################### laplacedot: Laplacian kernel function
kpc4 <- kpca(~.,data=as.data.frame(sin_curv), kernel = "laplacedot", kpar=list(sigma = 1))
pcv(kpc4)
plot(rotated(kpc4)[order(v),1:2],col=rainbow(n),
     xlab="1st Principal Component",ylab="2nd Principal Component",main="Plot of S_curv using Laplacian Kernel PCA - ") 

################ besseldot Bessel kernel function
kpc5 <- kpca(~.,data=as.data.frame(sin_curv), kernel = "besseldot", kpar=list(sigma = 1, order = 1, degree = 1))
pcv(kpc5)
plot(rotated(kpc5)[order(v),1:2],col=rainbow(n),
     xlab="1st Principal Component",ylab="2nd Principal Component",main="Plot of S_curv using Bessel Kernel PCA  ") 

############ anovadot ANOVA RBF kernel function
kpc6 <- kpca(~.,data=as.data.frame(sin_curv), kernel = "anovadot", kpar=list(sigma = 1, degree = 1))
pcv(kpc6)
plot(rotated(kpc6)[order(v),1:2],col=rainbow(n),
     xlab="1st Principal Component",ylab="2nd Principal Component",main="Plot of S_curv using Anova RBF Kernel PCA  ") 

########### splinedot pas de paramètre  Spline kernel
kpc7 <- kpca(~.,data=as.data.frame(sin_curv), kernel = "splinedot", kpar=list())
pcv(kpc7)
plot(rotated(kpc7)[order(v),1:2],col=rainbow(n),
     xlab="1st Principal Component",ylab="2nd Principal Component",main="Plot of S_curv using Spline Kernel PCA  ") 
plot(kpc7@eig,type='h') 
scatterplot3d(rotated(kpc7)[order(v),1:3], 
              col.axis="blue", col.grid="lightblue",
              main="S Curv - 3", pch=20,color = rainbow(n))