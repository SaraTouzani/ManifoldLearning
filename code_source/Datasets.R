library(scatterplot3d)

########### s_curv
n <- 2500 # Random position on the parameteric domain.
u <- matrix(runif(2 * n), ncol = 2)

v <- 3 * pi / 2 * (0.1 + 2 * u[, 1])

x <- -(v) 
y <- 20 * u[, 2]
z <- sin(v) 

sin_curv <- cbind(x, y , z)

scatterplot3d(sin_curv, highlight.3d=TRUE, col.axis="blue",
              col.grid="lightblue", main="scatterplot3d - 1", pch=20)


scatterplot3d(sin_curv[,3:1], highlight.3d=TRUE, col.axis="blue",
              col.grid="lightblue", main="scatterplot3d - 1", pch=20)

####################### hemisphere 

temp <- seq(-pi, 0, length = 50)

hemisphere<-data.frame(x=c(rep(1, 50) %*% t(cos(temp))),
                       y=c(cos(temp) %*% t(sin(temp))),
                       z=c(sin(temp) %*% t(sin(temp))))

scatterplot3d(hemisphere$x, hemisphere$y, hemisphere$z, highlight.3d=TRUE,
              col.axis="blue", col.grid="lightblue",
              main="scatterplot3d - 2", pch=20)


################## sphere

temp <- seq(-pi, pi, length = 100)

sphere<-data.frame(x=c(rep(1, 100) %*% t(cos(temp))),
                   y=c(cos(temp) %*% t(sin(temp))),
                   z=c(sin(temp) %*% t(sin(temp))))


scatterplot3d(sphere$z, sphere$y, sphere$x, highlight.3d=TRUE,
              col.axis="blue", col.grid="lightblue",
              main="scatterplot3d - 2", pch=20)



#################### twin peaks
n<-10000

x<-runif(n, -1, 1)
y<-runif(n, -1, 1)

twin_peaks<-data.frame(x = x,
                       y = y,
                       z = sin(pi * x) * tanh(3 * y))

scatterplot3d(twin_peaks$x, twin_peaks$y, twin_peaks$z, highlight.3d=TRUE,
              col.axis="blue", col.grid="lightblue",
              main="scatterplot3d - 2", pch=20)



########### fish boll
phi <- stats::runif(n, 0, 2 * pi)
psi <- acos(stats::runif(n, -1, 0.8))

fish_boll<-data.frame(x = cos(phi) * sin(psi),
                      y = sin(phi) * sin(psi),
                      z = cos(psi))

scatterplot3d(fish_boll$x, fish_boll$y, fish_boll$z, highlight.3d=TRUE,
              col.axis="blue", col.grid="lightblue",
              main="scatterplot3d - 2", pch=20)
z################ cube
n<-10000
sigma<-0.01
x = runif(n) + rnorm(n, sd = sigma)
y = runif(n) + rnorm(n, sd = sigma)
z = runif(n) + rnorm(n, sd = sigma)

cube<-data.frame(x = x,
                 y = y,
                 z = z)

scatterplot3d(cube$x, cube$y, cube$z, highlight.3d=TRUE,
              col.axis="blue", col.grid="lightblue",
              main="scatterplot3d - 2", pch=20)

