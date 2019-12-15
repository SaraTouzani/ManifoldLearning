fashion.mnist_test <- read.csv("/home/toutou/Téléchargements/fashion-mnist_test.csv")
View(fashion.mnist_test)
data<-as.matrix(fashion.mnist_test[,-1])
labels<-fashion.mnist_test[,1]


show_picture<-function(line){
  p<-sqrt(length(line))
  mat<-matrix(as.numeric(line), p)[,p:1]
  image(mat,col=grey(256:1/256))
  return(mat)
}

show_picture(data[1,])
labels[1]

Average_picture<-(colMeans(data)) 

pr.out = prcomp(data,center=TRUE, scale=FALSE)
#pr.out$sdev
pr.var=(pr.out$sdev)**2
#pve = pr.var/sum(pr.var) 

U = pr.out$rotation
Z = t(pr.out$x)
for (i in c(10,500))### on définit le nombre de dimention à garder 
{
  show_picture((U[,1:i]%*%Z[1:i,])[,1]) 
}

##### pour avoir encore l'image originale :

for (i in c(10,500))
{
  show_picture((U[,1:i]%*%Z[1:i,])[,1] + Average_picture) 
}
