####################################################
# Aggregation and Symmetry scores for matrices
# (c)2025 Anna van der Kaaden
####################################################
# Load packages
require(terra)
require(ggplot2)

# Matrix dimensions same as annotated cold-water coral images
dimI <- c(720,1280)
# Image divisions
divs <- c(2,4,5,8,10,16,20,40,80)
# The spatial scale belonging to each division
Scalesx <- 1/divs
Scalesy <- 1/divs

# Emtpy array to store the matrices
Matrices <- array(data=NA,dim=c(5,dimI))

# Create matrices
# 1. Small-box checkerboard matrix
# Empty matrix of correct dimensions
M <- matrix(0,nrow=2*dimI[1]/40,ncol=dimI[2]/40)
# Create small square with 1s
M[1:18,1:32] <- 1
# Paste the square next to a square with 0s of the same size
M <- cbind(M,-(M-1))
# Bind this configuration along columns and rows
M <- cbind(M,M,M,M,M,M,M,M,M,M,M,M,M,M,M,M,M,M,M,M)
M <- rbind(M,M,M,M,M,M,M,M,M,M,M,M,M,M,M,M,M,M,M,M)
# Plot matrix as image
image(M)
# Store matrix
Matrices[1,,] <- M

# 2. Matrix of two clusters of aggregated red and yellow species
# Emtpy matrix of correct dimensions
M <- matrix(NA,nrow=dimI[1],ncol=dimI[2])
# Four data frames with 10000 random x locations around 200, 250, and 500 and
# y locations around 350, 400 and 950
a <- data.frame( x=rnorm(50000, 200, 60), y=rnorm(50000, 400, 60) )
b <- data.frame( x=rnorm(50000, 250, 60), y=rnorm(50000, 350, 60) )
c <- data.frame( x=rnorm(50000, 500, 60), y=rnorm(50000, 980, 60) )
d <- data.frame( x=rnorm(50000, 495, 60), y=rnorm(50000, 950, 60))
# Round and bind locations together
data <- rbind(round(a),round(b),round(c),round(d))
# Make sure that the locations are within the dimensions
data$x[data$x>dimI[1]]<-dimI[1]
data$x[data$x<1]<-1
data$y[data$y>dimI[2]]<-dimI[2]
data$y[data$y<1]<-1
# Randomly assign species to locations
for (j in 1:dim(data)[1]){
  M[data$x[j],data$y[j]]<-sample(0:1,size=1)
}
# Plot matrix as image
image(M)
# Store matrix
Matrices[2,,]<-M

# 3. Two species segregated with noise at the edge
# Create a gradient from 0 to 14 with noise
M <- matrix(seq(from=0,to=14,length.out=dimI[1]*dimI[2])+rnorm(dimI[1]*dimI[2],mean=0,sd=.3),nrow=dimI[1],ncol=dimI[2])
# Everything below 7 becomes 0 and everything above 7 becomes 1
M[M<=7] <- 0
M[M>=7] <- 1
# Plot matrix as image
image(M)
# Store matrix
Matrices[3,,] <- M

# 4. A cluster of one species in habitat created by another species
# Empty matrix of correct dimensions
M <- matrix(0,dimI[1],dimI[2])
# Create 12000 x and y locations
a <- data.frame( x=rnorm(3000, 470, 50), y=rnorm(3000, 410, 75) )
b <- data.frame( x=rnorm(3000, 450, 75), y=rnorm(3000, 425, 60) )
c <- data.frame( x=rnorm(3000, 425, 30), y=rnorm(3000, 380, 95) )
d <- data.frame( x=rnorm(3000, 435, 55), y=rnorm(3000, 385, 55))
# Round and bind locations together
data <- rbind(round(a),round(b),round(c),round(d))
# Assign 1s to all locations
for (i in 1:12000){
  M[data$x[i],data$y[i]]<-1
}
# Store matrix
Matrices[4,,] <- M

# Calculate C- Aggregationi- and Symmetry- scores
# Empty array to store cscore data
CScore <- array(data=NA,dim=c(5,9))
# Emtpy array to store the null model of the C-score
NullModel <- array(data=NA,dim=c(5,9))
# Emtpy arrays to store the Aggregation and Symmetry scores
AScore <- array(data=NA,dim=c(5,9))
SScore <- array(data=NA,dim=c(5,9))

for (i in 1:4){
  M <- Matrices[i,,]
  # Calculate Cscore dividing images into divs along both dimensions
  for (j in 1:9){
    # Aggregate full images into larger blocks
    A <- array(0,dim(M))
    A[M==0] <- 1
    x <- aggregate(rast(A),c(dimI[1]/divs[j],dimI[2]/divs[j]),max)
    A <- as.matrix(x,dim(x)[1],dim(x)[2])
    B <- array(0,dim(M))
    B[M==1] <- 1
    x <- aggregate(rast(B),c(dimI[1]/divs[j],dimI[2]/divs[j]),max)
    B <- as.matrix(x,dim(x)[1],dim(x)[2])
    
    # ra = number of blocks containing species a and not b
    ra <- sum(A==1&B==0)
    # rb = number of blocks containing species b and not a
    rb <- sum(B==1&A==0)
    # Xab = number of blocks containing both species
    Xab <- sum(A==1&B==1)
    
    # Calculate Aggregation Score and Symmetry Score
    AScore[i,j] <- (Xab/(Xab+ra))*(Xab/(Xab+rb))
    SScore[i,j] <- min(c(ra,rb))/max(c(ra,rb))
    if (is.na((Xab/(Xab+ra))*(Xab/(Xab+rb)))){
      SScore[i,j] <- NaN
    }
    
    # Calculate Cscore
    # n = number of blocks containing species a, b, or both
    n <- sum(A==1 | B==1 | (A==1&B==1))
    N <- n*(n-1)/2
    CScore[i,j] <- (ra-Xab)*(rb-Xab)/N
    
    # Calculate null-model Cscore
    # randomize species presences 
    A <- A[sample(length(A),replace=FALSE)]
    B <- B[sample(length(B),replace=FALSE)]
    # ra = number of blocks containing species a and not b
    ra <- sum(A==1&B==0)
    # rb = number of blocks containing species b and not a
    rb <- sum(B==1&A==0)
    # Xab = number of blocks containing both species
    Xab <- sum(A==1&B==1)
    NullModel[i,j] <- (ra-Xab)*(rb-Xab)/N
  }
}

# plot Matrices, CScore and A- and S-score
quartz(width=1280/500,height=5*720/500)
par(mfrow=c(4,1))
for (i in 1:4){
  if (i==1){
    par(mar=c(.3,2,1,2))
    image(x=seq(0,1,length.out=1281),y=seq(0,1,length.out=721),z=t(Matrices[i,,]),ylab='',xlab='',col.axis='white',
          zlim=c(0,1))
    axis(2,col.axis='black',at=c(0,.2,.4,.6,.8,1.0),cex.axis=0.7)
  } else if (i<5){
    par(mar=c(.3,2,.3,2))
    image(x=seq(0,1,length.out=1280),y=seq(0,1,length.out=720),z=t(Matrices[i,,]),ylab='',xlab='',col.axis='white')
    axis(2,col.axis='black',at=c(0,.2,.4,.6,.8,1.0),cex.axis=0.7)
  } else {
    par(mar=c(2,2,.3,2))
    image(x=seq(0,1,length.out=1280),y=seq(0,1,length.out=720),z=t(Matrices[i,,]),ylab='',xlab='',cex.axis=0.7)
  }
}
quartz.save('Matrices.png',type="png",dpi=300)

quartz(height=3.5)
ggplot(data=data.frame(t(CScore)))+
  geom_point(shape=3,aes(x=Scalesx,y=CScore[1,]),col='blueviolet')+
  geom_line(linewidth=.5,aes(x=Scalesx,y=CScore[1,]),col='blueviolet')+
  geom_line(linewidth=.5,aes(x=Scalesx,y=NullModel[1,]),col='blueviolet',lty=2)+
  geom_point(shape=3,aes(x=Scalesx,y=CScore[2,]),col='cornflowerblue')+
  geom_line(linewidth=.5,aes(x=Scalesx,y=CScore[2,]),col='cornflowerblue')+
  geom_line(linewidth=.5,aes(x=Scalesx,y=NullModel[2,]),col='cornflowerblue',lty=2)+
  geom_point(shape=3,aes(x=Scalesx,y=CScore[3,]),col='aquamarine2')+
  geom_line(linewidth=.5,aes(x=Scalesx,y=CScore[3,]),col='aquamarine2')+
  geom_line(linewidth=.5,aes(x=Scalesx,y=NullModel[3,]),col='aquamarine2',lty=2)+
  geom_point(shape=3,aes(x=Scalesx,y=CScore[4,]),col='chartreuse3')+
  geom_line(linewidth=.5,aes(x=Scalesx,y=CScore[4,]),col='chartreuse3')+
  geom_line(linewidth=.5,aes(x=Scalesx,y=NullModel[4,]),col='chartreuse3',lty=2)+
  scale_y_continuous(labels=scales::number_format(accuracy=.01),limits=c(-.5,3.00))
quartz.save('Cscore_Matrices.png',type="png",dpi=300)

p <- ggplot(data=data.frame(t(AScore)))+
  geom_point(shape=4,aes(x=Scalesx,y=AScore[1,]),col='blueviolet')+
  geom_line(linewidth=.5,aes(x=Scalesx,y=AScore[1,]),col='blueviolet')+
  geom_point(shape=4,aes(x=Scalesx,y=AScore[2,]),col='cornflowerblue')+
  geom_line(linewidth=.5,aes(x=Scalesx,y=AScore[2,]),col='cornflowerblue')+
  geom_point(shape=4,aes(x=Scalesx,y=AScore[3,]),col='aquamarine2')+
  geom_line(linewidth=.5,aes(x=Scalesx,y=AScore[3,]),col='aquamarine2')+
  geom_point(shape=4,aes(x=Scalesx,y=AScore[4,]),col='chartreuse3')+
  geom_line(linewidth=.5,aes(x=Scalesx,y=AScore[4,]),col='chartreuse3')

text_labels1 <- data.frame(x=Scalesx[!is.na(SScore[1,])],y=AScore[1,!is.na(SScore[1,])],
                           label=paste(round(SScore[1,!is.na(SScore[1,])],digits=2)))
text_labels2 <- data.frame(x=Scalesx[!is.na(SScore[2,])],y=AScore[2,!is.na(SScore[2,])],
                           label=paste(round(SScore[2,!is.na(SScore[2,])],digits=2)))
text_labels3 <- data.frame(x=Scalesx[!is.na(SScore[3,])],y=AScore[3,!is.na(SScore[3,])],
                           label=paste(round(SScore[3,!is.na(SScore[3,])],digits=2)))
text_labels4 <- data.frame(x=Scalesx[!is.na(SScore[4,])],y=AScore[4,!is.na(SScore[4,])],
                           label=paste(round(SScore[4,!is.na(SScore[4,])],digits=2)))
text_labels5 <- data.frame(x=Scalesx[!is.na(SScore[5,])],y=AScore[5,!is.na(SScore[5,])],
                           label=paste(round(SScore[5,!is.na(SScore[5,])],digits=2)))
quartz(height=3.5)
p + geom_text(data=text_labels2, aes(x=x,y=y,label=label),vjust=c(1.5,1.5,1.5,1.5,2.3,2,1.9),size=3)+
  geom_text(data=text_labels3, aes(x=x,y=y,label=label),vjust=c(-.3,-.7,-.7,-.7,-.3,-.7,-.9,-.7),size=3)+
  geom_text(data=text_labels4, aes(x=x,y=y,label=label),vjust=-.7,size=3)+
  geom_text(data=text_labels1, aes(x=x,y=y,label=label),vjust= 1.5,size=3)
  
quartz.save('Scores_Matrices.png',type="png",dpi=300)

