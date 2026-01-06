#################################################################################
# This script simulates communities with interacting species
# (c)2025 A. van der Kaaden
############################################################################

# Imports necessary packages
library(ReacTran)
library(RColorBrewer)
require(ggplot2)

# The collection of functions used in the script
# Function to calculate ra, rb, and Xab
CalXrab <- function(Matrix,species,divs){
  require(terra)
  scale = length(divs)
  dimM = dim(Matrix)[1]
  Xrab<-array(data=NA,dim=c(3,scale))
  for (j in 1:scale){
    # For all species pairs
    for (a in 1){
      for (b in 2){
        if (a<b){
          # Aggregate full images into larger blocks
          A <- Matrix[,,a]
          A[A<.1] <- 0
          A <- as.array(A)
          x <- aggregate(rast(A),c(dimM/divs[j],dimM/divs[j]),max)
          A <- as.matrix(x,dim(x)[1],dim(x)[2])
          B <- Matrix[,,b]
          B[B<.1] <- 0
          B <- as.array(B)
          x <- aggregate(rast(B),c(dimM/divs[j],dimM/divs[j]),max)
          B <- as.matrix(x,dim(x)[1],dim(x)[2])
          
          # ra = number of blocks containing species a and not b
          ra <- sum(A>0&B==0)
          # rb = number of blocks containing species b and not a
          rb <- sum(B>0&A==0)
          # Xab = number of blocks containing both species
          Xab <- sum(A>0&B>0)
          
          Xrab[1,j]<-Xab
          Xrab[2,j]<-ra
          Xrab[3,j]<-rb
        }
      }
    }
  }
  return(Xrab)
}
# Function to initiate simulation grid
Initiate <- function(q,s,n,ID){
  N <- array(0,dim=c(q,q,s))
  for (i in 1:s){
    rx=runif(n=n,min=1,max=q)
    ry=runif(n=n,min=1,max=q)
    for (j in 1:n){
      N[(ry[j]-1):(ry[j]+1),(rx[j]-1):(rx[j]+1),i]=ID
    }
  }
  return(N)
}
# Function to perform the simulations
Simulate <- function(N, q, r, c, D, dx, dt, endT, vis) {
  numSpecies <- dim(N)[3]  # Get number of species
  # Initialize arrays for changes
  dNdt <- array(0, dim = dim(N))

  for (t in 1:endT) {
    # Stochastic mortality
    m <- array(rexp(q*q * numSpecies, rate = 12), dim=dim(N))  # Exponential random variables
    
    # Reaction
    for (i in 1:numSpecies) {
      reactionTerm <- r[i] * N[,,i]
      for (j in 1:numSpecies) {
        reactionTerm <- reactionTerm - c[i, j] * N[,,i] * N[,,j]
      }
      dNdt[,,i] <- reactionTerm - m[,,i] * N[,,i]
    }
    
    # Diffusion using ReacTran
    for (i in 1:numSpecies) {
      # Define the diffusion coefficient for each species
      tran <- tran.2D(N[,,i],D.x=D[i],dx=dx,dy=dx)
      
      # Update N with diffusion
      N[,,i] <- N[,,i] + (dNdt[,,i]+tran$dC) * dt
    }
    
    # Ensure populations are non-negative
    N[N < 1e-08] <- 0
    
    # Reintroduce species that went extinct
    for (i in 1:numSpecies){
      if (sum(N[,,i]>0)==0){
        rx=round(runif(n=1,min=6,max=(q-5)),0)
        ry=round(runif(n=1,min=6,max=(q-5)),0)
        N[(ry-5):(ry+5),(rx-5):(rx+5),i]<-100
      }
    }
    
    if (vis>0){
      # Visualization every 100 time steps
      if (t %% 100 == 0) {
        par(mfrow = c(1, numSpecies))
        for (i in 1:numSpecies) {
          image(1:q, 1:q, N[,,i], col = heat.colors(256), axes = FALSE)
          title(paste("Species", i))
          box()
        }
        text(q/5,q/10,paste("Time:", t))
        Sys.sleep(0.1)  # Pause to visualize
      }
    }
  }
  return(N)
}
# Function to create plots
CreatePlot <- function(Data,divs){
  QTLS_neutralF <- matrix(Data$c.QTLS_neutralF.,nrow=5,byrow=FALSE)
  QTLS_wcompF <- matrix(Data$c.QTLS_wcompF.,nrow=5,byrow=FALSE)
  QTLS_wfacF <- matrix(Data$c.QTLS_wfacF.,nrow=5,byrow=FALSE)
  QTLS_wdepF <- matrix(Data$c.QTLS_wdepF.,nrow=5,byrow=FALSE)
  QTLS_scompF <- matrix(Data$c.QTLS_scompF.,nrow=5,byrow=FALSE)
  QTLS_sfacF <- matrix(Data$c.QTLS_sfacF.,nrow=5,byrow=FALSE)
  QTLS_sdepF <- matrix(Data$c.QTLS_sdepF.,nrow=5,byrow=FALSE)
  QTLS_neutralD <- matrix(Data$c.QTLS_neutralD.,nrow=5,byrow=FALSE)
  QTLS_wcompD <- matrix(Data$c.QTLS_wcompD.,nrow=5,byrow=FALSE)
  QTLS_wfacD <- matrix(Data$c.QTLS_wfacD.,nrow=5,byrow=FALSE)
  QTLS_wdepD <- matrix(Data$c.QTLS_wdepD.,nrow=5,byrow=FALSE)
  QTLS_scompD <- matrix(Data$c.QTLS_scompD.,nrow=5,byrow=FALSE)
  QTLS_sfacD <- matrix(Data$c.QTLS_sfacD.,nrow=5,byrow=FALSE)
  QTLS_sdepD <- matrix(Data$c.QTLS_sdepD.,nrow=5,byrow=FALSE)
  
  # Weak interactions Aggregation scores
  p1 <- ggplot()+
    geom_line(aes(x=c(-1,0),y=c(0,0)),color='orange',lwd=6,alpha=0.5)+
    geom_line(aes(x=c(-1,0),y=c(0.25,0.25)),color='grey95',lwd=6,alpha=0.9)+
    geom_line(aes(x=c(-1,0),y=c(0.5,0.5)),color='yellow',lwd=6,alpha=0.5)+
    geom_line(aes(x=c(-1,0),y=c(1,1)),color='darkolivegreen1',lwd=6,alpha=0.5)+
    geom_line(aes(x=c(.5,1),y=c(0,0)),color='orange',lwd=6,alpha=0.5)+
    geom_line(aes(x=c(.5,1),y=c(0.25,0.25)),color='grey95',lwd=6,alpha=0.9)+
    geom_line(aes(x=c(.5,1),y=c(0.5,0.5)),color='yellow',lwd=6,alpha=0.5)+
    geom_line(aes(x=c(.5,1),y=c(1,1)),color='darkolivegreen1',lwd=6,alpha=0.5)+
    geom_polygon(aes(x=c(1/(divs),rev(1/(divs))),y=c(QTLS_neutralF[2,],rev(QTLS_neutralF[4,]))),fill='white')+
    geom_line(aes(x=1/(divs),y=QTLS_neutralF[3,]),color='grey65')+
    geom_point(aes(x=1/(divs),y=QTLS_neutralF[3,]),color='grey65')+
    geom_polygon(aes(x=c(1/(divs),rev(1/(divs))),y=c(QTLS_wcompF[2,],rev(QTLS_wcompF[4,]))),fill='orange',alpha=0.3)+
    geom_line(aes(x=1/(divs),y=QTLS_wcompF[3,]),color='darkorange')+
    geom_point(aes(x=1/(divs),y=QTLS_wcompF[3,]),color='darkorange')+
    geom_polygon(aes(x=c(1/(divs),rev(1/(divs))),y=c(QTLS_wfacF[2,],rev(QTLS_wfacF[4,]))),fill='darkolivegreen1',alpha=0.3)+
    geom_line(aes(x=1/(divs),y=QTLS_wfacF[3,]),color='darkolivegreen')+
    geom_point(aes(x=1/(divs),y=QTLS_wfacF[3,]),color='darkolivegreen')+
    geom_polygon(aes(x=c(1/(divs),rev(1/(divs))),y=c(QTLS_wdepF[2,],rev(QTLS_wdepF[4,]))),fill='gold',alpha=0.3)+
    geom_line(aes(x=1/(divs),y=QTLS_wdepF[3,]),color='gold3')+
    geom_point(aes(x=1/(divs),y=QTLS_wdepF[3,]),color='gold3')+
    coord_cartesian(xlim=c(0,0.5))+
    labs(x=NULL,y=NULL)
  
  # Strong interactions Aggregation scores
  p2 <- ggplot()+
    geom_line(aes(x=c(-1,0),y=c(0,0)),color='orange',lwd=6,alpha=0.5)+
    geom_line(aes(x=c(-1,0),y=c(0.25,0.25)),color='grey95',lwd=6,alpha=0.9)+
    geom_line(aes(x=c(-1,0),y=c(0.5,0.5)),color='yellow',lwd=6,alpha=0.5)+
    geom_line(aes(x=c(-1,0),y=c(1,1)),color='darkolivegreen1',lwd=6,alpha=0.5)+
    geom_line(aes(x=c(.5,1),y=c(0,0)),color='orange',lwd=6,alpha=0.5)+
    geom_line(aes(x=c(.5,1),y=c(0.25,0.25)),color='grey95',lwd=6,alpha=0.9)+
    geom_line(aes(x=c(.5,1),y=c(0.5,0.5)),color='yellow',lwd=6,alpha=0.5)+
    geom_line(aes(x=c(.5,1),y=c(1,1)),color='darkolivegreen1',lwd=6,alpha=0.5)+
    geom_polygon(aes(x=c(1/(divs),rev(1/(divs))),y=c(QTLS_neutralF[2,],rev(QTLS_neutralF[4,]))),fill='white')+
    geom_line(aes(x=1/(divs),y=QTLS_neutralF[3,]),color='grey65')+
    geom_point(aes(x=1/(divs),y=QTLS_neutralF[3,]),color='grey65')+
    geom_polygon(aes(x=c(1/(divs),rev(1/(divs))),y=c(QTLS_scompF[2,],rev(QTLS_scompF[4,]))),fill='orange',alpha=0.3)+
    geom_line(aes(x=1/(divs),y=QTLS_scompF[3,]),color='darkorange')+
    geom_point(aes(x=1/(divs),y=QTLS_scompF[3,]),color='darkorange')+
    geom_polygon(aes(x=c(1/(divs),rev(1/(divs))),y=c(QTLS_sfacF[2,],rev(QTLS_sfacF[4,]))),fill='darkolivegreen1',alpha=0.3)+
    geom_line(aes(x=1/(divs),y=QTLS_sfacF[3,]),color='darkolivegreen')+
    geom_point(aes(x=1/(divs),y=QTLS_sfacF[3,]),color='darkolivegreen')+
    geom_polygon(aes(x=c(1/(divs),rev(1/(divs))),y=c(QTLS_sdepF[2,],rev(QTLS_sdepF[4,]))),fill='gold',alpha=0.3)+
    geom_line(aes(x=1/(divs),y=QTLS_sdepF[3,]),color='gold3')+
    geom_point(aes(x=1/(divs),y=QTLS_sdepF[3,]),color='gold3')+
    coord_cartesian(xlim=c(0,0.5))+
    labs(x=NULL,y=NULL)
  
  # Weak interactions Symmetry scores
  p3 <- ggplot()+
    geom_line(aes(x=c(-1,0),y=c(0.5,0.5)),color='grey95',lwd=6,alpha=0.9)+
    geom_line(aes(x=c(-1,0),y=c(0,0)),color='yellow',lwd=6,alpha=0.5)+
    geom_line(aes(x=c(.5,1),y=c(0.5,0.5)),color='grey95',lwd=6,alpha=0.9)+
    geom_line(aes(x=c(.5,1),y=c(0,0)),color='yellow',lwd=6,alpha=0.5)+
    geom_polygon(aes(x=c(1/(divs),rev(1/(divs))),y=c(QTLS_neutralD[2,],rev(QTLS_neutralD[4,]))),fill='white')+
    geom_line(aes(x=1/(divs),y=QTLS_neutralD[3,]),color='grey65')+
    geom_point(aes(x=1/(divs),y=QTLS_neutralD[3,]),color='grey65')+
    geom_polygon(aes(x=c(1/(divs),rev(1/(divs))),y=c(QTLS_wcompD[2,],rev(QTLS_wcompD[4,]))),fill='orange',alpha=0.3)+
    geom_line(aes(x=1/(divs),y=QTLS_wcompD[3,]),color='darkorange')+
    geom_point(aes(x=1/(divs),y=QTLS_wcompD[3,]),color='darkorange')+
    geom_polygon(aes(x=c(1/(divs),rev(1/(divs))),y=c(QTLS_wfacD[2,],rev(QTLS_wfacD[4,]))),fill='darkolivegreen1',alpha=0.3)+
    geom_line(aes(x=1/(divs),y=QTLS_wfacD[3,]),color='darkolivegreen')+
    geom_point(aes(x=1/(divs),y=QTLS_wfacD[3,]),color='darkolivegreen')+
    geom_polygon(aes(x=c(1/(divs),rev(1/(divs))),y=c(QTLS_wdepD[2,],rev(QTLS_wdepD[4,]))),fill='gold',alpha=0.3)+
    geom_line(aes(x=1/(divs),y=QTLS_wdepD[3,]),color='gold3')+
    geom_point(aes(x=1/(divs),y=QTLS_wdepD[3,]),color='gold3')+
    coord_cartesian(xlim=c(0,0.5),ylim=c(0,1))+
    labs(x=NULL,y=NULL)
  
  # Strong interactions Symmetry scores
  p4 <- ggplot()+
    geom_line(aes(x=c(-1,0),y=c(0.5,0.5)),color='grey95',lwd=6,alpha=0.9)+
    geom_line(aes(x=c(-1,0),y=c(0,0)),color='yellow',lwd=6,alpha=0.5)+
    geom_line(aes(x=c(.5,1),y=c(0.5,0.5)),color='grey95',lwd=6,alpha=0.9)+
    geom_line(aes(x=c(.5,1),y=c(0,0)),color='yellow',lwd=6,alpha=0.5)+
    geom_polygon(aes(x=c(1/(divs),rev(1/(divs))),y=c(QTLS_neutralD[2,],rev(QTLS_neutralD[4,]))),fill='white')+
    geom_line(aes(x=1/(divs),y=QTLS_neutralD[3,]),color='grey65')+
    geom_point(aes(x=1/(divs),y=QTLS_neutralD[3,]),color='grey65')+
    geom_polygon(aes(x=c(1/(divs),rev(1/(divs))),y=c(QTLS_scompD[2,],rev(QTLS_scompD[4,]))),fill='orange',alpha=0.3)+
    geom_line(aes(x=1/(divs),y=QTLS_scompD[3,]),color='darkorange')+
    geom_point(aes(x=1/(divs),y=QTLS_scompD[3,]),color='darkorange')+
    geom_polygon(aes(x=c(1/(divs),rev(1/(divs))),y=c(QTLS_sfacD[2,],rev(QTLS_sfacD[4,]))),fill='darkolivegreen1',alpha=0.3)+
    geom_line(aes(x=1/(divs),y=QTLS_sfacD[3,]),color='darkolivegreen')+
    geom_point(aes(x=1/(divs),y=QTLS_sfacD[3,]),color='darkolivegreen')+
    geom_polygon(aes(x=c(1/(divs),rev(1/(divs))),y=c(QTLS_sdepD[2,],rev(QTLS_sdepD[4,]))),fill='gold',alpha=0.3)+
    geom_line(aes(x=1/(divs),y=QTLS_sdepD[3,]),color='gold3')+
    geom_point(aes(x=1/(divs),y=QTLS_sdepD[3,]),color='gold3')+
    coord_cartesian(xlim=c(0,0.5),ylim=c(0,1))+
    labs(x=NULL,y=NULL)
  
  Plots <- list()
  Plots[[1]] <- p1
  Plots[[2]] <- p2
  Plots[[3]] <- p3
  Plots[[4]] <- p4
  return(Plots)
}

# Simulation characteristics
q <- 180  # Grid size
# Grid divisions to calculate the scores for
divs <- c(2,3,4,5,6,9,10,12,15,18,20,30,36,45,60,90)
dx <- 1  # Spatial step
dt <- 1  # Time step
nns <- 100 # number of repetitions

# Two species
endT <- 1000  # End time
{
  numSpecies <- 2
  Abundance <- array(0,dim=c(1,numSpecies))
  r <- array(0.3,dim=c(numSpecies)) # Growth rates for each species
  D <- array(0.1,dim=c(numSpecies)) # Diffusion coefficients
# Competition array for neutral simulation
  c <- array(0.001,dim=c(numSpecies,numSpecies)) # competition array
# Initiate empty arrays to store scores
AScores <- array(data=NA,dim=c(length(divs),nns))
SScores <- array(data=NA,dim=c(length(divs),nns))
# Simulate for all repetitions
for (i in 1:nns){
  # Initiate grid for simulation
  N <- Initiate(q,numSpecies,3,300)
  # Simulate and store final grid state
  N_final <- Simulate(N,q,r,c,D,dx,dt,endT,0)
  # Calculate Xab, ra, and rb
  Xrab <- CalXrab(N_final,numSpecies,divs)
  Xab <- Xrab[1,]
  ra <- Xrab[2,]
  rb <- Xrab[3,]
  # Calculate the scores
  AScores[,i] <- (Xab/(Xab+ra))*(Xab/(Xab+rb))
  SScores[,i] <- apply(rbind(ra,rb),MARGIN=2,FUN=min)/apply(rbind(ra,rb),MARGIN=2,FUN=max)
  # Store abundances of all species
  Ab <- array(0,dim=c(1,numSpecies))
  for (j in 1:numSpecies){
    Ab[1,j] <- sum(N_final[,,j]>0)/(180*180)
  }
  Abundance <- rbind(Abundance,Ab)
}
# Calculate quantiles for plotting
QTLS_neutralD <- apply(X=SScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
QTLS_neutralF <- apply(X=AScores,MARGIN=1,FUN=quantile,na.rm=TRUE)

# Competition (weak)
c[1,2] <- 0.005
c[2,1] <- 0.005
AScores <- array(data=NA,dim=c(length(divs),nns))
SScores <- array(data=NA,dim=c(length(divs),nns))
for (i in 1:nns){
  N <- Initiate(q,numSpecies,3,300)
  N_final <- Simulate(N,q,r,c,D,dx,dt,endT,0)
  Xrab <- CalXrab(N_final,numSpecies,divs)
  Xab <- Xrab[1,]
  ra <- Xrab[2,]
  rb <- Xrab[3,]
  AScores[,i] <- (Xab/(Xab+ra))*(Xab/(Xab+rb))
  SScores[,i] <- apply(rbind(ra,rb),MARGIN=2,FUN=min)/apply(rbind(ra,rb),MARGIN=2,FUN=max)
  Ab <- array(0,dim=c(1,numSpecies))
  for (j in 1:numSpecies){
    Ab[1,j] <- sum(N_final[,,j]>0)/(180*180)
  }
  Abundance <- rbind(Abundance,Ab)
}
QTLS_wcompD <- apply(X=SScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
QTLS_wcompF <- apply(X=AScores,MARGIN=1,FUN=quantile,na.rm=TRUE)

# Competition (strong)
c[1,2] <- 0.01
c[2,1] <- 0.01
AScores <- array(data=NA,dim=c(length(divs),nns))
SScores <- array(data=NA,dim=c(length(divs),nns))
for (i in 1:nns){
  N <- Initiate(q,numSpecies,3,300)
  N_final <- Simulate(N,q,r,c,D,dx,dt,endT,0)
  Xrab <- CalXrab(N_final,numSpecies,divs)
  Xab <- Xrab[1,]
  ra <- Xrab[2,]
  rb <- Xrab[3,]
  AScores[,i] <- (Xab/(Xab+ra))*(Xab/(Xab+rb))
  SScores[,i] <- apply(rbind(ra,rb),MARGIN=2,FUN=min)/apply(rbind(ra,rb),MARGIN=2,FUN=max)
  Ab <- array(0,dim=c(1,numSpecies))
  for (j in 1:numSpecies){
    Ab[1,j] <- sum(N_final[,,j]>0)/(180*180)
  }
  Abundance <- rbind(Abundance,Ab)
}
QTLS_scompD <- apply(X=SScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
QTLS_scompF <- apply(X=AScores,MARGIN=1,FUN=quantile,na.rm=TRUE)

# Facilitation (weak)
c[1,2] <- -0.00005
c[2,1] <- -0.00005
r[1] <- r[2] <- 0.24
AScores <- array(data=NA,dim=c(length(divs),nns))
SScores <- array(data=NA,dim=c(length(divs),nns))
for (i in 1:100){
  N <- Initiate(q,numSpecies,3,300)
  N_final <- Simulate(N,q,r,c,D,dx,dt,endT,0)
  Xrab <- CalXrab(N_final,numSpecies,divs)
  Xab <- Xrab[1,]
  ra <- Xrab[2,]
  rb <- Xrab[3,]
  AScores[,i] <- (Xab/(Xab+ra))*(Xab/(Xab+rb))
  SScores[,i] <- apply(rbind(ra,rb),MARGIN=2,FUN=min)/apply(rbind(ra,rb),MARGIN=2,FUN=max)
  Ab <- array(0,dim=c(1,numSpecies))
  for (j in 1:numSpecies){
    Ab[1,j] <- sum(N_final[,,j]>0)/(180*180)
  }
  Abundance <- rbind(Abundance,Ab)
}
QTLS_wfacD <- apply(X=SScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
QTLS_wfacF <- apply(X=AScores,MARGIN=1,FUN=quantile,na.rm=TRUE)

# Facilitation (strong)
c[1,2] <- -0.0001
c[2,1] <- -0.0001
r[1] <- r[2] <- 0.24
AScores <- array(data=NA,dim=c(length(divs),nns))
SScores <- array(data=NA,dim=c(length(divs),nns))
for (i in 1:nns){
  N <- Initiate(q,numSpecies,3,300)
  N_final <- Simulate(N,q,r,c,D,dx,dt,endT,0)
  Xrab <- CalXrab(N_final,numSpecies,divs)
  Xab <- Xrab[1,]
  ra <- Xrab[2,]
  rb <- Xrab[3,]
  AScores[,i] <- (Xab/(Xab+ra))*(Xab/(Xab+rb))
  SScores[,i] <- apply(rbind(ra,rb),MARGIN=2,FUN=min)/apply(rbind(ra,rb),MARGIN=2,FUN=max)
  Ab <- array(0,dim=c(1,numSpecies))
  for (j in 1:numSpecies){
    Ab[1,j] <- sum(N_final[,,j]>0)/(180*180)
  }
  Abundance <- rbind(Abundance,Ab)
}
QTLS_sfacD <- apply(X=SScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
QTLS_sfacF <- apply(X=AScores,MARGIN=1,FUN=quantile,na.rm=TRUE)

# Dependence (weak)
c[1,2] <- 0.001
c[2,1] <- -0.00005
r[1] <- 0.3
r[2]<-0.08
AScores <- array(data=NA,dim=c(length(divs),nns))
SScores <- array(data=NA,dim=c(length(divs),nns))
for (i in 1:nns){
  N <- Initiate(q,numSpecies,3,300)
  N_final <- Simulate(N,q,r,c,D,dx,dt,endT,0)
  Xrab <- CalXrab(N_final,numSpecies,divs)
  Xab <- Xrab[1,]
  ra <- Xrab[2,]
  rb <- Xrab[3,]
  AScores[,i] <- (Xab/(Xab+ra))*(Xab/(Xab+rb))
  SScores[,i] <- apply(rbind(ra,rb),MARGIN=2,FUN=min)/apply(rbind(ra,rb),MARGIN=2,FUN=max)
  Ab <- array(0,dim=c(1,numSpecies))
  for (j in 1:numSpecies){
    Ab[1,j] <- sum(N_final[,,j]>0)/(180*180)
  }
  Abundance <- rbind(Abundance,Ab)
}
QTLS_wdepD <- apply(X=SScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
QTLS_wdepF <- apply(X=AScores,MARGIN=1,FUN=quantile,na.rm=TRUE)

# Dependence (strong)
c[1,2] <- 0.001
c[2,1] <- -0.0001
r[1] <- 0.3
r[2]<-0.08
AScores <- array(data=NA,dim=c(length(divs),nns))
SScores <- array(data=NA,dim=c(length(divs),nns))
for (i in 1:nns){
  N <- Initiate(q,numSpecies,3,300)
  N_final <- Simulate(N,q,r,c,D,dx,dt,endT,0)
  Xrab <- CalXrab(N_final,numSpecies,divs)
  Xab <- Xrab[1,]
  ra <- Xrab[2,]
  rb <- Xrab[3,]
  AScores[,i] <- (Xab/(Xab+ra))*(Xab/(Xab+rb))
  SScores[,i] <- apply(rbind(ra,rb),MARGIN=2,FUN=min)/apply(rbind(ra,rb),MARGIN=2,FUN=max)
  Ab <- array(0,dim=c(1,numSpecies))
  for (j in 1:numSpecies){
    Ab[1,j] <- sum(N_final[,,j]>0)/(180*180)
  }
  Abundance <- rbind(Abundance,Ab)
}
QTLS_sdepD <- apply(X=SScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
QTLS_sdepF <- apply(X=AScores,MARGIN=1,FUN=quantile,na.rm=TRUE)

Data<-data.frame(c(QTLS_neutralF),c(QTLS_scompF),c(QTLS_sdepF),c(QTLS_sfacF),c(QTLS_wcompF),c(QTLS_wdepF),c(QTLS_wfacF),
                 c(QTLS_neutralD),c(QTLS_scompD),c(QTLS_sdepD),c(QTLS_sfacD),c(QTLS_wcompD),c(QTLS_wdepD),c(QTLS_wfacD))
save(Data,file='Data_2sp.rdata')
write.csv2(Abundance,'Abundance_2sp.csv')
}

# Three species
endT <- 1000
{
  numSpecies <- 3
  Abundance <- array(0,dim=c(1,numSpecies))
  r <- array(0.3,dim=c(numSpecies)) # Growth rates for each species
  D <- array(0.1,dim=c(numSpecies)) # Diffusion coefficients
  # Competition array for neutral simulation
  c <- array(0.001,dim=c(numSpecies,numSpecies)) # competition array
  # Initiate empty arrays to store scores
  AScores <- array(data=NA,dim=c(length(divs),nns))
  SScores <- array(data=NA,dim=c(length(divs),nns))
  # Simulate for all repetitions
  for (i in 1:nns){
    # Initiate grid for simulation
    N <- Initiate(q,numSpecies,3,300)
    # Simulate and store final grid state
    N_final <- Simulate(N,q,r,c,D,dx,dt,endT,0)
    # Calculate Xab, ra, and rb
    Xrab <- CalXrab(N_final,numSpecies,divs)
    Xab <- Xrab[1,]
    ra <- Xrab[2,]
    rb <- Xrab[3,]
    # Calculate the scores
    AScores[,i] <- (Xab/(Xab+ra))*(Xab/(Xab+rb))
    SScores[,i] <- apply(rbind(ra,rb),MARGIN=2,FUN=min)/apply(rbind(ra,rb),MARGIN=2,FUN=max)
    # Store abundances of all species
    Ab <- array(0,dim=c(1,numSpecies))
    for (j in 1:numSpecies){
      Ab[1,j] <- sum(N_final[,,j]>0)/(180*180)
    }
    Abundance <- rbind(Abundance,Ab)
  }
  # Calculate quantiles for plotting
  QTLS_neutralD <- apply(X=SScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
  QTLS_neutralF <- apply(X=AScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
  
  # Competition (weak)
  c[1,2] <- 0.005
  c[2,1] <- 0.005
  AScores <- array(data=NA,dim=c(length(divs),nns))
  SScores <- array(data=NA,dim=c(length(divs),nns))
  for (i in 1:nns){
    N <- Initiate(q,numSpecies,3,300)
    N_final <- Simulate(N,q,r,c,D,dx,dt,endT,0)
    Xrab <- CalXrab(N_final,numSpecies,divs)
    Xab <- Xrab[1,]
    ra <- Xrab[2,]
    rb <- Xrab[3,]
    AScores[,i] <- (Xab/(Xab+ra))*(Xab/(Xab+rb))
    SScores[,i] <- apply(rbind(ra,rb),MARGIN=2,FUN=min)/apply(rbind(ra,rb),MARGIN=2,FUN=max)
    Ab <- array(0,dim=c(1,numSpecies))
    for (j in 1:numSpecies){
      Ab[1,j] <- sum(N_final[,,j]>0)/(180*180)
    }
    Abundance <- rbind(Abundance,Ab)
  }
  QTLS_wcompD <- apply(X=SScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
  QTLS_wcompF <- apply(X=AScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
  
  # Competition (strong)
  c[1,2] <- 0.01
  c[2,1] <- 0.01
  AScores <- array(data=NA,dim=c(length(divs),nns))
  SScores <- array(data=NA,dim=c(length(divs),nns))
  for (i in 1:nns){
    N <- Initiate(q,numSpecies,3,300)
    N_final <- Simulate(N,q,r,c,D,dx,dt,endT,0)
    Xrab <- CalXrab(N_final,numSpecies,divs)
    Xab <- Xrab[1,]
    ra <- Xrab[2,]
    rb <- Xrab[3,]
    AScores[,i] <- (Xab/(Xab+ra))*(Xab/(Xab+rb))
    SScores[,i] <- apply(rbind(ra,rb),MARGIN=2,FUN=min)/apply(rbind(ra,rb),MARGIN=2,FUN=max)
    Ab <- array(0,dim=c(1,numSpecies))
    for (j in 1:numSpecies){
      Ab[1,j] <- sum(N_final[,,j]>0)/(180*180)
    }
    Abundance <- rbind(Abundance,Ab)
  }
  QTLS_scompD <- apply(X=SScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
  QTLS_scompF <- apply(X=AScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
  
  # Facilitation (weak)
  c[1,2] <- -0.00005
  c[2,1] <- -0.00005
  r[1] <- r[2] <- 0.24
  AScores <- array(data=NA,dim=c(length(divs),nns))
  SScores <- array(data=NA,dim=c(length(divs),nns))
  for (i in 1:100){
    N <- Initiate(q,numSpecies,3,300)
    N_final <- Simulate(N,q,r,c,D,dx,dt,endT,0)
    Xrab <- CalXrab(N_final,numSpecies,divs)
    Xab <- Xrab[1,]
    ra <- Xrab[2,]
    rb <- Xrab[3,]
    AScores[,i] <- (Xab/(Xab+ra))*(Xab/(Xab+rb))
    SScores[,i] <- apply(rbind(ra,rb),MARGIN=2,FUN=min)/apply(rbind(ra,rb),MARGIN=2,FUN=max)
    Ab <- array(0,dim=c(1,numSpecies))
    for (j in 1:numSpecies){
      Ab[1,j] <- sum(N_final[,,j]>0)/(180*180)
    }
    Abundance <- rbind(Abundance,Ab)
  }
  QTLS_wfacD <- apply(X=SScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
  QTLS_wfacF <- apply(X=AScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
  
  # Facilitation (strong)
  c[1,2] <- -0.0001
  c[2,1] <- -0.0001
  r[1] <- r[2] <- 0.24
  AScores <- array(data=NA,dim=c(length(divs),nns))
  SScores <- array(data=NA,dim=c(length(divs),nns))
  for (i in 1:nns){
    N <- Initiate(q,numSpecies,3,300)
    N_final <- Simulate(N,q,r,c,D,dx,dt,endT,0)
    Xrab <- CalXrab(N_final,numSpecies,divs)
    Xab <- Xrab[1,]
    ra <- Xrab[2,]
    rb <- Xrab[3,]
    AScores[,i] <- (Xab/(Xab+ra))*(Xab/(Xab+rb))
    SScores[,i] <- apply(rbind(ra,rb),MARGIN=2,FUN=min)/apply(rbind(ra,rb),MARGIN=2,FUN=max)
    Ab <- array(0,dim=c(1,numSpecies))
    for (j in 1:numSpecies){
      Ab[1,j] <- sum(N_final[,,j]>0)/(180*180)
    }
    Abundance <- rbind(Abundance,Ab)
  }
  QTLS_sfacD <- apply(X=SScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
  QTLS_sfacF <- apply(X=AScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
  
  # Dependence (weak)
  c[1,2] <- 0.001
  c[2,1] <- -0.00005
  r[1] <- 0.3
  r[2]<-0.08
  AScores <- array(data=NA,dim=c(length(divs),nns))
  SScores <- array(data=NA,dim=c(length(divs),nns))
  for (i in 1:nns){
    N <- Initiate(q,numSpecies,3,300)
    N_final <- Simulate(N,q,r,c,D,dx,dt,endT,0)
    Xrab <- CalXrab(N_final,numSpecies,divs)
    Xab <- Xrab[1,]
    ra <- Xrab[2,]
    rb <- Xrab[3,]
    AScores[,i] <- (Xab/(Xab+ra))*(Xab/(Xab+rb))
    SScores[,i] <- apply(rbind(ra,rb),MARGIN=2,FUN=min)/apply(rbind(ra,rb),MARGIN=2,FUN=max)
    Ab <- array(0,dim=c(1,numSpecies))
    for (j in 1:numSpecies){
      Ab[1,j] <- sum(N_final[,,j]>0)/(180*180)
    }
    Abundance <- rbind(Abundance,Ab)
  }
  QTLS_wdepD <- apply(X=SScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
  QTLS_wdepF <- apply(X=AScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
  
  # Dependence (strong)
  c[1,2] <- 0.001
  c[2,1] <- -0.0001
  r[1] <- 0.3
  r[2]<-0.08
  AScores <- array(data=NA,dim=c(length(divs),nns))
  SScores <- array(data=NA,dim=c(length(divs),nns))
  for (i in 1:nns){
    N <- Initiate(q,numSpecies,3,300)
    N_final <- Simulate(N,q,r,c,D,dx,dt,endT,0)
    Xrab <- CalXrab(N_final,numSpecies,divs)
    Xab <- Xrab[1,]
    ra <- Xrab[2,]
    rb <- Xrab[3,]
    AScores[,i] <- (Xab/(Xab+ra))*(Xab/(Xab+rb))
    SScores[,i] <- apply(rbind(ra,rb),MARGIN=2,FUN=min)/apply(rbind(ra,rb),MARGIN=2,FUN=max)
    Ab <- array(0,dim=c(1,numSpecies))
    for (j in 1:numSpecies){
      Ab[1,j] <- sum(N_final[,,j]>0)/(180*180)
    }
    Abundance <- rbind(Abundance,Ab)
  }
  QTLS_sdepD <- apply(X=SScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
  QTLS_sdepF <- apply(X=AScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
  
  Data<-data.frame(c(QTLS_neutralF),c(QTLS_scompF),c(QTLS_sdepF),c(QTLS_sfacF),c(QTLS_wcompF),c(QTLS_wdepF),c(QTLS_wfacF),
                   c(QTLS_neutralD),c(QTLS_scompD),c(QTLS_sdepD),c(QTLS_sfacD),c(QTLS_wcompD),c(QTLS_wdepD),c(QTLS_wfacD))
  save(Data,file='Data_3sp.rdata')
  write.csv2(Abundance,'Abundance_3sp.csv')
}

# Ten species neutral
endT <- 1500
{
  numSpecies <- 10
  Abundance <- array(0,dim=c(1,numSpecies))
  r <- array(0.3,dim=c(numSpecies)) # Growth rates for each species
  D <- array(0.1,dim=c(numSpecies)) # Diffusion coefficients
  # Competition array for neutral simulation
  c <- array(0.001,dim=c(numSpecies,numSpecies)) # competition array
  # Initiate empty arrays to store scores
  AScores <- array(data=NA,dim=c(length(divs),nns))
  SScores <- array(data=NA,dim=c(length(divs),nns))
  # Simulate for all repetitions
  for (i in 1:nns){
    # Initiate grid for simulation
    N <- Initiate(q,numSpecies,3,300)
    # Simulate and store final grid state
    N_final <- Simulate(N,q,r,c,D,dx,dt,endT,0)
    # Calculate Xab, ra, and rb
    Xrab <- CalXrab(N_final,numSpecies,divs)
    Xab <- Xrab[1,]
    ra <- Xrab[2,]
    rb <- Xrab[3,]
    # Calculate the scores
    AScores[,i] <- (Xab/(Xab+ra))*(Xab/(Xab+rb))
    SScores[,i] <- apply(rbind(ra,rb),MARGIN=2,FUN=min)/apply(rbind(ra,rb),MARGIN=2,FUN=max)
    # Store abundances of all species
    Ab <- array(0,dim=c(1,numSpecies))
    for (j in 1:numSpecies){
      Ab[1,j] <- sum(N_final[,,j]>0)/(180*180)
    }
    Abundance <- rbind(Abundance,Ab)
  }
  # Calculate quantiles for plotting
  QTLS_neutralD <- apply(X=SScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
  QTLS_neutralF <- apply(X=AScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
  
  # Competition (weak)
  c[1,2] <- 0.005
  c[2,1] <- 0.005
  AScores <- array(data=NA,dim=c(length(divs),nns))
  SScores <- array(data=NA,dim=c(length(divs),nns))
  for (i in 1:nns){
    N <- Initiate(q,numSpecies,3,300)
    N_final <- Simulate(N,q,r,c,D,dx,dt,endT,0)
    Xrab <- CalXrab(N_final,numSpecies,divs)
    Xab <- Xrab[1,]
    ra <- Xrab[2,]
    rb <- Xrab[3,]
    AScores[,i] <- (Xab/(Xab+ra))*(Xab/(Xab+rb))
    SScores[,i] <- apply(rbind(ra,rb),MARGIN=2,FUN=min)/apply(rbind(ra,rb),MARGIN=2,FUN=max)
    Ab <- array(0,dim=c(1,numSpecies))
    for (j in 1:numSpecies){
      Ab[1,j] <- sum(N_final[,,j]>0)/(180*180)
    }
    Abundance <- rbind(Abundance,Ab)
  }
  QTLS_wcompD <- apply(X=SScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
  QTLS_wcompF <- apply(X=AScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
  
  # Competition (strong)
  c[1,2] <- 0.01
  c[2,1] <- 0.01
  AScores <- array(data=NA,dim=c(length(divs),nns))
  SScores <- array(data=NA,dim=c(length(divs),nns))
  for (i in 1:nns){
    N <- Initiate(q,numSpecies,3,300)
    N_final <- Simulate(N,q,r,c,D,dx,dt,endT,0)
    Xrab <- CalXrab(N_final,numSpecies,divs)
    Xab <- Xrab[1,]
    ra <- Xrab[2,]
    rb <- Xrab[3,]
    AScores[,i] <- (Xab/(Xab+ra))*(Xab/(Xab+rb))
    SScores[,i] <- apply(rbind(ra,rb),MARGIN=2,FUN=min)/apply(rbind(ra,rb),MARGIN=2,FUN=max)
    Ab <- array(0,dim=c(1,numSpecies))
    for (j in 1:numSpecies){
      Ab[1,j] <- sum(N_final[,,j]>0)/(180*180)
    }
    Abundance <- rbind(Abundance,Ab)
  }
  QTLS_scompD <- apply(X=SScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
  QTLS_scompF <- apply(X=AScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
  
  # Facilitation (weak)
  c[1,2] <- -0.00005
  c[2,1] <- -0.00005
  r[1] <- r[2] <- 0.24
  AScores <- array(data=NA,dim=c(length(divs),nns))
  SScores <- array(data=NA,dim=c(length(divs),nns))
  for (i in 1:100){
    N <- Initiate(q,numSpecies,3,300)
    N_final <- Simulate(N,q,r,c,D,dx,dt,endT,0)
    Xrab <- CalXrab(N_final,numSpecies,divs)
    Xab <- Xrab[1,]
    ra <- Xrab[2,]
    rb <- Xrab[3,]
    AScores[,i] <- (Xab/(Xab+ra))*(Xab/(Xab+rb))
    SScores[,i] <- apply(rbind(ra,rb),MARGIN=2,FUN=min)/apply(rbind(ra,rb),MARGIN=2,FUN=max)
    Ab <- array(0,dim=c(1,numSpecies))
    for (j in 1:numSpecies){
      Ab[1,j] <- sum(N_final[,,j]>0)/(180*180)
    }
    Abundance <- rbind(Abundance,Ab)
  }
  QTLS_wfacD <- apply(X=SScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
  QTLS_wfacF <- apply(X=AScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
  
  # Facilitation (strong)
  c[1,2] <- -0.0001
  c[2,1] <- -0.0001
  r[1] <- r[2] <- 0.24
  AScores <- array(data=NA,dim=c(length(divs),nns))
  SScores <- array(data=NA,dim=c(length(divs),nns))
  for (i in 1:nns){
    N <- Initiate(q,numSpecies,3,300)
    N_final <- Simulate(N,q,r,c,D,dx,dt,endT,0)
    Xrab <- CalXrab(N_final,numSpecies,divs)
    Xab <- Xrab[1,]
    ra <- Xrab[2,]
    rb <- Xrab[3,]
    AScores[,i] <- (Xab/(Xab+ra))*(Xab/(Xab+rb))
    SScores[,i] <- apply(rbind(ra,rb),MARGIN=2,FUN=min)/apply(rbind(ra,rb),MARGIN=2,FUN=max)
    Ab <- array(0,dim=c(1,numSpecies))
    for (j in 1:numSpecies){
      Ab[1,j] <- sum(N_final[,,j]>0)/(180*180)
    }
    Abundance <- rbind(Abundance,Ab)
  }
  QTLS_sfacD <- apply(X=SScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
  QTLS_sfacF <- apply(X=AScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
  
  # Dependence (weak)
  c[1,2] <- 0.001
  c[2,1] <- -0.00005
  r[1] <- 0.3
  r[2]<-0.08
  AScores <- array(data=NA,dim=c(length(divs),nns))
  SScores <- array(data=NA,dim=c(length(divs),nns))
  for (i in 1:nns){
    N <- Initiate(q,numSpecies,3,300)
    N_final <- Simulate(N,q,r,c,D,dx,dt,endT,0)
    Xrab <- CalXrab(N_final,numSpecies,divs)
    Xab <- Xrab[1,]
    ra <- Xrab[2,]
    rb <- Xrab[3,]
    AScores[,i] <- (Xab/(Xab+ra))*(Xab/(Xab+rb))
    SScores[,i] <- apply(rbind(ra,rb),MARGIN=2,FUN=min)/apply(rbind(ra,rb),MARGIN=2,FUN=max)
    Ab <- array(0,dim=c(1,numSpecies))
    for (j in 1:numSpecies){
      Ab[1,j] <- sum(N_final[,,j]>0)/(180*180)
    }
    Abundance <- rbind(Abundance,Ab)
  }
  QTLS_wdepD <- apply(X=SScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
  QTLS_wdepF <- apply(X=AScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
  
  # Dependence (strong)
  c[1,2] <- 0.001
  c[2,1] <- -0.0001
  r[1] <- 0.3
  r[2]<-0.08
  AScores <- array(data=NA,dim=c(length(divs),nns))
  SScores <- array(data=NA,dim=c(length(divs),nns))
  for (i in 1:nns){
    N <- Initiate(q,numSpecies,3,300)
    N_final <- Simulate(N,q,r,c,D,dx,dt,endT,0)
    Xrab <- CalXrab(N_final,numSpecies,divs)
    Xab <- Xrab[1,]
    ra <- Xrab[2,]
    rb <- Xrab[3,]
    AScores[,i] <- (Xab/(Xab+ra))*(Xab/(Xab+rb))
    SScores[,i] <- apply(rbind(ra,rb),MARGIN=2,FUN=min)/apply(rbind(ra,rb),MARGIN=2,FUN=max)
    Ab <- array(0,dim=c(1,numSpecies))
    for (j in 1:numSpecies){
      Ab[1,j] <- sum(N_final[,,j]>0)/(180*180)
    }
    Abundance <- rbind(Abundance,Ab)
  }
  QTLS_sdepD <- apply(X=SScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
  QTLS_sdepF <- apply(X=AScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
  
  Data<-data.frame(c(QTLS_neutralF),c(QTLS_scompF),c(QTLS_sdepF),c(QTLS_sfacF),c(QTLS_wcompF),c(QTLS_wdepF),c(QTLS_wfacF),
                   c(QTLS_neutralD),c(QTLS_scompD),c(QTLS_sdepD),c(QTLS_sfacD),c(QTLS_wcompD),c(QTLS_wdepD),c(QTLS_wfacD))
  save(Data,file='Data_10sp.rdata')
  write.csv2(Abundance,'Abundance_10sp.csv')
}

# 8 species non-neutral (fake cold-water coral community)
endT <- 1500
numSpecies <- 8
r <- array(0.3,dim=c(numSpecies)) # Growth rates for each species
c <- array(0.001,dim=c(numSpecies,numSpecies)) # Neutral competition array
## The non-neutral competition coefficients
# Species 1 and 2...
# ... have weak comp with live coral (3) and rubble (7)
c[1,3] <- c[3,1] <- c[2,3] <- c[3,2] <- c[1,7] <- c[7,1] <- c[2,7] <- c[7,2] <- 0.005
# ...strongly depend on LDF (4) and MDF (5)
c[1,4] <- c[1,5] <- c[2,4] <- c[2,5] <- -0.0001; r[1] <- 0.1; r[2] <- 0.1
c[4,1] <- c[5,1] <- c[4,2] <- c[5,2] <- 0.001
# ... have strong competition with soft corals (6) and sponges (8)
c[1,6] <- c[2,6] <- c[1,8] <- c[2,8] <- c[6,1] <- c[6,2] <- c[8,1] <- c[8,2] <- 0.01
## 3 = live framework-forming corals...
# ...have weak facilitation with LDF (4) and MDF (5)
c[3,4] <- c[4,3] <- c[3,5] <- c[5,3] <- -0.00005; r[3]<-0.15
# ... have strong competition with soft corals (6) and rubble (7) and weak competition with sponges (8)
c[3,6] <- c[6,3] <- c[3,7] <- c[7,3] <- 0.01
c[3,8] <- c[8,3] <- 0.005
## 4 = LDF...
# ... has strong facilitation with MDF (5)
c[4,5] <- c[5,4] <- -0.0001;r[4] <- r[5] <- 0.25
## 6 = Soft coral...
# ... has strong dependence on LDF (4) and MDF (5)
c[6,4] <- c[6,5] <- -0.0001; r[6] <- 0.1
c[4,6] <- c[5,6] <- 0.001
# ...has weak competition with coral rubble (7)
c[6,7] <- c[7,6] <- 0.005
## 8 = Sponge...
# ...has strong competition with soft corals (6)
c[6,8] <- c[8,6] <- 0.01
# ...has weak competition with rubble (7)
c[7,8] <- c[8,7] <- 0.005
# ... weakly depends on LDF (4) and MDF (5)
c[8,4] <- c[8,5] <- -0.00005; r[8] <- 0.1
c[4,8] <- c[5,8] <- 0.001
{
  Abundance <- array(0,dim=c(1,numSpecies))
  r <- array(0.3,dim=c(numSpecies)) # Growth rates for each species
  D <- array(0.1,dim=c(numSpecies)) # Diffusion coefficients
  # Competition array for neutral simulation
  c <- array(0.001,dim=c(numSpecies,numSpecies)) # competition array
  # Initiate empty arrays to store scores
  AScores <- array(data=NA,dim=c(length(divs),nns))
  SScores <- array(data=NA,dim=c(length(divs),nns))
  # Simulate for all repetitions
  for (i in 1:nns){
    # Initiate grid for simulation
    N <- Initiate(q,numSpecies,3,300)
    # Simulate and store final grid state
    N_final <- Simulate(N,q,r,c,D,dx,dt,endT,0)
    # Calculate Xab, ra, and rb
    Xrab <- CalXrab(N_final,numSpecies,divs)
    Xab <- Xrab[1,]
    ra <- Xrab[2,]
    rb <- Xrab[3,]
    # Calculate the scores
    AScores[,i] <- (Xab/(Xab+ra))*(Xab/(Xab+rb))
    SScores[,i] <- apply(rbind(ra,rb),MARGIN=2,FUN=min)/apply(rbind(ra,rb),MARGIN=2,FUN=max)
    # Store abundances of all species
    Ab <- array(0,dim=c(1,numSpecies))
    for (j in 1:numSpecies){
      Ab[1,j] <- sum(N_final[,,j]>0)/(180*180)
    }
    Abundance <- rbind(Abundance,Ab)
  }
  # Calculate quantiles for plotting
  QTLS_neutralD <- apply(X=SScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
  QTLS_neutralF <- apply(X=AScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
  
  # Competition (weak)
  c[1,2] <- 0.005
  c[2,1] <- 0.005
  AScores <- array(data=NA,dim=c(length(divs),nns))
  SScores <- array(data=NA,dim=c(length(divs),nns))
  for (i in 1:nns){
    N <- Initiate(q,numSpecies,3,300)
    N_final <- Simulate(N,q,r,c,D,dx,dt,endT,0)
    Xrab <- CalXrab(N_final,numSpecies,divs)
    Xab <- Xrab[1,]
    ra <- Xrab[2,]
    rb <- Xrab[3,]
    AScores[,i] <- (Xab/(Xab+ra))*(Xab/(Xab+rb))
    SScores[,i] <- apply(rbind(ra,rb),MARGIN=2,FUN=min)/apply(rbind(ra,rb),MARGIN=2,FUN=max)
    Ab <- array(0,dim=c(1,numSpecies))
    for (j in 1:numSpecies){
      Ab[1,j] <- sum(N_final[,,j]>0)/(180*180)
    }
    Abundance <- rbind(Abundance,Ab)
  }
  QTLS_wcompD <- apply(X=SScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
  QTLS_wcompF <- apply(X=AScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
  
  # Competition (strong)
  c[1,2] <- 0.01
  c[2,1] <- 0.01
  AScores <- array(data=NA,dim=c(length(divs),nns))
  SScores <- array(data=NA,dim=c(length(divs),nns))
  for (i in 1:nns){
    N <- Initiate(q,numSpecies,3,300)
    N_final <- Simulate(N,q,r,c,D,dx,dt,endT,0)
    Xrab <- CalXrab(N_final,numSpecies,divs)
    Xab <- Xrab[1,]
    ra <- Xrab[2,]
    rb <- Xrab[3,]
    AScores[,i] <- (Xab/(Xab+ra))*(Xab/(Xab+rb))
    SScores[,i] <- apply(rbind(ra,rb),MARGIN=2,FUN=min)/apply(rbind(ra,rb),MARGIN=2,FUN=max)
    Ab <- array(0,dim=c(1,numSpecies))
    for (j in 1:numSpecies){
      Ab[1,j] <- sum(N_final[,,j]>0)/(180*180)
    }
    Abundance <- rbind(Abundance,Ab)
  }
  QTLS_scompD <- apply(X=SScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
  QTLS_scompF <- apply(X=AScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
  
  # Facilitation (weak)
  c[1,2] <- -0.00005
  c[2,1] <- -0.00005
  AScores <- array(data=NA,dim=c(length(divs),nns))
  SScores <- array(data=NA,dim=c(length(divs),nns))
  for (i in 1:100){
    N <- Initiate(q,numSpecies,3,300)
    N_final <- Simulate(N,q,r,c,D,dx,dt,endT,0)
    Xrab <- CalXrab(N_final,numSpecies,divs)
    Xab <- Xrab[1,]
    ra <- Xrab[2,]
    rb <- Xrab[3,]
    AScores[,i] <- (Xab/(Xab+ra))*(Xab/(Xab+rb))
    SScores[,i] <- apply(rbind(ra,rb),MARGIN=2,FUN=min)/apply(rbind(ra,rb),MARGIN=2,FUN=max)
    Ab <- array(0,dim=c(1,numSpecies))
    for (j in 1:numSpecies){
      Ab[1,j] <- sum(N_final[,,j]>0)/(180*180)
    }
    Abundance <- rbind(Abundance,Ab)
  }
  QTLS_wfacD <- apply(X=SScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
  QTLS_wfacF <- apply(X=AScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
  
  # Facilitation (strong)
  c[1,2] <- -0.0001
  c[2,1] <- -0.0001
  AScores <- array(data=NA,dim=c(length(divs),nns))
  SScores <- array(data=NA,dim=c(length(divs),nns))
  for (i in 1:nns){
    N <- Initiate(q,numSpecies,3,300)
    N_final <- Simulate(N,q,r,c,D,dx,dt,endT,0)
    Xrab <- CalXrab(N_final,numSpecies,divs)
    Xab <- Xrab[1,]
    ra <- Xrab[2,]
    rb <- Xrab[3,]
    AScores[,i] <- (Xab/(Xab+ra))*(Xab/(Xab+rb))
    SScores[,i] <- apply(rbind(ra,rb),MARGIN=2,FUN=min)/apply(rbind(ra,rb),MARGIN=2,FUN=max)
    Ab <- array(0,dim=c(1,numSpecies))
    for (j in 1:numSpecies){
      Ab[1,j] <- sum(N_final[,,j]>0)/(180*180)
    }
    Abundance <- rbind(Abundance,Ab)
  }
  QTLS_sfacD <- apply(X=SScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
  QTLS_sfacF <- apply(X=AScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
  
  # Dependence (weak)
  c[1,2] <- 0.001
  c[2,1] <- -0.00005
  AScores <- array(data=NA,dim=c(length(divs),nns))
  SScores <- array(data=NA,dim=c(length(divs),nns))
  for (i in 1:nns){
    N <- Initiate(q,numSpecies,3,300)
    N_final <- Simulate(N,q,r,c,D,dx,dt,endT,0)
    Xrab <- CalXrab(N_final,numSpecies,divs)
    Xab <- Xrab[1,]
    ra <- Xrab[2,]
    rb <- Xrab[3,]
    AScores[,i] <- (Xab/(Xab+ra))*(Xab/(Xab+rb))
    SScores[,i] <- apply(rbind(ra,rb),MARGIN=2,FUN=min)/apply(rbind(ra,rb),MARGIN=2,FUN=max)
    Ab <- array(0,dim=c(1,numSpecies))
    for (j in 1:numSpecies){
      Ab[1,j] <- sum(N_final[,,j]>0)/(180*180)
    }
    Abundance <- rbind(Abundance,Ab)
  }
  QTLS_wdepD <- apply(X=SScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
  QTLS_wdepF <- apply(X=AScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
  
  # Dependence (strong)
  c[1,2] <- 0.001
  c[2,1] <- -0.0001
  AScores <- array(data=NA,dim=c(length(divs),nns))
  SScores <- array(data=NA,dim=c(length(divs),nns))
  for (i in 1:nns){
    N <- Initiate(q,numSpecies,3,300)
    N_final <- Simulate(N,q,r,c,D,dx,dt,endT,0)
    Xrab <- CalXrab(N_final,numSpecies,divs)
    Xab <- Xrab[1,]
    ra <- Xrab[2,]
    rb <- Xrab[3,]
    AScores[,i] <- (Xab/(Xab+ra))*(Xab/(Xab+rb))
    SScores[,i] <- apply(rbind(ra,rb),MARGIN=2,FUN=min)/apply(rbind(ra,rb),MARGIN=2,FUN=max)
    Ab <- array(0,dim=c(1,numSpecies))
    for (j in 1:numSpecies){
      Ab[1,j] <- sum(N_final[,,j]>0)/(180*180)
    }
    Abundance <- rbind(Abundance,Ab)
  }
  QTLS_sdepD <- apply(X=SScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
  QTLS_sdepF <- apply(X=AScores,MARGIN=1,FUN=quantile,na.rm=TRUE)
  
  Data<-data.frame(c(QTLS_neutralF),c(QTLS_scompF),c(QTLS_sdepF),c(QTLS_sfacF),c(QTLS_wcompF),c(QTLS_wdepF),c(QTLS_wfacF),
                   c(QTLS_neutralD),c(QTLS_scompD),c(QTLS_sdepD),c(QTLS_sfacD),c(QTLS_wcompD),c(QTLS_wdepD),c(QTLS_wfacD))
  save(Data,file='Data_8sp.rdata')
  write.csv2(Abundance,'Abundance_8sp.csv')
}

## Create plots
load(file='Data_2sp.rdata')
Plots <- CreatePlot(Data,divs)
quartz(width=6,height=5)
print(Plots[[1]])
quartz.save('2s_Fweak.png',type="png",dpi=300)
quartz(width=6,height=5)
print(Plots[[2]])
quartz.save('2s_Fstrong.png',type="png",dpi=300)
quartz(width=6,height=5)
print(Plots[[3]])
quartz.save('2s_Dweak.png',type="png",dpi=300)
quartz(width=6,height=5)
print(Plots[[4]])
quartz.save('2s_Dstrong.png',type="png",dpi=300)

load(file='Data_3sp.rdata')
Plots <- CreatePlot(Data,divs)
quartz(width=6,height=5)
print(Plots[[1]])
quartz.save('3s_Fweak.png',type="png",dpi=300)
quartz(width=6,height=5)
print(Plots[[2]])
quartz.save('3s_Fstrong.png',type="png",dpi=300)
quartz(width=6,height=5)
print(Plots[[3]])
quartz.save('3s_Dweak.png',type="png",dpi=300)
quartz(width=6,height=5)
print(Plots[[4]])
quartz.save('3s_Dstrong.png',type="png",dpi=300)

load(file='Data_10sp.rdata')
Plots <- CreatePlot(Data,divs)
quartz(width=6,height=5)
print(Plots[[1]])
quartz.save('10s_Fweak.png',type="png",dpi=300)
quartz(width=6,height=5)
print(Plots[[2]])
quartz.save('10s_Fstrong.png',type="png",dpi=300)
quartz(width=6,height=5)
print(Plots[[3]])
quartz.save('10s_Dweak.png',type="png",dpi=300)
quartz(width=6,height=5)
print(Plots[[4]])
quartz.save('10s_Dstrong.png',type="png",dpi=300)

load(file='Data_8sp.rdata')
Plots <- CreatePlot(Data,divs)
quartz(width=6,height=5)
print(Plots[[1]])
quartz.save('8s_Fweak.png',type="png",dpi=300)
quartz(width=6,height=5)
print(Plots[[2]])
quartz.save('8s_Fstrong.png',type="png",dpi=300)
quartz(width=6,height=5)
print(Plots[[3]])
quartz.save('8s_Dweak.png',type="png",dpi=300)
quartz(width=6,height=5)
print(Plots[[4]])
quartz.save('8s_Dstrong.png',type="png",dpi=300)

# Create plot with legend
ggplot()+
  geom_line(aes(x=Scales,y=QTLS_Fneutral[3,],color='Neutral simulation'))+
  geom_point(aes(x=Scales,y=QTLS_Fneutral[3,],color='Neutral simulation'))+
  geom_line(aes(x=Scales,y=QTLS_Fscomp[3,],color='Competition simulation'))+
  geom_point(aes(x=Scales,y=QTLS_Fscomp[3,],color='Competition simulation'))+
  geom_line(aes(x=Scales,y=QTLS_Fsfac[3,],color='Facilitation simulation'))+
  geom_point(aes(x=Scales,y=QTLS_Fsfac[3,],color='Facilitation simulation'))+
  geom_line(aes(x=Scales,y=QTLS_Fsdep[3,],color='Dependence simulation'))+
  geom_point(aes(x=Scales,y=QTLS_Fsdep[3,],color='Dependence simulation'))+
  geom_line(aes(x=c(-1,1),y=c(1,1),color='Facilitation'),lwd=6,alpha=0.3)+
  geom_line(aes(x=c(-1,1),y=c(0.5,0.5),color='Dependence'),lwd=6,alpha=0.3)+
  geom_line(aes(x=c(-1,1),y=c(0.25,0.25),color='Neutral'),lwd=6,alpha=0.7)+
  geom_line(aes(x=c(-1,1),y=c(0,0),color='Competition'),lwd=6,alpha=0.3)+
  scale_color_manual(name='Legend',
                     breaks=c('Neutral simulation','Competition simulation','Facilitation simulation',
                              'Dependence simulation','Facilitation','Dependence','Neutral','Competition'),
                     values=c('Neutral simulation'='gray65','Competition simulation'='darkorange',
                              'Facilitation simulation'='darkolivegreen','Dependence simulation'='gold3',
                              'Facilitation'='darkolivegreen1',
                              'Dependence'='yellow','Neutral'='grey95','Competition'='orange'))

########################################################################################
