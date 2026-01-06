#########################################################################################
# This script calculates the Aggregation and Symmetry scores for the annotated images
# (c)2025, Anna van der Kaaden
#########################################################################################

# Set working directory to folder location of image files
setwd("~/")

# Required packages
require(png)
require(readxl)
require(ggplot2)
require(terra)
require(RColorBrewer)

#########################################################################################
# The following code is used to calculate Aggregation and Symmetry scores across spatial resolutions
# The 2nd, 3rd, and 4th quantile are then plotted across resolutions, leading to figure 4 in the paper
#########################################################################################

# The list of RGB colors and corresponding species names, from Maier et al. (2021) https://zenodo.org/record/4076147
Labels <- t(matrix(c(255,0,255,0,0,255,0,255,0,251,175,93,255,255,0,255,255,255),
                 nrow=3,ncol=6))
LabNames <- c('Coral','LDF','MDF','Soft coral','Rubble','Sponge')

# Create a list of names of annotated images to loop over
Files <- list.files()

# The number of segments into which each image is divided along the x or y axis
# A division of 2 means that image is aggregated into 2x2=4 blocks
divs <- c(2,4,5,8,10,16,20,40,80)

# The size of the blocks is calculated assuming an average image size of 2.8m2,
# and an aspect ratio of 9/16, leading to a length and width of 2.23m and 1.26m
# We take the diameter of the aggregated blocks as the spatial length scale (m)
Scales <- sqrt((2.231/divs)^2 + (1.255/divs)^2)

# The arrays in which the data will be stored, 
# with 167 images, 9 length scales, and 6 species
# AScore will hold the Aggregation scores and SScore the Symmetry scores
AScore <- array(data=NA,dim=c(167,9,6,6))
SScore <- array(data=NA, dim=c(167,9,6,6))

# In this loop, we first load in the images, then convert them into species presence and absence
# Then we calculate the scores for each spatial scale and each species pair
# The spatial scales are created by aggregating the data into coarser resolutions
for (i in 1:167){
  # Read in an image from the list of files
  img <- round(255*readPNG(Files[i]))
  
  # Remove the alpha channel
  img <- img[,,-4]
  
  # Get image dimensions
  dimI <- dim(img)
  
  # Change to 2D image
  img <- matrix(img, ncol=3)
  
  # Create an empty species matrix to fil
  Species <- matrix(data=NA,nrow=dim(img)[1],ncol=1)
  
  #In this loop, we will in the species presence from the RGB labels
  for (j in 1:6){
    # Fill species matrix based on RGB codes and labels
    Species[apply(img, 1, identical, Labels[j,])]=j
  }
  
  # Convert species matrix into 2D
  Species <- matrix(Species, nrow=dimI[1], ncol=dimI[2])
  
  # In this loop, we calculate the scores 
  # over all spatial scales
  for (j in 1:9){
    # and for all species pairs
    # Note that we call the species a and b in the script and i and j in the paper
    for (a in 1:5){
      for (b in 2:6){
        if (a!=b){
          # Create new array for species a
          A <- array(0,dim(Species))
          
          # Fill in where species a is present, indicated by ones
          A[Species==a] <- 1
          
          # Aggregate the array into blocks as stipulated in 'divs'
          x <- aggregate(rast(A),c(dimI[1]/divs[j],dimI[2]/divs[j]),max)
          
          # Convert array to matrix
          A <- as.matrix(x,dim(x)[1],dim(x)[2])
          
          # Create a new array for species b
          B <- array(0,dim(Species))
          
          # Fill in where species b is present, indicated by ones
          B[Species==b] <- 1
          
          # Aggregate the array into blocks as stipulated in 'divs'
          x <- aggregate(rast(B),c(dimI[1]/divs[j],dimI[2]/divs[j]),max)
          
          # Convert array to matrix
          B <- as.matrix(x,dim(x)[1],dim(x)[2])
          
          # Here, we calculate the scores
          # Calculate ra = number of blocks containing species a and not b
          ra <- sum(A==1&B==0)
          
          # Calculate rb = number of blocks containing species b and not a
          rb <- sum(B==1&A==0)
          
          # Calculate Xab = number of blocks containing both species
          Xab <- sum(A==1&B==1)
          
          # Calculate A-score and S-score
          AScore[i,j,a,b] <- (Xab/(Xab+ra))*(Xab/(Xab+rb))
          SScore[i,j,a,b] <- min(c(ra,rb))/max(c(ra,rb))
          
          # If one of the species is not present, A-scores will be NaN
          # In which case the S-score should also be NaN
          if (is.na((Xab/(Xab+ra))*(Xab/(Xab+rb)))){
            SScore[i,j,a,b] <- NaN
          }
        }
      }
    }
  }
}

# Save the scores in a file
save(AScore,file='../../Scripts/AScore_CWC.rdata')
save(SScore,file='../../Scripts/SScore_CWC.rdata')

# Load the files to create the figures
load(file='../../Scripts/AScore_CWC.rdata')
load(file='../../Scripts/SScore_CWC.rdata')

# Here, for each species pair, we first calculate the quantiles of the scores
# Then we plot the second (25th), third (50th), and fourth (75th) quantile
# together with lines for the interactions belonging to score values in ideal situations
# as argued in the theoretical framework development (section 2) in the paper 
for (j in 1:5){
  for (k in 2:6){
    if (j<k){
      # Calculate the quantiles
      QTLS_F <- apply(X=AScore[,,j,k],MARGIN=2,FUN=quantile,na.rm=TRUE)
      QTLS_D <- apply(X=SScore[,,j,k],MARGIN=2,FUN=quantile,na.rm=TRUE)
      
      # Plot a polygon around the 2nd and 4th quantile, 
      # and plot the 3rd quantile (median) as a line,
      # S-scores in green and A-scores in blue
      # Plot the interaction types, 
      # orange=competition, grey=neutral, yellow=dependence, and green=facilitation
      p <- ggplot()+
        geom_line(aes(x=c(-1,2),y=c(0,0)),color='orange',lwd=6,alpha=0.3)+
        geom_line(aes(x=c(-1,2),y=c(0.25,0.25)),color='grey95',lwd=6,alpha=0.7)+
        geom_line(aes(x=c(-1,2),y=c(0.5,0.5)),color='yellow',lwd=6,alpha=0.3)+
        geom_line(aes(x=c(-1,2),y=c(1,1)),color='darkolivegreen1',lwd=6,alpha=0.3)+
        geom_polygon(aes(x=c(Scales,rev(Scales)),y=c(QTLS_Ts[2,],rev(QTLS_Ts[4,]))),fill='palegreen3',alpha=0.3)+
        geom_line(aes(x=Scales,y=QTLS_Ts[3,]),color='darkgreen')+
        geom_point(aes(x=Scales,y=QTLS_Ts[3,]),color='darkgreen')+
        geom_polygon(aes(x=c(Scales,rev(Scales)),y=c(QTLS_T[2,],rev(QTLS_T[4,]))),fill='lightsteelblue',alpha=0.3)+
        geom_line(aes(x=Scales,y=QTLS_T[3,]),color='blue')+
        geom_point(aes(x=Scales,y=QTLS_T[3,]),color='blue')+
        # This sets the axis limits
        coord_cartesian(ylim=c(-.25,1.25),xlim=c(0,1.3))+
        # The title indicates the species pair
        ggtitle(paste(as.character(j),as.character(k),sep=","))+
        # Axis labels
        labs(x="Spatial length-scale (m)",y="Scores")+
        # Create the legend
        scale_color_manual(name='Legend',
                         breaks=c('A-score','S-score','Facilitation','Dependence','Neutral','Competition'),
                         values=c('A-score'='blue','S-score'='darkgreen','Facilitation'='darkolivegreen1',
                                  'Dependence'='yellow','Neutral'='grey95','Competition'='orange'))
      
      # For mac OS, this creates a window for plotting
      quartz(width=5,height=5)
      # Print the figure
      print(p)
      # For mac OS, this saves the figure
      quartz.save(paste('Scores_CWC',as.character(j),as.character(k),'.png',sep=""),type="png",dpi=300)
    }
  }
}

#########################################################################################
# The following code calculate Aggregation and Symmetry scores for the entire images
# It plots the co-occurrence plots along with the scores (figure 5 in the paper)
#########################################################################################

# Read in the video data as provided by Maier et al. (2021), https://zenodo.org/record/4076147
video_data <- read_xlsx('../video_data.xlsx')
# Create empty arrays to store the Aggregation (F-) and Symmetry (D-) scores
AScore_whole <- array(data=NA,dim=c(6,6))
SScore_whole <- array(data=NA,dim=c(6,6))

# An array containing the percent cover of the 6 species of interest
Specs_whole <- array(c(video_data$coral,video_data$less_degr,video_data$more_degr,video_data$other_coral,
                       video_data$sed_rubble,video_data$erect_white_sponge),dim=c(167,6))

# Here, we loop over the species pairs to calculate the scores
for (a in 1:5){
  for (b in 2:6){
    if (a!=b){
      # Calculate the scores, note that here we refer to species a and b and in the paper to species i and j
      # Calculate ra = number of blocks containing species a and not b
      ra <- sum(Specs_whole[,a]>0 & Specs_whole[,b]==0)
      # Calculate rb = number of blocks containing species b and not a
      rb <- sum(Specs_whole[,a]==0 & Specs_whole[,b]>0)
      # Xab = number of blocks containing both species
      Xab <- sum(Specs_whole[,a]>0 & Specs_whole[,b]>0)
      # n = number of blocks containing species a, b, or both
      n <- sum(Specs_whole[,a]>0 | Specs_whole[,b]>0 | (Specs_whole[,a]>0 & Specs_whole[,b]>0))
      
      # Calculate the scores
      AScore_whole[a,b] <- (Xab/(Xab+ra))*(Xab/(Xab+rb))
      SScore_whole[a,b] <- min(c(ra,rb))/max(c(ra,rb))
      
      # If one of the species is not present, A-scores will be NaN
      # In which case the S-score should also be NaN
      if (is.na((Xab/(Xab+ra))*(Xab/(Xab+rb)))){
        SScore[i,j,a,b] <- NaN
      }
    }
  }
}

# Save the scores to a file
save(AScore_whole,file='../../Scripts/AScore_whole_CWC.rdata')
save(SScore_whole,file='../../Scripts/SScore_whole_CWC.rdata')

# Load the files to create the figures
load(file='../../Scripts/AScore_whole_CWC.rdata')
load(file='../../Scripts/SScore_whole_CWC.rdata')

# For all species pairs, plot co-occurrence plots with scores
# We color the window according to the interaction type
# Creat the color scale
Colorscale <- colorRampPalette(colors=c('orange','grey95','yellow','darkolivegreen1'))(16)
# Create the color scale with alpha-level
Colorscale_alpha <- sapply(Colorscale, function(col) rgb(red=col2rgb(col)[1]/255,
                                                     green=col2rgb(col)[2]/255,
                                                     blue=col2rgb(col)[3]/255,
                                                     alpha=0.3))
# Create the A-score values corresponding to the colors
ColValues <- c(seq(0,0.5,length.out=11),seq(0.6,1,length.out=5))

# Here, we plot the co-ocurrence plots for all species pairs, 
# as scatter plots using the % cover of both species from the video data
for (a in 1:5){
  for (b in 2:6){
    if (a<b){
      # For mac OS, create a quartz window for plotting
      quartz(width=4,height=4)
      
      # Create the ggplot scatter plot
      p<-ggplot(data=video_data,aes(x=Specs_whole[,b],y=Specs_whole[,a]))+
        geom_point(size=1,colour='black')+
        theme(panel.background = element_rect(fill=Colorscale_alpha[which.min(abs(ColValues-AScore_whole[a,b]))]),
              panel.grid.major=element_line(color='grey',size=.5))+
        labs(x="Cover (%)",y="Cover (%)")
      
      # Print the figure
      print(p)
      
      # Save the figure
      quartz.save(paste('CoOccurPlot_CWC',as.character(a),as.character(b),'.png',sep=""),type="png",dpi=300)
    }
  }
}

#########################################################################################
# Alternatively, the co-ocurrence plots can be plotted as bins instead of scatter plots
for (a in 1:5){
  for (b in 2:6){
    if (a<b){
      p<-ggplot(data=video_data,aes(x=Specs_whole[,b],y=Specs_whole[,a]))+
        geom_bin2d(binwidth=2,aes(fill=after_stat(ncount)))+
        lims(x=c(-10,100),y=c(-10,100))+
        ggtitle(paste(as.character(a),as.character(b),as.character(TScore_whole[a,b],sep=",")))+
        scale_fill_gradient(name="percent",labels=scales::percent,limits=c(0,1))
      print(p)
    }
  }
}
