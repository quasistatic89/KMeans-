#' K Means++ Algorithm for Movie Data using a subset of Covariates
#' to get optimal clustering

library(dplyr)
library(ggplot2)

D <- read.csv("Dataset_Revenue_v3.0.csv")
D1 <- subset(D,(!(is.na(Revenue))) & (!(is.na(Budget))) & Year > 1980) # Only considering movies after 1980

# Set of tunable paramters for the model
Dimension = 13
Num_Clust = 3
Num_Data_points = nrow(D1)

drops <- c("Revenue","MovieID","Movie",
           "MAge1","MAge7","FAge1","FAge7","MAge6",
           "FAge5","FAge6",
           "MAge2","MAge3", "MAge4", "MAge5", "FAge4","FAge2","FAge3") 

D2 <- D1[,!(names(D1) %in% drops)]  # Dropping the columns that will not be used in clustering or in the final analysis

#'Scaling is neccesary as different features have very different means and variances
Scaled_D <- scale(D2)   

#Initial_Center_using K means ++
Careful_Seeding <- function(){ 
    C <- matrix(NA,Num_Clust,Dimension)
    p=0
    while(p==0){
        t <- sample(1:Num_Data_points,1) #for 1st centre
        if (sum(is.na(Scaled_D[t,])) == 0){
            C[1,] <- Scaled_D[t,]      #Pick 1st Centre randomly
            p = 1
        }
    }
    for (i in 2:Num_Clust){
        Clusters <- Cluster(matrix(C[1:i-1,],ncol=Dimension)) # Find cluster that each point belongs to
        distance <- vector(length=Num_Data_points)
        for (j in 1:Num_Data_points){
            distance[j] <- Euc_Dist(Scaled_D[j,],C[Clusters[j],])^2
        }
        distance <- distance/(sum(distance)) # Normalizing the probabilities
        
        p=0
        while(p==0){
            t <- sample(1:Num_Data_points, size = 1, prob = distance) #for 1st centre
            if (sum(is.na(Scaled_D[t,])) == 0){
                C[i,] <- Scaled_D[t,]      #Next Centre
                p = 1
            }
                
        }
    }
    return(C)
}

# Euclidean Distance
Euc_Dist <- function(x1,x2){
    diff = sum(((x1-x2)^2),na.rm= TRUE)
    len = sqrt(length(na.omit(x1-x2)))
    return(sqrt(diff/len))
}


# Function to find the cluster for each data point
Cluster <- function(Centre){
    Clusters <- numeric(length = Num_Data_points)
    for (i in 1:Num_Data_points){
        d <- vector(length=nrow(Centre))
        for (m in 1:nrow(Centre)) {
            d[m] <- Euc_Dist(Scaled_D[i,],Centre[m,]) 
        }
        Clusters[i] <- which.min(d)
    }
    return (Clusters)
}

# Function to update Centres
Update_Centres <- function(Current_Center){
    # Using the Cluster function to calculate the closest center for each point in D
    Clusters <- Cluster(Current_Center)
    C <- matrix(NA,Num_Clust,Dimension)
    
    for (i in 1:Num_Clust){
        indices <- which(Clusters==i)
        for (j in 1:Dimension){
            C[i,j] <- mean(Scaled_D[indices,j],na.rm = TRUE)
            #C[i,2] <- mean(D[indices,2])
        }
    }
    return (C)
}

#Function to calculate the cost 
Cost_function <- function(centre){
    Clusters <- Cluster(centre)
    cost = 0
    for ( i in 1:Num_Data_points){
        cost =cost+Euc_Dist(Scaled_D[i,],centre[Clusters[i],])
    }
    return (cost)
}

# Lloyds function
Centres <- list() # storing centres for all iterations
Centres[[1]] <- Careful_Seeding() # Will change depending on whether we use kmeans or kmeans++
iterations = 0
repeat{
    iterations= iterations+1
    Centres[[iterations+1]] <- Update_Centres(Centres[[iterations]])
    if(sum(Centres[[iterations+1]] == Centres[[iterations]]) == (Num_Clust*Dimension)){
        break
    }
}

#Distance of centres at each iteration from final centres and the cost function
Cent_Dist <- matrix(NA,iterations,Num_Clust)
Cost <- vector(length=iterations)
for (i in 1:iterations){
    for (j in 1:Num_Clust){
        Cent_Dist[i,j] <- Euc_Dist(Centres[[i]][j,],Centres[[iterations+1]][j,])
    }
    Cost[i] = Cost_function(Centres[[i]])
}

#Plotting
par(mfrow=c(1,1))
matplot(1:iterations,Cent_Dist,type='l',lty=1,
        xlab="Number of Iterations",ylab="Distance",
        main="Dist of centres from final centres",col=1:Num_Clust,ylim=c(-0.5,4))
legend("topright", c("C1", "C2","C3"),pch = 5, col = 1:Num_Clust,cex=0.8)
plot(Cost,xlab="Number of Iterations",ylab="Cost",main="Cost function for k-means++",ylim=c(2000,2700),type='l',col='blue')

#Saving Final Centres and Cluster ID's for each point
Final_Centres = Centres[[iterations+1]]
Final_ClusterID = Cluster(Final_Centres)
Returns <- (D1$Revenue - D1$Budget)/D1$Budget
df <- data.frame( Returns = Returns, ClusterID= as.factor(Final_ClusterID))

#Plotting Box Plots for returns for movies (grouped by clusters)
Sum1 <- (summarise(group_by(df, ClusterID), N=n()))
par(mfrow=c(1,1))
boxplot(Returns~ClusterID,data=df, main="Performance of Different Clusters", 
        xlab="Cluster ID's", ylab="ROI (Profit/Budget)",outline=FALSE,ylim=c(-2,12))

# Summaries for each group
df1 <- D1
df1$ClusterID = df$ClusterID

Sum2 <- summarise(group_by(df1,ClusterID),N=n(),Perc_Actors = mean(Top.30.Actors),
                  Perc_Actress = mean(Top.30.Actresses), Perc_Directors = mean(Top.30.Directors),
                  Perc_Prod = mean(Major.Production.House),
                  Perc_Sequel = mean(Sequel),
                  Mean_Budget = mean(Budget),
                  Ratings = mean(MAge1+MAge2+ MAge3+ MAge4 + MAge5+ MAge6+MAge7+FAge1+
                                 FAge2 +FAge3+ FAge4+FAge5+FAge6+FAge7,na.rm = TRUE)/14,
                  Avg_Year = mean(Year))

par(mfrow=c(1,5))
barplot(Sum2$N, horiz=FALSE,ylab="No. of data points",ylim=c(0,800),xpd=FALSE,cex.lab=1.5,col='lightblue')
barplot(as.vector(Sum2$Perc_Actors), horiz=FALSE,ylim=c(0,2),ylab="% of Movies with Top Actors",cex.lab=1.5,col='lightblue')
barplot(as.vector(Sum2$Perc_Sequel), horiz=FALSE,ylim=c(0,.20),ylab="% of Sequel Movies ",cex.lab=1.5,col='lightblue')
barplot(as.vector(Sum2$Mean_Budget), horiz=FALSE,ylab="Mean Budget of Movies (Millions)",ylim=c(10000000,50000000),xpd=FALSE,cex.lab=1.5,col='lightblue')
barplot(as.vector(Sum2$Ratings),ylim=c(3,4),horiz=FALSE,ylab="Mean Rating of Movies ",xpd=FALSE,cex.lab=1.5,col='lightblue')
