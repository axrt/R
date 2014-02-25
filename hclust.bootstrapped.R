#Allows to cluster and bootstrap a given data.frame of eucledian coordinates.
#Param: data - data.frame of eucledian coordinates
#Param: method, that can be used for clustering (by the hclust() function), "average" by default
#Param: deg.freedom - allowed variation between the point coordinates, given as a percentage of the median distance between the point
#Param: num.bootstraps - a number of bootstrap iterations
#Param: create.plot - TRUE by default, plots the clustring and lays out the bootstrap values, switched off if FALSE
#Param: col - colour of the markings (pch 21 - default) that indicate the plot vertexes
#Param: pch - 21 by default, a pictogram that indicates the plot vertexes
#Returns a data.frame of bootstrap values, assigned to evely $merge field of an hclust structure
hclust.bootstrapped<-function(data, method="average", deg.freedom, num.bootstraps, create.plot=TRUE, col="red", pch=21){
  
  #FUNCTIONS
  #-----------------------------------------------------------------------------------------------------
  
  #Converts a hclust sctructure to a sparse matrix, 
  #where 1 indicates an existing connection and 0 - no connecton.
  #Note: Nodes are not connected to themselves and leaves 
  #are not connected to any other nodes.
  #Returns a matrix of connections between nodes
  convertHclustToMatrix<-function(hclust.data){
    
    if(!is(hclust.data, "hclust")){
      return(print(paste("ERR: Found:",class(hclust.data),"where hclust expected!",sep=" ")))
    }
    
    
    #Sort the $merge field accending to separeate the joins from leaves
    #Note: leaves in the hclust structure have negative values, 
    #nodes, which are not leaves, and indicate joins between nodes
    sorted<-sort(hclust.data$merge,decreasing=FALSE)
    
    #Separate the nodes from leaves
    nodes <- sorted[sorted>0 & !is.nan(sorted)]
    leaves <- sorted[sorted<0 & !is.nan(sorted)]
    
    #Create a blank matrix template to fit the tree ($merge lengt * $merge length)
    #!!!Unsure is this is a bad way to do things in R!!!
    matrix<-matrix(rep(0,length(hclust.data$merge)^2),nrow=length(hclust.data$merge))
    #assign the names to rows/columns in the matrix
    rownames(matrix)<-sorted
    colnames(matrix)<-sorted
    
    #Fill in the matrix by going through the node joins, 
    #finding 
    for(i in 1:length(nodes)){
      pair<-hclust.data$merge[nodes[i],]
      matchInMatrix<-match(pair,sorted) #This is to find the proper position in matrix
      index<-match(nodes[i],sorted)
      matrix[index,matchInMatrix[1]]<-matrix[index,matchInMatrix[1]]+1 
      matrix[matchInMatrix[1],index]<-matrix[matchInMatrix[1],index]+1
      matrix[index,matchInMatrix[2]]<-matrix[index,matchInMatrix[2]]+1
      matrix[matchInMatrix[2],index]<-matrix[matchInMatrix[2],index]+1
    }
    return(matrix)
  }
  #-----------------------------------------------------------------------------------------------------
  
  #Generates a given number of randomized dist matrices, primed by the given matrix.
  #Returns a list of newly created randomized matrices
  generateRandomizedDistance<-function(data,number,deg.freedom){
    
    #Prepare a list to hold the matrices
    #!!!Unsure is this is a bad way to do things in R!!!
    list<-list()
    
    #Fill in the list
    for(i in 1:number){
      list[[length(list)+1]]<-randomizeDistTable(data,deg.freedom)
      #!!!Could not use c(..,..) here because the output will not be a list... !!!
    }
    return (list)
  }
  #-----------------------------------------------------------------------------------------------------
  
  #Generates a randomized distance matrix.
  #Any point within the matrix can get as much randomization in a position as allowed by the deg.freedom parameter
  #allows.
  #Returns a randomized data.frame
  randomizeDistTable<-function(data, deg.freedom){
    
    #Copy the data.frame 
    new.data<-data
    for(i in 1:nrow(new.data)){
      for(j in 1:ncol(new.data)){
        new.data[i,j]<-new.data[i,j]+runif(1,-deg.freedom,deg.freedom)
        #!!!This generates tonns of warnings !!!
      }
    }
    return(data.frame(new.data))
  }
  #-----------------------------------------------------------------------------------------------------
  
  #Canculates a pair of x and y values to correctly indicate the vertex on the hclust plot
  #Returns a pair of x,y values of coordinates
  vertexMarkUpCoordinate<-function(x,y,hclust.struct){
    #If both values <0 then both are leaves and can be mapped directly
    if(x<0&&y<0){
      positionX=match(-x,hclust.struct$order)
      positionY=match(-y,hclust.struct$order)
      return(c(positionX,positionY))
      #Else if one of the values in a pair >0, then it it point to a join and the other is a leaf
    }else if(x<0&&y>0){
      positionX=match(-x,hclust.struct$order)
      #Find the proper position recursively
      positions<-vertexMarkUpCoordinate(hclust.struct$merge[y,1],hclust.struct$merge[y,2],hclust.struct)
      positionY<-(positions[1]+positions[2])/2
      return(c(positionX,positionY))
      #Else if one of the values in a pair >0, then it it point to a join and the other is a leaf
    }else if(x>0&&y<0){
      positionY=match(-y,hclust.struct$order)
      #Find the proper position recursively
      positions<-vertexMarkUpCoordinate(hclust.struct$merge[x,1],hclust.struct$merge[x,2],hclust.struct)
      positionX<-(positions[1]+positions[2])/2
      return(c(positionX,positionY))
      #If both are joins - find both coordinates recursively
    }else{
      #Find the proper position recursively
      positions<-vertexMarkUpCoordinate(hclust.struct$merge[x,1],hclust.struct$merge[x,2],hclust.struct)
      positionX<-(positions[1]+positions[2])/2
      #Find the proper position recursively
      positions<-vertexMarkUpCoordinate(hclust.struct$merge[y,1],hclust.struct$merge[y,2],hclust.struct)
      positionY<-(positions[1]+positions[2])/2
      return(c(positionX,positionY))
    }
  }
  #-----------------------------------------------------------------------------------------------------
  
  #Creates a data.frame, that has contains the initial coordinates data as well as the bootstrap values for the nodes
  #Returns the newly created data.frame
  hclustToDataFrameWithBootstraps<-function(hclust.struct, bootstrap.mtx){
    #Prepare a data.frame
    #!!!Unsure is this is a bad way to do things in R!!!
    frame<-data.frame(hclust.struct$merge)
    one<-list()
    two<-list()
    names1<-list()
    names2<-list()
    for(j in 1:nrow(frame)){
      one<-c(one,bootstrap.mtx[match(frame[j,1],row.names(bootstrap.mtx)),match(j,row.names(bootstrap.mtx))])
      two<-c(two,bootstrap.mtx[match(frame[j,2],row.names(bootstrap.mtx)),match(j,row.names(bootstrap.mtx))])
      if(frame[j,1]<0){
        names1<-c(names1,hclust.struct$labels[frame[j,1]*(-1)])
        
      }else{
        names1<-c(names1,frame[j,1])
      }
      if(frame[j,2]<0){
        names2<-c(names2,hclust.struct$labels[frame[j,2]*(-1)])
        
      }else{
        names2<-c(names2,frame[j,2])
      }
      
    }
    frame$bootstraps1<-c(one[1:length(one)-1],1)
    frame$bootstraps2<-c(two[1:length(two)-1],1)
    frame$names1<-names1
    frame$names2<-names2
    
    return(frame)
  }
  #-----------------------------------------------------------------------------------------------------
  #END FUNCTIONS
  
  
  #Check whether the dist.df is a data frame
  if(!is(data,"data.frame")){ 
    return(print(paste("ERR: Found:",class(data),"where data.frame expected!",sep=" ")))
  }
  #Check whether number of bootstraps is non-negative
  if(num.bootstraps<2){
    return(print("ERR: number of bootstraps must be >2"))
  }
  #Check whether the degree of freedom is within [0,1]
  if(deg.freedom<0||deg.freedom>1){
    return(print("ERR: The degree of freedom must be within [0:1]"))
  }else{
    #If a correct value was given - convert to a percent of the median-remote dots distance
    dist.data<-dist(data)
    deg.freedom<-median(dist.data)*deg.freedom
  }
  
  #Get the initial hclust structure
  print("Generating prime hclust")#TODO:comment out
  prime.hclust<-hclust(dist(x=data),method=method)
  
  #Convert the initial hclust structure to a sparse matrix of nodes
  print("Converting hclust to bootstrap matrix")#TODO:comment out
  prime.mtx<-convertHclustToMatrix(prime.hclust)
  
  #Randomize distance matrix allowing a given degree of freedom
  print("Randomizing the distance matrix")#TODO:comment out
  rand.mtxs<-generateRandomizedDistance(data,num.bootstraps,deg.freedom)
  
  #Cluster all the randomly generated distance data.frames and sum up in a new data.frame
  for(i in 1:length(rand.mtxs)){
    mtx<-convertHclustToMatrix(hclust(dist(rand.mtxs[[i]]),method=method))
    #If the matrix is bigger than the prime matrix - then sum up using the length
    #of the prime matrix,
    if(nrow(mtx)>=nrow(prime.mtx)){
      for(j in 1:ncol(prime.mtx)){
        for(k in 1:nrow(prime.mtx)){
          prime.mtx[j,k]<-prime.mtx[j,k]+mtx[j,k]
        }
      }
      #else - use the length of the newly created matrix
    }else{
      for(j in 1:ncol(mtx)){
        for(k in 1:nrow(mtx)){      
          prime.mtx[j,k]<-prime.mtx[j,k]+mtx[j,k]
        }
      }
    }
  }
  #Finally normalize by the number of bootstraps to get the bootsrap values
  for(i in 1:nrow(prime.mtx)){
    for(j in 1:ncol(prime.mtx)){
      prime.mtx[i,j]<-prime.mtx[i,j]/(num.bootstraps+1)
    }
  }
  
  #Compile a bootstrapped data.frame
  boot.df<-hclustToDataFrameWithBootstraps(prime.hclust,prime.mtx)
  
  #Finally plot the prime hclust and lay out the bootstrap values
  if(create.plot){
  plot(prime.hclust, hang=-1)
  for(i in length(prime.hclust$merge[,1]):1){
    x<-prime.hclust$merge[i,1]
    y<-prime.hclust$merge[i,2]
    positions<-vertexMarkUpCoordinate(x,y,prime.hclust)
    #Add ovals (pch 21) to indicate the vertex points
    points(positions[1],prime.hclust$height[i],col=col, pch=pch)
    points(positions[2],prime.hclust$height[i],col=col, pch=pch)
    text(positions[1],prime.hclust$height[i],labels=format(x=boot.df[i,3],digits=2),adj=c(0.5,-0.2), 
         col=rgb(255-255*as.numeric(boot.df[i,3]), 100, 255*as.numeric(boot.df[i,3]), 255, maxColorValue = 255),srt=90)
    text(positions[2],prime.hclust$height[i],labels=format(x=boot.df[i,4],digits=2),adj=c(0.5,1.2), 
         col=rgb(255-255*as.numeric(boot.df[i,4]), 100, 255*as.numeric(boot.df[i,4]), 255, maxColorValue = 255),srt=90)
  }
  }
  return(boot.df)
}