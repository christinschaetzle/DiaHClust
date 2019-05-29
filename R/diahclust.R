
#vnc function implements the ideas of the "Variability-based Neighbor Clustering" (VNC) algorithm developed by Gries and Hilpert (2008, 2012)
#d has to be a distance matrix
vnc=function(d, method=c("single","complete","average","median", "ward.D", "ward.D2", "mcquitty", "centroid")) 
{
  vnc.dist=distvnc(d)
  vnc.clust=hclust(vnc.dist, method=match.arg(method))
  vnc.clust$order=rev(-sort(vnc.clust$order)) #get diachronic order of labels
  return(vnc.clust)
}

#function which manipulates the distance matrix for VNC clustering via hclust()
distvnc=function(d)
{
  #manipulation of distance matrix to only allow for the clustering of temporally adjacent time periods
  if (!is.matrix(d)) d = as.matrix(d)
  
  N = nrow(d)
  C = ncol(d)
  z=max(d) #z is set to the value which will be assigned to values depicting the distance between temporally non-adjacent time periods
  
  #restrict possibilities for merge by assigning highest values to non-allowed combinations 
  
  for (k in seq(1,N))
  {
    l=k+2
    for (l in seq(l,C))
      if (l<N) d[k,l]=z
  }
  
  for (k in seq(1,N))
  {
    l=k-2
    
    
    if (l>0)
      for (l in seq(0,l))
        d[k,l]=z
  }
  
  for (k in seq(1,N-2))
  {
    d[k,C]=z
  }
  
  d = as.dist(d)   #distance matrix needed for hclust()
  return(d)
}


optimal_clust=function(x, y) #x has to be vnc() object and y a distance matrix, distance matrix needed for silhouettes
{
  K=nrow(x$merge) 
  
  outliers=c()
  sil_averages=c()
  
  y=distvnc(y)  
  #calculate silhouettes via silhouette() function for all clusters at each merging step
  for (n in seq(2,K)) {
    sil=silhouette(cutree(x, k=n), y)
    silsum=summary(sil)
    silavg=silsum$avg.width #average silhouette coefficients for all clusters at current merging step
    sil_averages=c(sil_averages,silavg) #collect average silhouette coefficients of each merge in a vector
  }
  
  #opt_clust corresponds to the number of clusters at the merging step with the highest average silhouette coefficient (+1 because previous iteration started at 2)
  opt_clust=which(sil_averages==max(sil_averages))+1
  opt_clust=min(opt_clust)    #if maximum value is shared by more than one clustering, choose the one with the lower number of clusters (arbitrary).
  
  print(paste0("Optimal number of clusters according to silhouette coefficient: ", opt_clust))
  silcoef=max(sil_averages)
  print(paste0("Average silhouette coefficient of clustering: ", silcoef))
  
  
  
  cutoff=opt_clust
  final_clust=cutree(x, k=cutoff)
  
  print("Cluster memberships:")
  print(final_clust)
  
  #outlier identification - current approach: individual texts that do not cluster with other texts, i.e., singletons. This is efficient when applying the diahclust() function.
  singletons=table(final_clust) #table gives an array, contigency table
  singles=which(singletons==1)
  
  print("Potential outliers:")
  
  for (j in seq(1,length(singles))) {
    find_singles=which(final_clust==singles[j]) #clusters with only one members
    
    #get only the ones that are still single texts
    if (length(find_singles)>0 && grepl("^X[0-9]+\\.[A-Z]+$", rownames(as.matrix(find_singles)))) {
      print(rownames(as.matrix(find_singles)))
    }
    
    
  }
  cat("\n")
  structure(list(opt_clust=opt_clust, final_clust=final_clust, silcoef=silcoef))
}



#aggregates data for iterative diahclust() according to outcome of optimal_clust()
aggregate_data=function(x,y)  #x has to be an optimal_clust object, and y is a data matrix
{
  aggregated_cluster_names=c()
  
  aggregated_data=cbind()
  
  #takes the optimal clustering as assessed via optimal_clust() and aggregates the data points in each cluster
  for (i in seq(1,x$opt_clust)) { 
    clusters=names(x$final_clust[x$final_clust==i])   #final_clust is a named vector
    C=length(clusters)
    
    cluster_vectors=cbind()
    cluster_names=c()
    for (j in seq(1,C)) {
      k=which(colnames(y)==clusters[j]) 
      
      cluster_vectors=cbind(cluster_vectors,y[,k])
      cluster_names=c(cluster_names,colnames(y)[k]) 
    }
    
    aggregated_cluster_name=paste(cluster_names, collapse=" + ")
    
    aggregated_vector=apply(cluster_vectors, 1, mean)
    aggregated_data=cbind(aggregated_data, aggregated_vector)
    
    aggregated_cluster_names=c(aggregated_cluster_names, aggregated_cluster_name)
    
    
    
  }
  
  colnames(aggregated_data)=aggregated_cluster_names
  aggregated_data=data.frame(aggregated_data)
  
  return(aggregated_data)
  
}


diahclust=function(x, y, method=c("single","complete","average","median", "ward.D", "ward.D2", "mcquitty", "centroid")) #x has to be optimal_clust object, y a data matrix
{
  count=1
  
  optimal_it=x
  data_it=y
  
  #continue optimal clustering until less than 10 clusters, aggregate data along the way
  
  while(optimal_it$opt_clust>9) {
    print(paste0("Iteration: ", count))
    data_it=aggregate_data(optimal_it, data_it)  #recursion needed
    aggregated.cor=cor(data_it) 
    aggregated.dist=dist(aggregated.cor)
    
    vnc_it=vnc(aggregated.dist, method=match.arg(method))
    
    #abbreviate label names for plotting
    J=ncol(data_it)
    abbreviated_names=c()
    
    
    for (i in seq(1,J))
    {
      last=stringr::str_extract(string = colnames(data_it)[i], pattern = "X[0-9]+\\.([A-Z]|[0-9])+$")
      first=stringr::str_extract(string = colnames(data_it)[i], pattern = "^X[0-9]+\\.([A-Z]|[0-9])+")
      
      abbreviated_name=paste(first, last, sep=" - ")
      if (first==last) {
        abbreviated_name=first
      }
      
      abbreviated_names=c(abbreviated_names, abbreviated_name)
    }
    
    
    plot(vnc_it, labels=abbreviated_names, hang=-1)
    
    optimal_it=optimal_clust(vnc_it, aggregated.dist)
    
    count=count+1
  }
  return(vnc_it)
}


####examples

icelandic=read.table("/../data/icelandic.txt", header=TRUE) #or use lazy load

#uncomment the following for outlier removal of, e.g., the text labeled "X1830.HELLISMENN"
#remove HELLISMEN
#icelandic[,"X1830.HELLISMENN"]=NULL

icelandic.cor=cor(icelandic[,-1])  #[,-1] because rows are labeled
icelandic.dist=dist(icelandic.cor)


icelandic.vnc=vnc(icelandic.dist, method="average")
plot(icelandic.vnc, hang=-1) 

optimal=optimal_clust(icelandic.vnc, icelandic.dist) 

icelandic.diahclust=diahclust(optimal, icelandic[,-1], method="average")

