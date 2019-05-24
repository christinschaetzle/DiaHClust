\name{vnc}
\alias{vnc}
\title{Variability-based Neighbor Clustering}
\description{
Implementation of the Variability-based Neighbor Clustering (VNC) algorithm developed by Gries and Hilpert (2008, 2012) for the identification of stages in language change which is based on \code{hclust} as part of the DiaHClust methodology (Schätzle and Booth 2019). 
}
\usage{
vnc(d, method = c("single", "complete", "average", "median", "ward.D", "ward.D2", "mcquitty", "centroid"))
}
\arguments{
  \item{d}{a distance matrix as produced by \code{dist}.}
  \item{method}{
the agglomeration method to be used. All methods from \code{hclust} are available for \code{vnc}, i.e., "single", "complete", "average", "median", "ward.D", "ward.D2", "mcquitty", and "centroid". 
}
}
\details{
The VNC approach is implemented by manipulating individual steps in the workflow behind R's standard agglomerative hierarchical clustering function \code{hclust}. In the vector-based approach to VNC by Gries and Hilpert (2008, 2012), a correlation statistic is calculated before clustering the data. This is also often done when clustering vectorial data with \code{hclust} (see, e.g., Baayen 2008). Thus, a correlation matrix should be calculated via \code{cor} which is then turned into a distance matrix via \code{dist} as input for \code{vnc}. For the analysis of syntactic change, the correlation matrix should be calculated based on a data matrix where each column represents a vector containing the changing syntactic features extracted from a text, see, e.g., the \code{icelandic} data set in the examples below. In the data matrix, the vectors have to be ordered from left to right according to the time stamp of the text. The time stamp should be encoded in the vector name, i.e., the name of the corresponding column in the data matrix. For the application of the full DiaHClust methodology, the vector name has to begin with a four digit year date followed by a dot and the text name, e.g., ``1250.STURLUNGA'', allowing one to easily identify individual texts in the clustering (This roughly corresponds to token IDs in Penn-style treebanks).  

\code{hclust} usually begins by clustering together the two most similar vectors, i.e., the data points with the smallest distance to one another, merging these two data points. This process continues until all data points have been clustered. The method chosen for clustering with \code{hclust} represents the method of agglomeration. For example, when \code{method="average"} is chosen for agglomeration, cluster similarity between two clusters is assessed based on the average of the data points in the clusters. Moreover, the two data points with the smallest distance are merged into a new data point by averaging the corresponding values after each iteration. This corresponds to the idea behind the amalgamation method in VNC (see Gries and Hilpert 2008, 2012). In general, all agglomeration methods available with \code{hclust} are available with \code{vnc}. It is recommended to use averages when applying \code{vnc} -- following Gries and Hilpert (2008, 2012) -- since, in quantitative corpus linguistics, (co-)occurence frequencies are usually assessed by averaging frequencies over texts/time periods.

The \code{vnc} function takes the same input as \code{hclust}, but manipulates the distance matrix so that only temporally adjacent data points (i.e. texts) are allowed to cluster with one another via the application of the \code{distvnc} function. \code{distvnc} sets all distance values which describe distances between non-temporally adjacent data points to the value which equals the maximum value of the distance matrix. As similarity is measured in terms of the minimum distance, it is highly unlikely that two data points which have these maximized distances to one another will be merged in the clustering process. This in turn allows for the application of \code{hclust} -- which is called inside \code{vnc} after the manipulation of the distance matrix -- for clustering according to the ideas of VNC, instead of having to implement a separate clustering algorithm. Thus, \code{vnc} returns an object of class \code{hclust}. Moreover, \code{vnc} adjusts the permutations of the data points which arise during the merging process in order to guarantee the diachronic ordering of data points for plotting via \code{plot}. For more details, see also Schätzle and Booth (2019).

}
\value{
	An object of class \code{hclust}, describing the tree produced by the clustering process. The object returns a list containing the following elements (information mostly taken from \code{hclust} help page):
	\item{merge}{an n-1 by 2 matrix. Row i of merge describes the merging of clusters at step i of the clustering. If an element j in the row is negative, then observation -j was merged at this stage. If j is positive then the merge was with the cluster formed at the (earlier) stage j of the algorithm. Thus negative entries in merge indicate agglomerations of singletons, and positive entries indicate agglomerations of non-singletons.}	
	\item{height}{a set of n-1 real values (non-decreasing for ultrametric trees). The clustering height: that is, the value of the criterion associated with the clustering method for the particular agglomeration.}	
	\item{order}{a vector giving the permutation of the original observations suitable for plotting, in the sense that a cluster plot using this ordering and matrix merge will not have crossings of the branches. The permutations of the observations are adjusted with \code{vnc} in order to maintain the diachronic ordering.}	
	\item{labels}{labels for each of the objects being clustered.}
	\item{call}{the call which produced the result.}
	\item{method}{the cluster method that has been used.}
	\item{dist.method}{the distance that has been used to create d (only returned if the distance object has a "method" attribute).}

}
\references{
R. Harald Baayen. 2008. Analyzing Linguistic Data. A Practical Introduction to Statistics Using R. Cambridge University Press, Cambridge. 
Stefan Th. Gries and Martin Hilpert. 2008. The identification of stages in diachronic data: variability-based neighbour clustering. Corpora, 3(1):59–81. 
Stefan Th. Gries and Martin Hilpert. 2012. Variability-based neighbor clustering: A bottom-up approach to periodization in historical linguistics. In Nevalainen Terttu and Elizabeth Closs Traugott, editors, The Oxford Handbook of the History of English, pages 134–144. Oxford University Press, Oxford. 
Christin Schätzle and Hannah Booth. 2019. DiaHClust: an iterative hierarchical clustering apprach for identifying stages in language change. to appear.
}
\author{Christin Schätzle}
\seealso{
	\code{\link{hclust}}
	\code{\link{diahclust}}
	\code{\link{distvnc}}
}
\examples{
icelandic=data(icelandic)
icelandic.cor=cor(icelandic[,-1])  #[,-1] because rows are labeled
icelandic.dist=dist(icelandic.cor)

icelandic.vnc=vnc(icelandic.dist, method="average")
plot(icelandic.vnc, hang=-1) #plotting the resulting dendrogram
}

