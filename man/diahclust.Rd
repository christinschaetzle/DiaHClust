\name{diahclust}
\alias{diahclust}
\title{
DiaHClust: an iterative hierarchical clustering approach for identifying stages in language change
}
\description{
This function implements the DiaHClust methodology (Schätzle and Booth 2019) for the identification of stages in historical linguistic change. DiaHClust is based on the iterative application of Variability-based Neighbor Clustering (VNC, Gries and Hilpert 2008, 2012) in combination with a cluster validation process using silhouette values in order to provide a multi-layered perspective on language change, from text-level to broader time stages, while also respecting outliers and genre effects.
}
\usage{
diahclust(x, y, method = c("single", "complete", "average", "median", "ward.D", "ward.D2", "mcquitty", "centroid"))
}
\arguments{
  \item{x}{an object of class optimal_clust.}
  \item{y}{a data matrix containing vectors of syntactic change (or other numerical data) as rows. The names of the vectors have to consist of a year date followed by a dot and a string of characters (e.g., the text name).}
    \item{method}{
the agglomeration method to be used. All methods from \code{hclust} are available for \code{diahclust}, i.e., "single", "complete", "average", "median", "ward.D", "ward.D2", "mcquitty", and "centroid". 
}
}
\details{
When the results of the \code{optimal_clust} function indicate that an initial clustering via \code{vnc} yields 10 or more clusters, the clustering process can be continued via the \code{diahclust} function. \code{diahclust} takes as input the original data matrix which is aggregated after each iteration and the output created by the \code{optimal_clust} function after the first clustering. Thus, before applying \code{diahclust}, clustering via \code{vnc} has to be performed and validated by \code{optimal_clust}. 

In \code{diahclust}, data points which belong to a single cluster are aggregated by averaging the corresponding syntactic vectors in the underlying dataset. This is done via \code{aggregate_data}. To keep track of the texts and time stages which form clusters across the iterations, the names of the new vectors consist of the sequence of the names of the aggregated vectors. The previously applied process of VNC with respect to the new dataset is repeated by calling \code{cor}, \code{dist}, and \code{vnc} in \code{diahclust}. Moreover, an agglomeration method has to be specified. When the agglomeration method chosen for VNC and DiaHClust is not ``average'', a different aggregation method, e.g., the minimum with single linkage clustering, should be applied.

\code{diahclust} automatically plots the clustering as a dendrogram. The labels on the dendrogram are abbreviated for better visibility, representing the range of previously aggregated vectors, with the oldest and the youngest text in the range connected via a hyphen.  The resulting clustering is again evaluated using the \code{optimal_clust} function, which returns the cluster memberships listing the full range of texts in the clusters. The application of this process is repeated until the final evaluation arrives at an optimal number of clusters less than 10. In this iterative process, the clusters, i.e., time stages, can be inspected at each step of the iteration, allowing one to track the composition of the clusters with respect to the individual texts from the first iteration onwards. 
}
\value{
	After the last iteration, \code{diahclust} returns object of class \code{hclust}, describing the tree produced by the final clustering process. The object returns a list containing the following elements (information mostly taken from \code{hclust} help page, see also \code{vnc}):
	\item{merge}{an n-1 by 2 matrix. Row i of merge describes the merging of clusters at step i of the clustering. If an element j in the row is negative, then observation -j was merged at this stage. If j is positive then the merge was with the cluster formed at the (earlier) stage j of the algorithm. Thus negative entries in merge indicate agglomerations of singletons, and positive entries indicate agglomerations of non-singletons.}	
	\item{height}{a set of n-1 real values (non-decreasing for ultrametric trees). The clustering height: that is, the value of the criterion associated with the clustering method for the particular agglomeration.}	
	\item{order}{a vector giving the permutation of the original observations suitable for plotting, in the sense that a cluster plot using this ordering and matrix merge will not have crossings of the branches. The permutations of the observations are adjusted with \code{vnc} in order to maintain the diachronic ordering.}	
	\item{labels}{labels for each of the objects being clustered.}
	\item{call}{the call which produced the result.}
	\item{method}{the cluster method that has been used.}
	\item{dist.method}{the distance that has been used to create d (only returned if the distance object has a "method" attribute).}
}
\references{
Stefan Th. Gries and Martin Hilpert. 2008. The identification of stages in diachronic data: variability-based neighbour clustering. Corpora, 3(1):59–81. 
Stefan Th. Gries and Martin Hilpert. 2012. Variability-based neighbor clustering: A bottom-up approach to periodization in historical linguistics. In Nevalainen Terttu and Elizabeth Closs Traugott, editors, The Oxford Handbook of the History of English, pages 134–144. Oxford University Press, Oxford. 
Christin Schätzle and Hannah Booth. 2019. DiaHClust: an iterative hierarchical clustering apprach for identifying stages in language change. to appear.
}
\author{Christin Schätzle}
\seealso{
	\code{\link{hclust}}
	\code{\link{diahclust}}
	\code{\link{distvnc}}
	\code{\link{optimal_clust}}
}
\examples{
icelandic=data(icelandic)
icelandic.cor=cor(icelandic[,-1])  #[,-1] because rows are labeled
icelandic.dist=dist(icelandic.cor)

icelandic.vnc=vnc(icelandic.dist, method="average")
plot(icelandic.vnc, hang=-1) #plotting the resulting dendrogram

#cluster validation for identification of time stages
optimal=optimal_clust(icelandic.vnc, icelandic.dist) 

#Iterative clustering with DiaHClust methodology
icelandic.diahclust=diahclust(optimal, icelandic[,-1], method="average")
}

