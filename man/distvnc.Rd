\name{distvnc}
\alias{distvnc}
\title{Manipulation of distance matrix for Variability-based Neighbor Clustering}
\description{
This function manipulates a distance matrix as created via \code{dist} for Variability-based Neighbor Clustering (VNC), see, e.g., Gries and Hilpert (2008, 2012), so that only temporally adjacent data points (i.e. texts) are allowed to cluster with one another. 
}
\usage{
distvnc(d)
}

\arguments{
  \item{d}{a distance matrix as produced by \code{dist}.}
}
\details{
The \code{vnc} function takes the same input as \code{hclust}, but manipulates the dissimilarity matrix so that only temporally adjacent data points (i.e. texts) are allowed to cluster with one another via the application of the \code{distvnc} function. \code{distvnc} sets all distance values which describe distances between non-temporally adjacent data points to the value which equals the maximum value of the dissimilarity matrix. As similarity is measured in terms of the minimum distance, it is highly unlikely that two data points which have these maximized distances to one another will be merged in the clustering process. This in turn allows for the application of \code{hclust} -- which is called inside \code{vnc} after the manipulation of the dissimilarity matrix -- for clustering according to the ideas of VNC, instead of having to implement a separate clustering algorithm. 
}
\value{
An object of class \code{dist}.
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
	\code{\link{dist}}
}
\examples{
icelandic=data(icelandic)
icelandic.cor=cor(icelandic[,-1])  #[,-1] because rows are labeled
icelandic.dist=dist(icelandic.cor)

##only manipulate distance matrix
icelandic.distvnc=distvnc(icelandic.dist)

###distvnc() implemented in vnc()
icelandic.vnc=vnc(icelandic.dist, method="average")
plot(icelandic.vnc, hang=-1) #plotting the resulting dendrogram
}

