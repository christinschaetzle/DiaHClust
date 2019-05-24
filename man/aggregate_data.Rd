\name{aggregate_data}
\alias{aggregate_data}
\title{Data Aggregation for DiaHClust}
\description{
This function aggregates vectorial data in the form of a matrix for the iterative application of Variability-based Neighbor Clustering (VNC, Gries and Hilpert 2008, 2012) and is part of the \code{diahclust} function.
}
\usage{
aggregate_data(x, y)
}
\arguments{
  \item{x}{
an object of class \code{optimal_clust}.
}
  \item{y}{
a data matrix where each column represents a vector containing the changing syntactic features extracted from a text.
}
}
\details{
 \code{aggregate_data} is called inside the \code{diahclust} function and is in charge of the aggregation of the data for the iterative application of the \code{vnc} function. Data points which belong to a single cluster are aggregated by averaging the corresponding syntactic vectors in the underlying dataset with \code{aggregate_data}. To keep track of the texts and time stages which form clusters across the iterations, the names of the new vectors consist of the sequence of the names of the aggregated vectors. 
}
\value{
A data matrix containing the aggregated vectors of change.
}
\references{
Stefan Th. Gries and Martin Hilpert. 2008. The identification of stages in diachronic data: variability-based neighbour clustering. Corpora, 3(1):59–81. 
Stefan Th. Gries and Martin Hilpert. 2012. Variability-based neighbor clustering: A bottom-up approach to periodization in historical linguistics. In Nevalainen Terttu and Elizabeth Closs Traugott, editors, The Oxford Handbook of the History of English, pages 134–144. Oxford University Press, Oxford. 
Christin Schätzle and Hannah Booth. 2019. DiaHClust: an iterative hierarchical clustering apprach for identifying stages in language change. to appear.
}
\author{Christin Schätzle}
\seealso{
	\code{\link{diahclust}}
}
\examples{
icelandic=data(icelandic)
icelandic.cor=cor(icelandic[,-1])  #[,-1] because rows are labeled
icelandic.dist=dist(icelandic.cor)

icelandic.vnc=vnc(icelandic.dist, method="average")
plot(icelandic.vnc, hang=-1) #plotting the resulting dendrogram

#cluster validation for identification of time stages
optimal=optimal_clust(icelandic.vnc, icelandic.dist) 

#aggregation of data for DiaHClust 
icelandic_aggregated=aggregate_data(optimal, icelandic[,-1])

#function is automatically called inside diahclust
icelandic.diahclust=diahclust(optimal, icelandic[,-1], method="average")
}

