\name{DiaHClust-package}
\alias{DiaHClust-package}
\alias{DiaHClust}
\docType{package}
\title{
DiaHClust: an iterative hierarchical clustering approach for identifying stages in language change
}
\description{
This package contain functions which implement the DiaHClust methodology (Schaetzle and Booth 2019) for the identification of stages in historical linguistic change. DiaHClust is based on the iterative application of Variability-based Neighbor Clustering (VNC, Gries and Hilpert 2008, 2012) in combination with a cluster validation process using silhouette values in order to provide a multi-layered perspective on language change, from text-level to broader time stages, while also respecting outliers and genre effects.
}
\details{
\tabular{ll}{
Package: \tab DiaHClust\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2019-04-30\cr
License: \tab GPL-3\cr
}
~~ An overview of how to use the package, including the most important functions ~~
}
\author{
Christin Schaetzle

Maintainer: Christin Schaetzle <christin.schaetzle@uni-konstanz.de>
}
\references{
Stefan Th. Gries and Martin Hilpert. 2008. The identification of stages in diachronic data: variability-based neighbour clustering. Corpora, 3(1):59–81. 
Stefan Th. Gries and Martin Hilpert. 2012. Variability-based neighbor clustering: A bottom-up approach to periodization in historical linguistics. In Nevalainen Terttu and Elizabeth Closs Traugott, editors, The Oxford Handbook of the History of English, pages 134–144. Oxford University Press, Oxford. 
Christin Schaetzle and Hannah Booth. 2019. DiaHClust: an iterative hierarchical clustering apprach for identifying stages in language change. to appear.
}
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
