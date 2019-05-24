# diaHClust

DiaHClust is an R package which implements the DiaHClust methodology (Schätzle and Booth 2019), 
    a new approach which can be used to identify stages in language change based on quantitative corpus-derived data. 
    DiaHClust is based on the hierarchical clustering approach for historical data from Hilpert and Gries (2008, 2012) 
    called 'Variability-based Neighbor Clustering' (VNC) and develops this further by adding an extra iterative approach 
    to the hierarchical clustering which results in a multi-layered perspective on change, from text-level to broader periods.

We plan to submit the DiaHClust package to CRAN in order to make it available in R directly. Thus, the package follows the typical R-package structure. The R code itself can be found in the R-subdirectory of DiaHClust including sample calls of the R commands. For usage in R, use the source()-function to import the code stored in diahclust.R for now. In man, you'll find a bunch of help pages describing each individual function in more detail. Example data is stored in the data-subdirectory. 
 
More information can be found in Schätzle and Booth (2019) which will be published soon. We'll update the references to this paper accordingly. 


 
References 

Stefan Th. Gries and Martin Hilpert. 2008. The identification of stages in diachronic data: variability-based neighbour clustering. Corpora, 3(1):59–81. 

Stefan Th. Gries and Martin Hilpert. 2012. Variability-based neighbor clustering: A bottom-up approach to periodization in historical linguistics. In Nevalainen Terttu and Elizabeth Closs Traugott, editors, The Oxford Handbook of the History of English, pages 134–144. Oxford University Press, Oxford. 

Christin Schaetzle and Hannah Booth. 2019. DiaHClust: an iterative hierarchical clustering apprach for identifying stages in language change. to appear.
