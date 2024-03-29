#' @importFrom AnnotationDbi mapIds
#' @importFrom DDRTree DDRTree
#' @importFrom DESeq2 DESeq DESeqDataSetFromMatrix results
#' @importFrom FastWilcoxTest CorMatrix NormalizeCells StatTest ZScore
#' @importFrom Matrix Matrix colSums drop0 rowSums sparseMatrix t
#' @importFrom R6 R6Class
#' @importFrom RDRToolbox Isomap LLE
#' @importFrom RSQLite dbClearResult dbConnect dbDisconnect dbFetch dbGetQuery dbHasCompleted dbSendQuery dbWriteTable fetch
#' @importFrom RSvgDevice devSVG
#' @importFrom Rtsne Rtsne
##' @importFrom Seurat CellCycleScoring CreateSeuratObject FindAllMarkers SetIdent
#' @importFrom dplyr progress_estimated
#' @importFrom gplots bluered
#' @importFrom grDevices pdf png rainbow x11
#' @importFrom graphics abline axis barplot boxplot hist image layout legend lines mtext par pie plot rect text title
#' @importFrom knitr kable
#' @importFrom limma toptable
#' @importFrom mclust hclass
#' @importFrom methods new setGeneric
#' @importFrom pcaMethods bpca
#' @importFrom randomForest randomForest
#' @importFrom rgl axis3d bg3d bgplot3d grid3d par3d
#' @importFrom stats anova cmdscale cor cutree density hclust kmeans predict quantile reorder runif sd var
#' @importFrom stringr str_replace str_replace_all str_split
#' @importFrom topGO GenTable runTest
#' @importFrom umap umap
#' @importFrom utils sessionInfo
NULL;


#' @name CellCycleGenes
#' @title A simple table containing the Human and mouse orthologe ensembl ids for the cell cycle genes used in PMID28263960 using BioMart and GRCh38.p10
#' @description The BioMart web interface was used to create the input data, which in turn was subselected using grep.
#' @docType data
#' @usage CellCycleGenes
#' @format data.frame
#' @keywords data
'CellCycleGenes'

#' @name CellCycle
#' @title A simple table containing the Human and mouse orthologe CellCycle genes from PMID17994010
#' "Genome-scale RNAi profiling of cell division in human tissue culture cells."
#' @description The data can be used by stating
#' useInbuiltGOIlists (cellexalObj, 'CellCycle' )
#' And it is used to visualize the cell cycle genes in the VR environment.
#' Only the genes also identifyable in mouse were used here.
#' @docType data
#' @usage CellCycle
#' @format data.frame
#' @keywords data
'CellCycle'

#' A list of mouse transcription factors.
#'
#' A list of mouse transcription factors. Used when generating TF networks
#'
#' @format A vector 1357 in length:
#' \describe{
#'   \item{tf}{transcription factor}
#'   ...
#' }
#' @source \url{http://bioinfo.life.hust.edu.cn/AnimalTFDB/}
#' @keywords data
"mouse.tfs"

#' A list of human transcription factors.
#'
#' A list of human transcription factors. Used when generating TF networks
#'
#' @format A vector 1468 in length:almdiR_Aurroa
#' 
#' \describe{
#'   \item{tf}{transcription factor}
#'   ...
#' }
#' @source \url{http://bioinfo.life.hust.edu.cn/AnimalTFDB/}
#' @keywords data
"human.tfs"

#' @name Epigenetic
#' @title A simple table containing the data from http://epifactors.autosome.ru/ as from 21st September 2017
#' @description This table can be used to create the epigenetics MDS objects.
#' @docType data
#' @usage Epigenetic
#' @format data.frame
#' @keywords data
'Epigenetic'

#' @name human.CellSurface
#' @title A simple list containing the human ENTREZ.gene.symbols from 'A Mass Spectrometric-Derived Cell Surface Protein Atlas' plos 2015
#' @description This list is used to crete a CellSurface object
#' @docType data
#' @usage human.CellSurface
#' @format vector
#' @keywords data
"human.CellSurface"

#' @name mouse.CellSurface
#' @title A simple list containing the mouse ENTREZ.gene.symbols from 'A Mass Spectrometric-Derived Cell Surface Protein Atlas' plos 2015
#' @description This list is used to crete a CellSurface object
#' @docType data
#' @usage mouse.CellSurface
#' @format vector
#' @keywords data
"mouse.CellSurface"

#' @name TestData
#' @title merged the two smallest Mouse Tabular muris datasets facs_Aorta_TabulaMuris and droplet_Heart_TabulaMuris
#' @description two TabulaMuris mouse datasets merged
#' @docType data
#' @usage TestData
#' @format BioData
#' @keywords data
"TestData"
