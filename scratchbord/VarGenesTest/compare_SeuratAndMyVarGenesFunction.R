## mainly follow the cell_cycle_vignette.html
# https://satijalab.org/seurat/v3.1/cell_cycle_vignette.html

library(Seurat)
packageVersion('Seurat')
#[1] ‘3.1.1’

exp.mat <- read.table(file = "data/nestorawa_forcellcycle_expressionMatrix.txt", header = TRUE,
    as.is = TRUE, row.names = 1)


s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

marrow <- CreateSeuratObject(counts = exp.mat)
marrow <- NormalizeData(marrow)
marrow <- FindVariableFeatures(marrow, selection.method = "vst")
marrow <- ScaleData(marrow, features = rownames(marrow))

marrow <- RunPCA(marrow, features = VariableFeatures(marrow), ndims.print = 6:10, nfeatures.print = 10)

marrow <- CellCycleScoring(marrow, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

grDevices::pdf( file="RidgePlot.pdf", width=10, height=10 )
RidgePlot(marrow, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
dev.off()


marrow <- RunPCA(marrow, features = VariableFeatures(marrow) )

grDevices::pdf( file="CellCylcle_first_2D_PCA_SeuratVarGenes.pdf", width=10, height=10 )
DimPlot(marrow)
dev.off()

marrow <- RunPCA(marrow, features = c(s.genes, g2m.genes))

grDevices::pdf( file="CellCylcle_first_2D_PCA_SeuratCellCycleGenes.pdf", width=10, height=10 )
DimPlot(marrow)
dev.off()

## right now skip to the alternative approach:

marrow$CC.Difference <- marrow$S.Score - marrow$G2M.Score

## wow - this step takes a lot of time!
## more than 20 min! with 774 cells and 24193 genes!
marrow <- ScaleData(marrow, vars.to.regress = "CC.Difference", features = rownames(marrow))

marrow <- RunPCA(marrow, features = VariableFeatures(marrow), nfeatures.print = 10)
#PC_ 1 
#Positive:  BLVRB, KLF1, ERMAP, FAM132A, CAR2, RHD, CES2G, SPHK1, AQP1, SLC38A5 
#Negative:  TMSB4X, CORO1A, PLAC8, H2AFY, LAPTM5, CD34, LCP1, TMEM176B, IGFBP4, EMB 
#PC_ 2 
#Positive:  APOE, GATA2, RAB37, ANGPT1, ADGRG1, MEIS1, MPL, F2R, PDZK1IP1, DAPP1 
#Negative:  CTSG, ELANE, LY6C2, HP, CLEC12A, ANXA3, IGSF6, TIFAB, SLPI, MPO 
#PC_ 3 
#Positive:  APOE, GATA2, NKG7, MUC13, ITGA2B, TUBA8, CPA3, RAB44, SLC18A2, CD9 
#Negative:  DNTT, FLT3, WFDC17, LSP1, MYL10, LAX1, GIMAP6, IGHM, CD24A, MN1 
#PC_ 4 
#Positive:  CSRP3, ST8SIA6, SCIN, LGALS1, APOE, ITGB7, MFSD2B, RGL1, DNTT, IGHV1-23 
#Negative:  MPL, MMRN1, PROCR, HLF, SERPINA3G, ESAM, PTGS1, D630039A03RIK, NDN, PPIC 
#PC_ 5 
#Positive:  HDC, LMO4, CSRP3, IFITM1, FCGR3, HLF, CPA3, PROCR, PGLYRP1, IKZF2 
#Negative:  GP1BB, PF4, SDPR, F2RL2, TREML1, RAB27B, SLC14A1, PBX1, PLEK, TUBA8 

grDevices::pdf( file="CellCylcle_first_2D_PCA_SeuratVarGenes_afterRegressingOutCellCycle.pdf", width=10, height=10 )
DimPlot(marrow)
dev.off()

## re obtain the unscaled data:
marrow.orig <- CreateSeuratObject(counts = exp.mat)
marrow.orig <- NormalizeData(marrow.orig)
marrow.orig <- FindVariableFeatures(marrow.orig, selection.method = "vst")
marrow.orig <- ScaleData(marrow.orig, features = rownames(marrow.orig))

marrow <- CellCycleScoring(marrow, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

BioData_VarGenes <- scan( 'BioData_varGenes_2000.txt', what=character())

marrow.orig <- RunPCA(marrow.orig, features = BioData_VarGenes, nfeatures.print = 10)

#PC_ 1
#Positive:  PCSK2, MRGPRX2, ZC3H12C, COL26A1, MAF, GATA3, TSHZ2, C530043K16RIK, TSKU, POU3F1
#Negative:  OPTN, DIAPH3, TRIB2, ATG2B, TULP4, USP32, TGFBR3, ADGRL2, HIF1AN, FAM179B
#PC_ 2
#Positive:  OPTN, TRIB2, RNF128, SYNGAP1, TSKU, SEMA6A, ABCA5, LRRC15, DCLK1, TMEM178B
#Negative:  CD34, CORO1A, H2AFY, PKM, TMSB4X, FXYD5, EMB, LCP1, ALOX5AP, SPI1
#PC_ 3
#Positive:  TXNIP, CTLA2A, RAB37, APOE, GATA2, ANGPT1, H2-Q6, F2R, S100A10, MEIS1
#Negative:  CTSG, ELANE, ANXA3, MPO, LY6C2, ATP8B4, PRTN3, ALAS1, TYROBP, HK3
#PC_ 4
#Positive:  FLT3, DNTT, WFDC17, MYL10, LSP1, 9030619P08RIK, LY6A, GIMAP6, GPR171, IGHM
#Negative:  APOE, NKG7, GATA2, CPA3, FAM46A, RAB44, HDC, FCGR3, CSRP3, MS4A3
#PC_ 5
#Positive:  MPL, SERPINA3G, CD63, CD63-PS, PROCR, F11R, TRPC6, PTGS1, GSTM2-PS1, ESAM
#Negative:  IRF8, S100A10, DNTT, CTSS, IGHM, TRIB2, SATB1, LGALS1, FLT3, DOCK10

grDevices::pdf( file="CellCylcle_first_2D_PCA_BioDataVarGenes.pdf", width=10, height=10 )
DimPlot(marrow.orig)
dev.off()

marrow <- RunPCA(marrow, features = BioData_VarGenes, ndims.print = 6:10, nfeatures.print = 10)

grDevices::pdf( file="CellCylcle_first_2D_PCA_BioDataVarGenes_afterRegressingOutCellCycle.pdf", width=10, height=10 )
DimPlot(marrow)
dev.off()

# What I have seen in Pavans data is, that using the SeuratVarGenes and intersecting that with the low expressed but in many cells gives an even better starting set for the analysis:

Intersect = names(which(table( c( BioData_VarGenes,  VariableFeatures(marrow.orig, selection.method = "vst") )) == 2))
length(Intersect)
#186 ## starting from 4000 genes!!
marrow <- RunPCA(marrow, features = Intersect, nfeatures.print = 10)
#PC_ 1 
#Positive:  APOE, F2R, FBXL13, ST8SIA6, CSRP3, SCIN, GATA2, PF4, PORCN, ITGB2 
#Negative:  CORO1A, CD34, H2AFY, TMSB4X, EMB, LCP1, PLAC8, CD53, PRTN3, MPO 
#PC_ 2 
#Positive:  ELANE, CTSG, ANXA3, LY6C2, ALAS1, ATP8B4, MPO, MS4A3, MTUS1, MS4A6C 
#Negative:  ANGPT1, MEIS1, MPL, RAB37, ZYX, DAPP1, GATA2, RGS1, MYCN, APOE 
#PC_ 3 
#Positive:  APOE, GATA2, NKG7, CPA3, SLC18A2, HDC, RAB44, CSRP3, MUC13, LMO4 
#Negative:  FLT3, DNTT, WFDC17, LSP1, IGHM, CTSS, MYL10, SATB1, IL12A, GPR171 
#PC_ 4 
#Positive:  CSRP3, SCIN, ST8SIA6, LGALS1, FOS, IER2, JUN, HSPA1A, DNTT, EGR1 
#Negative:  SERPINA3G, MPL, PTGS1, ESAM, PF4, HLF, GUCY1A3, CD9, GSTM1, MS4A3 
#PC_ 5 
#Positive:  CPA3, LMO4, HDC, HLF, MYL10, IFITM1, PROCR, FUT8, CDH1, LY6A 
#Negative:  PF4, PBX1, F2R, MEF2C, PLEK, IRF8, CSF1R, NDRG1, CCR2, IL6ST 

grDevices::pdf( file="CellCylcle_first_2D_PCA_INTERSECTGenes_afterRegressingOutCellCycle.pdf", width=10, height=10 )
DimPlot(marrow)
dev.off()


marrow.orig <- RunPCA(marrow.orig, features = Intersect, nfeatures.print = 10)
#PC_ 1 
#Positive:  APOE, F2R, FBXL13, ST8SIA6, CSRP3, SCIN, GATA2, PF4, PORCN, ITGB2 
#Negative:  CORO1A, CD34, H2AFY, EMB, TMSB4X, LCP1, PLAC8, CD53, PRTN3, MPO 
#PC_ 2 
#Positive:  ELANE, CTSG, LY6C2, ANXA3, ALAS1, ATP8B4, MPO, MS4A3, MTUS1, MS4A6C 
#Negative:  ANGPT1, MEIS1, RAB37, MPL, ZYX, DAPP1, GATA2, RGS1, APOE, MYCN 
#PC_ 3 
#Positive:  APOE, GATA2, NKG7, CPA3, SLC18A2, HDC, RAB44, CSRP3, MUC13, LMO4 
#Negative:  FLT3, DNTT, WFDC17, LSP1, IGHM, CTSS, MYL10, SATB1, IL12A, GPR171 
#PC_ 4 
#Positive:  CSRP3, SCIN, ST8SIA6, LGALS1, FOS, IER2, JUN, HSPA1A, DNTT, EGR1 
#Negative:  SERPINA3G, MPL, PTGS1, ESAM, PF4, HLF, GUCY1A3, CD9, GSTM1, MS4A3 
#PC_ 5 
#Positive:  CPA3, LMO4, HDC, HLF, MYL10, IFITM1, FUT8, CDH1, PROCR, IKZF2 
#Negative:  PF4, F2R, PBX1, MEF2C, PLEK, IRF8, CSF1R, NDRG1, CCR2, IL6ST 


grDevices::pdf( file="CellCylcle_first_2D_PCA_INTERSECTGenes.pdf", width=10, height=10 )
DimPlot(marrow.orig)
dev.off()


