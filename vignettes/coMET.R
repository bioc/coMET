### R code from vignette source 'vignettes/coMET/inst/doc/coMET.Rnw'

###################################################
### code chunk number 1: style
###################################################
## ----style, eval=TRUE, echo=FALSE,results="asis"----------------------------------------
#library("BiocStyle")
BiocStyle::latex()

###################################################
### code chunk number 2: style knitr
###################################################
## ----setup, include=FALSE, cache=FALSE, echo=F------------------------------------------
library(knitr)
opts_chunk$set(fig.path='figure/minimal-', fig.align='center', fig.show='hold')
options(replace.assign=TRUE,width=90)

## ----citatation-------------------------------------------------------------------------
citation(package='coMET')

###################################################
### code chunk number 3: install coMET
###################################################
## ----installPackages,eval=FALSE---------------------------------------------------------
#  source("http://bioconductor.org/biocLite.R")
#  biocLite("coMET")

## ----installPackagesdepend,eval=FALSE---------------------------------------------------
#  install.packages("psych")

## ----dd,echo=FALSE,eval=FALSE-----------------------------------------------------------
#  require("hash")
#  require("grid")
#  require("grDevices")
#  require("biomaRt")
#  require("Gviz")
#  require("ggbio")
#  require("rtracklayer")
#  require("GenomicRanges")
#  require("colortools")
#  require("gridExtra")
#  require("ggplot2")
#  require("trackViewer")
#  require("psych")
#  
#  rdir <- system.file("R", package="coMET",mustWork=TRUE)
#  source(file.path(rdir, "AnalyseFile.R"))
#  source(file.path(rdir, "BiofeatureGraphics.R"))
#  source(file.path(rdir, "comet.R"))
#  source(file.path(rdir, "cometWeb.R"))
#  source(file.path(rdir, "DrawPlot.R"))
#  source(file.path(rdir, "GeneralMethodComet.R"))
###################################################
### code chunk number 4: load coMET
###################################################
## ----loadPackages, eval=TRUE------------------------------------------------------------
library("coMET")

###################################################
### code chunk number 5: help
###################################################
## ----findHel, eval=FALSE----------------------------------------------------------------
#  ?comet
#  ?comet.web
#  ?comet.list

## ----installPackages_develbioc,eval=FALSE-----------------------------------------------
#  source("http://bioconductor.org/biocLite.R")
#  biocLite("coMET")

## ----installPackages_develversion,eval=FALSE--------------------------------------------
#  install.packages("YourPath/coMET_YourVersion.tar.gz",repos=NULL,type="source")
#  ##This is an example
#  install.packages("YourPath/coMET_0.99.9.tar.gz",repos=NULL,type="source")

###################################################
### code chunk number 6: read info file (site without association)
###################################################
## ----Filedata---------------------------------------------------------------------------
extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
infofile <- file.path(extdata, "cyp1b1_infofile.txt")

data_info <-read.csv(infofile, header = TRUE,
                     sep = "\t", quote = "")

head(data_info)

###################################################
### code chunk number 7: read supplementary file (region with association)
###################################################

## ----expMatrix--------------------------------------------------------------------------
extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
infoexp <- file.path(extdata, "cyp1b1_infofile_exprGene_region.txt")

data_infoexp <-read.csv(infoexp, header = TRUE, sep = "\t", quote = "")

head(data_infoexp)

###################################################
### code chunk number 8: read correlation matrix
###################################################

## ----corrMatrix-------------------------------------------------------------------------
extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
corfile <- file.path(extdata, "cyp1b1_res37_rawMatrix.txt")

data_cor <-read.csv(corfile, header = TRUE,
                    sep = "\t", quote = "")
data_cor[1:6,1:6]

###################################################
### code chunk number 9: read configuration file
###################################################
## ----configMatrix-----------------------------------------------------------------------
extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
configfile <- file.path(extdata, "config_cyp1b1_zoom_4webserver_Grch38.txt")

data_config <-read.csv(configfile, quote = "", sep="\t", header=FALSE)
data_config

###################################################
### code chunk number 10: Creation plot with comet.web
###################################################

## ----cometwebPlotText, eval=FALSE, fig.keep='last'--------------------------------------
#  extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
#  myinfofile <- file.path(extdata, "cyp1b1_infofile_Grch38.txt")
#  myexpressfile <- file.path(extdata, "cyp1b1_infofile_exprGene_region_Grch38.txt")
#  mycorrelation <- file.path(extdata, "cyp1b1_res37_rawMatrix.txt")
#  configfile <- file.path(extdata, "config_cyp1b1_zoom_4webserver_Grch38.txt")
#  comet.web(config.file=configfile, mydata.file=myinfofile,
#            cormatrix.file=mycorrelation ,mydata.large.file=myexpressfile,
#            print.image=FALSE,verbose=FALSE)

## ----cometwebPlot, echo=FALSE, fig.keep='last'------------------------------------------
extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
myinfofile <- file.path(extdata, "cyp1b1_infofile_Grch38.txt")
myexpressfile <- file.path(extdata, "cyp1b1_infofile_exprGene_region_Grch38.txt")
mycorrelation <- file.path(extdata, "cyp1b1_res37_rawMatrix.txt")
configfile <- file.path(extdata, "config_cyp1b1_zoom_4webserver_Grch38.txt")
comet.web(config.file=configfile, mydata.file=myinfofile,
          cormatrix.file=mycorrelation ,mydata.large.file=myexpressfile,
          print.image=FALSE,verbose=FALSE)

###################################################
### code chunk number 11: Creation plot with comet
###################################################
## ----cometPlotText, eval=FALSE, fig.keep='last'-----------------------------------------
#  extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
#  configfile <- file.path(extdata, "config_cyp1b1_zoom_4comet.txt")
#  myinfofile <- file.path(extdata, "cyp1b1_infofile.txt")
#  myexpressfile <- file.path(extdata, "cyp1b1_infofile_exprGene_region.txt")
#  mycorrelation <- file.path(extdata, "cyp1b1_res37_rawMatrix.txt")
#  
#  chrom <- "chr2"
#  start <- 38290160
#  end <- 38303219
#  gen <- "hg19"
#  strand <- "*"
#  
#  BROWSER.SESSION="UCSC"
#  mySession <- browserSession(BROWSER.SESSION)
#  genome(mySession) <- gen
#  
#  genetrack <-genes_ENSEMBL(gen,chrom,start,end,showId=TRUE)
#  snptrack <- snpBiomart_ENSEMBL(gen,chrom, start, end, dataset="hsapiens_snp_som",showId=FALSE)
#  cpgIstrack <- cpgIslands_UCSC(gen,chrom,start,end)
#  
#  prombedFilePath <- file.path(extdata, "/RoadMap/regions_prom_E063.bed")
#  promRMtrackE063<- DNaseI_RoadMap(gen,chrom,start, end, prombedFilePath,
#                               featureDisplay='promotor', stacking_type="squish")
#  
#  bedFilePath <- file.path(extdata, "RoadMap/E063_15_coreMarks_mnemonics.bed")
#  chromHMM_RoadMapAllE063 <- chromHMM_RoadMap(gen,chrom,start, end, bedFilePath, featureDisplay = "all", colorcase='roadmap15' )
#  
#  listgviz <- list(genetrack,snptrack,cpgIstrack,promRMtrackE063,chromHMM_RoadMapAllE063)
#  
#  comet(config.file=configfile, mydata.file=myinfofile, mydata.type="file",
#        cormatrix.file=mycorrelation,  cormatrix.type="listfile",
#        mydata.large.file=myexpressfile, mydata.large.type="listfile",
#        tracks.gviz=listgviz, verbose=FALSE, print.image=FALSE)

## ----cometPlotfile, echo=FALSE, fig.keep='last'-----------------------------------------
extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
configfile <- file.path(extdata, "config_cyp1b1_zoom_4comet.txt")
myinfofile <- file.path(extdata, "cyp1b1_infofile.txt")
myexpressfile <- file.path(extdata, "cyp1b1_infofile_exprGene_region.txt")
mycorrelation <- file.path(extdata, "cyp1b1_res37_rawMatrix.txt")

chrom <- "chr2"
start <- 38290160
end <- 38303219
gen <- "hg19"
strand <- "*"

data(geneENSEMBLtrack)
data(snpBiomarttrack)
data(cpgIslandtrack)
data(promRMtrackE063)
data(chromHMM_RoadMapAllE063)

listgviz <- list(genetrack,snptrack,cpgIstrack,promRMtrackE063,chromHMM_RoadMapAllE063)

comet(config.file=configfile, mydata.file=myinfofile, mydata.type="file",
      cormatrix.file=mycorrelation,  cormatrix.type="listfile",
      mydata.large.file=myexpressfile, mydata.large.type="listfile",
      tracks.gviz=listgviz, verbose=FALSE, print.image=FALSE)

###################################################
### code chunk number 11: Creation plot with comet without pvalue plot
###################################################
## ----cometPlotmatrix, eval=FALSE, fig.keep='last'---------------------------------------
#  extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
#  configfile <- file.path(extdata, "config_cyp1b1_zoom_4comet.txt")
#  myinfofile <- file.path(extdata, "cyp1b1_infofile.txt")
#  myexpressfile <- file.path(extdata, "cyp1b1_infofile_exprGene_region.txt")
#  mycorrelation <- file.path(extdata, "cyp1b1_res37_rawMatrix.txt")
#  
#  chrom <- "chr2"
#  start <- 38290160
#  end <- 38303219
#  gen <- "hg19"
#  strand <- "*"
#  
#  BROWSER.SESSION="UCSC"
#  mySession <- browserSession(BROWSER.SESSION)
#  genome(mySession) <- gen
#  
#  genetrack <-genes_ENSEMBL(gen,chrom,start,end,showId=TRUE)
#  snptrack <- snpBiomart_ENSEMBL(chrom, start, end, dataset="hsapiens_snp_som",showId=FALSE)
#  iscatrack <-ISCA_UCSC(gen,chrom,start,end,mySession, table="iscaPathogenic")
#  
#  listgviz <- list(genetrack,snptrack,iscatrack)
#  
#  matrix.dnamethylation <- read.delim(myinfofile, header=TRUE, sep="\t", as.is=TRUE,
#                                      blank.lines.skip = TRUE, fill=TRUE)
#  matrix.expression <- read.delim(myexpressfile, header=TRUE, sep="\t", as.is=TRUE,
#                                  blank.lines.skip = TRUE, fill=TRUE)
#  cormatrix.data.raw <- read.delim(mycorrelation, sep="\t", header=TRUE, as.is=TRUE,
#                                   blank.lines.skip = TRUE, fill=TRUE)
#  
#  listmatrix.expression <- list(matrix.expression)
#  listcormatrix.data.raw <- list(cormatrix.data.raw)
#  comet(config.file=configfile, mydata.file=matrix.dnamethylation,
#        mydata.type="dataframe",cormatrix.file=listcormatrix.data.raw,
#        cormatrix.type="listdataframe",cormatrix.sig.level=0.05,
#        cormatrix.conf.level=0.05, cormatrix.adjust="BH",
#        mydata.large.file=listmatrix.expression, mydata.large.type="listdataframe",
#        fontsize.gviz =12,
#        tracks.gviz=listgviz,verbose=FALSE, print.image=FALSE)

## ----cometPlotMatrix, echo=FALSE, fig.keep='last'---------------------------------------
extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
configfile <- file.path(extdata, "config_cyp1b1_zoom_4comet.txt")
myinfofile <- file.path(extdata, "cyp1b1_infofile.txt")
myexpressfile <- file.path(extdata, "cyp1b1_infofile_exprGene_region.txt")
mycorrelation <- file.path(extdata, "cyp1b1_res37_rawMatrix.txt")

#configfile <- "../inst/extdata/config_cyp1b1_zoom_4comet.txt" 
chrom <- "chr2"
start <- 38290160
end <- 38303219
gen <- "hg19"
strand <- "*"

data(geneENSEMBLtrack)
data(snpBiomarttrack)
data(ISCAtrack)

listgviz <- list(genetrack,snptrack,iscatrack)

matrix.dnamethylation <- read.delim(myinfofile, header=TRUE, sep="\t", as.is=TRUE, blank.lines.skip = TRUE, fill=TRUE)
matrix.expression <- read.delim(myexpressfile, header=TRUE, sep="\t", as.is=TRUE, blank.lines.skip = TRUE, fill=TRUE)
cormatrix.data.raw <- read.delim(mycorrelation, sep="\t", header=TRUE, as.is=TRUE, blank.lines.skip = TRUE, fill=TRUE)

listmatrix.expression <- list(matrix.expression)
listcormatrix.data.raw <- list(cormatrix.data.raw)
comet(config.file=configfile, mydata.file=matrix.dnamethylation, 
      mydata.type="dataframe",cormatrix.file=listcormatrix.data.raw,
      cormatrix.type="listdataframe",cormatrix.sig.level=0.05,
      cormatrix.conf.level=0.05, cormatrix.adjust="BH",
      mydata.large.file=listmatrix.expression, mydata.large.type="listdataframe",
      fontsize.gviz =12,
      tracks.gviz=listgviz,verbose=FALSE, print.image=FALSE)

## ----cometPlotNopvalText, eval=FALSE, fig.keep='last'-----------------------------------
#  extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
#  configfile <- file.path(extdata, "config_cyp1b1_zoom_4cometnopval.txt")
#  myinfofile <- file.path(extdata, "cyp1b1_infofile.txt")
#  mycorrelation <- file.path(extdata, "cyp1b1_res37_rawMatrix.txt")
#  
#  chrom <- "chr2"
#  start <- 38290160
#  end <- 38303219
#  gen <- "hg19"
#  strand <- "*"
#  
#  genetrack <-genes_ENSEMBL(gen,chrom,start,end,showId=FALSE)
#  snptrack <- snpBiomart_ENSEMBL(chrom, start, end,
#                         dataset="hsapiens_snp_som",showId=FALSE)
#  strutrack <- structureBiomart_ENSEMBL(chrom, start, end,
#                                strand, dataset="hsapiens_structvar_som")
#  clinVariant<-ClinVarMain_UCSC(gen,chrom,start,end)
#  clinCNV<-ClinVarCnv_UCSC(gen,chrom,start,end)
#  gwastrack <-GWAScatalog_UCSC(gen,chrom,start,end)
#  geneRtrack <-GeneReviews_UCSC(gen,chrom,start,end)
#  
#  listgviz <- list(genetrack,snptrack,strutrack,clinVariant,
#                   clinCNV,gwastrack,geneRtrack)
#  comet(config.file=configfile, mydata.file=myinfofile, mydata.type="file",
#        cormatrix.file=mycorrelation, cormatrix.type="listfile",
#        fontsize.gviz =12,
#        tracks.gviz=listgviz, verbose=FALSE, print.image=FALSE,disp.pvalueplot=FALSE)

## ----cometPlot_nopval, echo=FALSE, fig.keep='last'--------------------------------------
extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
configfile <- file.path(extdata, "config_cyp1b1_zoom_4cometnopval.txt")
#configfile <- "../inst/extdata/config_cyp1b1_zoom_4comet.txt" 
myinfofile <- file.path(extdata, "cyp1b1_infofile.txt")
mycorrelation <- file.path(extdata, "cyp1b1_res37_rawMatrix.txt")

chrom <- "chr2"
start <- 38290160
end <- 38303219
gen <- "hg19"
strand <- "*"

data(geneENSEMBLtrack)
data(snpBiomarttrack)
data(strucBiomarttrack)
data(ClinVarCnvTrack)
data(clinVarMaintrack)
data(GWASTrack)
data(GeneReviewTrack)

listgviz <- list(genetrack,snptrack,strutrack,clinVariant,
                 clinCNV,gwastrack,geneRtrack)
comet(config.file=configfile, mydata.file=myinfofile, mydata.type="file",
      cormatrix.file=mycorrelation, cormatrix.type="listfile",
      fontsize.gviz =12,
      tracks.gviz=listgviz, verbose=FALSE, print.image=FALSE,disp.pvalueplot=FALSE)

## ----cometPlotNomatrixtext, eval=FALSE, fig.keep='last'---------------------------------
#  extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
#  configfile <- file.path(extdata, "config_cyp1b1_zoom_4nomatrix.txt")
#  myinfofile <- file.path(extdata, "cyp1b1_infofile.txt")
#  mycorrelation <- file.path(extdata, "cyp1b1_res37_rawMatrix.txt")
#  
#  chrom <- "chr2"
#  start <- 38290160
#  end <- 38303219
#  gen <- "hg19"
#  strand <- "*"
#  
#  genetrack <-genes_ENSEMBL(gen,chrom,start,end,showId=FALSE)
#  snptrack <- snpBiomart_ENSEMBL(chrom, start, end,
#                         dataset="hsapiens_snp_som",showId=FALSE)
#  strutrack <- structureBiomart_ENSEMBL(chrom, start, end,
#                                strand, dataset="hsapiens_structvar_som")
#  clinVariant<-ClinVarMain_UCSC(gen,chrom,start,end)
#  clinCNV<-ClinVarCnv_UCSC(gen,chrom,start,end)
#  gwastrack <-GWAScatalog_UCSC(gen,chrom,start,end)
#  geneRtrack <-GeneReviews_UCSC(gen,chrom,start,end)
#  
#  listgviz <- list(genetrack,snptrack,strutrack,clinVariant,
#                   clinCNV,gwastrack,geneRtrack)
#  comet(config.file=configfile, mydata.file=myinfofile, mydata.type="file",
#        cormatrix.file=mycorrelation, cormatrix.type="listfile",
#        fontsize.gviz =12, font.factor=3,
#        tracks.gviz=listgviz, verbose=FALSE, print.image=FALSE)

## ----cometPlot_nomatrix, echo=FALSE, fig.keep='last'------------------------------------
extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
configfile <- file.path(extdata, "config_cyp1b1_zoom_4nomatrix.txt")
myinfofile <- file.path(extdata, "cyp1b1_infofile.txt")
mycorrelation <- file.path(extdata, "cyp1b1_res37_rawMatrix.txt")

chrom <- "chr2"
start <- 38290160
end <- 38303219
gen <- "hg19"
strand <- "*"

data(geneENSEMBLtrack)
data(snpBiomarttrack)
data(strucBiomarttrack)
data(ClinVarCnvTrack)
data(clinVarMaintrack)
data(GWASTrack)
data(GeneReviewTrack)

listgviz <- list(genetrack,snptrack,strutrack,clinVariant,
                 clinCNV,gwastrack,geneRtrack)

comet(config.file=configfile, mydata.file=myinfofile, mydata.type="file",
      cormatrix.file=mycorrelation, cormatrix.type="listfile",
      fontsize.gviz =12, font.factor=3,
      tracks.gviz=listgviz, verbose=FALSE, print.image=FALSE)

## ----cometlist--------------------------------------------------------------------------
extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
mycorrelation <- file.path(extdata, "cyp1b1_res37_rawMatrix.txt")
myoutput <- file.path(extdata, "cyp1b1_res37_cormatrix_list_BH05.txt")


comet.list(cormatrix.file=mycorrelation,cormatrix.method = "spearman",
           cormatrix.format= "raw", cormatrix.conf.level=0.05,
           cormatrix.sig.level= 0.05, cormatrix.adjust="BH",
           cormatrix.type = "listfile", cormatrix.output=myoutput, 
           verbose=FALSE)

listcorr <- read.csv(myoutput, header = TRUE,
                     sep = "\t", quote = "")
dim(listcorr)
head(listcorr)

## ----GenesAndTranscript,eval=FALSE------------------------------------------------------
#  gen <- "hg38"
#  chr <- "chr15"
#  start <- 75011669
#  end <- 75019876
#  interestfeatures <- rbind(c("75011883","75013394","bad"),c("75013932","75014410","good"))
#  interestcolor <- list("bad"="red", "good"="green")
#  
#  interestgenesENSMBLtrack<-interestGenes_ENSEMBL(gen,chr,start,end,interestfeatures,
#                                                 interestcolor,showId=TRUE)
#  plotTracks(interestgenesENSMBLtrack, from=start, to=end)

## ----GenesAndTranscript2,echo=FALSE-----------------------------------------------------
gen <- "hg38"
chr <- "chr15"
start <- 75011669
end <- 75019876

data(interestgenesENSMBLtrack)
plotTracks(interestgenesENSMBLtrack, from=start, to=end)

## ----Roadmap_example1e,eval=FALSE-------------------------------------------------------
#  chr<-"chr2"
#  start <- 38290160
#  end <- 38303219
#  gen<-"hg19"
#  
#  extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
#  prombedFilePath <- file.path(extdata, "/RoadMap/regions_prom_E001.bed")
#  
#  promRMtrack<- DNaseI_RoadMap(gen,chr,start, end, prombedFilePath,
#                               featureDisplay='promotor', type_stacking="squish")
#  
#  enhbedFilePath <- file.path(extdata, "/RoadMap/regions_enh_E001.bed")
#  
#  enhRMtrack<- DNaseI_RoadMap(gen,chr,start, end, enhbedFilePath,
#                              featureDisplay='enhancer', type_stacking="squish")
#  
#  dyabedFilePath <- file.path(extdata, "/RoadMap/regions_dyadic_E001.bed")
#  
#  dyaRMtrack<- DNaseI_RoadMap(gen,chr,start, end, dyabedFilePath,
#                              featureDisplay='dyadic', type_stacking="squish")
#  
#  genetrack <-genes_ENSEMBL(gen,chr,start,end,showId=TRUE)
#  
#  listRoadMap <- list(genetrack,promRMtrack,enhRMtrack,dyaRMtrack)
#  plotTracks(listRoadMap, chromosome=chr,from=start,to=end)

## ----Roadmap_example1e2,echo=FALSE------------------------------------------------------
chr<-"chr2"
start <- 38290160
end <- 38303219
gen<-"hg19"

data(promRMtrack)
data(enhRMtrack)
data(dyaRMtrack)
data(genetrack4RM)

listRoadMap <- list(genetrack,promRMtrack,enhRMtrack,dyaRMtrack)
plotTracks(listRoadMap, chromosome=chr,from=start,to=end)

## ----encore_genes,eval=FALSE------------------------------------------------------------
#  #Genes from GENCODE
#  chr<-3
#  start <- 132239976
#  end <- 132541303
#  gen<-"hg19"
#  
#  extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
#  gtfFilePath <- file.path(extdata, "/GTEX/gencode.v19.genes.patched_contigs.gtf")
#  options(ucscChromosomeNames=FALSE)
#  grtrack <- GeneRegionTrack(range=gtfFilePath ,chromosome = chr, start= start,
#                             end= end, name = "Gencode V19",
#                             collapseTranscripts=TRUE, showId=TRUE,shape="arrow")
#  plotTracks(grtrack, chromosome=chr,from=start,to=end)

## ----encore_genes2,echo=FALSE-----------------------------------------------------------
#Genes from GENCODE
chr<-3
start <- 132239976
end <- 132541303
gen<-"hg19"

data(genesGencodetrack) 

plotTracks(grtrack, chromosome=chr,from=start,to=end)

## ----encode_example,eval=FALSE----------------------------------------------------------
#  #TF Chip-seq data
#  gen <- "hg19"
#  chr<-"chr1"
#  start <- 1000
#  end <- 329000
#  extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
#  bedFilePath <- file.path(extdata, "ENCODE/motifs1000_matches_ENCODE.txt")
#  motif_color <- file.path(extdata, "ENCODE/TFmotifs_colors.csv")
#  
#  chipTFtrack <- ChIPTF_ENCODE(gen,chr,start, end, bedFilePath,
#                               featureDisplay=c("AHR::ARNT::HIF1A_1",
#                                                "AIRE_1","AIRE_2","AHR::ARNT_1"),
#                               motif_color,type_stacking="squish",showId=TRUE)
#  
#  plotTracks(chipTFtrack, chromosome=chr,from=start,to=end)

## ----encode_example2,echo=FALSE---------------------------------------------------------
#TF Chip-seq data
gen <- "hg19"
chr<-"chr1"
start <- 1000
end <- 329000

data(chipTFtrack)
plotTracks(chipTFtrack, from = start, to = end)

## ----GTEX_example1,eval=FALSE-----------------------------------------------------------
#  ## eQTL data
#  chr<-"chr3"
#  start <- 132239976
#  end <- 132541303
#  gen<-"hg19"
#  
#  extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
#  bedFilePath <- file.path(extdata, "/GTEX/eQTL_Uterus_Analysis_extract100.snpgenes")
#  
#  eGTex<- eQTL_GTEx(gen,chr, start, end, bedFilePath, featureDisplay = 'all',
#                    showId=TRUE, type_stacking="squish", just_group="left" )
#  
#  eGTex_SNP<- eQTL_GTEx(gen,chr, start, end, bedFilePath,
#                        featureDisplay = 'SNP', showId=FALSE,
#                        type_stacking="dense", just_group="left")
#  
#  #Genes from
#  gtfFilePath <- file.path(extdata, "/GTEX/gencode.v19.genes.patched_contigs.gtf")
#  options(ucscChromosomeNames=FALSE)
#  grtrack <- GeneRegionTrack(genome="hg19",range=gtfFilePath ,chromosome = chr,
#                             start= start, end= end, name = "Gencode V19",
#                             collapseTranscripts=TRUE, showId=TRUE,shape="arrow")
#  eGTexTracklist <- list(grtrack,eGTexTrackSNP)
#  plotTracks(eGTexTracklist, chromosome=chr,from=start,to=end)

## ----GTEX_example12,echo=FALSE----------------------------------------------------------
## eQTL data
chr<-"chr3"
start <- 132239976
end <- 132541303
gen<-"hg19"

data(eGTexTrackSNP)
data(eGTexTrackall)
data(grtrack4eGTex)
#Genes from

eGTexTracklist <- list(grtrack,eGTexTrackSNP)
plotTracks(eGTexTracklist, chromosome=chr,from=start,to=end)

## ----GTEX_example4e,eval=FALSE----------------------------------------------------------
#  ### psiQTL
#  chr<-"chr13"
#  start <- 52713837
#  end <- 52715894
#  gen<-"hg19"
#  
#  extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
#  psiQTLFilePath <- file.path(extdata, "/GTEX/psiQTL_Assoc-total.AdiposeTissue.txt")
#  
#  psiGTex<- psiQTL_GTEx(gen,chr,start, end, psiQTLFilePath, featureDisplay = 'all',
#                        showId=TRUE, type_stacking="squish",just_group="above" )
#  
#  genetrack <-genes_ENSEMBL(gen,chr,start,end,showId=TRUE)
#  
#  psiTrack <- list(genetrack,psiGTex)
#  plotTracks(psiTrack, chromosome=chr,from=start,to=end)

## ----GTEX_example4e2,echo=FALSE---------------------------------------------------------
### psiQTL 
chr<-"chr13"
start <- 52713837
end <- 52715894
gen<-"hg19"

data(psiGTexTrackall)
data(genetrack4psiGTEX)

psiTrack <- list(genetrack,psiGTexTrackall)

plotTracks(psiTrack, chromosome=chr,from=start,to=end)

## ----TEX_example5e1---------------------------------------------------------------------
data(imprintedGenesGTEx)
as.character(unique(imprintedGenesGTEx$Tissue.Name))
as.character(unique(imprintedGenesGTEx$Classification))

## ----GTEX_example5e,eval=FALSE----------------------------------------------------------
#  ### inprinted genes
#  chr<- "chr1"
#  start <- 7895752
#  end <- 7914572
#  gen<-"hg19"
#  
#  genesTrack <- genes_ENSEMBL(gen,chr,start,end,showId=TRUE)
#  
#  allIG <- imprintedGenes_GTEx(gen,chr,start, end, tissues="all",
#                               classification="imprinted",showId=TRUE)
#  
#  allimprintedIG <- imprintedGenes_GTEx(chr,start, end, tissues="all",
#                                        classification="imprinted",showId=TRUE)
#  
#  StomachIG <-imprintedGenes_GTEx(gen,chr,start, end, tissues="Stomach",
#                                  classification="all",showId=TRUE)
#  
#  PancreasIG <- imprintedGenes_GTEx(gen,chr,start, end,
#                                    tissues="Pancreas",
#                                    classification="all",showId=TRUE)
#  PancreasimprintedIG <- imprintedGenes_GTEx(gen,chr,start, end, tissues="Pancreas",
#                                             classification="imprinted",showId=TRUE)
#  
#  plotTracks(list(genesTrack, allIG, allimprintedIG,
#                  StomachIG,PancreasIG,PancreasimprintedIG),
#             chromosome=chr, from=start, to=end)
#  

## ----GTEX_example5e2,echo=FALSE---------------------------------------------------------
### inprinted genes
chr<- "chr1"
start <- 7895752
end <- 7914572
gen<-"hg19"

genesTrack <- genes_ENSEMBL(gen,chr,start,end,showId=TRUE)

data(allIGtrack)
data(allimprintedIGtrack)
data(StomachIGtrack)
data(PancreasIGtrack)
data(PancreasimprintedIGtrack)
imprintinglist <- list(genesTrack,allIGtrack,allimprintedIGtrack,StomachIGtrack,
                       PancreasIGtrack,PancreasimprintedIGtrack)
plotTracks(imprintinglist, from = start, to = end)


## ----GTEX_example3,eval=FALSE,echo=FALSE------------------------------------------------
#  ## sQTL data
#  chr<-3
#  start <- 132239976
#  end <- 132541303
#  gen<-"hg19"
#  
#  extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
#  
#  ### sQTL from Altran methods
#  sQTLAFilePath <- file.path(extdata, "/GTEX/sQTL_HeartLeftVentricle.Altrans.FDR05.bestPerLink")
#  
#  sAGTex<- sQTL_Altrans_GTEx(gen,chr,start, end,
#                             sQTLAFilePath, featureDisplay = 'all', showId=TRUE,
#                             type_stacking="squish",just_group="left" )

## ----extract_data_1kb,eval=FALSE--------------------------------------------------------
#  library('corrplot')
#  #Hi-C data
#  gen <- "hg19"
#  chr<-"chr1"
#  start <- 5000000
#  end <- 9000000
#  extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
#  bedFilePath <- file.path(extdata, "HiC/chr1_1mb.RAWobserved")
#  
#  matrix_HiC <- HiCdata2matrix(chr,start, end, bedFilePath)
#  cor_matrix_HiC <- cor(matrix_HiC)
#  diag(cor_matrix_HiC)<-1
#  corrplot(cor_matrix_HiC, method = "circle")

## ----extract_data_1kb2,echo=FALSE-------------------------------------------------------
library('corrplot')
#Hi-C data
gen <- "hg19"
chr<-"chr1"
start <- 5000000
end <- 9000000

data(matrix_HiC_Rao)
cor_matrix_HiC <- cor(matrix_HiC_Rao)
diag(cor_matrix_HiC)<-1
corrplot(cor_matrix_HiC, method = "circle")

## ----extract_data,eval=FALSE------------------------------------------------------------
#  extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
#  info_HiC <- file.path(extdata, "Human_IMR90_Fibroblast_topological_domains.txt")
#  data_info_HiC <-read.csv(info_HiC, header = FALSE, sep = "\t", quote = "")
#  
#  intrachr_HiC <- file.path(extdata, "Human_IMR90_Fibroblast_Normalized_Matrices.txt")
#  data_intrachr_HiC <- read.csv(intrachr_HiC, header = TRUE, sep = "\t", quote = "")
#  
#  chr_interest <- "chr2"
#  start_interest <- "1"
#  end_interest <- "160000"
#  list_bins <- which(data_info_HiC[,1] == chr_interest &
#                       data_info_HiC[,2] >= start_interest &
#                       data_info_HiC[,2] <= end_interest )
#  
#  subdata_info_Hic <- data_info_HiC[list_bins,]
#  subdata_intrachr_HiC <- data_intrachr_HiC[list_bins,list_bins]
#  

## ----FANTOM_example,eval=FALSE----------------------------------------------------------
#  gen <- "hg19"
#  chr<- "chr1"
#  start <- 6000000
#  end <- 6500000
#  
#  extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
#  
#  ##Enhancer
#  enhFantomFile <- file.path(extdata, "/FANTOM/human_permissive_enhancers_phase_1_and_2.bed")
#  enhFANTOMtrack <-DNaseI_FANTOM(gen,chr,start, end, enhFantomFile, featureDisplay='enhancer')
#  
#  ### TFBS motif
#  AP1FantomFile <- file.path(extdata, "/FANTOM/Fantom_hg19.AP1_MA0099.2.sites.txt")
#  tfbsFANTOMtrack <- TFBS_FANTOM(gen,chr,start, end, AP1FantomFile)

## ----FANTOM_example2,echo=FALSE---------------------------------------------------------
gen <- "hg19"
chr<- "chr1"
start <- 6000000
end <- 6500000

extdata <- system.file("extdata", package="coMET",mustWork=TRUE)

##Enhancer
data(enhFANTOMtrack)

### TFBS motif
data(tfbsFANTOMtrack)

Fantom5list <- list(enhFANTOMtrack,tfbsFANTOMtrack)
plotTracks(Fantom5list, from = start, to = end)

## ----faq0,eval=FALSE--------------------------------------------------------------------
#  genetrack <-genesENSEMBL(gen,chrom,start,end,showId=TRUE)
#  
#  plotTracks(genetrack)
#  
#  str(genetrack)

## ----faq1,eval=FALSE--------------------------------------------------------------------
#  genetrack <-genesENSEMBL(gen,chrom,start,end,showId=TRUE)
#  
#  displayPars(genetrack)

## ----faq2,eval=FALSE--------------------------------------------------------------------
#  comet(config.file = configfile, mydata.file = myinfofile, mydata.format = "file",
#        cormatrix.file = mycorrelation, cormatrix.type = "listfile",
#        mydata.large.file = mylargedata,mydata.large.type = "listfile",
#        tracks.gviz = listGviz, verbose = TRUE,
#        print.image=TRUE,fontsize.gviz=10)

## ----faq3,eval=FALSE--------------------------------------------------------------------
#  geneTrack <- refGenesUCSC(gen, chr, start, end, IdType ="name", showId = TRUE)
#  
#  geneTrackShow <- geneTrack[gene(geneTrack) %in% c("AHRR")]

## ----session-info, results="asis"-------------------------------------------------------
toLatex(sessionInfo())

## ----Development, eval=FALSE,echo=FALSE-------------------------------------------------
#  #Need to have the last version of R associated with Bioconductor in parallele of your original version
#  #install from source  http://bioconductor.org/developers/how-to/useDevel/
#  mkdir /home/tiphaine/Rdevel2.2/
#  cd Rdevel2.2
#  mkdir tar
#  #Extract R source in this latter folder
#  ./configure -prefix=/home/tiphaine/Rdevel2.2/
#  
#  sudo make
#  
#  sudo make install
#  
#  #Need to have the last version of Bioconductor and BiocCheck
#  #Need to update under Rdev
#  remove.packages("BiocInstaller")
#  source("http://bioconductor.org/biocLite.R")
#  biocLite("BiocInstaller")
#  
#  #Need to update different packages associated with coMET
#  biocLite("coMET")
#  
#  #Go to the parent folder of the package
#  # To create the manual documentation, need to run
#  Rdevel CMD Rd2pdf --pdf coMET
#  
#  #Need to build the  new package
#  R-devel CMD build coMET --resave-data --no-build-vignettes
#  
#  #Need to check if the new package follow the rules of R
#  R-devel CMD check coMET
#  
#  #To check if the new package follow the rules from bioconductor
#  R-devel CMD BiocCheckdev3.1 coMET

