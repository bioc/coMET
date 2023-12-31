\name{eQTL_GTEx}
\alias{eQTL_GTEx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates a eQTL track from GTEx
}
\description{
Creates a track of eQTL from GTEx using the Gviz bioconductor package. A complete list of features and their associated colours can be found in the user guide.
}
\usage{
eQTL_GTEx(gen="hg19",chr,start, end, bedFilePath, featureDisplay = 'all',
showId=FALSE, type_stacking="squish",just_group="above",title="eQTL GTEX")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{gen}{
the name of the genome. Default value=hg19
}
  \item{chr}{
The chromosome of interest
}
  \item{start}{
The starting position in the region of interest (the smallest value)
}
  \item{end}{
The end position in the region of interest (the largest value)
}
 \item{bedFilePath}{
The path of the BED file from Kheradpour and Kellis, 2014.
}
  \item{featureDisplay}{
A vector of regulatory features to be displayed, such as Predicted heterochomatin. Spelling and capitalisation of features must be identical to those in the user guide. There are three possibilities. First, the visualisation of only one feature (e.g. featureDisplay <- "Predicted heterochomatin"), only the name of the specific feature is required. Second, visualisation of a set of features, for this a vector of features is required (e.g. featureDisplay <- c("Predicted low activity","Predicted heterochomatin")). Finally, visualison all features in the genomic region, achived by using the word "all" (e.g. featureDisplay <- "all"), "all" is set by default. You can find the complete list of features and their associated colours in the user guide.
}
\item{showId}{
logical. say if we write the name of group
}
  \item{type_stacking}{
Object of class"character", the stacking type of overlapping items on the final plot.One in c(hide, dense, squish, pack,full). More information cf the option "stacking" in Gviz
}
\item{just_group}{
position. say where we write the name of group (choice in c("above","righ","left"))
}
  \item{title}{
 The name of the annotation track
}
}
%% \details{
%%  ~~ If necessary, more details than the description above ~~
%% }
\value{
An AnnotationTrack object of Gviz
}
\references{
http://bioconductor.org/packages/release/bioc/html/Gviz.html

Got to BindingMotifsBiomart binding motif biomart
}
\author{
Tiphaine Martin
}
%% \note{
%%  ~~further notes~~
%% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

%% \seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%% }
\examples{
library("Gviz")
gen <- "hg19"
chr<-"chr3"
start <- 132423172
end <- 132442807
extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
bedFilePath <- file.path(extdata, "/GTEX/eQTL_Uterus_Analysis_extract100.snpgenes")

%all
if(interactive()){
 eGTexTrackall <- eQTL_GTEx(gen,chr,start, end, bedFilePath,
 featureDisplay="all", showId=TRUE,just_group="left")
  plotTracks(eGTexTrackall, from = start, to = end,
                                   fontfamily="sans",fontfamily.title="sans")
} else {
  data(eGTexTrackall)
  plotTracks(eGTexTrackall, from = start, to = end,
                                   fontfamily="sans",fontfamily.title="sans")
}

%SNP
if(interactive()){
 eGTexTrackSNP <- eQTL_GTEx(gen,chr,start, end, bedFilePath,
 featureDisplay="SNP", showId=TRUE,just_group="left")
  plotTracks(eGTexTrackSNP, from = start, to = end,
                                   fontfamily="sans",fontfamily.title="sans")
} else {
  data(eGTexTrackSNP)
  plotTracks(eGTexTrackSNP, from = start, to = end,
                                   fontfamily="sans",fontfamily.title="sans")
}

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
