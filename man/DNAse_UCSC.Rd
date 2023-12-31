\name{DNAse_UCSC}
\alias{DNAse_UCSC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creation of an UCSC's DNase clusters track - obselete function
}
\description{
Creation of DNase cluster track from a connection to UCSC genome browser in using the GViz bioconductor package. Obselete function 
}
\usage{
DNAse_UCSC(gen, chr, start, end, mySession, title="DNA cluster",
track.name = "DNase Clusters", table.name = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gen}{
the name of the genome. Data is not currently available for GRCh38 (hg38).
}
  \item{chr}{
the chromosome of interest
}
  \item{start}{
the first position in the region of interest (the smallest value)
}
  \item{end}{
the last position in the region of interest (the largest value)
}
  \item{mySession}{
the object session from the function browserSession of rtracklayer
}
    \item{title}{
    Name of tracks
}
  \item{track.name}{
the name of the track DNAse_UCSC. "DNase Clusters"(default)
}
  \item{table.name}{
the name of the table from the track
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
An AnnotationTrack object of Gviz
}
\references{
http://bioconductor.org/packages/release/bioc/html/Gviz.html

http://genome-euro.ucsc.edu/cgi-bin/hgTrackUi?hgsid=202839739_2hYQ1BAOuBMAR620GjrtdrFAy6dn&c=chr6&g=wgEncodeDNAseSuper
}
\author{
Tiphaine Martin
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{
# library("Gviz")
# library("rtracklayer")

# gen <- "hg19"
# chr <- "chr7"
# start <- 38290160
# end <- 38303219
# if(interactive()){
#   BROWSER.SESSION="UCSC"
#   mySession <- browserSession(BROWSER.SESSION)
#   genome(mySession) <- gen
#   track.name="Broad ChromHMM"
#   tablestrack<-tableNames(ucscTableQuery(mySession, track=track.name))
#   table.name<-tablestrack[1]
#   dnasetrack<-DNAse_UCSC(gen,chr,start,end,mySession)
#   plotTracks(dnasetrack, from = start, to =end,
#                                   fontfamily="sans",fontfamily.title="sans")
# }else {
#     data(dnasetrack)
#    plotTracks(dnasetrack, from = start, to =end,
#                                   fontfamily="sans",fontfamily.title="sans")
# }

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }

