\name{cpgIslands_UCSC}
\alias{cpgIslands_UCSC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
create track CpG Island from UCSC
}
\description{
create track CpG Island from UCSC using the Gviz bioconductor package
}
\usage{
cpgIslands_UCSC(gen, chr, start, end, title="CpG Islands UCSC")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{gen}{
the name of the genome
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
    \item{title}{
    Name of tracks
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
An UcscTrack object of Gviz
}
\references{
http://bioconductor.org/packages/release/bioc/html/Gviz.html

http://genome-euro.ucsc.edu/cgi-bin/hgTrackUi?hgsid=202839739_2hYQ1BAOuBMAR620GjrtdrFAy6dn&c=chr6&g=cpgIslandSuper
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
library("Gviz")
chrom <- "chr2"
start <- 100000
end <- 1000000
gen <- "hg38"

if(interactive()) {
    cpgIstrack<-cpgIslands_UCSC(gen, chrom, start, end)
    plotTracks(cpgIstrack, from = start, to =end,
    fontfamily="sans",fontfamily.title="sans")
}else {
    data(cpgIslandtrack)
    plotTracks(cpgIstrack, from = start, to =end,
    fontfamily="sans",fontfamily.title="sans")
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
