#' Import genomic signals
#' @description
#' Import genomic signals for given range.
#' @param paths Vector of character. File path
#' @param range GRanges object. The coordinates.
#' @param cols The colors for each signal.
#' @param format The format of the file.
#' @return A list of track object
#' @export
#' @importFrom methods is
#' @importFrom trackViewer importScore setTrackStyleParam
#' @examples
#' # example code
#' 

importGenomicSigs <- function(paths, range, cols, format='BigWig'){
  stopifnot(is.character(paths) || is.character(unlist(paths)))
  stopifnot(is(range, 'GRanges'))
  mapply(paths, cols, format, FUN=function(p, color, f){
    dat <- importScore(p, format=f, ranges=range)
    setTrackStyleParam(dat, 'color', color)
    dat
  })
}

#' @importFrom BiocFileCache BiocFileCache getBFCOption bfcrpath
downloadFiles <- function(urls){
  stopifnot(is.character(urls))
  bfc <- BiocFileCache(cache = getBFCOption("CACHE"), ask = interactive())
  cpath <- lapply(urls, bfcrpath, x=bfc)
  return(cpath)
}

#' Import genomic signals from URL
#' @description
#' Download files and import signals.
#' @param urls Vector of character. URLs.
#' @param range GRanges object. The coordinates.
#' @param cols The colors for each signal.
#' @param format The format of the file.
#' @return A list of track object
#' @importFrom methods is
#' @export
#' @examples
#' # example code
#' 
importSignalFromUrl <- function(urls, range, cols, format='BigWig'){
  stopifnot(is.character(urls))
  stopifnot(is(range, 'GRanges'))
  cpath <- downloadFiles(urls)
  sigs <- importGenomicSigs(cpath, range = range, cols = cols, format = format)
}

#' Import FLAMINGO output
#' @description
#' Import FLAMINGO output to a GRangesList Object
#' @param filenames Vector of character. The file names of output of FLAMINGO.
#' @return An object of GRangesList
#' @importFrom GenomicRanges GRanges GRangesList
#' @importFrom IRanges IRanges
#' @importFrom utils read.delim
#' @export
#' @examples
#' # example code
#' 
importFLAMINGO <- function(filenames){
  dat <- lapply(filenames, function(.ele) read.delim(.ele, header = FALSE))
  gr <- lapply(dat, function(.ele){
    if(.ele[1, 1]=='chr'){
      .ele <- .ele[-1, ]
    }
    with(.ele, GRanges(V1, IRanges(as.numeric(V2)+1, as.numeric(V3)), 
                       x=as.numeric(V4), y=as.numeric(V5), z=as.numeric(V6)))
  })
  GRangesList(gr)
}

#' Import genomic interactions from URL
#' @description
#' Download files and import interactions
#' @param urls Vector of character. URLs.
#' @param range GRanges object. The coordinates.
#' @param resolution The resolution of the matrix.
#' @param format The format of the file.
#' @param normalization Normalization matrix in the file.
#' @return A list of GInteractions object
#' @export
#' @importFrom trackViewer importGInteractions
#' @importFrom utils read.delim
#' @examples
#' # example code
#' 
importGInteractionsFromUrl <- function(urls, resolution, range,
                                       format='cool', normalization='balanced'){
  stopifnot(is.character(urls))
  stopifnot(is(range, 'GRanges'))
  stopifnot(is.numeric(resolution))
  cpath <- downloadFiles(urls)
  mapply(cpath, format, resolution, normalization, FUN=function(p, f, r, n){
    dat <- importGInteractions(p, format=f, resolution = r,
                               ranges=range, normalization=n,
                               out = 'GInteractions')
    dat
  })
}