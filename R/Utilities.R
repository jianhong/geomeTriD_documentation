#' Import genomic signals
#' @description
#' Import genomic signals for given range.
#' @param paths Vector of character. File path
#' @param range \link[GenomicRanges:GRanges-class]{GRanges} object. The coordinates.
#' @param cols The colors for each signal.
#' @param format The format of the file.
#' @return A list of track object
#' @export
#' @importFrom methods is
#' @importFrom trackViewer importScore setTrackStyleParam
#' @importFrom R.utils gunzip
#' @examples
#' # example code
#' 

importGenomicSigs <- function(paths, range, cols, format='BigWig'){
  stopifnot(is.character(paths) || is.character(unlist(paths)))
  stopifnot(is(range, 'GRanges'))
  mapply(paths, cols, format, FUN=function(p, color, f){
    if(grepl('.gz$', p)){
      tmp <- sub('.gz$', '', p)
      tryCatch(gunzip(p, destname=tmp, remove=FALSE), error=function(.e){
        message(.e, 'File may be decompressed already.')
      })
      p <- tmp
    }
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
#' @param range \link[GenomicRanges:GRanges-class]{GRanges} object. The coordinates.
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

#' Import Tensor-FLAMINGO output
#' @description
#' Import Tensor-FLAMINGO output to a GRangesList Object
#' @param filenames Vector of character. The file names of output of Tensor-FLAMINGO.
#' @param binsize Bin size of the interaction. It will be the width of 
#' GRanges. 
#' @param chr Chromosome name
#' @return An object of GRangesList
#' @importFrom GenomicRanges GRanges GRangesList
#' @importFrom IRanges IRanges
#' @importFrom utils read.delim
#' @export
#' @examples
#' # example code
#' 
importTensorFLAMINGO <- function(filenames, binsize, chr){
  stopifnot(is.numeric(binsize))
  stopifnot(is.character(chr))
  dat <- lapply(filenames, function(.ele) read.delim(.ele, header = FALSE))
  gr <- lapply(dat, function(.ele){
    with(.ele, GRanges(chr, IRanges(as.numeric(V1)*binsize+1, width = binsize), 
                       x=as.numeric(V2), y=as.numeric(V3), z=as.numeric(V4)))
  })
  GRangesList(gr)
}

#' Import SuperRec output
#' @description
#' Import SuperRec output to a GRangesList Object
#' @param filenames Vector of character. The file names of output of SuperRec
#' @param superRecInputFilenames Vector of character.
#'  The filenames of input of SuperRec.
#' @param binsize Bin size of the interaction. It will be the width of 
#' GRanges. 
#' @param chr Chromosome name
#' @return An object of GRangesList
#' @importFrom GenomicRanges GRanges GRangesList
#' @importFrom IRanges IRanges
#' @importFrom utils read.delim read.table
#' @export
#' @examples
#' # example code
#' 
importSuperRec <- function(filenames, superRecInputFilenames, binsize, chr){
  stopifnot(is.numeric(binsize))
  stopifnot(is.character(chr))
  stopifnot(length(filenames)==length(superRecInputFilenames))
  dat <- lapply(filenames, function(.ele) read.delim(.ele, header = FALSE))
  input <- lapply(superRecInputFilenames, function(.ele) read.table(.ele, header = FALSE))
  gr <- mapply(dat, input, FUN=function(.ele, .loci){
    .loci <- sort(unique(c(.loci[, 1], .loci[, 2])))
    with(.ele, GRanges(chr, IRanges(.loci[seq.int(nrow(.ele))]*binsize+1, width=binsize), 
                       x=as.numeric(V1), y=as.numeric(V2), z=as.numeric(V3)))
  }, SIMPLIFY = FALSE)
  GRangesList(gr)
}

#' Import Hickit output
#' @description
#' Import 3dg file from Hickit output to a GRangesList Object
#' @param filenames Vector of character. The 3dg file names of output of Hickit
#' @param comment.char,... parameters pass to \link[utils:read.table]{rea.delim}.
#' @param parental_postfix the postfix of chromosome names. Default is 
#' c("(pat)", "(mat)"). 
#' For hickit, it maybe c('a', 'b').
#' @return An object of GRangesList
#' @importFrom GenomicRanges GRanges GRangesList
#' @importFrom IRanges IRanges
#' @importFrom utils read.delim
#' @export
#' @examples
#' f3dg <- system.file('extdata', 'GSE162511',
#'  'GSM4382149_cortex-p001-cb_001.20k.1.clean.3dg.txt.gz',
#'   package='geomeTriD.documentation')
#' xyz <- import3dg(f3dg)
import3dg <- function(filenames, comment.char="#", ...,
                      parental_postfix=c("(pat)", "(mat)")){
  stopifnot(length(parental_postfix)==2)
  stopifnot(is.character(parental_postfix))
  dat <- lapply(filenames, function(.ele) read.delim(.ele, header = FALSE,
                                                     comment.char=comment.char,
                                                     ...))
  gr <- lapply(dat, function(.ele){
    .ele <- split(.ele, .ele$V1)
    .ele <- lapply(.ele, function(.e){
      w <- .e$V2[-1] - .e$V2[-length(.e$V2)]
      w <- c(w, w[length(w)])
      .e$wid <- w
      .e
    })
    .ele <- do.call(rbind, .ele)
    
    parental1 <- grepl(parental_postfix[1], .ele[, 1])
    parental2 <- grepl(parental_postfix[2], .ele[, 1])
    if(any(parental1) | any(parental2)){
      .ele$parental[parental1] <- parental_postfix[1]
      .ele$parental[parental2] <- parental_postfix[2]
      .ele$parental[!(parental1|parental2)] <- NA
      parental_postfix <- gsub('\\(|\\)', '.', parental_postfix)
      .ele$V1 <- sub(paste0(parental_postfix[2], '$'), '',
                     sub(paste0(parental_postfix[1], '$'), '',
                           .ele$V1))
    }else{
      .ele$parental <- NA
    }
    
    with(.ele, GRanges(V1, IRanges(as.numeric(V2)+1, width=wid), 
                       x=as.numeric(V3), y=as.numeric(V4), z=as.numeric(V5),
                       parental=parental))
  })
  GRangesList(gr)
}

#' Import genomic interactions from URL
#' @description
#' Download files and import interactions
#' @param urls Vector of character. URLs.
#' @param range \link[GenomicRanges:GRanges-class]{GRanges} object. The coordinates.
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

#' Create annotation features
#' @description
#' Generate annotation features (feature.gr for \link[geomeTriD]{view3dStructure})
#' from TxDb and Org object.
#' @param txdb An \link[GenomicFeatures:TxDb-class]{TxDb} object.
#' @param org An \link[AnnotationDbi:AnnotationDb-class]{OrgDb} object.
#' @param range \link[GenomicRanges:GRanges-class]{GRanges} object. The coordinates.
#' @param geneSymbolColumn,keytype The column names in the OrgDb for key of TxDb
#' and the gene symnbols.
#' @param cols The colors for each gene.
#' @return A \link[GenomicRanges:GRanges-class]{GRanges} object
#' @export
#' @importFrom AnnotationDbi select
#' @importFrom GenomicFeatures genes
#' @importFrom IRanges subsetByOverlaps
#' @examples
#' # example code
#' 
getFeatureGR <- function(txdb, org, range,
                         keytype='ENTREZID',
                         geneSymbolColumn = 'SYMBOL',
                         cols=seq.int(7)){
  stopifnot(is(txdb, "TxDb"))
  stopifnot(is(org, 'OrgDb'))
  stopifnot(is(range, 'GRanges'))
  #### get all genes
  feature.gr <- genes(txdb)
  #### subset the data by target viewer region
  feature.gr <- subsetByOverlaps(feature.gr, range)
  if(length(feature.gr)){
    #### assign symbols for each gene
    symbols <- select(org, 
                      feature.gr$gene_id,
                      columns= geneSymbolColumn,
                      keytype = keytype)
    feature.gr$label <- symbols[match(feature.gr$gene_id, symbols[, keytype]),
                                geneSymbolColumn]
    #### assign colors for each gene
    feature.gr$col <- sample(cols, length(feature.gr), replace = TRUE)
    feature.gr$type <- 'gene'
  }
  return(feature.gr)
}

#' Show a pair of geometries
#' @description
#' Show a pair of geometries side by side.
#' @param a,b A list of \link[geomeTriD:threeJsGeometry-class]{threeJsGeometry} object.
#' @param height The height of the widgets.
#' @param ... The parameter for \link[geomeTriD:threeJsViewer]{threeJsViewer}.
#' @return A \link[GenomicRanges:GRanges-class]{GRanges} object
#' @export
#' @importFrom geomeTriD threeJsViewer
#' @examples
#' # example code
#' 
showPairs <- function(a, b, height = '95vh', ...){
  b <- lapply(b, function(.ele) {
    .ele$side = 'right'
    .ele
  })
  threeJsViewer(a, b, height=height, ...)
}

#' Padding GRanges List
#' @description
#' Make all GRanges object in a list same length by filling with NA
#' @param xyzs A list of GRanges.
#' @return A list with element 'gr' represent the genomic ranges and
#' 'xyzs' represent the x,y,z positions in a data.frame.
#' @export
#' @importFrom GenomicRanges GRangesList
#' @importFrom IRanges findOverlaps disjoin slidingWindows
#' @importFrom S4Vectors queryHits subjectHits mcols mcols<-
#' 
paddingGRangesList <- function(xyzs){
  if(!is(xyzs, 'GRangesList')) xyzs <- GRangesList(xyzs)
  gr <- unlist(xyzs)
  gr <- sort(disjoin(unique(gr)))
  width <- table(width(gr))
  width <- sort(width, decreasing = TRUE)
  width <- as.integer(names(width)[1])
  if(any(width(gr)<width)) stop('can not handle unfixed bin size.')
  gr.m <- gr[width(gr)==width]
  gr.l <- gr[width(gr)>width]
  gr.l <- slidingWindows(gr.l, width = width, step = width)
  gr <- sort(disjoin(c(gr.m, unlist(gr.l))))
  colnames(mcols(gr)) <- tolower(colnames(mcols(gr)))
  mcols(gr)$x <- NA
  mcols(gr)$y <- NA
  mcols(gr)$z <- NA
  names(gr) <- NULL
  xyzs <- lapply(xyzs, function(.ele){
    colnames(mcols(.ele)) <- tolower(colnames(mcols(.ele)))
    ol <- findOverlaps(gr, .ele, minoverlap = width)
    newEle <- gr
    mcols(newEle[queryHits(ol)]) <- mcols(.ele[subjectHits(ol)])
    as.data.frame(mcols(newEle))
  })
  return(list(gr=gr, xyzs=xyzs))
}

#' Aggregate xyz values by xyz group
#' @description
#' Aggregate the corresponding xyz values of xyz values for each xyz group
#' by given function.
#' @param xyz.list A list of list with xyz values.
#' @param FUN The function for aggregate the values.
#' @param na.rm Remove the NA values or not. If TRUE, NA values will not be 
#' considered. Otherwise, the NA values will be filled with the mean of
#' nearby values.
#' @param ... Other parameters for the FUN.
#' @return A list of xyz.
#' @importFrom geomeTriD fill_NA
#' @export
#' 
aggregateXYZs <- function(xyz.list, FUN=mean, na.rm=FALSE, ...){
  stopifnot(is.list(xyz.list))
  cn <- c('x', 'y', 'z')
  xyz.list <- lapply(xyz.list, function(xyzs){
    lapply(xyzs, function(.ele){
      if(!(is.matrix(.ele) || is.data.frame(.ele))){
        stop("The elements for each group in xyz.list must be a matrix or data.frame.")
      }
      colnames(.ele) <- tolower(colnames(.ele))
      if(!all(cn %in% colnames(.ele))){
        stop('The elements for each group in xyz.list must contain colnames "x", "y" and "z".')
      }
      .ele
    })
  })
  if(!na.rm){
    xyz.list <- lapply(xyz.list, function(.ele) lapply(.ele, fill_NA))
  }
  
  cg_xyz <- lapply(xyz.list, function(.ele){
    xyz <- lapply(cn, function(coln){
      apply(do.call(cbind, lapply(.ele, function(.e){
        .e[, coln]
      })), 1, FUN, na.rm = na.rm, ...)
    })
    .ele <- as.data.frame(.ele[[1]])
    .ele$x <- xyz[[1]]
    .ele$y <- xyz[[2]]
    .ele$z <- xyz[[3]]
    .ele
  })
}

#' Find the cluster medoid sample
#' @description
#' Find the medoid sample of each cluster in a given distance matrix
#' @param dist A dist object
#' @param cluster The group information for the elements in the dist.
#' @param N The number of medoid sample to be returned.
#' @param na.rm A logical evaluating to TRUE or FALSE indicating whether NA 
#' values should be stripped before the computation proceeds. 
#' @return The medoid element labels for each cluster.
#' @export
#' @importFrom utils head
#' @examples
#' x <- matrix(rnorm(100), nrow = 5)
#' dist <- dist(x)
#' group <- c('A', 'A', 'B', 'B', 'B')
#' getMedoidId(dist, group)
getMedoidId <- function(dist, cluster, N=1L, na.rm=TRUE){
  stopifnot(is(dist, 'dist'))
  stopifnot(is.numeric(N))
  n <- attr(dist, 'Size')
  stopifnot(n==length(cluster))
  dist <- split(as.data.frame(as.matrix(dist)), cluster)
  idx <- lapply(dist, function(.ele){
    rs <- rowSums(.ele[, rownames(.ele), drop=FALSE], na.rm = na.rm)
    rs <- sort(rs, decreasing = FALSE)
    names(head(rs, n=N))
  })
  return(idx)
}