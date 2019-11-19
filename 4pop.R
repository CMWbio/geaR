#' 4 pop outgroup statistics
#'
#' @description Imports genotypes across a locus into matrix
#'
#' @details Authours: Chris Ward & Alastair Ludington
#' Uses a GRanges locus to import genotypes, either nucleotide or RAW, from a GDS file
#'
#'
#' @param GDS \code{GDS} object with variant data to import genotypes from
#' @param locus \code{GRanges} Locus to import genotypes for
#' @param minSites \code{numeric} minimum number of sites as a proportion of loci length
#' @param nucleotide \code{logical} Import RAW genotypes or nucleotides
#' @param ploidy \code{numeric} ploidy of sample
#' @param pops \code{data_frame} populaiton dataFrame
#' @param removeIndels removes indels
#'
#'
#' @return A \code{matrix} of genotypes
#'
#'
#' @examples
#'
#'
#'
#' @rdname getFStats






getFStats <- function(GDS, loci, pops, tests, nCores, minSites){


if(length(pops) == 0){
  samples <- seqGetData(gdsfile = GDS, var.name = "sample.id")
  pops <- data_frame(Sample = samples, Population = samples)
}

popList <- split(pops, pops$Population)

stopifnot(minSites < 1 & minSites > 0)

plan(multiprocess, workers = nCores)

fstats <- future_map(seq(length(loci)), function(y){

  locus <- loci[[y]]

  testD <- map(tests, function(x){

    if(length(x) == 4) {

    }
    if(length(x) == 3) fun <- "threePop"

    f <- .fourPop(GDS, locus, pops, x)

    scaf <- locus@seqnames[1]
    st <- locus@ranges@start[1]
    end <- st + sum(locus@ranges@width)

    bind_cols(tibble(scaf = as.character(scaf), start = st, end, snpMid = median(pos), nSites = length(pos)), f)
  })

  bind_rows(testD)

})

}


.fourPop <- function(GDS, locus, pops, x){

  p <- filter(pops, pops$Population %in% x)

  AF <- getAF(GDS, locus, pops = p, minSites = minSites)

  AF <- as.matrix(AF[complete.cases(AF),])

  pos <- as.numeric(AF[,1])




  f4Stats <- apply(AF, 1, function(z){

    #### get base AF
    p1 <- as.numeric(z[4])
    p2 <- as.numeric(z[5])
    p3 <- as.numeric(z[6])
    p4 <- as.numeric(z[7])

    ## get pd for fd
    f4 <- .f4fun(p1,p2,p3,p4)

    fd <- .f4fun(p1,pd,pd,p4)

    fhom <- .f4fun(p1,p3,p3,p4)


    tibble(f4,fhom,fd)

  })

  f4Stats <- bind_rows(f4Stats)

  fd <- abs(sum(f4Stats$f4)/sum(f4Stats$fd))

  f4 <- mean(f4Stats$f4)

  fhom <- abs(sum(f4Stats$f4)/sum(f4Stats$fhom))

  tibble(f4, fhom, fd)

}

.f4fun <- function(p1,p2,p3,p4){
  q1 <- 1 - p1
  q2 <- 1 - p2
  q3 <- 1 - p3
  q4 <- 1 - p4
  f4p <- (p1 - p2) * (p3 - p4)
  f4q <- (q1 - q2) * (q3 - q4)
  f4Corr <- -1*((f4p + f4q)/2)

}
