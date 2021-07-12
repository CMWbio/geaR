### three and four pop admixture helper functions


## calculate f3 
.threePop <- function(AF, locus, pops, x){
  
  AF <- AF[complete.cases(AF),]
  AF <- AF[AF[[paste0(x[[3]], "_AF")]] %in% c(1,0),]
  AF <- AF[AF$alt != "N",]
  pos <- AF$pos
  
  AF <- AF[,paste0(x, "_AF")]
  
  AF <- abs(as.matrix(AF) - AF[[paste0(x[[3]], "_AF")]])
  nSites <- length(pos)
  snpMid <- median(pos)
  AF <- AF[apply(AF[,-1], 1, function(x) !all(x==0)),]
  
  f3Stats <- apply(AF, 1, function(z){
    
    #### get base AF
    p1 <- as.numeric(z[paste0(x[[1]], "_AF")])
    p2 <- as.numeric(z[paste0(x[[2]], "_AF")])
    p3 <- as.numeric(z[paste0(x[[3]], "_AF")])
    
    f3 <- .f3fun(p1, p2, p3)
    
    tibble::tibble(f3)
  })
  
  f3Stats <- dplyr::bind_rows(f3Stats)
  mean(f3Stats$f3)
  
}



.f3fun <- function(p1,p2,p3){
  q1 <- 1 - p1
  q2 <- 1 - p2
  q3 <- 1 - p3
  
  
  f3p <- (p3 - p1) * (p3 - p2)
  f3q <- (q3 - q1) * (q3 - q2)
  f3Corr <- ((f3p + f3q)/2)
  f3Corr
}


## calculate fd and f4 stat at a locus
.fourPop <- function(AF, locus, pops, x){

  AF <- AF[complete.cases(AF),]
  AF <- AF[AF[[paste0(x[[4]], "_AF")]] %in% c(1,0),]
  AF <- AF[AF$alt != "N",]
  pos <- AF$pos

  AF <- AF[,paste0(x, "_AF")]

  AF <- abs(as.matrix(AF) - AF[[paste0(x[[4]], "_AF")]])
  nSites <- length(pos)
  snpMid <- median(pos)
  AF <- AF[apply(AF[,-1], 1, function(x) !all(x==0)),]
  
  ### should I use data.table?
  # DT <- data.table(AF2)
  # DT[, ..I := .I]
  # system.time(DT[, f4Stats(z = .SD), by = ..I])

  f4Stats <- apply(AF, 1, function(z){

    #### get base AF
    p1 <- as.numeric(z[paste0(x[[1]], "_AF")])
    p2 <- as.numeric(z[paste0(x[[2]], "_AF")])
    p3 <- as.numeric(z[paste0(x[[3]], "_AF")])
    p4 <- as.numeric(z[paste0(x[[4]], "_AF")])

    # p1 <- tryoni_AF
    # p2 <- hybrid_AF
    # p3 <- dorsalis_AF
    # p4 <- oleae_AF

    pd <- (p2>= p3)*1 + (p2<p3)*1

    ## get pd for fd
    f4 <- .f4fun(p1,p2,p3,p4)

    fd <- .f4fun(p1,pd,pd,p4)

    fhom <- .f4fun(p1,p3,p3,p4)


    tibble(f4,fhom,fd)

  })

  f4Stats <- bind_rows(f4Stats)

  fd <- abs(sum(f4Stats$f4)/sum(f4Stats$fd))

  f4 <- mean(f4Stats$f4)

  df <- tibble(snpMid, nSites, f4, fd)

  colnames(df) <- c("snpMid", "nSites", paste0(c("f4(", "fd("), paste(x[1:3], collapse = ","), ";", x[4], ")"))

  
  df
}

### f4 base function
.f4fun <- function(p1,p2,p3,p4){
  q1 <- 1 - p1
  q2 <- 1 - p2
  q3 <- 1 - p3
  q4 <- 1 - p4
  
  
  f4p <- (p1 - p2) * (p3 - p4)
  f4q <- (q1 - q2) * (q3 - q4)
  f4Corr <- -1*((f4p + f4q)/2)

}



