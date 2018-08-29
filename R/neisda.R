neisDa <- function(dxy, pi){

  da <- lapply(seq(from = 1, to = length(dxy), by = 2), function(x){
    #get dxy pairwise comparison name from colnames of dxy
    dxyName <- colnames(dxy)[x]
    #remove the dxy from the end
    compName <- gsub("_dxy", "", dxyName)
    #get sample names
    samples <- unlist(strsplit(compName, "v"))

    #get pi for population x
    xPi <- pi[[paste0(samples[1], "_pi")]]
    #pi for population y
    yPi <- pi[[paste0(samples[2], "_pi")]]

    #carry out da calculation from Nei 1987
    da <- dxy[[dxyName]] - ((xPi + yPi)/2)
    #make tibble
    da <- data_frame(da)
    #colnames
    colnames(da) <- paste0(compName, "_da")
    da

  })

  bind_cols(da)


}
