require(stringr)
remove.excess.space.tabs <- function(iStr){
    oStr <- gsub("\t", "", oStrTmp <- iStr, fixed = T)
    while(oStr != oStrTmp){
        oStrTmp <- oStr
        oStr <- gsub("  ", " ", oStrTmp, fixed = T)
    }
    oStr
}
get.input.names <- function(inDir, iStrVar, iIntEquals, nr = NULL)
{
    if(is.null(nr))
        nr <- as.integer(system(paste("ls", inDir, "| grep -n \"NNI_\" | wc -l"), intern = T))
    pp <- readLines(paste0(inDir, "mod_inputs.txt"))
    intpp <- grep(paste(iStrVar,"= "), pp)
    regionsAll <- word(pp[intpp], iIntEquals + 1, sep = "=")
    ## Have we got all the variable information
    while(length(grep(";", regionsAll)) == 0) { ## There is more info to add
        intpp <- intpp + 1
        regionsAll <- paste(regionsAll, pp[intpp])
    }
    ## Get the region names in vectorised form
    regionsNew <- remove.excess.space.tabs(regionsAll)
    intRegs <- 1:nr
    regions <- word(regionsNew, intRegs, intRegs, sep = " ")
    while(regions[1] == ""){
        intRegs <- intRegs + 1
        regions <- word(regionsNew, intRegs, intRegs, sep = " ")
    }
    regions <- word(regions, 1, 1, sep = ",")
    word(regions, 1, 1, sep = ";")
}
get.variable.value <- function(inDir, iStrVar, iFile = "mod_inputs.txt", eq.delim = TRUE)
{
    pp <- readLines(paste0(inDir, iFile))
    intpp <- grep(iStrVar, pp)
    intpp <- intpp[-grep("#", pp[intpp])]
    oStrVar <- remove.excess.space.tabs(pp[intpp])
    ind <- iFile == "mod_pars.txt"
    iWhichWord <- 2 + (ind)
    delimSymb <- ifelse(eq.delim, "=", " ")
    outVal <- word(oStrVar, iWhichWord, iWhichWord, delimSymb)
    if(eq.delim){
        outTmp <- NULL
        i <- 1
        repeat{
            outTmp <- c(outTmp, word(outVal, i, i, ","))
            i <- i+1
            if(length(grep(";", outTmp)) > 0)
                break
        }
        outVal <- outTmp
    }
    as.numeric(str_extract(outVal, "\\-*\\d+\\.*\\d*"))
}

