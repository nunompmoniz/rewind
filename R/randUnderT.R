#' Temporal-Biased Random Undersampling
#'
#' @param form a model formula
#' @param data the original training set (with the unbalanced distribution)
#' @param rel is the relevance determined automatically (default: "auto") or provided by the user through a matrix. See examples.
#' @param thr.rel is the relevance threshold above which a case is considered as an extreme value
#' @param C.perc is a list containing the under-sampling percentage/s to apply to all/each "class" obtained with the relevance threshold. This percentage represents the percentage of examples that is maintained in each "class". Examples are randomly removed in each "class". Moreover, different percentages may be provided for each "class". Alternatively, it may be "balance" or "extreme", cases where the under-sampling percentages are automatically estimated.
#' @param repl is it allowed to perform sampling with replacement (bootstrapping)
#'
#' @return
#' @export
#'
#' @examples
#' library(rewind)
#' data(temp)
#' ds <- create.data(temp,10)
#' C.perc <- list(0.5)
#' underT <- randUnderT(V10 ~ ., ds, C.perc=C.perc)
#' underT.Bal <- randUnderT(V10 ~ ., ds, C.perc="balance")
#' underT.Ext <- randUnderT(V10 ~ ., ds, C.perc="extreme")
#'
randUnderT <- function(form, data, rel="auto", thr.rel=0.5, C.perc="balance", repl=FALSE) {

  suppressWarnings(suppressPackageStartupMessages(library('uba')))

  if(any(is.na(data))){
    stop("The data set provided contains NA values!")
  }

  # the column where the target variable is
  tgt <- which(names(data) == as.character(form[[2]]))

  y <- data[,tgt]
  attr(y,"names") <- rownames(data)
  s.y <- sort(y)

  if (is.matrix(rel)){
    pc <- uba::phi.control(y, method="range", control.pts=rel)
  }else if(is.list(rel)){
    pc <- rel
  }else if(rel=="auto"){
    pc <- uba::phi.control(y, method="extremes")
  }else{# TODO: handle other relevance functions and not using the threshold!
    stop("future work!")
  }

  temp <- y.relev <- phi(s.y,pc)
  if(!length(which(temp<1)))stop("All the points have relevance 1. Please, redefine your relevance function!")
  if(!length(which(temp>0)))stop("All the points have relevance 0. Please, redefine your relevance function!")

  temp[which(y.relev>thr.rel)] <- -temp[which(y.relev>thr.rel)]
  bumps <- c()
  for(i in 1:(length(y)-1)){if(temp[i]*temp[i+1]<0) bumps <- c(bumps,i)}
  nbump <- length(bumps)+1 # number of different classes

  # collect the indexes in each "class"
  nbumps <- length(bumps) +1
  obs.ind <- as.list(rep(NA, nbump))
  last <- 1
  for(i in 1:length(bumps)) {
    obs.ind[[i]] <- s.y[last:bumps[i]]
    last <- bumps[i] +1
  }
  obs.ind[[nbumps]] <- s.y[last:length(s.y)]

  imp <- sapply(obs.ind, function(x)mean(phi(x,pc)))

  und <- which(imp<thr.rel)
  ove <- which(imp>thr.rel)

  newdata <- NULL
  for(j in 1:length(ove)){
    newdata <- rbind(newdata,data[names(obs.ind[[ove[j]]]),]) # start with the examples from the minority "classes"
  }

  # set the undersampling percentages
  if(is.list(C.perc)){
    if(length(und) > 1 & length(C.perc)==1){ # the same under-sampling percentage is applied to all the "classes"
      C.perc <- rep(C.perc[1],length(und))
    } else if(length(und)>length(C.perc) & length(C.perc)>1){
      stop("The number of under-sampling percentages must be equal to the number of bumps below the threshold defined!")
    }else if(length(und)< length(C.perc)){
      stop("the number of under-sampling percentages must be at most the number of bumps below the threshold defined!")
    }# each class has its predefined over-sampling percentage
  }else if(C.perc == "balance"){
    B <- sum(sapply(obs.ind[ove],length))
    obj <- B/length(und)
    C.perc <- as.list(round(obj/sapply(obs.ind[und],length),5))
  }else if(C.perc== "extreme"){
    Bove <- sum(sapply(obs.ind[ove],length))/length(ove)
    obj <- Bove^2/sapply(obs.ind[und],length)
    C.perc <- as.list(round(obj/sapply(obs.ind[und], length),5))
  }

  for(j in 1:length(und)){
    chdata <- data[names(obs.ind[[und[j]]]),] # select examples in the bump to undersample
    sdata <- chdata[order(names(obs.ind[[und[j]]])),]
    r <- nrow(sdata)
    probs <- c()
    for(k in 1:r){
      probs <- c(probs,k/r)
    }

    sel <- c()
    if(length(rownames(sdata))==1) {
      sel <- rownames(sdata)
    } else {
      if((C.perc[[j]]*r)>r) {
        sel <- rownames(sdata)
      } else {
        sel <- sample(rownames(sdata),C.perc[[j]]*r, replace=repl, prob=probs)
      }

    }

    newdata <- rbind(newdata, data[sel,])
  }

  newdata <- newdata[order(rownames(newdata)),]
  newdata
}
