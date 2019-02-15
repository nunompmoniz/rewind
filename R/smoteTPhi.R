#' Temporal and Relevance Biased SMOTE
#'
#' @param form a model formula
#' @param data the original training set (with the unbalanced distribution)
#' @param rel is the relevance determined automatically (default: "auto") or provided by the user through a matrix. See examples.
#' @param thr.rel is the relevance threshold above which a case is considered as an extreme value
#' @param C.perc is a list containing the over-sampling percentage/s to apply to
#' all/each "class" obtained with the relevance threshold. The percentage represents the percentage of replicas that are added. Replicas of the examples are added randomly in each "class". Moreover, different percentages may be provided for each "class". Alternatively, it may be "balance" (the default) or "extreme", cases where the over-sampling percentages are automatically estimated.
#' @param k is the number of neighbours to consider as the pool from where the new generated examples are generated
#' @param repl is it allowed to perform sampling with replacement (bootstrapping).
#' @param dist is the distance measure to be used (defaults to "Euclidean"). Use "HEOM" if there are nominal and numerical predictors
#' @param p is a parameter used when a p-norm is computed
#' @param delta for cases where data is extremely unbalanced and the majority of sampling probabilities are close zero, this is added to such probabilities
#'
#' @return a new training data set resulting from the application of the resampling strategy
#' @export
#'
#' @examples
#' library(rewind)
#' data(temp)
#' ds <- create.data(temp,10)
#' C.perc <- list(4,0.5,4)
#' smoteTPhi <- smoteTPhi(V10 ~ ., ds, C.perc=C.perc)
#' smoteTPhi.Bal <- smoteTPhi(V10 ~ ., ds, C.perc="balance")
#' smoteTPhi.Ext <- smoteTPhi(V10 ~ ., ds, C.perc="extreme")
#'
smoteTPhi <- function(form, data, rel="auto", thr.rel=0.5, C.perc="balance",
  k=5, repl=FALSE, dist="Euclidean", p=2, delta=0.01) {

  suppressWarnings(suppressPackageStartupMessages(library('uba')))

  if(any(is.na(data))){
    stop("The data set provided contains NA values!")
  }

  # the column where the target variable is
  tgt <- which(names(data) == as.character(form[[2]]))

  if (tgt < ncol(data)) {
    orig.order <- colnames(data)
    cols <- 1:ncol(data)
    cols[c(tgt,ncol(data))] <- cols[c(ncol(data),tgt)]
    data <-  data[,cols]
  }
  if(is.na(thr.rel)){
    stop("Future work!")
  }


  y <- resp(form,data)
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

  newdata <- data.frame()

  if(is.list(C.perc)){
    if(length(C.perc)!= nbump) stop(paste0("The percentages provided must be the same length as the number of bumps (",nbump,")!"))
  }else if(C.perc=="balance"){ # estimate the percentages of over/under sampling
    B <- round(nrow(data)/nbump,0)
    C.perc <- B/sapply(obs.ind, length)
  } else if(C.perc == "extreme"){
    B <- round(nrow(data)/nbump,0)
    rescale <- nbump*B/sum(B^2/sapply(obs.ind,length))
    obj <- round((B^2/sapply(obs.ind, length))*rescale,2)
    C.perc <- round(obj/sapply(obs.ind, length),1)
  }

  for(i in 1:nbump){

    if(length(names(obs.ind[[i]]))==1) {
      newdata <- rbind(newdata, data[names(obs.ind[[i]]),])
    } else {
      if(C.perc[[i]]==1){
        newdata <- rbind(newdata, data[names(obs.ind[[i]]),])
      }else if(C.perc[[i]]>1){

        if(length(names(obs.ind[[i]]))<=k) {
          #newdata <- rbind(newdata, data[names(obs.ind[[i]]),])
          newExs <- smote.exsTPhi(data[names(obs.ind[[i]]),],
            ncol(data),
            C.perc[[i]],
            (length(names(obs.ind[[i]]))-1),
            dist,
            p,
            pc)
          # add original rare examples and synthetic generated examples
          newdata <- rbind(newdata, newExs, data[names(obs.ind[[i]]),])
        } else {
          newExs <- smote.exsTPhi(data[names(obs.ind[[i]]),],
            ncol(data),
            C.perc[[i]],
            k,
            dist,
            p,
            pc)
          # add original rare examples and synthetic generated examples
          newdata <- rbind(newdata, newExs, data[names(obs.ind[[i]]),])
        }

      }else if(C.perc[[i]]<1){

        chdata <- data[names(obs.ind[[i]]),] # select examples in the bump to undersample
        sdata <- chdata[order(names(obs.ind[[i]])),]
        r <- nrow(sdata)
        #    timeRange <- length(unique(as.POSIXct(rownames(sdata))))
        probs <- c()
        for(pr in 1:r){
          probs <- c(probs,pr/r)
        }
        s.rel <- phi(sdata[,tgt],pc)
        probsU <- probs*s.rel

        sel = tryCatch({
          sample(rownames(sdata),C.perc[[i]]*r, replace=repl, prob=probsU)
        }, error = function(e) {
          sample(rownames(sdata),C.perc[[i]]*r, replace=repl, prob=probsU+delta)
        })

        newdata <- rbind(newdata, data[sel,])

      }
    }


  }

  if (tgt < ncol(data)) {
    newdata <- newdata[,cols]
    data <- data[,cols]
  }

  newdata
}

#' Generate synthetic cases for smoteTPhi
#'
#' The result of the function is a (N-1)*nrow(data) set of synthetically generated examples with rare values on the target
#'
#' @param data the rare cases (the minority "class" cases)
#' @param tgt the column nr of the target variable
#' @param N the percentage of over-sampling to carry out;
#' @param k the number of nearest neighours
#' @param dist the distance function used for the neighours computation
#' @param p an integer used when a "p-norm" distance is selected
#' @param pc relevance function
#'
#' @return
#'
smote.exsTPhi <- function(data, tgt, N, k, dist, p, pc) {

  nomatr <- c()
  T <- matrix(nrow=dim(data)[1],ncol=dim(data)[2])
  for(col in seq.int(dim(T)[2]))
    if (class(data[,col]) %in% c('factor','character')) {
      T[,col] <- as.integer(data[,col])
      nomatr <- c(nomatr,col)
    } else T[,col] <- data[,col]

  nC <- dim(T)[2]
  nT <- dim(T)[1]


  ranges <- rep(1,nC)
  if(length(nomatr)){
    for(x in (1:nC)[-c(nomatr)]) ranges[x] <- max(T[,x]) - min(T[,x])
  } else{
    for(x in (1:nC)) ranges[x] <- max(T[,x]) - min(T[,x])
  }

  # order the data by time so that the most recent examples have the higher line numbers
  data <- data[order(rownames(data)),]

  kNNs <-neighbours(tgt, data, dist, p, k)

  nexs <-  as.integer(N-1) # nr of examples to generate for each rare case
  extra <- as.integer(nT*(N-1-nexs)) # the extra examples to generate
  idx <- sample(1:nT, extra)
  new <- matrix(nrow=nexs*nT+extra,ncol=nC)    # the new cases

  if(nexs){
    for(i in 1:nT) {
      for(n in 1:nexs) {
        # select the nearest neighbour more recent and with higher phi
        y.rel <- phi(data[kNNs[i,],tgt], pc)
        pos.eval <- (kNNs[i,]/max(kNNs[i,]))*y.rel
        neig <- match(max(pos.eval), pos.eval)
        # the attribute values of the generated case
        difs <- T[kNNs[i,neig],-tgt]-T[i,-tgt]
        new[(i-1)*nexs+n,-tgt] <- T[i,-tgt]+runif(1)*difs
        for(a in nomatr) # nominal attributes are randomly selected among the existing values of seed and the selected neighbour
          new[(i-1)*nexs+n,a] <- c(T[kNNs[i,neig],a],T[i,a])[1+round(runif(1),0)]

        # now the target value (weighted (by inverse distance) average)
        d1 <- d2 <- 0
        for(x in (1:nC)[-c(nomatr, tgt)]) {
          d1 <- abs(T[i,x] - new[(i-1)*nexs+n,x])/ranges[x]
          d2 <- abs(T[kNNs[i,neig],x] - new[(i-1)*nexs+n,x])/ranges[x]
        }
        if (length(nomatr)) {
          d1 <- d1 + sum(T[i,nomatr] != new[(i-1)*nexs+n,nomatr])
          d2 <- d2 + sum(T[kNNs[i,neig],nomatr] != new[(i-1)*nexs+n,nomatr])
        }
        # (d2+d1-d1 = d2 and d2+d1-d2 = d1) the more distant the less weight
        new[(i-1)*nexs+n,tgt] <- if (d1 == d2) (T[i,tgt]+T[kNNs[i,neig],tgt])/2 else (d2*T[i,tgt]+d1*T[kNNs[i,neig],tgt])/(d1+d2)

      }
    }
  }
  if(extra){
    count<-1
    for (i in idx){

      # select the nearest neighbour more recent and with higher phi
      y.rel <- phi(data[kNNs[i,],tgt], pc)
      pos.eval <- (kNNs[i,]/max(kNNs[i,]))*y.rel
      neig <- match(max(pos.eval), pos.eval)

      # the attribute values of the generated case
      difs <- T[kNNs[i,neig],-tgt]-T[i,-tgt]
      new[nexs*nT+count,-tgt] <- T[i,-tgt]+runif(1)*difs
      for(a in nomatr)
        new[nexs*nT+count,a] <- c(T[kNNs[i,neig],a],T[i,a])[1+round(runif(1),0)]


      # now the target value (weighted (by inverse distance) average)
      d1 <- d2 <- 0
      for(x in (1:nC)[-c(nomatr,tgt)]) {
        d1 <- abs(T[i,x] - new[nexs*nT+count,x])/ranges[x]
        d2 <- abs(T[kNNs[i,neig],x] - new[nexs*nT+count,x])/ranges[x]
      }
      if (length(nomatr)) {
        d1 <- d1 + sum(T[i,nomatr] != new[(i-1)*nexs+n,nomatr])
        d2 <- d2 + sum(T[kNNs[i,neig],nomatr] != new[(i-1)*nexs+n,nomatr])
      }
      # (d2+d1-d1 = d2 and d2+d1-d2 = d1) the more distant the less weight
      new[nexs*nT+count,tgt] <- if (d1 == d2) (T[i,tgt]+T[kNNs[i,neig],tgt])/2 else (d2*T[i,tgt]+d1*T[kNNs[i,neig],tgt])/(d1+d2)

      count <- count+1
    }
  }

  newCases <- data.frame(new)

  for(a in nomatr)
    newCases[,a] <- factor(newCases[,a],levels=1:nlevels(data[,a]),labels=levels(data[,a]))

  colnames(newCases) <- colnames(data)
  newCases

}
