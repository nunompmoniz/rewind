#' Temporal Biased Random Oversampling
#'
#' @param form a model formula
#' @param data the original training set (with the unbalanced distribution)
#' @param rel is the relevance determined automatically (default: "auto") or provided by the user through a matrix. See examples.
#' @param thr.rel is the relevance threshold above which a case is considered as an extreme value
#' @param C.perc is a list containing the over-sampling percentage/s to apply to
#' all/each "class" obtained with the relevance threshold. The percentage represents the percentage of replicas that are added. Replicas of the examples are added randomly in each "class". Moreover, different percentages may be provided for each "class". Alternatively, it may be "balance" (the default) or "extreme", cases where the over-sampling percentages are automatically estimated.
#' @param repl is it allowed to perform sampling with replacement (bootstrapping). Defaults to TRUE because if the over-sampling percentage is >2 this is necessary.
#'
#' @return a new training data set resulting from the application of the resampling strategy
#' @export
#'
#' @examples
#' library(rewind)
#' data(temp)
#' ds <- create.data(temp,10)
#' C.perc <- list(4)
#' overT <- randOverT(V10 ~ ., ds, C.perc=C.perc)
#' overT.Bal <- randOverT(V10 ~ ., ds, C.perc="balance")
#' overT.Ext <- randOverT(V10 ~ ., ds, C.perc="extreme")
#'
randOverT <- function(form, dat, rel = "auto", thr.rel = 0.5,
  C.perc = "balance", repl = TRUE) {

  if(is.list(C.perc) & any(unlist(C.perc)<1)){
    stop("The over-sampling percentages provided in parameter C.perc
      can not be lower than 1!")
  }

  # the column where the target variable is
  tgt <- which(names(dat) == as.character(form[[2]]))

  y <- dat[, tgt]
  attr(y, "names") <- rownames(dat)

  s.y <- sort(y)
  if (is.matrix(rel)) {
    pc <- phi.control(y, method = "range", control.pts = rel)
  } else if (is.list(rel)) {
    pc <- rel
  } else if (rel == "auto") {
    pc <- phi.control(y, method = "extremes")
  } else {# TODO: handle other relevance functions and not using the threshold!
    stop("future work!")
  }

  temp <- y.relev <- phi(s.y, pc)
  if (!length(which(temp < 1))) {
    stop("All the points have relevance 1. Please, redefine your relevance
      function!")
  }
  if (!length(which(temp > 0))) {
    stop("All the points have relevance 0. Please, redefine your relevance
      function!")
  }

  bumps <- c()
  for (i in 1:(length(y) - 1)) {

    if ((temp[i] >= thr.rel && temp[i+1] < thr.rel) ||
        (temp[i] < thr.rel && temp[i+1] >= thr.rel)) {
      bumps <- c(bumps, i)
    }

  }
  nbump <- length(bumps) + 1 # number of different "classes"

  # collect the indexes in each "class"
  obs.ind <- as.list(rep(NA, nbump))
  last <- 1
  for (i in 1:length(bumps)) {
    obs.ind[[i]] <- s.y[last:bumps[i]]
    last <- bumps[i] + 1
  }
  obs.ind[[nbump]] <- s.y[last:length(s.y)]

  imp <- sapply(obs.ind, function(x) mean(phi(x, pc)))

  ove <- which(imp >= thr.rel)
  und <- which(imp < thr.rel)

  # set the over-sampling percentages
  if (is.list(C.perc)) {
    if (length(ove) > 1 & length(C.perc) == 1) {
      # only one percentage to apply to all the "classes"
      C.perc <- rep(C.perc[1], length(ove))
    } else if (length(ove) > length(C.perc) & length(C.perc) > 1) {
      stop("The number of over-sampling percentages must be equal to the
        number of bumps above the threshold defined!")
    } else if (length(ove) < length(C.perc)) {
      stop("The number of over-sampling percentages must be at most the
        number of bumps above the threshold defined!")
    }
    } else if (C.perc == "balance") {
      B <- sum(sapply(obs.ind[und], length))
      obj <- B/length(ove)
      C.perc <- as.list(round(obj/sapply(obs.ind[ove], length), 5))
  } else if (C.perc == "extreme") {
    Bund <- sum(sapply(obs.ind[und], length))/length(und)
    obj <- Bund^2/sapply(obs.ind[ove], length)
    C.perc <- as.list(round(obj/sapply(obs.ind[ove], length), 5))
  }

  newdata <- dat
  for (j in 1:length(ove)) {
    chdata <- dat[names(obs.ind[[ove[j]]]),] # select examples in the bump to undersample
    sdata <- chdata[order(names(obs.ind[[ove[j]]])),]
    r <- nrow(sdata)
    probs <- c()
    for(k in 1:r){
      probs <- c(probs,k/r)
    }
    sel <- sample(rownames(sdata),C.perc[[j]]*r, replace=repl, prob=probs)
    newdata <- rbind(newdata, dat[sel,])
  }

  newdata
  }
