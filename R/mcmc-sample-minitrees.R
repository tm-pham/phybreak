### functions to simulate mini-trees ###


sample_coaltimes <- function(tiptimes, inftime, parameters, historyhost = FALSE) {
  ### tests
  if(!historyhost) if(min(tiptimes) < inftime) stop("sample_coaltimes with negative tip times")

  ### function body
  if(length(tiptimes) < 2) return(c())
  
  if(historyhost) {
    # transform times so that fixed rate 1 can be used
    #ttrans <- sort(tiptimes/parameters$wh.history, decreasing = TRUE, na.last = TRUE)
    ttrans <- tiptimes/parameters$wh.history
    ttrans <- ttrans[order(ttrans, decreasing = TRUE)]
    tnodetrans <- .sctwh3(ttrans)
    
    #res <- sort(parameters$wh.history * tnodetrans, na.last = TRUE)
    res <- (parameters$wh.history * tnodetrans)[order(parameters$wh.history * tnodetrans)]
    
    return(res)
  }
  
  switch(
    parameters$wh.model,
    #coalescence at transmission
    single = head(sort(tiptimes, na.last = TRUE), -1),
    #coalescence at infection
    infinite = inftime + 0*tiptimes[-1],
    #linear increase
    linear = {
      # transform times so that fixed rate 1 can be used
      #ttrans <- sort(log(parameters$wh.level/parameters$wh.slope + tiptimes - inftime)/(parameters$wh.slope), decreasing = TRUE, na.last = TRUE)
      ttrans <- log(parameters$wh.level/parameters$wh.slope + tiptimes - inftime)/(parameters$wh.slope)
      ttrans <- ttrans[order(ttrans, decreasing = TRUE)]
      tnodetrans <- .sctwh3(ttrans)
      
      #res <- sort(exp(parameters$wh.slope * tnodetrans), na.last = TRUE)
      res <- (parameters$wh.slope * tnodetrans)[order(parameters$wh.slope * tnodetrans)]
      res <- apply(cbind(res,
                         min(10^-5, (parameters$wh.level/parameters$wh.slope + tiptimes - inftime)/length(tiptimes))),
                   1, max)
      res <- res + inftime - parameters$wh.level/parameters$wh.slope
      
      return(res)
    },
    #exponential increase
    exponential = {
      # transform times so that fixed rate 1 can be used
      ttrans <- sort(-exp(-parameters$wh.exponent * (tiptimes - inftime))/(parameters$wh.exponent * parameters$wh.level), 
                     decreasing = TRUE, na.last = TRUE)
      tnodetrans <- .sctwh3(ttrans)
      
      res <- inftime + sort(-log(-parameters$wh.exponent * parameters$wh.level * tnodetrans)/(parameters$wh.exponent), na.last = TRUE)

      return(res)
    },
    #constant level
    constant = {
      # transform times so that fixed rate 1 can be used
      ttrans <- sort((tiptimes - inftime)/parameters$wh.level, decreasing = TRUE, na.last = TRUE)
      tnodetrans <- .sctwh3(ttrans)
      
      res <- inftime + sort(parameters$wh.level * tnodetrans, na.last = TRUE)
      
      return(res)
    }
  )
}

sample_topology <- function(nodeIDs, nodetimes, nodetypes, infectornodes) {
  res <- rep(-1, length(nodeIDs))
  tochoose <- infectornodes
  for(i in 1:length(nodeIDs)) {
    res[i] <- tochoose[1]
    if(nodetypes[i] == "c") {
      tochoose <- sample(c(tochoose[-1], nodeIDs[i], nodeIDs[i]))
    } else {
      tochoose <- tochoose[-1]
    }
  }
  return(res)
}

sample_singlecoaltime <- function(oldtiptimes, oldcoaltimes, newtiptime, inftime, parameters, historyhost = FALSE) {
  ### return -Inf if there is no existing tree
  if(length(oldtiptimes) + length((oldcoaltimes)) == 0) return(-Inf)
  
  ### add 0s to prevent problems with empty vectors
  oldtiptimes <- c(inftime, oldtiptimes)
  oldcoaltimes <- c(inftime, oldcoaltimes)
  
  
  ### function body
  if(historyhost) {
    # transform times so that fixed rate 1 can be used
    transtiptimes <- (oldtiptimes)/parameters$wh.history
    transcoaltimes <- (oldcoaltimes)/parameters$wh.history
    transcurtime <- (newtiptime)/parameters$wh.history
  } else switch(
    parameters$wh.model,
    #coalescence at transmission
    single = return(min(newtiptime, max(oldtiptimes))),
    #coalescence at infection
    infinite = return(inftime),
    linear = {
      # transform times so that fixed rate 1 can be used
      transtiptimes <- log(parameters$wh.level/parameters$wh.slope + oldtiptimes - inftime)/parameters$wh.slope
      transcoaltimes <- log(parameters$wh.level/parameters$wh.slope + oldcoaltimes - inftime)/parameters$wh.slope
      transcurtime <- log(parameters$wh.level/parameters$wh.slope + newtiptime - inftime)/parameters$wh.slope
    },
    exponential = {
      transtiptimes <- -exp(-parameters$wh.exponent * (oldtiptimes - inftime))/(parameters$wh.exponent * parameters$wh.level)
      transcoaltimes <- -exp(-parameters$wh.exponent * (oldcoaltimes - inftime))/(parameters$wh.exponent * parameters$wh.level)
      transcurtime <- -exp(-parameters$wh.exponent * (newtiptime - inftime))/(parameters$wh.exponent * parameters$wh.level)
    },
    constant = {
      transtiptimes <- (oldtiptimes - inftime)/parameters$wh.level
transcoaltimes <- (oldcoaltimes - inftime)/parameters$wh.level
      transcurtime <- (newtiptime - inftime)/parameters$wh.level
    }
  )

  # start with number of edges at current time and time of next node (backwards next)
  curnedge <- sum(transtiptimes >= transcurtime) - sum(transcoaltimes >= transcurtime)
  transnexttime <- max(c(-Inf, transtiptimes[transtiptimes < transcurtime],  
                         transcoaltimes[transcoaltimes < transcurtime]))
  
  # traverse minitree node by node, subtracting coalescence rate until cumulative rate reaches 0
  #frailty <- stats::rexp(1)
  frailty <- rexp(1)
  while(curnedge * (transcurtime - transnexttime) < frailty) {
    transcurtime <- transnexttime
    curnedge <- sum(transtiptimes >= transcurtime) - sum(transcoaltimes >= transcurtime)
    transnexttime <- max(c(-Inf, transtiptimes[transtiptimes < transcurtime],  
                           transcoaltimes[transcoaltimes < transcurtime]))
    #frailty <- stats::rexp(1)
    frailty <- rexp(1)
  }
  
  # calculate transformed node time
  transreturntime <- transcurtime - frailty / curnedge
  
  # transform to real time
  if(historyhost) {
    parameters$wh.history * transreturntime
  } else switch(
    parameters$wh.model, single =, infinite = ,
    linear = {
      res <- min(c(10^-5, (oldcoaltimes[-1] - inftime)/2))
      res <- max(res, exp(parameters$wh.slope * transreturntime))
      res <- inftime + res - parameters$wh.level/parameters$wh.slope
      return(res)
      },
    exponential = inftime - log(-parameters$wh.exponent * parameters$wh.level * transreturntime)/(parameters$wh.exponent),
    constant = inftime + parameters$wh.level * transreturntime
  )
}

sample_coaltimes_history <- function(init = FALSE) {
  ### First, sample coalescentimes
  # edges entering hostID, with endtimes
  edgesin <- which(pbe1$v$nodehosts == 0 & pbe1$v$nodetypes != "c")
  edgeintimes <- pbe1$v$nodetimes[edgesin]
  
  # all coalescent nodes in new infector and hostID
  coalnodes <- which(pbe1$v$nodehosts == 0 & pbe1$v$nodetypes == "c")
  newcoaltimes <- rev(sample_coaltimes(edgeintimes, NA, pbe1$p, TRUE))
  
  if (init) newcoaltimes <- newcoaltimes - abs(min(pbe1$v$inftimes) - max(newcoaltimes))
  
  ### starting in the history tips, go backwards and give new time to 
  # sampled coalescent node
  tips <- which(pbe1$v$nodetypes == "t" & pbe1$v$nodehosts == 0)
  npaths <- lapply(tips, function(x) return(.ptr(pbe1$v$nodeparents, x)[-1]))
  
  break.proposal <- FALSE
  
  for (i in seq_along(newcoaltimes)){
    freq <- table(sapply(npaths, function(x) x[[1]]))
    cnode <- names(freq)[freq==2]
    free_cnodes <- sapply(cnode, function(x) {
      if (all(pbe1$v$nodetimes[pbe1$v$nodeparents == x] > newcoaltimes[i]))
        return(x)
      else
        return(NA)
    })
    
    if (!all(is.na(free_cnodes))) {
      cnode <- sample(free_cnodes[!is.na(free_cnodes)], size = 1)
    } else {
      break.proposal <- TRUE
      break
    }
    cnode <- as.numeric(cnode)
    pbe1$v$nodetimes[cnode] <- newcoaltimes[i]
    
    remove <- c()
    for (j in seq_along(npaths)){
      if (cnode %in% npaths[[j]] & is.null(remove)) {
        remove <- j
      } else if (cnode %in% npaths[[j]]){
        npaths[[j]] <- tail(npaths[[j]], -1)
      }
    }
    npaths[[remove]] <- NULL
  }
  
  ### if not all coalescent times assigned correctly, stop proposal
  if(break.proposal) return(TRUE)
  else return(FALSE)
  
}

sample_singlechildnode <- function(nodeIDs, nodeparents, nodetimes, newnodetime) {
  candidatenodes <- nodeIDs[nodetimes >= newnodetime & 
                           c(-Inf, nodetimes)[1 + match(nodeparents, nodeIDs, nomatch = 0)] < newnodetime]
  if(length(candidatenodes) == 0) candidatenodes <- nodeIDs[nodetimes == max(nodetimes)]
  return(sample(rep(candidatenodes, 2), 1))
}


### Should I ever consider weighted topology sampling, the topology should be sampled backwards in time instead of forwards.
### This function does exactly that and is slightly less efficient than the now-used forward sampling
# sample_topology_backwards <- function(nodeIDs, nodetypes, infectornodes) {
#   res <- rep(-1, length(nodeIDs))
#   tochoose <- c()
#   for(i in length(nodeIDs):1) {
#     if(nodetypes[i] == "c") {
#       tojoin <- sample(length(tochoose), 2)
#       res[tochoose[tojoin]] <- nodeIDs[i]
#       tochoose <- c(tochoose[-tojoin], i)
#     } else {
#       tochoose <- c(tochoose, i)
#     }
#   }
#   res[res == -1] <- infectornodes
#   return(res)
# }
