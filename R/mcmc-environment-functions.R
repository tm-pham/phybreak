# The environments pbe0 and pbe1 will contain phybreak objects (excluding mcmc samples) + likelihood
# calculations.  pbe0 is consistent with the current state of updating phybreak.object; pbe1 is used to
# propose states and calculate likelihoods for these proposals.

# likarray is an array with likelihoods per nucleotide per SNP per node (sampling and coalescence). They are stored as
# vectors of length Nnodes*nSNPs*4, with implicit dim = c(4,nSNPs,Nnodes).

# logLikgen and logLiksam are the log-likelihood values of generation and sampling times with means mG and mS. logLikcoal is
# the log-likelihood value of the coalescent model.

# The function 'build_pbe' is used to initialise pbe0 at the start of mcmc-sampling. The function 'destroy_pbe' 
# is used at the end to return a phybreak object with current state. The function 'prepare_pbe' is used to 
# prepare pbe1 for another proposal. The proposal itself is then made in update-functions. The function 'propose_pbe' is
# used to calculate the likelihoods in pbe1 for the new proposed state.  The function
# 'accept_pbe' is to change pbe0 by copying pbe1.


# The environments are used during MCMC-updating, and in sim_phybreak and phybreak to simulate the phylogenetic tree.
pbe0 <- new.env()
pbe1 <- new.env()
userenv <- new.env()

# Copy functions to phybreak environments
copy2pbe0 <- function(var, env) {
  assign(var, get(var, env), pbe0)
}
copy2pbe0_2 <- function(var, env) {
  assign(var, get(var, env), pbe0_2)
}
copy2pbe1 <- function(var, env) {
  assign(var, get(var, env), pbe1)
}
copy2pbe1_2 <- function(var, env) {
  assign(var, get(var, env), pbe1_2)
}
copy2userenv <- function(var, env) {
  assign(var, get(var, env), userenv)
}

### build the pbe0 at the start of an mcmc chain by copying the fixed parameters and phybreak object, and by
### calculating likarray and the log-likelihoods 
build_pbe <- function(phybreak.obj) {
  ### Making everything available within the function
  le <- environment()
  d <- phybreak.obj$d
  h <- phybreak.obj$h
  v <- phybreak.obj$v
  p <- phybreak.obj$p
  lik_func <- phybreak.obj$likelihoods
  SNP <- t(matrix(unlist(d$sequences), ncol = d$nsamples))
  SNPfr <- attr(d$sequences, "weight")
  
  
  ### Jukes-Cantor: reduce diversity by naming each nucleotide by its frequency order, and grouping SNPs by same pattern across
  ### hosts
  # interpreting the phangorn sequence codes
  codematrix <- t(matrix(c(1,0,0,0,0,1,1,1,0,0,0,1,1,1,0,1,1,1,
                           0,1,0,0,0,1,0,0,1,1,0,1,1,0,1,1,1,1,
                           0,0,1,0,0,0,1,0,1,0,1,1,0,1,1,1,1,1,
                           0,0,0,1,1,0,0,1,0,1,1,0,1,1,1,1,1,1),
                         ncol = 4))
  # function to rename nucleotides by their frequency order (most frequent is 1, second most is 2, etc)
  fn <- function(snpvector) {
    bincodes <- codematrix[,snpvector]
    bincodes <- bincodes[order(rowSums(bincodes), decreasing = TRUE),]
    numcodes <- colSums(bincodes * c(1,2,4,8))
    snpcodes <- match(numcodes, colSums(codematrix * c(1,2,4,8)))
    return(snpcodes)
  }
  # rename the nucleotides by frequency order
  snpreduced <- apply(SNP, 2, fn)
  snpfrreduced <- SNPfr
  # remove identical SNP patterns
  if (ncol(SNP) > 1) {
    # create progress bar
    #pb <- txtProgressBar(min = 0, max = (ncol(SNP)-1), style = 3)
    for (i in (ncol(SNP) - 1):1) {
    #  setTxtProgressBar(pb, i)
      for (j in length(snpfrreduced):(i + 1)) {
        if (all(snpreduced[, i] == snpreduced[, j])) {
          snpfrreduced[i] <- snpfrreduced[i] + snpfrreduced[j]
          snpfrreduced <- snpfrreduced[-j]
          snpreduced <- snpreduced[, -j, drop = FALSE]
        }
      }
    }
  #  close(pb)
  }
  likarrayfreq <- snpfrreduced
  
  
  ### initialize all dimensions of likarray
  likarray <- array(1, dim = c(4, length(snpfrreduced), 2 * d$nsamples - 1))
  ### initialize likarray with observations on sampling nodes: 0 or 1
  likarray[cbind(1, rep(1:length(snpfrreduced), each = d$nsamples), rep(1:d$nsamples, 
                                                                        length(snpfrreduced)))] <- 1 * (snpreduced %in% c(1, 6, 7, 8, 12, 13, 14, 16))
  likarray[cbind(2, rep(1:length(snpfrreduced), each = d$nsamples), rep(1:d$nsamples, 
                                                                        length(snpfrreduced)))] <- 1 * (snpreduced %in% c(2, 6, 9, 10, 12, 13, 15, 16))
  likarray[cbind(3, rep(1:length(snpfrreduced), each = d$nsamples), rep(1:d$nsamples, 
                                                                        length(snpfrreduced)))] <- 1 * (snpreduced %in% c(3, 7, 9, 11, 12, 14, 15, 16))
  likarray[cbind(4, rep(1:length(snpfrreduced), each = d$nsamples), rep(1:d$nsamples, 
                                                                        length(snpfrreduced)))] <- 1 * (snpreduced %in% c(4, 5, 8, 10, 11, 13, 14, 15, 16))
  
  ### change the variables slot to an environmental variables slot (with transmission nodes in the tree)
  v <- phybreak2environment(v)
  
  ### complete likarray and calculate log-likelihood of sequences
  .likseqenv(le, (d$nsamples + 1):(2 * d$nsamples - 1), 1:d$nsamples)
  
  ### initialize all dimensions of contactarray
  if (inherits(d$contact, "matrix")){
    contactarray <- array(1, dim = c(ncol(d$contact), nrow(d$contact), 1))
  } else {
    contactarray <- array(1, dim = c(length(unique(d$hostnames)), 
                                   length(unique(d$hostnames)),
                                   length(d$contact)))
  }

  if (p$contact){
    ### initialize contact array with initial transmission tree
    for (i in seq_len(dim(contactarray)[3])){
      for(host_i in seq_len(ncol(contactarray[,,i]))){
        for(host_j in seq_len(ncol(contactarray[,,i]))){
          contactarray[host_i,host_j,i] <- phybreak:::get_contact_probability(d$contact, host_i, host_j, le)
        }
      }
    }
  }
  copy2pbe0("contactarray", le)


  ### calculate the other log-likelihoods
  logLiksam <- lik_sampletimes(p$obs, p$sample.shape, p$sample.mean, v$nodetimes, v$inftimes)
  logLikgen <- lik_gentimes(le)
  logLikcoal <- lik_coaltimes(le)

  if(length(lik_func) > 0){
    for(n in names(lik_func)){
      assign(n, lik_func[[n]](le))
      copy2pbe0(n, le)
    }
  }
  
  # logLikdist <- lik_distances(p$dist.model, p$dist.exponent, p$dist.scale, p$dist.mean, 
  #                             v$infectors, d$distances, d$area)
  # logLikcontact <- lik_contact(v$infectors, d$contact.matrix, p$cnt.invest.trans, p$cnt.invest.nontrans,
  #                              p$cnt.rep, p$cnt.rep.false)
  
  ### copy everything into pbe0
  copy2pbe0("d", le)
  copy2pbe0("h", le)
  copy2pbe0("v", le)
  copy2pbe0("p", le)
  copy2pbe0("lik_func", le)
  copy2pbe0("likarrayfreq", le)
  copy2pbe0("likarray", le)
  copy2pbe0("logLikseq", le)
  copy2pbe0("logLiksam", le)
  copy2pbe0("logLikgen", le)
  copy2pbe0("logLikcoal", le)
}


### take the elements d, v, p, and h from pbe0, and s from the function arguments, and make a new phybreak-object. Then
### empty the environments and return the new object.  
destroy_pbe <- function(phybreak.obj.samples) {
  #remove_history()
  res <- list(d = pbe0$d, v = environment2phybreak(pbe0$v), p = pbe0$p, h = pbe0$h, s = phybreak.obj.samples,
              hist = pbe0$v$inftimes[1])
  class(res) <- c("phybreak", "list")
  rm(list = ls(pbe0), envir = pbe0)
  rm(list = ls(pbe1), envir = pbe1)
  return(res)
}


### copy the elements from pbe0 into pbe1 to prepare for a proposal
prepare_pbe <- function() {
  copy2pbe1("d", pbe0)
  copy2pbe1("v", pbe0)
  copy2pbe1("p", pbe0)
  copy2pbe1("h", pbe0)
  pbe1$likarray <- pbe0$likarray + 0  #make a true copy, not a pointer
  copy2pbe1("likarrayfreq", pbe0)
  pbe1$contactarray <- pbe0$contactarray + 0 # make a true copy, not a pointer
  pbe1$logLikseq <- pbe0$logLikseq + 0 #make a true copy, not a pointer
  pbe1$logLiktoporatio <- 0
}


### calculate the new log-likelihoods where necessary and adjust likarray. Argument f indicates which type of function it is
### called from 
propose_pbe <- function(f) {
  ### Making variables and parameters available within the function
  le <- environment()
  d <- pbe0$d
  v <- pbe1$v
  p <- pbe1$p
  h <- pbe0$h
  contactarray <- pbe0$contactarray
  lik_func <- pbe0$lik_func
  hostID <- pbe1$hostID
  
  if (f == "phylotrans" || f == "withinhost") {
    # identify changed nodes
    chnodes <- c(which((v$nodeparents != pbe0$v$nodeparents[1:length(v$nodeparents)]) | 
                         (v$nodetimes != pbe0$v$nodetimes[1:length(v$nodetimes)])),
                 which(is.na(pbe0$v$nodeparents[1:length(v$nodeparents)])))
    chnodes <- unique(unlist(sapply(chnodes, .ptr, pars = v$nodeparents)))
    chnodes <- chnodes[chnodes > d$nsamples & chnodes < 2 * d$nsamples]
    # identify nodetips
    nodetips <- which(v$nodeparents %in% chnodes & v$nodetypes != "0")
    while(any(nodetips >= 2 * d$nsamples)) {
      nodetips[nodetips >= 2 * d$nsamples] <- which(v$nodeparents %in% nodetips[nodetips >= 2 * d$nsamples])
    }
    nodetips <- setdiff(nodetips, chnodes)
  } else if (f == "mu") {
    chnodes <- (d$nsamples + 1):(2 * d$nsamples - 1)
    nodetips <- 1:d$nsamples
  } else {
    chnodes <- NULL
  }
    
  if (!is.null(chnodes)) {
    .likseqenv(pbe1, chnodes, nodetips)

    if (p$contact == TRUE){
      chhosts <- which(v$infectors != pbe0$v$infectors)
      if (length(chhosts) > 0){
        if (length(chhosts) > 1){
          host.pairs <- combn(chhosts[chhosts <= ncol(d$contact)], m=2)
        } else if (length(chhosts) == 1){
          host.pairs <- t(t(c(v$infectors[chhosts], chhosts)))
        }
        if (!(0 %in% host.pairs)){
          for (i in seq_len(ncol(host.pairs))){
            for (array_i in seq_len(dim(contactarray)[3])){
              contactarray[host.pairs[1,i],host.pairs[2,i],array_i] <- 
              get_contact_probability(contactarray[,,array_i], host.pairs[1,i], host.pairs[2,i], le)
            }
          }
        }
      }
      copy2pbe1("contactarray", le)
    }
  }
  
  if (f == "phylotrans" || f == "trans" || f == "mG" || f == "ir" || f == "R") {
    logLikgen <- lik_gentimes(le)
    copy2pbe1("logLikgen", le)
  }
  
  if (f == "phylotrans" || f == "trans" || f == "mS") {
    logLiksam <- lik_sampletimes(p$obs, p$sample.shape, p$sample.mean, v$nodetimes, v$inftimes)
    copy2pbe1("logLiksam", le)
  }
  
  if (f == "trans" || (f == "mS" && p$wh.bottleneck == "wide") || f == "wh.slope" || f == "wh.exponent" || f == "wh.level" || f == "wh.history") {
    logLikcoal <- lik_coaltimes(le)
    copy2pbe1("logLikcoal", le)
  }
  
  if (f == "contact"){
    for (array_i in seq_len(dim(contactarray)[3])){
      for(host_i in seq_len(dim(contactarray)[1])){
        for(host_j in seq_len(dim(contactarray)[2])){
          if (dim(contactarray)[3] == 1) contactarray[host_i,host_j,array_i] <- get_contact_probability(d$contact, host_i, host_j, le)
          else contactarray[host_i,host_j,array_i] <- get_contact_probability(d$contact[[array_i]], host_i, host_j, le)
        }
      }
    }
    copy2pbe1("contactarray", le)
  }

  logLiks <- names(pbe0)[grepl("logLik", names(pbe0))]
  if (length(setdiff(logLiks, c("logLikseq", "logLikcoal", "logLiksam", "logLikgen"))) > 0){
    for(n in names(lik_func)){
      assign(n, lik_func[[n]](le))
      copy2pbe1(n, le)
    }
  }
  
  # if (f == "trans" || f == "phylotrans" || f == "dist.exponent" || f == "dist.scale" || f == "dist.mean") {
  #   logLikdist <- lik_distances(p$dist.model, p$dist.exponent, p$dist.scale, p$dist.mean,
  #                               v$infectors, d$distances, d$area)
  #   copy2pbe1("logLikdist", le)
  # }
  # 
  # if (f == ){}
  #   logLikcontact <- lik_contact(p)
  #   copy2pbe1("logLikcontact")
  # }
  
  # if (f == "phylotrans" || f == "hist.mean"){
  #   logLikintro <- lik_introductions(p$hist.mean, sum(v$infectors == 1), 
  #                                    max(d$sample.times) - min(v$inftimes[-1]))
  #   copy2pbe1("logLikintro", le)
  # }
  
  if (f == "mS" && p$wh.bottleneck == "complete") {
    copy2pbe1("logLikcoal", pbe0)
  }
}


### copy the elements from pbe1 into pbe0 upon acceptance of a proposal 
accept_pbe <- function(f) {
  if(f == "phylotrans" || f == "withinhost" || f == "trans") {
    copy2pbe0("v", pbe1)
  }
  
  if(f == "mG" || f == "mS" || f == "mu" || f == "ir" || 
     f == "wh.slope" || f == "wh.exponent" || f == "wh.level" || f == "wh.history" ||
     f == "dist.exponent" || f == "dist.scale" || f == "dist.mean") {
    copy2pbe0("p", pbe1)
  }
  
  if(f == "phylotrans" || f == "withinhost" || f == "mu") {
    copy2pbe0("likarray", pbe1)
    copy2pbe0("logLikseq", pbe1)
  }
  
  if(f == "phylotrans" || f == "trans" || f == "mS") {
    copy2pbe0("logLiksam", pbe1)
  }
  
  if(f == "phylotrans" || f == "trans" || f == "mG" || f == "ir") {
    copy2pbe0("logLikgen", pbe1)
  }
  
  # if (f == "trans" || f == "phylotrans" || f == "dist.exponent" || f == "dist.scale" || f == "dist.mean") {
  #   copy2pbe0("logLikdist", pbe1)
  # }
  
  if(f == "trans" || (f == "mS" && pbe0$p$wh.bottleneck == "wide") || f == "wh.slope" || f == "wh.exponent" || f == "wh.level" || f == "wh.history") {
    copy2pbe0("logLikcoal", pbe1)
  }
  
  if(f == "withinhost" || f == "phylotrans") {
    logLikcoal <- lik_coaltimes(pbe1)
    copy2pbe0("logLikcoal", environment())
  }
  
  logLiks <- names(pbe0)[grepl("logLik", names(pbe0))]
  if (length(setdiff(logLiks, c("logLikseq", "logLikcoal", "logLiksam", "logLikgen"))) > 0){
    logLiks <- setdiff(logLiks, c("logLikseq", "logLikcoal", "logLiksam", "logLikgen"))
    for(n in logLiks){
      copy2pbe0(n, pbe1)
      if (n == "logLikcontact"){
        copy2pbe0("contactarray", pbe1)
        copy2pbe0("p", pbe1)
      }
    }
  }
  
}

get_contact_probability <- function(contact, host_i, host_j, le){
  v = le$v
  p = le$p
  if (host_i != host_j){
    if (contact[host_i, host_j] == 0) return(p$cnt.eta * (1-p$cnt.epsilon))
    else if (contact[host_i, host_j] == 1) return(p$cnt.eta * (1-p$cnt.zeta))

    # #tau <- ifelse(dim(contactarray)[3] == 1, 1, floor(v$inftimes[host_j]))
    # if (v$infectors[host_j] == host_i || v$infectors[host_i] == host_j){
    #   if (contact[host_i,host_j] == 1){
    #     return((1-p$cnt.eta)*p$cnt.zeta + p$cnt.eta * p$cnt.epsilon)
    #   } else{
    #     return((1-p$cnt.eta)*(1-p$cnt.zeta) + p$cnt.eta*(1-p$cnt.epsilon))
    #   }
    # # } else if (v$infectors[host_i] == host_j){
    # #   tau <- ifelse(dim(contactarray)[3] == 1, 1, floor(v$inftimes[host_j]))
    # #   if (d$contact[[tau]][j,i] == 1){
    # #     lik <- c(lik, ((1-p$cnt.eta)*p$cnt.zeta + p$cnt.eta * p$cnt.epsilon) * contact.probs[categories[i], categories[j]])
    # #   } else{
    # #     lik <- c(lik, (1-p$cnt.eta)*(1-p$cnt.zeta) + p$cnt.eta*(1-p$cnt.epsilon) * (1-contact.probs[categories[i], categories[j]]))
    # #   }
    # } else {
    #   if (contact[host_i,host_j] == 1){
    #     return((1-p$cnt.lambda)*p$cnt.zeta + p$cnt.lambda * p$cnt.epsilon)
    #   } else {
    #     return((1-p$cnt.lambda)*(1-p$cnt.zeta) + p$cnt.lambda * (1-p$cnt.epsilon))
    #   }
    # }
  } else {
    return(1)
  }
}