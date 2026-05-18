#' Log-likelihood of a phybreak-object.
#' 
#' The likelihood of a \code{phybreak}-object is calculated, with the option to include or exclude parts of the 
#' likelihood for genetic data, phylogenetic tree (within-host model), sampling times and generation times.
#' 
#' The sequence likelihood is calculated by Felsenstein's pruning algorithm, assuming a prior probability of 0.25 
#' for each nucleotide. The within-host likelihood is the likelihood of coalescence times given the within-host model 
#' and slope. The generation interval and sampling interval likelihood are log-densities of the gamma distributions 
#' for these variables.
#' 
#' @param object An object of class \code{phybreak}.
#' @param genetic Whether to include the likelihood of the mutation model.
#' @param withinhost Whether to include the likelihood of within-host (coalescent) model.
#' @param sampling Whether to include the likelihood of the sampling model (sampling intervals).
#' @param generation Whether to include the likelihood of the transmission model (generation intervals).
#' @param ... Some methods for this generic require additional arguments. None are used in this method.
#' @return The log-likelihood as an object of class logLik.
#' @author Don Klinkenberg \email{don@@xs4all.nl}
#' @references \href{http://dx.doi.org/10.1371/journal.pcbi.1005495}{Klinkenberg et al. (2017)} Simultaneous 
#'   inference of phylogenetic and transmission trees in infectious disease outbreaks. 
#'   \emph{PLoS Comput Biol}, \strong{13}(5): e1005495.
#' @examples 
#' #First build a phybreak-object containing samples.
#' simulation <- sim_phybreak(obsize = 5)
#' MCMCstate <- phybreak(dataset = simulation)
#' logLik(MCMCstate)
#' 
#' MCMCstate <- burnin_phybreak(MCMCstate, ncycles = 20)
#' logLik(MCMCstate)
#' 
#' tree0 <- get_phylo(MCMCstate)
#' seqdata <- get_data(MCMCstate)$sequences
#' phangorn::pml(tree0, seqdata, rate = 0.75*get_parameters(MCMCstate, "mu") 
#' logLik(MCMCstate, genetic = TRUE, withinhost = FALSE, 
#'        sampling = FALSE, generation = FALSE) #should give the same result as 'pml'
#' @export
logLik.phybreak <- function(object, genetic = TRUE, withinhost = TRUE, sampling = TRUE,
                            generation = TRUE, 
                            distance = TRUE, ...) {
  
  res <- 0
  if (genetic) {
    res <- res + with(object, .likseq(matrix(unlist(d$sequences), ncol = d$nsamples), 
                                      attr(d$sequences, "weight"), 
                                      v$nodeparents, v$nodetimes, p$mu, d$nsamples))
  }
  if (generation) {
    res <- res + with(object, lik_gentimes(list(p = p, v = v)))
  }
  if (sampling) {
    res <- res + with(object, lik_sampletimes(p$obs, p$sample.shape, p$sample.mean, v$nodetimes, v$inftimes, d$last.negative))
  }
  if (withinhost) {
    objectenv <- object
    objectenv$v <- phybreak2environment(objectenv$v)
    res <- res + with(object, lik_coaltimes(objectenv))
  }
  if(distance && !is.null(object$d$distances)) {
    res <- res + with(object, lik_distances(p$dist.model, p$dist.exponent, p$dist.scale, p$dist.mean, 
                                            v$infectors, d$distances))
  }
  attributes(res) <- list(
    nobs = object$p$obs,
    df = 1 + object$h$est.mG + object$h$est.mS + object$h$est.wh.s + object$h$est.wh.e + object$h$est.wh.0 +
      object$h$est.wh.h + object$h$est.dist.e + object$h$est.dist.s + object$h$est.dist.m,
    genetic = genetic, withinhost = withinhost, sampling = sampling, generation = generation, distance = distance
  )
  class(res) <- "logLik"
  return(res)
}


### calculate the log-likelihood of generation intervals 
# lik_gentimes <- function(shapeG, meanG, sampleScale, cullingScale, transModel, inftimes, infectors,
#                          nodetimes, cultimes) {
lik_gentimes <- function(le){
  p <- le$p
  v <- le$v
  indices <- v$infectors == 0
  othercases <- v$infectors > 0
  
  intro.rate <- ifelse(is.null(p$intro.rate), 1, p$intro.rate)
  R <- ifelse(is.null(p$R), 1, p$R)
  
  L <- 0 + # force of infection from external source
    (- intro.rate * (max(v$nodetimes) - min(v$inftimes))) +  
    log(intro.rate) * sum(indices) +
    - R * length(v$infectors)
  
  if(p$infectivity && sum(othercases) == 0)
    return(L)
  else
    return( L +
            sum(log(R) + infect_distribution(time = v$inftimes[othercases],
                                             inftimes = v$inftimes[v$infectors[othercases]],
                                             nodetimes = v$nodetimes[v$nodetypes=="s"][v$infectors[othercases]],
                                             le = le, log = TRUE)))
}

### calculate the log-likelihood of sampling intervals
###
### When last.negative[i] is set, host i's sampling interval D_i = nodetime_i - tinf_i
### is constrained to D_i < M_i where M_i = nodetime_i - last.negative_i (perfect-test
### assumption: tinf_i > last.negative_i). The likelihood for such hosts is the
### Gamma density truncated to (0, M_i):
###
###   log f(D_i) - log F(M_i, shape, scale)   if D_i < M_i
###   -Inf                                    otherwise
###
### The -log F(M_i) normalizer is required for unbiased inference of shapeS/meanS:
### without it, the conditioning event D_i < M_i is absorbed into the data without
### the corresponding likelihood penalty for parameters that put little mass in (0, M_i).
lik_sampletimes <- function(obs, shapeS, meanS, nodetimes, inftimes, last.negative = NULL) {
  scaleS <- meanS / shapeS
  D <- nodetimes[1:obs] - inftimes[1:obs]
  log_dens <- dgamma(D, shape = shapeS, scale = scaleS, log = TRUE)

  if (length(last.negative) > 0) {
    ln <- last.negative[1:obs]
    has_lastneg <- !is.na(ln)
    if (any(has_lastneg)) {
      M <- nodetimes[1:obs][has_lastneg] - ln[has_lastneg]
      if (any(D[has_lastneg] >= M)) return(-Inf)
      log_trunc_norm <- pgamma(M, shape = shapeS, scale = scaleS, log.p = TRUE)
      return(sum(log_dens) - sum(log_trunc_norm))
    }
  }

  sum(log_dens)
}


### calculate the log-likelihood of distances 
# lik_distances <- function(dist.model, dist.exponent, dist.scale, dist.mean, infectors, distances, area) {
#   dist.model <- ifelse(is.null(dist.model), "none", dist.model)
#   if(dist.model == "none") return(0)
#   distancevector <- distances[cbind(which(infectors!=0), infectors[infectors!=0])]
#   sum((infectors == 0) * dunif(1, min = 0, max = area, log = TRUE)) + 
#   switch(dist.model,
#          power = sum(log(
#            dist.exponent * sin(pi/dist.exponent) / 
#              (dist.scale * pi * (1 + (distancevector/dist.scale)^dist.exponent))
#            )),
#          exponential = sum(
#            log(dist.exponent) - dist.exponent * distancevector
#          ),
#          poisson = sum(
#            -dist.mean + distancevector * log(dist.mean) - lgamma(1 + distancevector)
#          )
#   )
# }

### calculate the log-likelihood of contacts
# lik_contact <- function(infectors, cnt.matrix, cnt.invest.trans, cnt.invest.nontrans,
#                         cnt.rep, cnt.rep.false) {
#   if(is.null(cnt.matrix)) return(0)
#   
#   lik <- c()
#   for(i in 1:ncol(cnt.matrix)){
#     for(j in 1:ncol(cnt.matrix)){
#       if (i != j){
#         if (infectors[j] == i){
#           if (cnt.matrix[i,j] == 1){
#             lik <- c(lik, (1-cnt.rep.false)*cnt.invest.trans + cnt.invest.trans*cnt.rep)
#           } else{
#             lik <- c(lik, (1-cnt.rep.false)*(1-cnt.invest.trans) + cnt.invest.trans*(1-cnt.rep))
#           }
#         } else {
#           if (cnt.matrix[i,j] == 1){
#             lik <- c(lik, (1-cnt.invest.nontrans)*cnt.rep.false + cnt.invest.nontrans*cnt.rep)
#           } else {
#             lik <- c(lik, (1-cnt.invest.nontrans)*(1-cnt.rep.false) + cnt.invest.nontrans*(1-cnt.rep))
#           }
#         }
#       }
#     }
#   }
#   
#   return(sum(log(lik)))
# }

### calculate the log-likelihood of coalescent intervals 
lik_coaltimes <- function(phybreakenv) {
  if (phybreakenv$p$wh.model %in% c(1, 2, "single", "infinite")) 
    return(0)
  
  if(phybreakenv$p$wh.model == "linear" && phybreakenv$p$wh.bottleneck == "wide") {
    if(min(phybreakenv$v$inftimes) - min(phybreakenv$v$nodetimes[phybreakenv$v$nodetypes == "c"]) > 
       phybreakenv$p$sample.mean + phybreakenv$p$wh.level/phybreakenv$p$wh.slope) return(-Inf)
  }
  
  mult.intro <- ifelse(is.null(phybreakenv$p$mult.intro), FALSE, phybreakenv$p$mult.intro)
  
  remove0nodes <- phybreakenv$v$nodetypes != "0"
  nodetypes <- phybreakenv$v$nodetypes[remove0nodes]
  nodehosts <- phybreakenv$v$nodehosts[remove0nodes]
  nodetimes <- phybreakenv$v$nodetimes[remove0nodes]
  inftimes <- c(min(phybreakenv$v$inftimes) - phybreakenv$p$sample.mean, phybreakenv$v$inftimes)
  
  coalnodes <- nodetypes == "c"
  orderednodes <- order(nodehosts, nodetimes)
  
  coalnodes <- coalnodes[orderednodes]
  orderedhosts <- nodehosts[orderednodes]
  
  bottlenecks <- sapply(0:phybreakenv$p$obs, function(i) sum((orderedhosts == i) * (1 - 2 * coalnodes))) - 1
  dlineage <- 2 * c(FALSE, head(coalnodes, -1)) - 1
  dlineage[!duplicated(orderedhosts)] <- bottlenecks
  nrlineages <- 1 + cumsum(dlineage)
  
  whtimes <- nodetimes[orderednodes] - inftimes[orderedhosts + 1]
  if(mult.intro) whtimes[orderedhosts == 0] <- whtimes[orderedhosts == 0] - min(whtimes[orderedhosts == 0])
  whtimes[c(!duplicated(orderedhosts)[-1], FALSE)] <- 0
  
  logcoalrates <- switch(phybreakenv$p$wh.model, single =, infinite =,
                         linear = -log(phybreakenv$p$wh.level + phybreakenv$p$wh.slope * whtimes[coalnodes]),
                         exponential = 
                           -log(phybreakenv$p$wh.level * 
                                  exp(phybreakenv$p$wh.exponent * 
                                        whtimes[coalnodes])),
                         constant = -log(phybreakenv$p$wh.level) * coalnodes)
  cumcoalrates <- switch(phybreakenv$p$wh.model, single =, infinite=,
                       linear = log(whtimes + phybreakenv$p$wh.level/phybreakenv$p$wh.slope + 
                                      ((whtimes + phybreakenv$p$wh.level/phybreakenv$p$wh.slope) == 0)) / phybreakenv$p$wh.slope,
                       exponential  = -1/(phybreakenv$p$wh.level * phybreakenv$p$wh.exponent * 
                                            exp(phybreakenv$p$wh.exponent * whtimes)),
                       constant = whtimes/phybreakenv$p$wh.level)
  if (mult.intro) {
    logcoalrates[orderedhosts[coalnodes] == 0] <- -log(phybreakenv$p$wh.history)
    cumcoalrates[orderedhosts == 0] <- whtimes[orderedhosts == 0]/phybreakenv$p$wh.history
  }
  
  coalratediffs <- cumcoalrates - c(0, head(cumcoalrates, -1))
  logcoalescapes <- -coalratediffs * choose(nrlineages, 2)
  
  return(sum(logcoalrates) + sum(logcoalescapes))
}

### replacement proposed for lik_coaltimes
# to be used with original style phybreak object, ie without history host in v$inftimes or v$infectors
# to be used with original phybreak2environment (called in logLik.phybreak)
lik_coaltimes_new <- function(phybreakenv) {
  if (phybreakenv$p$wh.model %in% c(1, 2, "single", "infinite")) 
    return(0)
  
  if(phybreakenv$p$wh.model == "linear" && phybreakenv$p$wh.bottleneck == "wide") {
    if(min(phybreakenv$v$inftimes) - min(phybreakenv$v$nodetimes[phybreakenv$v$nodetypes == "c"]) > 
       phybreakenv$p$sample.mean + phybreakenv$p$wh.level/phybreakenv$p$wh.slope) return(-Inf)
  }
  
  remove0nodes <- phybreakenv$v$nodetypes != "0"
  nodetypes <- phybreakenv$v$nodetypes[remove0nodes]
  nodehosts <- phybreakenv$v$nodehosts[remove0nodes]
  nodetimes <- phybreakenv$v$nodetimes[remove0nodes]
  inftimes <- c(min(phybreakenv$v$inftimes) - phybreakenv$p$sample.mean, phybreakenv$v$inftimes)
  
  coalnodes <- nodetypes == "c"
  orderednodes <- order(nodehosts, nodetimes)
  
  coalnodes <- coalnodes[orderednodes]
  orderedhosts <- nodehosts[orderednodes]
  
  bottlenecks <- sapply(0:phybreakenv$p$obs, function(i) sum((orderedhosts == i) * (1 - 2 * coalnodes))) - 1
  dlineage <- 2 * c(FALSE, head(coalnodes, -1)) - 1
  dlineage[!duplicated(orderedhosts)] <- bottlenecks
  nrlineages <- 1 + cumsum(dlineage)
  
  whtimes <- nodetimes[orderednodes] - inftimes[orderedhosts + 1]
  if(phybreakenv$p$hist) whtimes[orderedhosts == 0] <- whtimes[orderedhosts == 0] - min(whtimes[orderedhosts == 0])
  whtimes[c(!duplicated(orderedhosts)[-1], FALSE)] <- 0
  
  logcoalrates <- switch(phybreakenv$p$wh.model, single =, infinite =,
                         linear = -log(phybreakenv$p$wh.level + phybreakenv$p$wh.slope * whtimes[coalnodes]),
                         exponential = 
                           -log(phybreakenv$p$wh.level * 
                                  exp(phybreakenv$p$wh.exponent * 
                                        whtimes[coalnodes])),
                         constant = -log(phybreakenv$p$wh.level) * coalnodes)
  cumcoalrates <- switch(phybreakenv$p$wh.model, single =, infinite=,
                         linear = log(whtimes + phybreakenv$p$wh.level/phybreakenv$p$wh.slope + 
                                        ((whtimes + phybreakenv$p$wh.level/phybreakenv$p$wh.slope) == 0)) / phybreakenv$p$wh.slope,
                         exponential  = -1/(phybreakenv$p$wh.level * phybreakenv$p$wh.exponent * 
                                              exp(phybreakenv$p$wh.exponent * whtimes)),
                         constant = whtimes/phybreakenv$p$wh.level)
  if(phybreakenv$p$hist) {
    logcoalrates[orderedhosts[coalnodes] == 0] <- -log(phybreakenv$p$wh.history) 
    cumcoalrates[orderedhosts == 0] <- whtimes[orderedhosts == 0]/phybreakenv$p$wh.history
  }
  
  
  coalratediffs <- cumcoalrates - c(0, head(cumcoalrates, -1))
  logcoalescapes <- -coalratediffs * choose(nrlineages, 2)
  
  return(sum(logcoalrates) + sum(logcoalescapes))
}


### calculate the log-likelihood of coalescent intervals in a single host
lik_coaltimes_host <- function(phybreakenv, hostID) {
  if (phybreakenv$p$wh.model %in% c(1, 2, "single", "infinite")) 
    return(0)
  
  if(phybreakenv$p$wh.model == "linear" && phybreakenv$p$wh.bottleneck == "wide") {
    if(min(phybreakenv$v$inftimes) - min(phybreakenv$v$nodetimes[phybreakenv$v$nodetypes == "c"]) > 
       phybreakenv$p$sample.mean + phybreakenv$p$wh.level/phybreakenv$p$wh.slope) return(-Inf)
  }
  
  selecthostnodes <- phybreakenv$v$nodehosts == hostID
  nodetypes <- phybreakenv$v$nodetypes[selecthostnodes]
  nodetimes <- phybreakenv$v$nodetimes[selecthostnodes]
  inftime <- phybreakenv$v$inftimes[hostID]
  
  coalnodes <- nodetypes == "c"
  orderednodes <- order(nodetimes)
  
  coalnodes <- coalnodes[orderednodes]

  bottlenecks <- sum(1 - 2 * coalnodes) - 1
  dlineage <- 2 * c(FALSE, head(coalnodes, -1)) - 1
  dlineage[1] <- bottlenecks
  nrlineages <- 1 + cumsum(dlineage)
  
  whtimes <- nodetimes[orderednodes] - inftime

  logcoalrates <- switch(phybreakenv$p$wh.model, single =, infinite =,
                         linear = -log(phybreakenv$p$wh.level + phybreakenv$p$wh.slope * whtimes[coalnodes]),
                         exponential = 
                           -log(phybreakenv$p$wh.level * 
                                  exp(phybreakenv$p$wh.exponent * 
                                        whtimes[coalnodes])),
                         constant = -log(phybreakenv$p$wh.level) * coalnodes)
  cumcoalrates <- switch(phybreakenv$p$wh.model, single =, infinite=,
                         linear = log(whtimes + phybreakenv$p$wh.level/phybreakenv$p$wh.slope + 
                                        ((whtimes + phybreakenv$p$wh.level/phybreakenv$p$wh.slope) == 0)) / phybreakenv$p$wh.slope,
                         exponential  = -1/(phybreakenv$p$wh.level * phybreakenv$p$wh.exponent * 
                                              exp(phybreakenv$p$wh.exponent * whtimes)),
                         constant = whtimes/phybreakenv$p$wh.level)
  coalratediffs <- cumcoalrates - c(0, head(cumcoalrates, -1))
  logcoalescapes <- -coalratediffs * choose(nrlineages, 2)
  
  return(sum(logcoalrates) + sum(logcoalescapes))
}


### calculate the log-likelihood of the within-host topology
lik_topology_host <- function(phybreakenv, hostID) {
  if (phybreakenv$p$wh.model %in% c(1, 2, "single", "infinite")) 
    return(0)
  
  selecthostnodes <- phybreakenv$v$nodehosts == hostID
  nodetypes <- phybreakenv$v$nodetypes[selecthostnodes]
  nodetimes <- phybreakenv$v$nodetimes[selecthostnodes]
  inftime <- phybreakenv$v$inftimes[hostID]
  
  coalnodes <- nodetypes == "c"
  orderednodes <- order(nodetimes)
  
  coalnodes <- coalnodes[orderednodes]
  
  bottlenecks <- sum(1 - 2 * coalnodes) - 1
  #dlineage <- 2 * c(FALSE, head(coalnodes, -1)) - 1
  dlineage <- 2 * c(FALSE, coalnodes[1:(length(coalnodes)-1)]) - 1
  dlineage[1] <- bottlenecks
  nrlineages <- 1 + cumsum(dlineage)
  
  #logcoalprobabilities <- -log(choose(nrlineages[c(FALSE, head(coalnodes, -1))], 2))
  logcoalprobabilities <- -log(choose(nrlineages[c(FALSE, coalnodes[1:(length(coalnodes)-1)])], 2))
  
  return(sum(logcoalprobabilities))
}
