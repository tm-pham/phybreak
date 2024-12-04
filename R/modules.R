#' Create a phybreak-object from data and prior distributions.
#' @export
add_modules_to_phybreak <- function(le, 
                                    multiple.introductions = TRUE, 
                                    spatial = FALSE, 
                                    contact = FALSE, 
                                    infectivity = FALSE, 
                                    phyb.obj = TRUE, ...){
  
  extras <- list(...)
  extras.introductions <- extras[names(extras) %in% names(as.list(args(introductions_parameters)))]
  extras.spatial <- extras[names(extras) %in% names(as.list(args(spatial_parameters)))]
  extras.contact <- extras[names(extras) %in% names(as.list(args(contact_parameters)))]
  extras.infectivity <- extras[names(extras) %in% names(as.list(args(infectivity_parameters)))]
  
  if(phyb.obj == TRUE){
    if(multiple.introductions) le <- do.call(introductions_parameters, c(le, extras.introductions))
    else {
      le$parameterslot$mult.intro = FALSE
      le$parameterslot$wh.history = 1
    }
    
    if(spatial) le <- do.call(spatial_parameters, c(le, extras.spatial))
    else le$parameterslot$spatial = FALSE
    
    if(contact) le <- do.call(contact_parameters, c(le, extras.contact))
    else le$parameterslot$contact = FALSE
    
    if(infectivity | !is.null(extras$infectivity_file) | !is.null(extras$removal.times)){ 
      cat("Test: infectivity is called.\n")
      le <- do.call(infectivity_parameters, c(le, extras.infectivity))
    }else le$parameterslot$infectivity = FALSE
    
  } else {
    le[["likelihoods"]] <- list()
    le[["updaters"]] <- list()
    if(le$p$mult.intro) le <- introductions_functions(le)
    if(le$p$spatial) le <- spatial_functions(le)
    if(le$p$contact) le <- contact_functions(le)
    if(le$p$infectivity) le <- infectivity_functions(le)
  }
  
  return(le)
}

#####
### Multiple introductions ###

# Add data and parameters to phybreak.object
#' @export
introductions_parameters <- function(le, introductions = 1, 
    wh.history = 1, intro.rate = 1, reproduction.rate = 1,
    est.intro.rate = TRUE, prior.intro.rate.mean = 1, prior.intro.rate.shape = 1,
    est.reproduction = TRUE, prior.reproduction.rate = 0.1,
    est.wh.history = TRUE, prior.wh.history.shape = 1, prior.wh.history.mean = 100,
    use.NJtree = TRUE){
  
  # Parameterslot
  le$parameterslot <- c(le$parameterslot, list(
    introductions = introductions,
    intro.rate = intro.rate,
    R = reproduction.rate,
    mult.intro = TRUE,
    use.NJtree = use.NJtree,
    wh.history = wh.history))
  
  # Helperslot
  le$helperslot <- c(le$helperslot, list(
    si.ir = 2.38*sqrt(trigamma(introductions)),
    est.ir = est.intro.rate,
    est.R = est.reproduction,
    est.wh.h = est.wh.history,
    ir.av = prior.intro.rate.mean,
    ir.sh = prior.intro.rate.shape,
    R.rate = prior.reproduction.rate,
    wh.h.sh = prior.wh.history.shape,
    wh.h.av = prior.wh.history.mean))
  
  le$sampleslot <- c(le$sampleslot, list(
    introductions = c(),
    intro.rate = c(),
    wh.history = c(),
    R = c()))
  
  return(le)
}

# Add likelihood functions and parameter functions to object
introductions_functions <- function(le){
  le$updaters[["update_wh_history"]] <- function(){
    ### create an up-to-date proposal-environment
    prepare_pbe()
    
    ### making variables and parameters available within the function
    le <- environment()
    h <- pbe0$h
    p <- pbe1$p
    v <- pbe1$v
    
    ### check whether to estimate
    if(!h$est.wh.h) return()
    
    ### change to proposal state
    p$wh.history <- exp(log(p$wh.history) + rnorm(1, 0, h$si.wh))
    #if (p$wh.history > 1) return()
    
    ### update proposal environment
    copy2pbe1("p", le)
    
    ### calculate proposalratio
    logproposalratio <- log(p$wh.history) - log(pbe0$p$wh.history)
    
    ### calculate likelihood
    propose_pbe("wh.history")
    
    ### calculate acceptance probability
    logaccprob <- pbe1$logLikcoal - pbe0$logLikcoal + logproposalratio + 
      dgamma(pbe1$p$wh.history, shape = h$wh.h.sh, scale = h$wh.h.av/h$wh.h.sh, log = TRUE) - 
      dgamma(pbe0$p$wh.history, shape = h$wh.h.sh, scale = h$wh.h.av/h$wh.h.sh, log = TRUE)
    
    ### accept or reject
    if (runif(1) < exp(logaccprob)) {
      accept_pbe("wh.history")
    }
  }
  le$updaters[["update_ir"]] <- function() {
    ### create an up-to-date proposal-environment
    prepare_pbe()
    
    ### making variables and parameters available within the function
    le <- environment()
    h <- pbe0$h
    p <- pbe1$p
    v <- pbe1$v
    
    ### check whether to estimate
    if(!h$est.ir) return()
    
    p$intro.rate <- exp(log(p$intro.rate) + rnorm(1, 0, h$si.ir))
    
    ### update proposal environment
    copy2pbe1("p", le)
    
    ### calculate proposalratio
    logproposalratio <- log(p$intro.rate) - log(pbe0$p$intro.rate)
    
    ### calculate likelihood
    propose_pbe("ir")
    
    ### calculate acceptance probability
    logaccprob <- pbe1$logLikgen - pbe0$logLikgen + logproposalratio + 
      dgamma(pbe1$p$intro.rate, shape = h$ir.sh, scale = h$ir.av/h$ir.sh, log = TRUE) - 
      dgamma(pbe0$p$intro.rate, shape = h$ir.sh, scale = h$ir.av/h$ir.sh, log = TRUE)
    
    ### accept
    if (runif(1) < exp(logaccprob)) {
      accept_pbe("ir")
    }
  }

  le$updaters[["update_R"]] <- function(){
     ### create an up-to-date proposal-environment
    prepare_pbe()
    
    ### making variables and parameters available within the function
    le <- environment()
    h <- pbe0$h
    p <- pbe1$p
    v <- pbe1$v
    
    ### check whether to estimate
    if(!h$est.ir) return()
    
    p$R <- exp(log(p$R) + rnorm(1, 0, h$si.ir))
    
    ### update proposal environment
    copy2pbe1("p", le)
    
    ### calculate proposalratio
    logproposalratio <- log(p$R) - log(pbe0$p$R)
    
    ### calculate likelihood
    propose_pbe("R")
    
    ### calculate acceptance probability
    logaccprob <- pbe1$logLikgen - pbe0$logLikgen + logproposalratio
    
    ### accept
    if (runif(1) < exp(logaccprob)) {
      accept_pbe("R")
    }
  }
  
  return(le)
}

#####
### Spatial distance ###
spatial_parameters <- function(le,
    dist.model = "power", dist.exponent = 2, dist.scale = 1, dist.mean = 1,
    est.dist.exponent = TRUE, prior.dist.exponent.shape = 1, prior.dist.exponent.mean = 1,
    est.dist.scale = TRUE, prior.dist.scale.shape = 1, prior.dist.scale.mean = 1,
    est.dist.mean = TRUE, prior.dist.mean.shape = 1, prior.dist.mean.mean = 1){
  
  dist.model <- choose_distmodel(dist.model, le$dataset$distances)
  
  # Dataslot
  le$dataslot <- c(le$dataslot, list(
   locations = le$dataset$locations,
   distances = le$dataset$distances))
  
  if(is.null(le$dataslot$locations)){
    if(!is.null(le$dataslot$distances)) {
      l <- cmdscale(le$dataslot$distances)
    }
  } else {
    l <- le$dataslot$locations
  }
  if(exists("l")) {
    le$dataslot$area <- (max(l[,1])-min(l[,2])) * (max(l[,2]) - min(l[,2]))
  }
  
  #Parameterslot
  le$parameterslot <- c(le$parameterslot, list(
    spatial = TRUE,
    dist.model = dist.model,
    dist.exponent = dist.exponent,
    dist.scale = dist.scale,
    dist.mean = dist.mean))
  
  #Helperslot
  le$helperslot <- c(le$helperslot, list(
    si.dist = 2.38*sqrt(trigamma(le$parameterslot$obs - 1)),
    est.dist.e = est.dist.exponent && dist.model %in% c("power", "exponential"),
    est.dist.s = est.dist.scale && dist.model == "power",
    est.dist.m = est.dist.mean && dist.model == "poisson",
    dist.e.sh = prior.dist.exponent.shape,
    dist.e.av = prior.dist.exponent.mean,
    dist.s.sh = prior.dist.scale.shape,
    dist.s.av = prior.dist.scale.mean,
    dist.m.sh = prior.dist.mean.shape,
    dist.m.av = prior.dist.mean.mean))
  
  # Sampleslot
  le$sampleslot <- c(le$sampleslot, list(
    dist.exponent = c(),
    dist.scale = c(),
    dist.mean = c()))
  
  return(le)
}

spatial_functions <- function(le){
  
  # calculate the log-likelihood of distances
  le$likelihoods[["logLikdist"]] <- function(le){#dist.model, dist.exponent, dist.scale, dist.mean, infectors, distances, area) {
    dist.model = le$p$dist.model
    dist.scale = le$p$dist.scale
    dist.exponent = le$p$dist.exponent
    dist.mean = le$p$dist.mean
    infectors = le$v$infectors
    distances = le$d$distances
    area = le$d$area
    
    dist.model <- ifelse(is.null(dist.model), "none", dist.model)
    if(dist.model == "none") return(0)
    distancevector <- distances[cbind(which(infectors!=0), infectors[infectors!=0])]
    sum((infectors == 0) * dunif(1, min = 0, max = area, log = TRUE)) + 
      switch(dist.model,
             power = sum(log(
               dist.exponent * sin(pi/dist.exponent) / 
                 (dist.scale * pi * (1 + (distancevector/dist.scale)^dist.exponent))
             )),
             exponential = sum(
               log(dist.exponent) - dist.exponent * distancevector
             ),
             poisson = sum(
               -dist.mean + distancevector * log(dist.mean) - lgamma(1 + distancevector)
             )
      )
  }
  
  # update parameters of loglik of distances
  le$updaters[["update_dist_exponent"]] <- function() {
    ### create an up-to-date proposal-environment
    prepare_pbe()
    
    ### making variables and parameters available within the function
    le <- environment()
    h <- pbe0$h
    p <- pbe1$p
    v <- pbe1$v
    
    ### change to proposal state
    p$dist.exponent <- 1 + exp(log(p$dist.exponent - 1) + rnorm(1, 0, h$si.dist))
    
    ### update proposal environment
    copy2pbe1("p", le)
    
    ### calculate proposalratio
    logproposalratio <- log(p$dist.exponent - 1) - log(pbe0$p$dist.exponent - 1)
    
    ### calculate likelihood
    propose_pbe("dist.exponent")
    
    ### calculate acceptance probability
    logaccprob <- pbe1$logLikdist - pbe0$logLikdist + logproposalratio + 
      dgamma(pbe1$p$dist.exponent - 1, shape = h$dist.e.sh, scale = h$dist.e.av/h$dist.e.sh, log = TRUE) - 
      dgamma(pbe0$p$dist.exponent - 1, shape = h$dist.e.sh, scale = h$dist.e.av/h$dist.e.sh, log = TRUE)
    
    ### accept or reject
    if (runif(1) < exp(logaccprob)) {
      accept_pbe("dist.exponent")
    }
  }
  le$updaters[["update_dist_scale"]] <- function() {
    ### create an up-to-date proposal-environment
    prepare_pbe()
    
    ### making variables and parameters available within the function
    le <- environment()
    h <- pbe0$h
    p <- pbe1$p
    v <- pbe1$v
    
    ### change to proposal state
    p$dist.scale <- exp(log(p$dist.scale) + rnorm(1, 0, h$si.dist))
    
    ### update proposal environment
    copy2pbe1("p", le)
    
    ### calculate proposalratio
    logproposalratio <- log(p$dist.scale) - log(pbe0$p$dist.scale)
    
    ### calculate likelihood
    propose_pbe("dist.scale")
    
    ### calculate acceptance probability
    logaccprob <- pbe1$logLikdist - pbe0$logLikdist + logproposalratio + 
      dgamma(pbe1$p$dist.scale, shape = h$dist.s.sh, scale = h$dist.s.av/h$dist.s.sh, log = TRUE) - 
      dgamma(pbe0$p$dist.scale, shape = h$dist.s.sh, scale = h$dist.s.av/h$dist.s.sh, log = TRUE)
    
    ### accept or reject
    if (runif(1) < exp(logaccprob)) {
      accept_pbe("dist.scale")
    }
  }
  le$updaters[["update_dist_mean"]] <- function() {
    ### create an up-to-date proposal-environment
    prepare_pbe()
    
    ### making variables and parameters available within the function
    le <- environment()
    h <- pbe0$h
    p <- pbe1$p
    v <- pbe1$v
    
    ### change to proposal state
    p$dist.mean <- exp(log(p$dist.mean) + rnorm(1, 0, h$si.dist))
    
    ### update proposal environment
    copy2pbe1("p", le)
    
    ### calculate proposalratio
    logproposalratio <- log(p$dist.mean) - log(pbe0$p$dist.mean)
    
    ### calculate likelihood
    propose_pbe("dist.mean")
    
    ### calculate acceptance probability
    logaccprob <- pbe1$logLikdist - pbe0$logLikdist + logproposalratio + 
      dgamma(pbe1$p$dist.mean, shape = h$dist.m.sh, scale = h$dist.m.av/h$dist.m.sh, log = TRUE) - 
      dgamma(pbe0$p$dist.mean, shape = h$dist.m.sh, scale = h$dist.m.av/h$dist.m.sh, log = TRUE)
    
    ### accept or reject
    if (runif(1) < exp(logaccprob)) {
      accept_pbe("dist.mean")
    }
  }
  
  return(le)
}

#####
### Contact
#' Module for contact data
#' 
#' @param cnt.invest.trans      Probability for transmission pair that contact took place (\eqn{\eta})
#' @param cnt.invest.nontrans   Probability for non-transmission pair that contact took place (\eqn{\lambda})
#' @param cnt.report            Reporting coverage - probability that contact was reported (\eqn{\epsilon})
#' @param cnt.report.false      False reporting probability - probability that report was false (\eqn{\zeta})
#' 
#' @param est.cnt.invest.trans      Estimate \eqn{\eta}, default = TRUE
#' @param est.cnt.invest.nontrans   Estimate \eqn{\lambda}, default = TRUE
#' @param est.cnt.report            Estimate \eqn{\epsilon}, default = TRUE
#' @param prior.cnt.eta.alpha       Beta distributed prior for \eqn{\eta}, \eqn{\alpha} parameter
#' @param prior.cnt.eta.beta        Beta distributed prior for \eqn{\eta}, \eqn{\beta} parameter
#' @param prior.cnt.epsilon.alpha   Beta distributed prior for \eqn{\epsilon}, \eqn{\alpha} parameter
#' @param prior.cnt.epsilon.beta    Beta distributed prior for \eqn{\epsilon}, \eqn{\beta} parameter
#' @param prior.cnt.lambda.alpha    Beta distributed prior for \eqn{\lambda}, \eqn{\alpha} parameter
#' @param prior.cnt.lambda.beta     Beta distributed prior for \eqn{\lambda}, \eqn{\beta} parameter
#'  
#' @export
contact_parameters <- function(le,
    cnt.invest.trans = 1, cnt.invest.nontrans = 0.5, cnt.report = 1, cnt.report.false = 0,
    cnt.types = c(0.5, 0.5), est.cnt.types = F,
    est.cnt.invest.trans = T, est.cnt.invest.nontrans = T, est.cnt.report = T,
    prior.cnt.eta.alpha = 1, prior.cnt.eta.beta = 1,
    prior.cnt.lambda.alpha = 1, prior.cnt.lambda.beta = 1,
    prior.cnt.epsilon.alpha = 1, prior.cnt.epsilon.beta = 1){
  
  # Dataslot
  le$dataslot$contact <- le$dataset$contact.matrix
  
  # Parameterslot
  le$parameterslot <- c(le$parameterslot, list(
    contact = TRUE,
    cnt.eta = cnt.invest.trans,
    cnt.lambda = cnt.invest.nontrans,
    cnt.epsilon = cnt.report,
    cnt.zeta = cnt.report.false,
    cnt.types = cnt.types))
 
  # Helperslot
  le$helperslot <- c(le$helperslot, list(
    est.cnt.eta = est.cnt.invest.trans,
    est.cnt.lambda = est.cnt.invest.nontrans,
    est.cnt.epsilon = est.cnt.report,
    cnt.eta.alpha = prior.cnt.eta.alpha,
    cnt.eta.beta = prior.cnt.eta.beta,
    cnt.l.alpha = prior.cnt.lambda.alpha,
    cnt.l.beta = prior.cnt.lambda.beta,
    cnt.e.alpha = prior.cnt.epsilon.alpha,
    cnt.e.beta = prior.cnt.epsilon.beta,
    est.cnt.types = est.cnt.types
  ))
  
  # Sampleslot
  le$sampleslot <- c(le$sampleslot, list(
    cnt.eta = c(),
    cnt.lambda = c(),
    cnt.epsilon = c(),
    cnt.zeta = c(),
    cnt.types = matrix(NA, nrow = 2, ncol = 1)))
  
  return(le)
}

contact_functions <- function(le){
  
  # calculate the log-likelihood of contacts
  # le$likelihoods[["logLikcontact"]] <- function(le){#infectors, cnt.matrix, cnt.invest.trans, cnt.invest.nontrans,
  #                                         #cnt.rep, cnt.rep.false) {

  #   lik <- with(le, {
  #   #   if(is.null(d$contact)) return(0)
  #     lik <- unlist(lapply(seq_len(dim(contactarray)[3]), function(i){
  #       mat <- contactarray[,,i]
  #       diag(mat) <- NA
  #       mat <- mat #* p$cnt.types[i]
  #       return(mat)
  #     }))
  #     lik <- lik[!is.na(lik)]

  #     # lik <- c()
  #     # for(i in seq_len(ncol(d$contact))){
  #     #   for(j in seq_len(ncol(d$contact))){
  #     #     if (i != j){
  #     #       if (v$infectors[j] == i){
  #     #         #tau <- floor(v$inftimes[j])
  #     #         if (d$contact[i,j] == 1){
  #     #           lik <- c(lik, ((1-p$cnt.eta)*p$cnt.zeta + p$cnt.eta * p$cnt.epsilon)) # * contact.probs[categories[i], categories[j]])
  #     #         } else{
  #     #           lik <- c(lik, (1-p$cnt.eta)*(1-p$cnt.zeta) + p$cnt.eta*(1-p$cnt.epsilon)) #* (1-contact.probs[categories[i], categories[j]]))
  #     #         }
  #     #       } else if (v$infectors[i] == j){
  #     #         #tau <- floor(v$inftimes[i])
  #     #         if (d$contact[j,i] == 1){
  #     #           lik <- c(lik, ((1-p$cnt.eta)*p$cnt.zeta + p$cnt.eta * p$cnt.epsilon)) #* contact.probs[categories[i], categories[j]])
  #     #         } else{
  #     #           lik <- c(lik, (1-p$cnt.eta)*(1-p$cnt.zeta) + p$cnt.eta*(1-p$cnt.epsilon)) #* (1-contact.probs[categories[i], categories[j]]))
  #     #         }
  #     #       } else {
  #     #         if (d$contact[i,j] == 1){
  #     #           lik <- c(lik, (1-p$cnt.lambda)*p$cnt.zeta + p$cnt.lambda * p$cnt.epsilon) #* contact.probs[categories[i], categories[j]])
  #     #         } else {
  #     #           lik <- c(lik, (1-p$cnt.lambda)*(1-p$cnt.zeta) + p$cnt.lambda * (1-p$cnt.epsilon)) #* (1-contact.probs[categories[i], categories[j]]))
  #     #         }
  #     #       }
  #     #     }
  #     #   }
  #     # }
  #     #if (any(lik == 0)) stop("zeros in likelihood")
  #     #print(ifelse(is.infinite(sum(log(lik))), -1e5, sum(log(lik))))
  #     return(ifelse(is.infinite(sum(log(lik))), -1e5, sum(log(lik))))
  #   })
  #   return(lik)
  # }

  le$updaters[["update_cnt_eta"]] <- function(){
    ### create an up-to-date proposal-environment
    prepare_pbe()
    
    ### making variables and parameters available within the function
    le <- environment()
    h <- pbe0$h
    p <- pbe1$p
    v <- pbe1$v
    
    ### check whether to estimate
    if(!h$est.cnt.eta) return()
    
    ### change to proposal state
    p$cnt.eta <- rnorm(1, mean = p$cnt.eta, sd = 0.05)
    
    if (p$cnt.eta < 0 | p$cnt.eta > 1) return()
    
    ### update proposal environment
    copy2pbe1("p", le)
    
    ### calculate proposalratio
    logproposalratio <- log(p$cnt.eta) - log(pbe0$p$cnt.eta)
    
    ### calculate likelihood
    propose_pbe("contact")
    
    ### calculate acceptance probability
    logaccprob <- pbe1$logLikcontact - pbe0$logLikcontact + logproposalratio +
      dbeta(pbe1$p$cnt.eta, shape1 = h$cnt.eta.alpha, shape2 = h$cnt.eta.beta, log = TRUE) - 
      dbeta(pbe0$p$cnt.eta, shape1 = h$cnt.eta.alpha, shape2 = h$cnt.eta.beta, log = TRUE)
      
    ### accept or reject
    if (runif(1) < exp(logaccprob)) {
      accept_pbe("contact")
    }
  }

  # update parameters of loglik of contacts
  le$updaters[["update_cnt_lambda"]] <- function(){
    ### create an up-to-date proposal-environment
    prepare_pbe()
    
    ### making variables and parameters available within the function
    le <- environment()
    h <- pbe0$h
    p <- pbe1$p
    v <- pbe1$v
    
    ### check whether to estimate
    if(!h$est.cnt.lambda) return()
    
    ### change to proposal state
    p$cnt.lambda <- rnorm(1, mean = p$cnt.lambda, sd = 0.05)
    
    if (p$cnt.lambda < 0 | p$cnt.lambda > 1) return()
    
    ### update proposal environment
    copy2pbe1("p", le)
    
    ### calculate proposalratio
    logproposalratio <- log(p$cnt.lambda) - log(pbe0$p$cnt.lambda)
    
    ### calculate likelihood
    propose_pbe("contact")
    
    ### calculate acceptance probability
    logaccprob <- pbe1$logLikcontact - pbe0$logLikcontact + logproposalratio +
      dbeta(pbe1$p$cnt.lambda, shape1 = h$cnt.l.alpha, shape2 = h$cnt.l.beta, log = TRUE) - 
      dbeta(pbe0$p$cnt.lambda, shape1 = h$cnt.l.alpha, shape2 = h$cnt.l.beta, log = TRUE)

    ### accept or reject
    if (runif(1) < exp(logaccprob)) {
      accept_pbe("contact")
    }
  }
  
  le$updaters[["update_cnt_epsilon"]] <- function(){
    ### create an up-to-date proposal-environment
    prepare_pbe()
    
    ### making variables and parameters available within the function
    le <- environment()
    h <- pbe0$h
    p <- pbe1$p
    v <- pbe1$v
    
    ### check whether to estimate
    if(!h$est.cnt.epsilon) return()
    
    ### change to proposal state
    p$cnt.epsilon <- rnorm(1, mean = p$cnt.epsilon, sd = 0.05)
    
    if (p$cnt.epsilon < 0 | p$cnt.epsilon > 1) return()
    
    ### update proposal environment
    copy2pbe1("p", le)
    
    ### calculate proposalratio
    logproposalratio <- log(p$cnt.epsilon) - log(pbe0$p$cnt.epsilon)
    
    ### calculate likelihood
    propose_pbe("contact")
    
    ### calculate acceptance probability
    logaccprob <- pbe1$logLikcontact - pbe0$logLikcontact + logproposalratio +
      dbeta(pbe1$p$cnt.epsilon, shape1 = h$cnt.e.alpha, shape2 = h$cnt.e.beta, log = TRUE) - 
      dbeta(pbe0$p$cnt.epsilon, shape1 = h$cnt.e.alpha, shape2 = h$cnt.e.beta, log = TRUE)
    
    
    ### accept or reject
    if (runif(1) < exp(logaccprob)) {
      accept_pbe("contact")
    }
  }

  # le$updaters[["update_cnt_types"]] <- function(){
  #   prepare_pbe()

  #   ### making variables and parameters available within the function
  #   le <- environment()
  #   h <- pbe0$h
  #   p <- pbe1$p
  #   v <- pbe1$v

  #   ### change to proposal state
  #   p$cnt.type.1 <- rnorm(1, mean = p$cnt.types[1], sd = 0.05)
    
  #   if (p$cnt.epsilon < 0 | p$cnt.epsilon > 1){ 
  #     return()
  #   } else {
  #     p$cnt.types <- c(p$cnt.type.1, 1 - p$cnt.type.1)
  #   }

  #   ### update proposal environment
  #   copy2pbe1("p", le)
    
  #   ### calculate proposalratio
  #   logproposalratio <- sum(log(p$cnt.types)) - sum(log(pbe0$p$cnt.types))
    
  #   ### calculate likelihood
  #   propose_pbe("contact")
    
  #   ### calculate acceptance probability
  #   logaccprob <- pbe1$logLikcontact - pbe0$logLikcontact + logproposalratio
    
  #   ### accept or reject
  #   if (runif(1) < exp(logaccprob)) {
  #     accept_pbe("contact")
  #   }
  # }

  return(le)
}

#####
### Generation time
infectivity_parameters <- function(le, admission.times = NULL, removal.times = NULL,
    trans.model = "gamma", trans.growth = 1,
    infectivity_file = NULL){
  
  if(!is.null(infectivity_file)) trans.model = "user"
  # Use the Gamma distribution
  if(trans.model == "gamma"){
    # If no removal times, use standard Gamma distribution
    if (is.null(removal.times)) return(le)
    
    # If removal times are present, use adjusted Gamma distribution
    # Dataslot
    if (is.null(removal.times)) {
      if (!is.null(le$dataset$removal.times)) {
        removal.times = le$dataset$removal.times
      }
    }

    le$dataslot <- c(le$dataslot, list(
      admission.times = admission.times,
      removal.times = removal.times
    ))
    
    # Parameterslot
    le$parameterslot <- c(le$parameterslot, list(
      infectivity = TRUE,
      removal.rate = 5
    ))
    
    # Helperslot
    le$helperslot <- c(le$helperslot, list(
      est.removal.rate = FALSE 
    ))
    
    # Sampleslot
    le$sampleslot <- c(le$sampleslot, list(
      removal.rate = c()
    ))
    cat("trans model gamma: infectivity_parameters is called. Creating inf_function. \n")
    le$parameterslot[["inf_function"]] <- function(time, inftimes, le, nodetimes, 
                                                   host, log = FALSE,
                                                   test.arguments = FALSE){
      
      p <- le$p
      v <- le$v
      
      # Calculate normalization factor by calculating mean AUC of infectiousness function
      AUCs <- unlist(lapply(1:length(v$inftimes), function(i){
        remtime = as.numeric(le$d$removal.times)[i] - v$inftimes[i]
        remtime.prob = dgamma(remtime,
                              shape = le$p$gen.shape,
                              scale = le$p$gen.mean/le$p$gen.shape)
        probs <- sum(pgamma(remtime,
                            shape = le$p$gen.shape,
                            scale = le$p$gen.mean/le$p$gen.shape,
                            log = FALSE),
                     (1 - exp(-5*le$p$removal.rate) * (remtime.prob/le$p$removal.rate) ))
        return(probs)
      }))
      norm_factor <- 1/mean(AUCs)
      
      remtimes <- as.numeric(le$d$removal.times[match(inftimes, v$inftimes)] - inftimes)
      admission.times <- if (is.null(le$d$admission.times)) rep(0, length(v$inftimes)) else le$d$admission.times[match(inftimes,v$inftimes)] - inftimes
      hosttimes <- as.numeric(time - inftimes)
      
      C = p$removal.rate 
      
      if (length(hosttimes) == 0){
        probs = 1
      } else {
        
        if(is.null(host)){
          if(length(hosttimes) != length(nodetimes)){
            probs <- 0.1
            j <- 1
          } else {
            probs <- c()
            j <- 0
          }
          for (i in 1:length(remtimes)){
            if(hosttimes[i+j] < admission.times[i])
              probs <- c(probs, 0)
            else if(hosttimes[i+j] <= remtimes[i])
              probs <- c(probs, dgamma(hosttimes[i+j],
                                       shape = p$gen.shape,
                                       scale = p$gen.mean/p$gen.shape))
            else if(hosttimes[i+j] >= remtimes[i])
              probs <- c(probs, dgamma(hosttimes[i+j], 
                                       shape = p$gen.shape,
                                       rate = p$gen.mean/p$gen.shape) * exp(-C*(hosttimes[i+j]-remtimes[i])))
            # else 
            #   probs <- c(probs, 0)
          }
        }
      }
          
      if(log)
        return(log(probs*norm_factor))
      else 
        return(probs*norm_factor)
    }
    return(le)
  }
  
  else {
    # Check for removal times
    if (is.null(removal.times)){
      if (!is.null(le$dataset$removal.times)) {
        removal.times = le$dataset$removal.times
      } else {
        stop("Provide removal times in same order as hosts")
      }
    }
    
    # Load user-defined infectivity
    # if(trans.model == "user") {
    #   if(is.null(infectivity_file))
    #     stop("Please provide a R file stating the infectivity function")
    #   else 
    #     source(infectivity_file, local = userenv)
    # } else {
    #   datas <- NULL
    #   parameters <- NULL
    #   helpers <- NULL
    #   samplers <- NULL
    # }

    #Add data for user-defined function
    le$dataslot <- c(le$dataslot, list(
      removal.times = removal.times
    ))
    
    # Parameterslot
    le$parameterslot <- c(le$parameterslot, list(
      trans.init = 1e-4,
      trans.removal = 5,
      trans.growth = trans.growth,
      trans.sample = 1,
      trans.model = "user",
      infectivity = TRUE
    ))
    
    # Helperslot
    le$helperslot <- c(le$helperslot, userenv$helperslot)
    
    # Sampleslot
    le$sampleslot <- c(le$sampleslot, userenv$sampleslot)
    
    # Infectivity function
    # Changed inf_function instead of infect_function
    cat("trans model user: infectivity_parameters is called. Creating inf_function. \n")
    le$parameterslot[["inf_function"]] <- function(time, inftimes, le, nodetimes, 
                                                                      host, log = FALSE,
                                                                      test.arguments = FALSE){
      
      d <- le$d
      p <- le$p
      v <- le$v
      
      if (is.null(d$removal.times)) {
        stop("removal times of hosts must be provided")
      } else {
        removal.times <- d$removal.times
        if(class(removal.times[1]) == "Date"){
          removal.times <- as.numeric(removal.times - d$reference.date)
        }
      }
      
      if(test.arguments) return()
      
      if(is.null(p$trans.init))
        stop("initial fraction infected is missing")
      if(is.null(p$trans.growth))
        stop("growth factor of infectiousness is missing")
      if(is.null(p$trans.sample))
        stop("reduction factor after first positive sample is missing")
      if(is.null(p$trans.removal))
        stop("decay factor after removal is missing")
      
      a <- (1-p$trans.init)/p$trans.init
      r <- p$trans.growth
      S <- p$trans.sample
      C <- p$trans.removal
      
      # Calculate normalization factor by calculating mean AUC of infectiousness function
      AUCs <- unlist(lapply(1:length(v$inftimes), function(i){
        samtime = v$nodetimes[i] - v$inftimes[i]
        cultime = removal.times[i] - v$inftimes[i]
        if (r*samtime < 100){
          probs = sum((log(a+exp(r*samtime)) - log(a+1)) / r,
                      S * ( log(a+exp(r*cultime)) - log(a+exp(r*samtime)) ) / r,
                      (S / (1 + a*exp(-r*cultime))) / C)
        } else {
          probs = sum((r*samtime - log(a+1)) / r,
                      S * ( r*(cultime - samtime) ) / r,
                      S / C)
        }
        return(probs)
      }))
      norm_factor <- 1/mean(AUCs)
      
      # Use removal times of infectors in rest of calculations
      cultimes <- removal.times[match(inftimes, v$inftimes)]
      samtimes <- as.numeric(nodetimes - inftimes)
      cultimes <- as.numeric(cultimes - inftimes)
      hosttimes <- as.numeric(time - inftimes)
      
      if (length(hosttimes) == 0){
        probs = 1
      }else{
        
        if(is.null(host)){
          if(length(hosttimes) != length(nodetimes)){
            probs <- 0.1
            j <- 1
          } else {
            probs <- c()
            j <- 0
          }
          for (i in 1:length(samtimes)){
            if(hosttimes[i+j] < 0)
              probs <- c(probs, 0)
            else if(hosttimes[i+j] < samtimes[i])
              probs <- c(probs, 1/(1+a*exp(-r*hosttimes[i+j])))
            else if(hosttimes[i+j] >= samtimes[i] & hosttimes[i+j] < cultimes[i])
              probs <- c(probs, S/(1+a*exp(-r*hosttimes[i+j])))
            else if(hosttimes[i+j] >= cultimes[i] & hosttimes[i+j] < cultimes[i] + 5)
              probs <- c(probs, S/(1+a*exp(-r*cultimes[i])) * exp(-C*(hosttimes[i+j]-cultimes[i])))
            else 
              probs <- c(probs, 0)
          }
        }
      }
      
      if(log){
        return(log(probs*norm_factor))
      }else{
        return(probs*norm_factor)
      }

    }
    
    return(le)
  }
}

infectivity_functions <- function(le){
  le$updaters[["removal.rate"]] <- function(){
    
  }
  
  return(le)
}
