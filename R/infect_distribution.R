#' Infectiousness distribution.
#' 
#' This function calculates the probability of transmission given the infectiousness 
#' distribution for the hosts.
#' 
#' @param time The proposed infection time of the current host
#' @param vars The list of variables describing the nodes in the tree. This list should 
#'   at least include the infection times of the hosts. For \code{p$trans.model = "sample+culling"}
#'   also nodetimes of the sample nodes and the culling times of the hosts are needed as input.
#' @param p The list with parameters for the model.  
#' 
#' @return The transmission probabilities of all hosts infected before the current host. For each host
#'   with an infection time before the proposed infection time, the probability of being an infector
#'   of the current host is calculated according to the given infectiousness function.
#' 
#' @export
infect_distribution <- function(time, inftimes, le, 
                                nodetimes = NULL,  
                                host = NULL, log = FALSE){
  
  ### Gamma distributed ####
  trans.model <- ifelse(is.null(le$p$trans.model), "gamma", le$p$trans.model)
  
  if(trans.model == "gamma") {
    if(is.null(le$d$removal.times)){
    
      if(log)
        return(dgamma(time - inftimes, 
                      shape = le$p$gen.shape, 
                      scale = le$p$gen.mean/le$p$gen.shape,
                      log = TRUE))
      else 
        return(dgamma(time - inftimes, 
                      shape = le$p$gen.shape, 
                      scale = le$p$gen.mean/le$p$gen.shape,
                      log = FALSE))
    } else {
      return(le$p$inf_function(time, inftimes, le,
                               nodetimes, host, log))
    }
    

  ### User-defined generation distribution ###
  } else if(trans.model =="user") {
    #print(list(time, inftimes, le, nodetimes, host, log))
    return(le$p$inf_function(time, inftimes, le, 
                           nodetimes, host, log))
  }
}