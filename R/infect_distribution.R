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
        prob <- dgamma(time - inftimes, 
                      shape = le$p$gen.shape, 
                      scale = le$p$gen.mean/le$p$gen.shape,
                      log = TRUE)
      else 
        prob <- dgamma(time - inftimes, 
                      shape = le$p$gen.shape, 
                      scale = le$p$gen.mean/le$p$gen.shape,
                      log = FALSE)
    } else {
      prob <- le$p$inf_function(time, inftimes, le,
                                 nodetimes, host, log)
    }
    

  ### User-defined generation distribution ###
  } else if(trans.model =="user") {
    #print(list(time, inftimes, le, nodetimes, host, log))
    prob <- le$p$inf_function(time, inftimes, le, 
                              nodetimes, host, log)
    return(prob)
  }

  if (le$p$contact){
    if (length(time) == 1){
      if (length(pbe0) == 0){
        # If in phybreak function
        host_i <- which(inftimes == time)
        prob.cnt <- sapply(1:length(inftimes), function(host_j){
          get_contact_probability(le$d$contact.matrix, host_i, host_j, le)
        })
        return(prob * prob.cnt)
      }

      # If calculating proposed ID
      if (log) return(prob + pbe0$contactarray[pbe1$hostID,,1])
      else return(prob * pbe0$contactarray[pbe1$hostID,,1])
    } else {
      # If calculating likelihood
      prob.cnt <- do.call(c, lapply(seq_along(le$v$infectors), function(i){
        if (le$v$infectors[i] > 0) pbe0$contactarray[le$v$infectors[i], i, 1]
        else return(1)
      }))
      if (log) return(prob + log(prob.cnt))
      else return(prob * prob.cnt)
    }
  }
}