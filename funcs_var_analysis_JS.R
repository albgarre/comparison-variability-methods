
## Function for simulations

simulate_spread <- function(n_strains, n_replicates, n_exp,
                            var_S, var_R, var_E, mean_logD = 0) {
    
    strain_error <- rnorm(n_strains, 0, sqrt(var_S))
    bio_error <- rnorm(n_strains*n_replicates, 0, sqrt(var_R))
    exp_error <- rnorm(n_replicates*n_strains*n_exp, 0, sqrt(var_E))
    
    data_frame(
        strain_e = rep(strain_error, rep(n_replicates*n_exp, n_strains)),
        biological_e = rep(bio_error, rep(n_exp, n_strains*n_replicates)),
        rep_e = exp_error
    ) %>%
        mutate(
            strain = rep(1:n_strains, rep(n_replicates*n_exp, n_strains)),
            bio_rep = rep(rep(1:n_replicates, rep(n_exp, n_replicates)), n_strains),
            tech_rep = rep(1:n_exp, n_replicates*n_strains)
        ) %>%
        mutate(total_e = strain_e + biological_e + rep_e) %>%
        mutate(bio_rep_c = paste(strain, bio_rep, sep = "-"),
               tech_rep_c = paste(bio_rep_c, tech_rep, sep = "-")) %>%
        mutate(strain = factor(strain), 
               bio_rep = factor(bio_rep),
               tech_rep = factor(tech_rep)) %>%
        mutate(log_D = total_e + mean_logD)
    
}

group_variance <- function(my_e) {
    
    my_e %>%
        group_by(bio_rep_c) %>%
        mutate(X_BS = mean(log_D)) %>%
        ungroup() %>%
        group_by(strain) %>%
        mutate(X_S = mean(log_D)) %>%
        ungroup() %>%
        mutate(SS_exp = (log_D - X_BS)^2,
               SS_bio = (X_BS - X_S)^2,
               SS_strain = (mean(log_D) - X_S)^2
        )
    
}

extract_estimates <- function(my_e, dof_strain, dof_bio, dof_exp) {
    
    strain_var <- my_e %>%
        ungroup() %>%
        group_by(strain, SS_strain) %>%
        summarize() %>%
        .$SS_strain %>%
        sum() 
    
    bio_var <- my_e %>%
        ungroup() %>%
        group_by(strain, SS_bio) %>%
        summarize() %>%
        .$SS_bio %>%
        sum() 
    
    exp_var <- my_e %>%
        ungroup() %>%
        group_by(strain, SS_exp) %>%
        summarize() %>%
        .$SS_exp %>%
        sum() 
    
    total_var <- var(my_e$log_D)
    
    tibble(total = total_var,
           strain = strain_var/dof_strain,
           bio = bio_var/dof_bio,
           exp = exp_var/dof_exp)
    
}


group_variance_mcmc_Dvalues_Listeria <- function(my_e, nstrain, nbio, nexp, eps_e, eps_b, eps_s, niter){
  
  x.ebs <- array(NA, dim=c(nstrain, nbio, nexp))

  tec_repv <- unlist(my_e[,"tec_rep"])
  log_Dv <- unlist(my_e[,"log_D"])
  
  for(k.str in 1:nstrain){
    # strn <- unique(my_e$strain)[k.str]
    strn <- paste0(unique(my_e$strain)[k.str],"-")
    for(k.bio in 1:nbio){
      for(k.exp in 1:nexp){
        
        x.ebs[k.str,k.bio,k.exp] <- log_Dv[intersect(grep(strn,tec_repv),
                                                     grep(paste0(LETTERS[k.bio],"-",k.exp),tec_repv))] # for ARYANI DATA
        # x.ebs[k.str,k.bio,k.exp] <- log_Dv[intersect(grep(strn,tec_repv),
        #                                              grep(paste0("-",k.bio,"-",k.exp),tec_repv))]  # for ALL DATA
        
      }
    }
  }
  
  data <- list("nstrain"=nstrain,"nbio"=nbio,"nexp"=nexp,"x.ebs"=x.ebs)
  
  modelstring=paste0("
  model {
  for(k.str in 1:nstrain){
    x.s[k.str] ~ dnorm(icpt,prec.str)
    for(k.bio in 1:nbio){
      x.bs[k.str,k.bio] ~ dnorm(x.s[k.str],prec.bio)
      for(k.exp in 1:nexp){
        x.ebs[k.str,k.bio,k.exp] ~ dnorm(x.bs[k.str,k.bio],prec.exp)
    }
    }
  }
  
  icpt ~ dnorm(0, .0001)
  prec.str ~ dgamma(",eps_s,",",eps_s,")
  var.str <- 1/prec.str
  prec.bio ~ dgamma(",eps_b,",",eps_b,")
  var.bio <- 1/prec.bio
  prec.exp ~ dgamma(",eps_e,",",eps_e,")
  var.exp <- 1/prec.exp
  
  }
  ")
  
  # MCMC
  model=jags.model(textConnection(modelstring), data=data)
  update(model,n.iter=100)
  output=coda.samples(model=model,variable.names=c("icpt","var.str","var.bio","var.exp"), n.iter=niter, thin=1)
  # plot(output)
  # return(summary(output))
  return(output)
  
}

group_variance_mcmc_Dvalues_Listeria_exponential <- function(my_e, nstrain, nbio, nexp, eps_e, eps_b, eps_s, niter){
  
  x.ebs <- array(NA, dim=c(nstrain, nbio, nexp))
  
  tec_repv <- unlist(my_e[,"tec_rep"])
  log_Dv <- unlist(my_e[,"log_D"])
  
  for(k.str in 1:nstrain){
    # strn <- unique(my_e$strain)[k.str]
    strn <- paste0(unique(my_e$strain)[k.str],"-")
    for(k.bio in 1:nbio){
      for(k.exp in 1:nexp){
        
        x.ebs[k.str,k.bio,k.exp] <- log_Dv[intersect(grep(strn,tec_repv),
                                                     grep(paste0(LETTERS[k.bio],"-",k.exp),tec_repv))] # for ARYANI DATA
        # x.ebs[k.str,k.bio,k.exp] <- log_Dv[intersect(grep(strn,tec_repv),
        #                                              grep(paste0("-",k.bio,"-",k.exp),tec_repv))]  # for ALL DATA
        
      }
    }
  }
  
  data <- list("nstrain"=nstrain,"nbio"=nbio,"nexp"=nexp,"x.ebs"=x.ebs)
  
  modelstring=paste0("
  model {
  for(k.str in 1:nstrain){
    x.s[k.str] ~ dnorm(icpt,prec.str)
    for(k.bio in 1:nbio){
      x.bs[k.str,k.bio] ~ dnorm(x.s[k.str],prec.bio)
      for(k.exp in 1:nexp){
        x.ebs[k.str,k.bio,k.exp] ~ dnorm(x.bs[k.str,k.bio],prec.exp)
    }
    }
  }
  
  icpt ~ dnorm(0, .0001)
  var.str ~ dexp(",eps_s,")
  prec.str <- 1/var.str
  var.bio ~ dexp(",eps_b,")
  prec.bio <- 1/var.bio
  var.exp ~ dexp(",eps_e,")
  prec.exp <- 1/var.exp
  
  }
  ")
  
  # MCMC
  model=jags.model(textConnection(modelstring), data=data)
  update(model,n.iter=100)
  output=coda.samples(model=model,variable.names=c("icpt","var.str","var.bio","var.exp"), n.iter=niter, thin=1)
  # plot(output)
  # return(summary(output))
  return(output)
  
}

group_variance_mcmc <- function(my_e, nstrain, nbio, nexp, eps_e, eps_b, eps_s){
  
  x.ebs <- array(NA, dim=c(nstrain, nbio, nexp))
  
  tec_repv <- unlist(my_e[,"tec_rep"])
  log_Dv <- unlist(my_e[,"log_D"])
  
  for(k.str in 1:nstrain){
    strn <- unique(my_e$strain)[k.str]
    for(k.bio in 1:nbio){
      for(k.exp in 1:nexp){
        
        # x.ebs[k.str,k.bio,k.exp] <- log_Dv[intersect(grep(strn,tec_repv),
        #                                              grep(paste0(LETTERS[k.bio],"-",k.exp),tec_repv))] # for ARYANI DATA
        x.ebs[k.str,k.bio,k.exp] <- log_Dv[intersect(grep(strn,tec_repv),
                                                     grep(paste0("-",k.bio,"-",k.exp),tec_repv))]  # for ALL DATA
        
      }
    }
  }
  
  data <- list("nstrain"=nstrain,"nbio"=nbio,"nexp"=nexp,"x.ebs"=x.ebs)
  
  modelstring=paste0("
  model {
  for(k.str in 1:nstrain){
  x.s[k.str] ~ dnorm(icpt,prec.str)
  for(k.bio in 1:nbio){
  x.bs[k.str,k.bio] ~ dnorm(x.s[k.str],prec.bio)
  for(k.exp in 1:nexp){
  x.ebs[k.str,k.bio,k.exp] ~ dnorm(x.bs[k.str,k.bio],prec.exp)
  }
  }
  }
  
  icpt ~ dnorm(0, .0001)
  prec.str ~ dgamma(",eps_s,",",eps_s,")
  var.str <- 1/prec.str
  prec.bio ~ dgamma(",eps_b,",",eps_b,")
  var.bio <- 1/prec.bio
  prec.exp ~ dgamma(",eps_e,",",eps_e,")
  var.exp <- 1/prec.exp
  
  }
  ")
  
  # MCMC
  model=jags.model(textConnection(modelstring), data=data)
  update(model,n.iter=100)
  output=coda.samples(model=model,variable.names=c("icpt","var.str","var.bio","var.exp"), n.iter=10000, thin=1)
  return(summary(output))
  
}