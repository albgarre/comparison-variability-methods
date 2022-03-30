
## Libraries

library(tidyverse)
library(readxl)
library(rjags)

## Load data

aryani_data <- read_xlsx("./data/D-values.xlsx", skip = 2)

aryani_data <- aryani_data  %>%
    set_names(., paste(names(.), .[1,], .[2,], sep = "-")) %>% # add ref to biological replicate (1-3) and to experimental replicate (1-2) in colname 
    rename(., strain = "...1-NA-NA") %>% # rename first column
    slice(., 3:n()) %>% # select rows by position i.e. drop first 2 rows
    gather(condition, D_val, -strain) %>% # collapse columns into key-value pairs. "condition" is new column name. "D-val" contains corresponding values
    separate(condition, into = c("temp", "bio_rep", "tec_rep"), sep = "-") %>% # split condition column in 3.
    separate(temp, into = c("temp", "foo"), sep = 2) %>%
    select(-foo) %>%
    mutate(temp = as.numeric(temp)) %>% # change format of temp column
    mutate(bio_rep = LETTERS[as.numeric(bio_rep)]) %>%
    mutate(bio_rep = paste(strain, temp, bio_rep, sep = "-")) %>%
    mutate(tec_rep = paste(bio_rep, tec_rep, sep = "-")) %>%
    filter(!is.na(D_val)) # deletes rows with D-val not available

## Prepare the data

n_strain <- 20
n_bio <- 3
n_exp <- 2

my_e <- aryani_data %>%
    filter(temp == 55) %>%
    mutate(log_D = log10(D_val)) 

x.ebs <- array(NA, dim=c(n_strain, n_bio, n_exp))

tec_repv <- unlist(my_e[,"tec_rep"])
log_Dv <- unlist(my_e[,"log_D"])

for(k.str in 1:n_strain){
    # strn <- unique(my_e$strain)[k.str]
    strn <- paste0(unique(my_e$strain)[k.str],"-")
    for(k.bio in 1:n_bio){
        for(k.exp in 1:n_exp){
            
            x.ebs[k.str,k.bio,k.exp] <- log_Dv[intersect(grep(strn,tec_repv),
                                                         grep(paste0(LETTERS[k.bio],"-",k.exp),tec_repv))] # for ARYANI DATA
            # x.ebs[k.str,k.bio,k.exp] <- log_Dv[intersect(grep(strn,tec_repv),
            #                                              grep(paste0("-",k.bio,"-",k.exp),tec_repv))]  # for ALL DATA
            
        }
    }
}

## Functions for fitting the MCMC

#' MCMC model with priors defined on the precision
#' 
mcmc_Dvalues_prior_precision <- function(x.ebs, nstrain, nbio, nexp, 
                                         prior_e, prior_b, prior_s,
                                         prior_icpt,
                                         niter){
    
    
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
          ",
          prior_icpt,
          prior_e,
          "\n var.str <- 1/prec.str\n",
          prior_b,
          "\nvar.bio <- 1/prec.bio\n",
          prior_s,
          "\nvar.exp <- 1/prec.exp\n",
          "}")
    
    # MCMC
    model=jags.model(textConnection(modelstring), data=data)
    update(model,n.iter=100)
    output=coda.samples(model=model,variable.names=c("icpt","var.str","var.bio","var.exp"), n.iter=niter, thin=1)
    # plot(output)
    # return(summary(output))
    return(output)
    
}

## Some plots

ggplot(my_e) +
    stat_bin(aes(x = log_D, y = ..count../sum(..count..))) +
    # geom_histogram(aes(log_D)) +
    geom_line(aes(x, y),
              data = tibble(x = seq(-1, 3, length = 100), y = dnorm(x, 0, 2)),
              colour = "red",
              inherit.aes = FALSE)
sample_prior <- function(n_sims, eps_icpt, eps_exp, eps_bio, eps_str) {
    
    tibble(sim = rep(1:n_sims, each = n_exp*n_bio*n_strain),
           strain = rep( rep(1:n_strain, each = n_bio*n_exp), n_sims),
           bio = rep( rep(rep(1:n_bio, each = n_exp), n_strain), n_sims),
           exp = rep( rep(1:n_exp, n_bio*n_strain), n_sims)
    ) %>%
        
        mutate( # Priors
            icpt = rep(rnorm(n_sims, 0, 1/eps_icpt), each = n_strain*n_bio*n_exp), 
            prec.exp = rep(rgamma(n_sims, eps_exp, eps_exp), each = n_strain*n_bio*n_exp),
            prec.bio = rep(rgamma(n_sims, eps_bio, eps_bio), each = n_strain*n_bio*n_exp),
            prec.str = rep(rgamma(n_sims, eps_str, eps_str), each = n_strain*n_bio*n_exp)
        ) %>% 
        mutate( # Variances
            var.exp = 1/prec.exp,
            var.bio = 1/prec.bio,
            var.str = 1/prec.str
        ) %>%
        group_by(sim, strain) %>%
        mutate(logD_str = rnorm(1, icpt, var.str)) %>% 
        ungroup() %>%
        group_by(sim, strain, bio) %>%
        mutate(logD_bio = rnorm(1, logD_str, var.bio)) %>%
        ungroup() %>%
        group_by(sim, strain, bio, exp) %>%
        mutate(logD_exp = rnorm(1, logD_bio, var.exp))
    
}

aa <- sample_prior(1000, 10, 10, 10, 10)

ggplot(aa) +
    geom_density(aes(logD_exp)) +
    stat_bin(aes(x = log_D, y = ..count../sum(..count..)), data = my_e, binwidth = .1)

bb <- aa %>%
    group_by(sim, strain, bio) %>%
    mutate(X_BS = mean(logD_exp)) %>%
    ungroup() %>%
    group_by(sim, strain) %>%
    mutate(X_S = mean(logD_exp)) %>%
    ungroup() %>%
    group_by(sim) %>%
    mutate(SS_exp = (logD_exp - X_BS)^2,
           SS_bio = (X_BS - X_S)^2,
           SS_strain = (mean(logD_exp) - X_S)^2
           )


##


mcmc_Dvalues_prior_precision(x.ebs, n_strain, n_bio, n_exp, 
                             "prec.exp ~ dgamma(0.001, 0.001)",
                             "prec.bio ~ dgamma(0.001, 0.001)",
                             "prec.str ~ dgamma(0.001, 0.001)",
                             "icpt ~ dnorm(0, .0001)\n",
                             niter = 1000) %>% plot()





























