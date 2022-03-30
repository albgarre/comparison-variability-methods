
## Load libraries

library(tidyverse)
library(readxl)
library(lme4)
library(rjags)
library(coda)

source("funcs_var_analysis_JS.R")

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

## Fit the method by Aryani

n_strain <- 20
n_bio <- 3
n_exp <- 2

dof_strain <- n_strain - 1
dof_bio <- n_strain*n_bio - n_strain
dof_exp <- n_strain*n_bio*n_exp - n_strain*n_bio

aryani_fit <- aryani_data %>% 
    filter(temp %in% c(55, 60, 65)) %>% # find rows where temp is 55/60/65
    mutate(log_D = log10(D_val)) %>% # 
    rename(bio_rep_c = bio_rep) %>%
    split(.$temp) %>% # splits tibble in 3 tibble for different temperatures
    map(group_variance) %>% # use functions as defined by Aryani
    map(., ~ extract_estimates(., dof_strain, dof_bio, dof_exp)) %>% # extract estimates
    # map(., ~ set_names(., paste0("var_", names(.)))) %>%
    imap_dfr(., ~ mutate(.x, temp = as.numeric(.y))) %>% # Apply a function to each element of a vector, and its index
    mutate(method = "Aryani") %>%
    gather(source, variance, -method, -temp)

## Fit the mixed-effects model

mixed_fit <- aryani_data %>% 
    filter(temp %in% c(55, 60, 65)) %>%
    mutate(log_D = log10(D_val)) %>%
    split(.$temp) %>%
    map(., ~ lmer(log_D ~ (1|strain)+ (1|bio_rep), data = ., REML = TRUE))

## Prepare the data for the Bayesian model

n_strain <- 20
n_bio <- 3
n_exp <- 2

prepare_data <- function(my_e) {
    
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
    
    x.ebs
    
}

my_x <- aryani_data %>%
    mutate(log_D = log10(D_val)) %>%
    filter(temp %in% c(55, 60, 65)) %>%
    split(.$temp) %>%
    map(prepare_data)

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

#' MCMC model with priors defined on the variance
#' 
mcmc_Dvalues_prior_exponential <- function(x.ebs, nstrain, nbio, nexp, 
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
                       "\n prec.str <- 1/var.str\n",
                       prior_b,
                       "\n prec.bio <- 1/var.bio\n",
                       prior_s,
                       "\n prec.exp <- 1/var.exp\n",
                       "}")
    
    # MCMC
    model=jags.model(textConnection(modelstring), data=data)
    update(model,n.iter=100)
    output=coda.samples(model=model,variable.names=c("icpt","var.str","var.bio","var.exp"), n.iter=niter, thin=1)
    # plot(output)
    # return(summary(output))
    return(output)
    
}

## Fitting the bayesian model with invgamma priors

set.seed(1212)

bayesian_prec_0p001 <- my_x %>%
    map(., 
        ~ mcmc_Dvalues_prior_precision(., n_strain, n_bio, n_exp, 
                                       "prec.exp ~ dgamma(0.001, 0.001)",
                                       "prec.bio ~ dgamma(0.001, 0.001)",
                                       "prec.str ~ dgamma(0.001, 0.001)",
                                       "icpt ~ dnorm(0, .0001)\n",
                                       niter = 1000)
    )

bayesian_prec_0p01 <- my_x %>%
    map(., 
        ~ mcmc_Dvalues_prior_precision(., n_strain, n_bio, n_exp, 
                                       "prec.exp ~ dgamma(0.01, 0.01)",
                                       "prec.bio ~ dgamma(0.01, 0.01)",
                                       "prec.str ~ dgamma(0.01, 0.01)",
                                       "icpt ~ dnorm(0, .0001)\n",
                                       niter = 1000)
    )

bayesian_prec_0p0001 <- my_x %>%
    map(., 
        ~ mcmc_Dvalues_prior_precision(., n_strain, n_bio, n_exp, 
                                       "prec.exp ~ dgamma(0.0001, 0.0001)",
                                       "prec.bio ~ dgamma(0.0001, 0.0001)",
                                       "prec.str ~ dgamma(0.0001, 0.0001)",
                                       "icpt ~ dnorm(0, .0001)\n",
                                       niter = 1000)
    )

aa <- bayesian_prec_0p001 %>%
    map(., ~ summary(.)[]$statistics) %>%
    map(as_tibble, rownames = "par") %>%
    imap_dfr(., ~ mutate(.x, temp = .y)) %>%
    mutate(prior = "invGamma", prior_val = 0.001)

bb <- bayesian_prec_0p01 %>%
    map(., ~ summary(.)[]$statistics) %>%
    map(as_tibble, rownames = "par") %>%
    imap_dfr(., ~ mutate(.x, temp = .y)) %>%
    mutate(prior = "invGamma", prior_val = 0.01)

fits_bayesian_gamma <- bayesian_prec_0p0001 %>%
    map(., ~ summary(.)[]$statistics) %>%
    map(as_tibble, rownames = "par") %>%
    imap_dfr(., ~ mutate(.x, temp = .y)) %>%
    mutate(prior = "invGamma", prior_val = 0.0001) %>%
    bind_rows(., aa, bb)

aa <- bayesian_prec_0p001 %>%
    map(., ~ summary(.)[]$quantiles) %>%
    map(as_tibble, rownames = "par") %>%
    imap_dfr(., ~ mutate(.x, temp = .y)) %>%
    mutate(prior = "invGamma", prior_val = 0.001)

bb <- bayesian_prec_0p01 %>%
    map(., ~ summary(.)[]$quantiles) %>%
    map(as_tibble, rownames = "par") %>%
    imap_dfr(., ~ mutate(.x, temp = .y)) %>%
    mutate(prior = "invGamma", prior_val = 0.01)

quantiles_bayesian_gamma <- bayesian_prec_0p0001 %>%
    map(., ~ summary(.)[]$quantiles) %>%
    map(as_tibble, rownames = "par") %>%
    imap_dfr(., ~ mutate(.x, temp = .y)) %>%
    mutate(prior = "invGamma", prior_val = 0.0001) %>%
    bind_rows(., aa, bb)

## Fitting the bayesian model with exponential priors

set.seed(12412)

bayesian_exp_0p001 <- my_x %>%
    map(.,
        ~ mcmc_Dvalues_prior_exponential(., n_strain, n_bio, n_exp, 
                                         "var.exp ~ dexp(0.001)",
                                         "var.bio ~ dexp(0.001)",
                                         "var.str ~ dexp(0.001)",
                                         "icpt ~ dnorm(0, .0001)\n",
                                         niter = 1000)
        )

bayesian_exp_0p1 <- my_x %>%
    map(.,
        ~ mcmc_Dvalues_prior_exponential(., n_strain, n_bio, n_exp, 
                                         "var.exp ~ dexp(0.1)",
                                         "var.bio ~ dexp(0.1)",
                                         "var.str ~ dexp(0.1)",
                                         "icpt ~ dnorm(0, .0001)\n",
                                         niter = 1000)
    )

bayesian_exp_0p01 <- my_x %>%
    map(.,
        ~ mcmc_Dvalues_prior_exponential(., n_strain, n_bio, n_exp, 
                                         "var.exp ~ dexp(0.01)",
                                         "var.bio ~ dexp(0.01)",
                                         "var.str ~ dexp(0.01)",
                                         "icpt ~ dnorm(0, .0001)\n",
                                         niter = 1000)
    )

aa <- bayesian_exp_0p01 %>%
    map(., ~ summary(.)[]$statistics) %>%
    map(as_tibble, rownames = "par") %>%
    imap_dfr(., ~ mutate(.x, temp = .y)) %>%
    mutate(prior = "exponential", prior_val = 0.01)

bb <- bayesian_exp_0p1 %>%
    map(., ~ summary(.)[]$statistics) %>%
    map(as_tibble, rownames = "par") %>%
    imap_dfr(., ~ mutate(.x, temp = .y)) %>%
    mutate(prior = "exponential", prior_val = 0.1)

fits_bayesian_exp <- bayesian_exp_0p001 %>%
    map(., ~ summary(.)[]$statistics) %>%
    map(as_tibble, rownames = "par") %>%
    imap_dfr(., ~ mutate(.x, temp = .y)) %>%
    mutate(prior = "exponential", prior_val = 0.001) %>%
    bind_rows(., aa, bb)

aa <- bayesian_exp_0p01 %>%
    map(., ~ summary(.)[]$quantiles) %>%
    map(as_tibble, rownames = "par") %>%
    imap_dfr(., ~ mutate(.x, temp = .y)) %>%
    mutate(prior = "exponential", prior_val = 0.01)

bb <- bayesian_exp_0p1 %>%
    map(., ~ summary(.)[]$quantiles) %>%
    map(as_tibble, rownames = "par") %>%
    imap_dfr(., ~ mutate(.x, temp = .y)) %>%
    mutate(prior = "exponential", prior_val = 0.1)

quant_bayesian_exp <- bayesian_exp_0p001 %>%
    map(., ~ summary(.)[]$quantiles) %>%
    map(as_tibble, rownames = "par") %>%
    imap_dfr(., ~ mutate(.x, temp = .y)) %>%
    mutate(prior = "exponential", prior_val = 0.001) %>%
    bind_rows(., aa, bb)

## Assess the impact of the priors

bind_rows(fits_bayesian_gamma, fits_bayesian_exp) %>%
    filter(par != "icpt") %>%
    mutate(full_prior = paste(prior, prior_val, sep = "-")) %>%
    ggplot(aes(x = full_prior, y = Mean, colour = factor(prior))) +
        geom_point(position = position_dodge(width = 1)) +
        geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                      position = position_dodge(width = 1)) +
        facet_wrap(par~temp, scales = "free") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Compare the Bayesian and other fits

my_vars <- mixed_fit %>%
    map(summary)  %>%
    map(., ~ .$varcor) %>%
    map(as.data.frame) %>%
    imap_dfr(., ~ mutate(.x, temp = as.numeric(.y))) %>%
    mutate(method = "mixed") %>%
    mutate(source = ifelse(grp == "bio_rep", "bio",
                           ifelse(grp == "Residual", "exp",
                                  "strain"))) %>%
    select(temp, method, source, variance = vcov) %>%
    bind_rows(., aryani_fit)

dof_strain <- n_strain - 1
dof_bio <- n_strain*n_bio - n_strain
dof_exp <- n_strain*n_bio*n_exp - n_strain*n_bio

dof_map <- tibble(source = rep(c("exp", "bio", "strain"), 3),
                  method = rep(c("Aryani", "mixed", "mixed_allT"), each = 3),
                  dof = c(rep(c(dof_exp, dof_bio, dof_strain), 2),
                          n_strain*n_bio*n_exp*3 - n_strain*n_bio*3, n_strain*n_bio*3 - n_strain*3, n_strain*3-1)
)

aa <- my_vars %>%
    filter(source != "total") %>%
    full_join(dof_map) %>%
    mutate(chi0p05 = qchisq(0.05, dof),
           chi0p95 = qchisq(0.95, dof)) %>%
    mutate(var_left = variance/chi0p95*dof,
           var_right = variance/chi0p05*dof) %>%
    mutate(source = factor(source, levels = c("exp", "bio", "strain")))

## Figure 4

# bind_rows(quant_bayesian_exp, quantiles_bayesian_gamma) %>%
#     select(prior, prior_val, temp, par,`50%`, `2.5%`, `97.5%`) %>%
#     mutate(model = "bayesian") %>%
#     unite("fullPrior", c("prior", "prior_val")) %>%
#     # .$fullPrior %>% unique()
#     filter(fullPrior %in% c("invGamma_1e-04", "exponential_0.1")) %>% 
#     separate(fullPrior, c("prior", "aa"), sep = "_") %>%
#     unite("method", c("model", "prior"), sep = "-") %>%
#     filter(par != "icpt") %>%
#     rename(var_left = `2.5%`, var_right =`97.5%`, variance = `50%`) %>%
#     select(-aa) %>%
#     separate("par", c("foo", "aa")) %>%
#     mutate(source = ifelse(aa == "str", "strain", aa)) %>%
#     select(-aa, -foo) %>%
#     mutate(temp = as.numeric(temp)) %>%
#     bind_rows(., aa) %>% 
#     filter(method != "mixed_allT") %>%
#     mutate(temp = paste0(temp, "ºC")) %>%
#     mutate(source = ifelse(source == "bio", "Within-strain",
#                            ifelse(source == "exp", "Experimental",
#                                   "Between-strain"))) %>%
#     mutate(source = factor(source, 
#                            levels = c("Experimental", "Within-strain", "Between-strain"))) %>%
#     ggplot(aes(x = source, y = variance, colour = method)) +
#     geom_point(position = position_dodge(1)) +
#     geom_errorbar(aes(ymin = var_left, ymax = var_right),
#                   position = position_dodge(1)) +
#     facet_wrap("temp", scales="free") +
#     theme_bw() +
#     xlab("") +
#     scale_y_continuous(name = "Estimated variance", 
#                        sec.axis = sec_axis(name = "Standard error", ~ sqrt(.))) +
#     theme(legend.position = "top",
#           legend.title = element_blank(),
#           axis.text = element_text(size = 14),
#           axis.title = element_text(size = 16),
#           legend.text = element_text(size = 14, face="plain"),
#           strip.text = element_text(size = 14))

bind_rows(quant_bayesian_exp, quantiles_bayesian_gamma) %>%
    select(prior, prior_val, temp, par,`50%`, `2.5%`, `97.5%`) %>%
    mutate(model = "Bayesian") %>%
    unite("fullPrior", c("prior", "prior_val")) %>%
    # .$fullPrior %>% unique()
    filter(fullPrior %in% c("invGamma_1e-04", "exponential_0.1")) %>% 
    separate(fullPrior, c("prior", "aa"), sep = "_") %>%
    unite("method", c("model", "prior"), sep = ": ") %>%
    filter(par != "icpt") %>%
    rename(var_left = `2.5%`, var_right =`97.5%`, variance = `50%`) %>%
    select(-aa) %>%
    separate("par", c("foo", "aa")) %>%
    mutate(source = ifelse(aa == "str", "strain", aa)) %>%
    select(-aa, -foo) %>%
    mutate(temp = as.numeric(temp)) %>%
    bind_rows(., aa) %>% 
    filter(method != "mixed_allT") %>%
    filter(temp == 55) %>%
    mutate(source = ifelse(source == "bio", "Within-strain",
                           ifelse(source == "exp", "Experimental",
                                  "Between-strain"))) %>%
    mutate(source = factor(source, 
                           levels = c("Between-strain", "Within-strain", "Experimental"))) %>%
    mutate(method = ifelse(method == "mixed", "Mixed-effects", method)) %>%
    mutate(method = factor(method, levels = c("Aryani", "Mixed-effects", "Bayesian: invGamma", "Bayesian: exponential"))) %>%
    ggplot(aes(x = method, y = variance, colour = source)) +
    geom_point(position = position_dodge(.5)) +
    geom_errorbar(aes(ymin = var_left, ymax = var_right),
                  position = position_dodge(.5), width = .5) +
    # theme_bw() +
    xlab("") +
    scale_y_continuous(name = expression(Estimated~variance~(log~min^2)), 
                       sec.axis = sec_axis(name = expression(Standard~deviation~(log~min)), ~ sqrt(.))) +
    theme(legend.position = "top",
          legend.title = element_blank(),
          # axis.text = element_text(size = 14),
          # axis.title = element_text(size = 16),
          legend.text = element_text(size = 14, face="plain"),
          # strip.text = element_text(size = 14)
          )

bind_rows(quant_bayesian_exp, quantiles_bayesian_gamma) %>%
    select(prior, prior_val, temp, par,`50%`, `2.5%`, `97.5%`) %>%
    mutate(model = "bayesian") %>%
    unite("fullPrior", c("prior", "prior_val")) %>%
    # .$fullPrior %>% unique()
    filter(fullPrior %in% c("invGamma_1e-04", "exponential_0.1")) %>% 
    separate(fullPrior, c("prior", "aa"), sep = "_") %>%
    unite("method", c("model", "prior"), sep = "-") %>%
    filter(par != "icpt") %>%
    rename(var_left = `2.5%`, var_right =`97.5%`, variance = `50%`) %>%
    select(-aa) %>%
    separate("par", c("foo", "aa")) %>%
    mutate(source = ifelse(aa == "str", "strain", aa)) %>%
    select(-aa, -foo) %>%
    mutate(temp = as.numeric(temp)) %>%
    bind_rows(., aa) %>% 
    filter(method != "mixed_allT") %>%
    filter(temp == 55) %>%
    mutate(source = ifelse(source == "bio", "Within-strain",
                           ifelse(source == "exp", "Experimental",
                                  "Between-strain"))) %>%
    mutate(source = factor(source, 
                           levels = c("Experimental", "Within-strain", "Between-strain"))) %>%
    group_by(method, source) %>%
    summarize(mean(variance), sqrt(mean(variance)))


















