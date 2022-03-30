
## Load libraries

library(tidyverse)
library(readxl)
library(lme4)
library(rjags)
library(coda)
library(cowplot)

source("funcs_var_analysis.R")

## Load the data

data_LA <- c("0mM", "3mM", "4mM") %>%
    set_names(., .) %>%
    map(., ~ read_excel("./data/Result_combined_LA_var_fitting.xlsx", sheet = .)) %>%
    imap_dfr(., ~ mutate(.x, cond = .y)) %>%
    mutate(bio_rep = paste(strain, cond, day, sep = "-"),
           tec_rep = paste(strain, cond, day, rep, sep = "-"),
           bio_rep_c = bio_rep,
           tec_rep_c = tec_rep,
           log_D = mu
    )

data_aw <- c("aw_10", "aw_9", "aw_7p5", "aw_5", "aw_2p5", "aw_0") %>%
    set_names(., .) %>%
    map(., ~ read_excel("./data/result combined_Aw_variability.xlsx", sheet = .)) %>%
    imap_dfr(., ~ mutate(.x, cond = .y)) %>%
    mutate(bio_rep = paste(strain, cond, day, sep = "-"),
           tec_rep = paste(strain, cond, day, repetition, sep = "-"),
           bio_rep_c = bio_rep,
           tec_rep_c = tec_rep,
           log_D = mu
    ) %>%
    rename(rep = repetition)


data_pH <- c("pH7", "pH6p5", "pH6", "pH5p5", "pH5") %>%
    set_names(., .) %>%
    map(., ~ read_excel("./data/result combined_pHvariability.xlsx", sheet = .)) %>%
    imap_dfr(., ~ mutate(.x, cond = .y)) %>%
    mutate(bio_rep = paste(strain, cond, day, sep = "-"),
           tec_rep = paste(strain, cond, day, rep, sep = "-"),
           bio_rep_c = bio_rep,
           tec_rep_c = tec_rep,
           log_D = mu
    )

data_T <- c("T5", "T10", "T20", "T30pH", "T30aw") %>%
    set_names(., .) %>%
    map(., ~ read_excel("./data/result combined_T_var.xlsx", sheet = .)) %>%
    imap_dfr(., ~ mutate(.x, cond = .y)) %>%
    mutate(bio_rep = paste(strain, cond, day, sep = "-"),
           tec_rep = paste(strain, cond, day, rep, sep = "-"),
           bio_rep_c = bio_rep,
           tec_rep_c = tec_rep,
           log_D = mu
    )

all_data <- bind_rows(data_LA, data_aw, data_pH, data_T)

## Fit the mixed-effects model

mixed_models <- all_data %>%
    split(.$cond) %>%
    map(., ~ lmer(mu ~ (1|strain)+ (1|bio_rep), data = ., REML = TRUE))

## Aryani model when bio == 3

n_strain <- 20
n_bio <- 3
n_exp <- 2

dof_strain <- n_strain - 1
dof_bio <- n_strain*n_bio - n_strain
dof_exp <- n_strain*n_bio*n_exp - n_strain*n_bio

Aryani_bio3 <- all_data %>%
    group_by(cond) %>%
    mutate(n_bio = max(day)) %>%
    filter(n_bio == 3) %>%
    na.omit() %>%
    split(.$cond) %>%
    map(., ~ group_variance(.)) %>%
    map(., ~ extract_estimates(., dof_strain, dof_bio, dof_exp)) %>%
    map(., ~ mutate_all(., sqrt)) %>%
    map(., ~ gather(., source, std)) %>%
    imap_dfr(., ~ mutate(.x, cond = .y)) %>%
    mutate(model = "Aryani")

# mixed_bio_3 <- all_data %>%
#     group_by(cond) %>%
#     mutate(n_bio = max(day)) %>%
#     filter(n_bio == 3) %>%
#     split(.$cond) %>%
#     map(., ~ lmer(mu ~ (1|strain)+ (1|bio_rep), data = ., REML = TRUE))

## Aryani model when bio == 4

n_strain <- 20
n_bio <- 4
n_exp <- 2

dof_strain <- n_strain - 1
dof_bio <- n_strain*n_bio - n_strain
dof_exp <- n_strain*n_bio*n_exp - n_strain*n_bio

Aryani_bio4 <- all_data %>%
    group_by(cond) %>%
    mutate(n_bio = max(day)) %>%
    filter(n_bio == 4) %>%
    na.omit() %>%
    split(.$cond) %>%
    map(., ~ group_variance(.)) %>%
    map(., ~ extract_estimates(., dof_strain, dof_bio, dof_exp)) %>%
    map(., ~ mutate_all(., sqrt)) %>%
    map(., ~ gather(., source, std)) %>%
    imap_dfr(., ~ mutate(.x, cond = .y)) %>%
    mutate(model = "Aryani")

# mixed_bio_4 <- all_data %>%
#     group_by(cond) %>%
#     mutate(n_bio = max(day)) %>%
#     filter(n_bio == 4) %>%
#     split(.$cond) %>%
#     map(., ~ lmer(mu ~ (1|strain)+ (1|bio_rep), data = ., REML = TRUE))


## Aryani model when bio == 6

n_strain <- 20
n_bio <- 6
n_exp <- 2

dof_strain <- n_strain - 1
dof_bio <- n_strain*n_bio - n_strain
dof_exp <- n_strain*n_bio*n_exp - n_strain*n_bio

Aryani_bio6 <- all_data %>%
    group_by(cond) %>%
    mutate(n_bio = max(day)) %>%
    filter(n_bio == 6) %>%
    na.omit() %>%
    split(.$cond) %>%
    map(., ~ group_variance(.)) %>%
    map(., ~ extract_estimates(., dof_strain, dof_bio, dof_exp)) %>%
    map(., ~ mutate_all(., sqrt)) %>%
    map(., ~ gather(., source, std)) %>%
    imap_dfr(., ~ mutate(.x, cond = .y)) %>%
    mutate(model = "Aryani")

# mixed_bio_6 <- all_data %>%
#     group_by(cond) %>%
#     mutate(n_bio = max(day)) %>%
#     filter(n_bio == 6) %>%
#     split(.$cond) %>%
#     map(., ~ lmer(mu ~ (1|strain)+ (1|bio_rep), data = ., REML = TRUE))

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
                
                # print(k.str)
                # print(k.bio)
                # print(k.exp)
                # print("*")
                # x.ebs[k.str,k.bio,k.exp] <- log_Dv[intersect(grep(strn,tec_repv),
                #                                              grep(paste0(LETTERS[k.bio],"-",k.exp),tec_repv))] # for ARYANI DATA
                x.ebs[k.str,k.bio,k.exp] <- log_Dv[intersect(grep(strn,tec_repv),
                                                             grep(paste0("-",k.bio,"-",k.exp),tec_repv))]  # for ALL DATA
                
            }
        }
    }
    
    x.ebs
    
}

x_aw10 <- all_data %>%
    filter(cond == "aw_10") %>%
    prepare_data()

x_aw5 <- all_data %>%
    filter(cond == "aw_5") %>%
    prepare_data()

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

## Bayesian model with inverse gamma priors

set.seed(97192)

pred_0p001 <- mcmc_Dvalues_prior_precision(x_aw5, n_strain, n_bio, n_exp, 
                             "prec.exp ~ dgamma(0.001, 0.001)",
                             "prec.bio ~ dgamma(0.001, 0.001)",
                             "prec.str ~ dgamma(0.001, 0.001)",
                             "icpt ~ dnorm(0, .0001)\n",
                             niter = 1000)

pred_0p01 <- mcmc_Dvalues_prior_precision(x_aw5, n_strain, n_bio, n_exp, 
                                           "prec.exp ~ dgamma(0.01, 0.01)",
                                           "prec.bio ~ dgamma(0.01, 0.01)",
                                           "prec.str ~ dgamma(0.01, 0.01)",
                                           "icpt ~ dnorm(0, .0001)\n",
                                           niter = 1000)

pred_0p0001 <- mcmc_Dvalues_prior_precision(x_aw5, n_strain, n_bio, n_exp, 
                                           "prec.exp ~ dgamma(0.0001, 0.0001)",
                                           "prec.bio ~ dgamma(0.0001, 0.0001)",
                                           "prec.str ~ dgamma(0.0001, 0.0001)",
                                           "icpt ~ dnorm(0, .0001)\n",
                                           niter = 1000)


pars_prec <- list(pred_0p0001,
     pred_0p001,
     pred_0p01) %>%
    map(., ~ summary(.)[]$statistics) %>%
    map(as_tibble, rownames = "par") %>%
    map2_dfr(., list(.0001, 0.001, .01),
         ~ mutate(.x, prior_val = .y)) %>%
    mutate(prior = "invGamma") 

quant_prec <- list(pred_0p0001,
     pred_0p001,
     pred_0p01) %>%
    map(., ~ summary(.)[]$quantiles) %>%
    map(as_tibble, rownames = "par") %>%
    map2_dfr(., list(.0001, 0.001, .01),
             ~ mutate(.x, prior_val = .y)) %>%
    mutate(prior = "invGamma") 

pred_0p001_aw10 <- mcmc_Dvalues_prior_precision(x_aw10, n_strain, n_bio, n_exp, 
                                           "prec.exp ~ dgamma(0.001, 0.001)",
                                           "prec.bio ~ dgamma(0.001, 0.001)",
                                           "prec.str ~ dgamma(0.001, 0.001)",
                                           "icpt ~ dnorm(0, .0001)\n",
                                           niter = 1000)

pred_0p01_aw10 <- mcmc_Dvalues_prior_precision(x_aw10, n_strain, n_bio, n_exp, 
                                          "prec.exp ~ dgamma(0.01, 0.01)",
                                          "prec.bio ~ dgamma(0.01, 0.01)",
                                          "prec.str ~ dgamma(0.01, 0.01)",
                                          "icpt ~ dnorm(0, .0001)\n",
                                          niter = 1000)

pred_0p0001_aw10 <- mcmc_Dvalues_prior_precision(x_aw10, n_strain, n_bio, n_exp, 
                                            "prec.exp ~ dgamma(0.0001, 0.0001)",
                                            "prec.bio ~ dgamma(0.0001, 0.0001)",
                                            "prec.str ~ dgamma(0.0001, 0.0001)",
                                            "icpt ~ dnorm(0, .0001)\n",
                                            niter = 1000)


quant_prec_aw10 <- list(pred_0p0001_aw10,
                   pred_0p001_aw10,
                   pred_0p01_aw10) %>%
    map(., ~ summary(.)[]$quantiles) %>%
    map(as_tibble, rownames = "par") %>%
    map2_dfr(., list(.0001, 0.001, .01),
             ~ mutate(.x, prior_val = .y)) %>%
    mutate(prior = "invGamma") 

## Fitting the bayesian model with exponential priors

set.seed(12412)

bayesian_exp_0p001 <- mcmc_Dvalues_prior_exponential(x_aw5, n_strain, n_bio, n_exp, 
                                                     "var.exp ~ dexp(0.001)",
                                                     "var.bio ~ dexp(0.001)",
                                                     "var.str ~ dexp(0.001)",
                                                     "icpt ~ dnorm(0, .0001)\n",
                                                     niter = 1000)

bayesian_exp_0p1 <- mcmc_Dvalues_prior_exponential(x_aw5, n_strain, n_bio, n_exp, 
                                                   "var.exp ~ dexp(0.1)",
                                                   "var.bio ~ dexp(0.1)",
                                                   "var.str ~ dexp(0.1)",
                                                   "icpt ~ dnorm(0, .0001)\n",
                                                   niter = 1000)

bayesian_exp_0p01 <- mcmc_Dvalues_prior_exponential(x_aw5, n_strain, n_bio, n_exp, 
                                                    "var.exp ~ dexp(0.01)",
                                                    "var.bio ~ dexp(0.01)",
                                                    "var.str ~ dexp(0.01)",
                                                    "icpt ~ dnorm(0, .0001)\n",
                                                    niter = 1000)




pars_exp <- list(bayesian_exp_0p001,
                 bayesian_exp_0p1,
                 bayesian_exp_0p01) %>%
    map(., ~ summary(.)[]$statistics) %>%
    map(as_tibble, rownames = "par") %>%
    map2_dfr(., list(.001, 0.1, .01),
             ~ mutate(.x, prior_val = .y)) %>%
    mutate(prior = "exponential") 

quant_exp <- list(bayesian_exp_0p001,
                 bayesian_exp_0p1,
                 bayesian_exp_0p01) %>%
    map(., ~ summary(.)[]$quantiles) %>%
    map(as_tibble, rownames = "par") %>%
    map2_dfr(., list(.001, 0.1, .01),
             ~ mutate(.x, prior_val = .y)) %>%
    mutate(prior = "exponential") 

set.seed(12412)

bayesian_exp_0p001_aw10 <- mcmc_Dvalues_prior_exponential(x_aw10, n_strain, n_bio, n_exp, 
                                                     "var.exp ~ dexp(0.001)",
                                                     "var.bio ~ dexp(0.001)",
                                                     "var.str ~ dexp(0.001)",
                                                     "icpt ~ dnorm(0, .0001)\n",
                                                     niter = 1000)

bayesian_exp_0p1_aw10 <- mcmc_Dvalues_prior_exponential(x_aw10, n_strain, n_bio, n_exp, 
                                                   "var.exp ~ dexp(0.1)",
                                                   "var.bio ~ dexp(0.1)",
                                                   "var.str ~ dexp(0.1)",
                                                   "icpt ~ dnorm(0, .0001)\n",
                                                   niter = 1000)

bayesian_exp_0p01_aw10 <- mcmc_Dvalues_prior_exponential(x_aw10, n_strain, n_bio, n_exp, 
                                                    "var.exp ~ dexp(0.01)",
                                                    "var.bio ~ dexp(0.01)",
                                                    "var.str ~ dexp(0.01)",
                                                    "icpt ~ dnorm(0, .0001)\n",
                                                    niter = 1000)

quant_exp_aw10 <- list(bayesian_exp_0p001_aw10,
                  bayesian_exp_0p1_aw10,
                  bayesian_exp_0p01_aw10) %>%
    map(., ~ summary(.)[]$quantiles) %>%
    map(as_tibble, rownames = "par") %>%
    map2_dfr(., list(.001, 0.1, .01),
             ~ mutate(.x, prior_val = .y)) %>%
    mutate(prior = "exponential") 

## Impact of the priors

bind_rows(pars_exp, pars_prec) %>%
    filter(par != "icpt") %>%
    unite(full_prior, c("prior", "prior_val")) %>%
    ggplot(aes(x = par, y = Mean, colour = factor(full_prior))) +
    geom_point(position = position_dodge(1)) +
    geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                  position = position_dodge(1))

## Compare the parameter estimates

n_strain <- 20
n_exp <- 2

dof_map <- all_data %>%
    group_by(cond) %>%
    summarize(n_bio = max(day)) %>%
    mutate(
        dof_strain = n_strain - 1,
        dof_bio = n_strain*n_bio - n_strain,
        dof_exp = n_strain*n_bio*n_exp - n_strain*n_bio
    ) %>%
    select(-n_bio) %>%
    gather(source, dof, -cond) %>%
    mutate(source = gsub("dof_", "", source))

my_pars <- mixed_models %>%
    map(summary) %>%
    map(., ~ .$varcor) %>%
    map(as.data.frame) %>%
    imap_dfr(., ~ mutate(.x, cond = .y)) %>%
    mutate(model = "mixed") %>%
    mutate(source = ifelse(grp == "bio_rep", "bio",
                           ifelse(grp == "Residual", "exp",
                                  "strain"))) %>%
    select(source, model, std = sdcor, cond) %>%
    bind_rows(., Aryani_bio3, Aryani_bio4, Aryani_bio6) %>%
    filter(source != "total") %>%
    full_join(., dof_map) %>%
    mutate(variance = std^2) %>%
    mutate(chi0p05 = qchisq(0.05, dof),
           chi0p95 = qchisq(0.95, dof)) %>%
    mutate(var_left = variance/chi0p95*dof,
           var_right = variance/chi0p05*dof) %>%
    mutate(source = factor(source, levels = c("exp", "bio", "strain"))) %>%
    select(model, source, cond, variance, var_left, var_right)
    
p1 <- bind_rows(quant_exp, quant_prec) %>%
    filter(par != "icpt") %>%
    unite("fullPrior", c("prior", "prior_val")) %>%
    filter(fullPrior %in% c("exponential_0.1", "invGamma_1e-04")) %>%
    separate(fullPrior, c("model", "foo"), sep = "_") %>%
    separate(par, c("foo", "source"), sep = "\\.") %>%
    mutate(model = paste("Bayesian:", model)) %>%
    select(model, source, variance = `50%`, var_left = `2.5%`, var_right = `97.5%`) %>%
    mutate(cond = "aw_5") %>%
    bind_rows(., my_pars) %>% 
    filter(cond == "aw_5") %>%
    mutate(source = ifelse(source == "strain", "str", source)) %>%
    mutate(source = ifelse(source == "str", "Between-strain",
                           ifelse(source == "exp", "Experimental",
                                  "Within-strain"))) %>%
    mutate(source = factor(source,
                           levels = c("Between-strain", "Within-strain", "Experimental"))) %>%
    mutate(model = ifelse(model == "mixed", "Mixed-effects", model)) %>%
    mutate(model = factor(model, levels = c("Aryani", "Mixed-effects", "Bayesian: invGamma", "Bayesian: exponential"))) %>%
    ggplot(aes(x = model, y = variance, colour = factor(source))) +
    geom_point(position = position_dodge(.5)) +
    geom_errorbar(aes(ymin = var_left, ymax = var_right),
                  position = position_dodge(.5), width = .5) +
    scale_y_continuous(name = expression(Estimated~variance~(h^-2)),
                       sec.axis = sec_axis(name = expression(Standard~deviation~(h^-1)), ~sqrt(.))) +
    xlab("") +
    theme(legend.position = "top",
          legend.title = element_blank())

# p2 <- bind_rows(quant_exp_aw10, quant_prec_aw10) %>%
#     filter(par != "icpt") %>%
#     unite("fullPrior", c("prior", "prior_val")) %>%
#     filter(fullPrior %in% c("exponential_0.1", "invGamma_1e-04")) %>%
#     separate(fullPrior, c("model", "foo"), sep = "_") %>%
#     separate(par, c("foo", "source"), sep = "\\.") %>%
#     mutate(model = paste("Bayesian:", model)) %>%
#     select(model, source, variance = `50%`, var_left = `2.5%`, var_right = `97.5%`) %>%
#     mutate(cond = "aw_10") %>%
#     bind_rows(., my_pars) %>% 
#     filter(cond == "aw_10") %>%
#     mutate(source = ifelse(source == "strain", "str", source)) %>%
#     mutate(source = ifelse(source == "str", "Between-strain",
#                            ifelse(source == "exp", "Experimental",
#                                   "Within-strain"))) %>%
#     mutate(source = factor(source,
#                            levels = c("Between-strain", "Within-strain", "Experimental"))) %>%
#     mutate(model = factor(model, levels = c("Aryani", "Mixed-effects", "Bayesian:invGamma", "Bayesian:exponential"))) %>%
#     ggplot(aes(x = model, y = variance, colour = factor(source))) +
#     geom_point(position = position_dodge(.5)) +
#     geom_errorbar(aes(ymin = var_left, ymax = var_right),
#                   position = position_dodge(.5), width = .5) +
#     scale_y_continuous(name = "Estimated variance",
#                        sec.axis = sec_axis(name = "Standard deviation", ~sqrt(.))) +
#     xlab("") +
#     theme(legend.position = "top",
#           legend.title = element_blank())
#     
# 
# plot_grid(p1, p2, labels = "AUTO", nrow = 2)

p1

####

bind_rows(quant_exp, quant_prec) %>%
    filter(par != "icpt") %>%
    unite("fullPrior", c("prior", "prior_val")) %>%
    filter(fullPrior %in% c("exponential_0.1", "invGamma_1e-04")) %>%
    separate(fullPrior, c("model", "foo"), sep = "_") %>%
    separate(par, c("foo", "source"), sep = "\\.") %>%
    mutate(model = paste("Bayesian:", model)) %>%
    select(model, source, variance = `50%`, var_left = `2.5%`, var_right = `97.5%`) %>%
    mutate(cond = "aw_5") %>%
    bind_rows(., my_pars) %>% 
    filter(cond == "aw_5") %>%
    mutate(source = ifelse(source == "strain", "str", source)) %>%
    mutate(source = ifelse(source == "str", "Between-strain",
                           ifelse(source == "exp", "Experimental",
                                  "Within-strain"))) %>%
    mutate(source = factor(source,
                           levels = c("Between-strain", "Within-strain", "Experimental"))) %>%
    group_by(model, source) %>%
    summarize(mean(variance), sqrt(mean(variance)))

## Supp Table 1

my_pars %>%
    filter(model == "mixed") %>%
    mutate(source = ifelse(source == "strain", "str", source)) %>%
    mutate(source = ifelse(source == "str", "Between-strain",
                           ifelse(source == "exp", "Experimental",
                                  "Within-strain"))) %>%
    write_excel_csv(., "supp_Table_1.csv")


# ## Compare the results
# 
# mixed_models %>%
#     map(summary) %>%
#     map(., ~ .$varcor) %>%
#     map(as.data.frame) %>%
#     imap_dfr(., ~ mutate(.x, cond = .y)) %>%
#     mutate(model = "mixed") %>%
#     mutate(source = ifelse(grp == "bio_rep", "bio",
#                            ifelse(grp == "Residual", "exp",
#                                   "strain"))) %>%
#     select(source, model, std = sdcor, cond) %>%
#     bind_rows(., Aryani_bio3, Aryani_bio4, Aryani_bio6) %>%
#     filter(source != "total") %>%
#     mutate(source = factor(source, levels = c("exp", "bio", "strain"))) %>%
#     ggplot() + 
#     geom_bar(aes(x = cond, y = std, fill = model), stat = "identity", position = "dodge") +
#     xlab("") + ylab("Estimate of the contribution to the standard deviation of mu") +
#     facet_wrap("source")
# 
# mixed_models %>%
#     map(summary) %>%
#     map(., ~ .$varcor) %>%
#     map(as.data.frame) %>%
#     imap_dfr(., ~ mutate(.x, cond = .y)) %>%
#     mutate(model = "mixed") %>%
#     mutate(source = ifelse(grp == "bio_rep", "bio",
#                            ifelse(grp == "Residual", "exp",
#                                   "strain"))) %>%
#     select(source, model, std = sdcor, cond) %>%
#     bind_rows(., Aryani_bio3, Aryani_bio4, Aryani_bio6) %>%
#     filter(source != "total") %>%
#     mutate(source = factor(source, levels = c("exp", "bio", "strain"))) %>%
#     ggplot() + 
#     geom_bar(aes(x = source, y = std, fill = model), stat = "identity", position = "dodge") +
#     xlab("") + ylab("Estimate of the contribution to the standard deviation of mu") +
#     facet_wrap("cond")
# 
# 
# ## Add the variances
# 
# n_strain <- 20
# n_exp <- 2
# 
# dof_map <- all_data %>%
#     group_by(cond) %>%
#     summarize(n_bio = max(day)) %>%
#     mutate(
#         dof_strain = n_strain - 1,
#         dof_bio = n_strain*n_bio - n_strain,
#         dof_exp = n_strain*n_bio*n_exp - n_strain*n_bio
#     ) %>%
#     select(-n_bio) %>%
#     gather(source, dof, -cond) %>%
#     mutate(source = gsub("dof_", "", source))
# 
# mixed_models %>%
#     map(summary) %>%
#     map(., ~ .$varcor) %>%
#     map(as.data.frame) %>%
#     imap_dfr(., ~ mutate(.x, cond = .y)) %>%
#     mutate(model = "mixed") %>%
#     mutate(source = ifelse(grp == "bio_rep", "bio",
#                            ifelse(grp == "Residual", "exp",
#                                   "strain"))) %>%
#     select(source, model, std = sdcor, cond) %>%
#     bind_rows(., Aryani_bio3, Aryani_bio4, Aryani_bio6) %>%
#     filter(source != "total") %>%
#     full_join(., dof_map) %>% 
#     mutate(variance = std^2) %>%
#     mutate(chi0p05 = qchisq(0.05, dof),
#            chi0p95 = qchisq(0.95, dof)) %>%
#     mutate(var_left = variance/chi0p95*dof,
#            var_right = variance/chi0p05*dof) %>%
#     mutate(source = factor(source, levels = c("exp", "bio", "strain"))) %>%
#     filter(model == "mixed") %>%
#     ggplot(aes(x = cond, y = variance)) +
#     geom_point() +
#     geom_errorbar(aes(ymin = var_left, ymax = var_right)) +
#     facet_wrap("source") +
#     scale_y_sqrt(name = "Variance",
#                  sec.axis = sec_axis(~sqrt(.), name = "Standard deviation")) +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#     xlab("")
# 
# mixed_models %>%
#     map(summary) %>%
#     map(., ~ .$varcor) %>%
#     map(as.data.frame) %>%
#     imap_dfr(., ~ mutate(.x, cond = .y)) %>%
#     mutate(model = "mixed") %>%
#     mutate(source = ifelse(grp == "bio_rep", "bio",
#                            ifelse(grp == "Residual", "exp",
#                                   "strain"))) %>%
#     select(source, model, std = sdcor, cond) %>%
#     bind_rows(., Aryani_bio3, Aryani_bio4, Aryani_bio6) %>%
#     filter(source != "total") %>%
#     full_join(., dof_map) %>% 
#     mutate(variance = std^2) %>%
#     mutate(chi0p05 = qchisq(0.05, dof),
#            chi0p95 = qchisq(0.95, dof)) %>%
#     mutate(var_left = variance/chi0p95*dof,
#            var_right = variance/chi0p05*dof) %>%
#     mutate(source = factor(source, levels = c("exp", "bio", "strain"))) %>%
#     ggplot(aes(x = source, y = variance, colour = model)) +
#     geom_point() +
#     geom_errorbar(aes(ymin = var_left, ymax = var_right)) +
#     facet_wrap("cond") +
#     scale_y_sqrt(name = "Variance",
#                  sec.axis = sec_axis(~sqrt(.), name = "Standard deviation")) +
#     xlab("")
# 
# ## Relative standard deviations
# 
# mean_mus <- all_data %>%
#     group_by(cond) %>%
#     summarize(m_mu = mean(mu, na.rm = TRUE))
# 
# c(mixed_bio_3, mixed_bio_4, mixed_bio_6) %>%
#     map(summary) %>%
#     map(., ~ .$varcor) %>%
#     map(as.data.frame) %>%
#     imap_dfr(., ~ mutate(.x, cond = .y)) %>%
#     mutate(model = "mixed") %>%
#     mutate(source = ifelse(grp == "bio_rep", "bio",
#                            ifelse(grp == "Residual", "exp",
#                                   "strain"))) %>%
#     select(source, model, std = sdcor, cond) %>%
#     left_join(., mean_mus) %>%
#     mutate(rel_std = std/m_mu) %>%
#     ggplot() +
#         geom_col(aes(x = cond, y = rel_std, fill = source), position = "dodge") +
#         ylab("Relative standard deviation") + xlab("")
# 
# c(mixed_bio_3, mixed_bio_4, mixed_bio_6) %>%
#     map(summary) %>%
#     map(., ~ .$varcor) %>%
#     map(as.data.frame) %>%
#     imap_dfr(., ~ mutate(.x, cond = .y)) %>%
#     mutate(model = "mixed") %>%
#     mutate(source = ifelse(grp == "bio_rep", "bio",
#                            ifelse(grp == "Residual", "exp",
#                                   "strain"))) %>%
#     select(source, model, std = sdcor, cond) %>%
#     left_join(., mean_mus) %>%
#     mutate(rel_std = std/m_mu) %>%
#     ggplot() +
#     geom_col(aes(x = cond, y = std, fill = source), position = "dodge") +
#     xlab("") + ylab("Standard deviation")










