
eps_e=1/1000
eps_b=1/1000
eps_s=1/1000

## Load libraries

library(tidyverse)
library(readxl)
library(sofa)
library(lme4)
library(rjags)
library(ggplot2)
library(coda)

source("funcs_var_analysis_JS.R")

## Data of variability (strain/rep) on L. monocytogenes (Figure 3.4)
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

## Some plot

# aryani_data %>%
#     ggplot(aes(x = temp, y = log10(D_val), colour = strain)) +
#         geom_point() +
#         geom_smooth(method = "lm", se = FALSE, linetype = 2)

## Analysis of variance

n_strain <- 20
n_bio <- 3
n_exp <- 2

dof_strain <- n_strain - 1
dof_bio <- n_strain*n_bio - n_strain
dof_exp <- n_strain*n_bio*n_exp - n_strain*n_bio

## Fit the models

all_fix <- aryani_data %>% 
    filter(temp %in% c(55, 60, 65)) %>% # find rows where temp is 55/60/65
    mutate(log_D = log10(D_val)) %>% # 
    rename(bio_rep_c = bio_rep) %>%
    split(.$temp) %>% # splits tibble in 3 tibble for different temperatures
    map(group_variance) %>% # use functions as defined by Aryani
    map(., ~ extract_estimates(., dof_strain, dof_bio, dof_exp)) %>% # extract estimates
    # map(., ~ set_names(., paste0("var_", names(.)))) %>%
    imap_dfr(., ~ mutate(.x, temp = as.numeric(.y))) %>% # Apply a function to each element of a vector, and its index
    mutate(model = "fixed") %>%
    gather(source, variance, -model, -temp)

all_mixed <- aryani_data %>% 
    filter(temp %in% c(55, 60, 65)) %>%
    mutate(log_D = log10(D_val)) %>%
    split(.$temp) %>%
    map(., ~ lmer(log_D ~ (1|strain)+ (1|bio_rep), data = ., REML = TRUE)) # DOES THIS ACCOUNT FOR NESTING?? https://stats.stackexchange.com/questions/228800/crossed-vs-nested-random-effects-how-do-they-differ-and-how-are-they-specified
# only crossed effects??? YES< this is okay because the nesting is explicit:
# xtabs(~strain+bio_rep, aryani_data)

mixed_all_temps <- aryani_data %>%
    filter(temp %in% c(55, 60 ,65)) %>%
    mutate(log_D = log10(D_val)) %>%
    lmer(log_D ~ 1 + temp + (1|strain)+ (1|bio_rep), data = ., REML = TRUE)

all_mcmc <- aryani_data %>% 
  filter(temp %in% c(55, 60, 65)) %>%
  mutate(log_D = log10(D_val)) %>%
  split(.$temp) %>%
  map(group_variance_mcmc_Dvalues_Listeria, nstrain=n_strain, nbio=n_bio, 
      nexp=n_exp, eps_e=eps_e, eps_b=eps_b, eps_s=eps_s, niter = 1e4)


all_mcmc %>% map(plot)
all_mcmc <- all_mcmc %>% map(summary)

all_mcmc_exp <- aryani_data %>% 
  filter(temp %in% c(55, 60, 65)) %>%
  mutate(log_D = log10(D_val)) %>%
  split(.$temp) %>%
  map(group_variance_mcmc_Dvalues_Listeria_exponential, nstrain=n_strain, nbio=n_bio, 
      nexp=n_exp, eps_e=1, eps_b=1, eps_s=1, niter = 1e4)

all_mcmc_exp %>% map(plot)
all_mcmc_exp <- all_mcmc_exp %>% map(summary)

# all_mcmc_all_temps <- aryani_data %>% 
#   filter(temp %in% c(55, 60, 65)) %>%
#   mutate(log_D = log10(D_val)) %>%
#   map(group_variance_mcmc_Dvalues_Listeria) # TO DO: WE SHOULD DEFINE ANOTHER FUNCTION. HOW DOES TEMPERATURE AFFECT VARIABILITY?


## Compare them
my_vars_mcmc <- lapply(all_mcmc, function(l) l[[1]])
my_vars_mcmc_quant <- lapply(all_mcmc, function(l) l[[2]])

my_vars_mcmc <- my_vars_mcmc %>% 
  map(as.data.frame) %>%
  imap_dfr(., ~ slice(.x, 2:n())) %>%
  add_column(source=rep(c("bio", "exp", "strain"),3), .before = 1) %>%
  add_column(model=rep("mcmc",9), .before = 1) %>%
  add_column(temp=rep(c(55, 60, 65),each=3), .before = 1) %>% 
  rename(variance = Mean) %>%
  select(temp, model, source, variance)

my_vars_mcmc_quant <- my_vars_mcmc_quant %>% 
  map(as.data.frame) %>%
  imap_dfr(., ~ slice(.x, 2:n())) %>%
  add_column(source=rep(c("bio", "exp", "strain"),3), .before = 1) %>%
  add_column(model=rep("mcmc",9), .before = 1) %>%
  add_column(temp=rep(c(55, 60, 65),each=3), .before = 1)

my_vars <- all_mixed %>%
    map(summary)  %>%
    map(., ~ .$varcor) %>%
    map(as.data.frame) %>%
    imap_dfr(., ~ mutate(.x, temp = as.numeric(.y))) %>%
    mutate(model = "mixed") %>%
    mutate(source = ifelse(grp == "bio_rep", "bio",
                           ifelse(grp == "Residual", "exp",
                                  "strain"))) %>%
    select(temp, model, source, variance = vcov) %>%
    bind_rows(., all_fix) %>%
   bind_rows(., my_vars_mcmc)

# pdf(file=paste0("output/sdplot_DvalListeria_prior",eps_e,"_",eps_b,"_",eps_s,".pdf"))
# mixed_all_temps %>%
#     summary() %>%
#     .$varcor %>%
#     as.data.frame() %>%
#     mutate(temp = NA) %>%
#     mutate(model = "mixed_allT") %>%
#     mutate(source = ifelse(grp == "bio_rep", "bio",ifelse(grp == "Residual", "exp","strain"))) %>%
#     select(temp, model, source, variance = vcov) %>%
#     bind_rows(., my_vars) %>%
#     filter(source != "total") %>%
#     mutate(source = factor(source, levels = c("exp", "bio", "strain"))) %>%
#     ggplot() +
#     geom_bar(aes(x = factor(temp), y = sqrt(variance), fill = model),
#              position = "dodge", stat = "identity") +
#     facet_grid(~source) +
#     xlab("Temperature (ºC)") + ylab("Attributed standard deviation")
# dev.off()

## Add the uncertainties

dof_strain <- n_strain - 1
dof_bio <- n_strain*n_bio - n_strain
dof_exp <- n_strain*n_bio*n_exp - n_strain*n_bio

dof_map <- tibble(source = rep(c("exp", "bio", "strain"), 3),
                  model = rep(c("fixed", "mixed", "mixed_allT"), each = 3),
                  dof = c(rep(c(dof_exp, dof_bio, dof_strain), 2),
                          n_strain*n_bio*n_exp*3 - n_strain*n_bio*3, n_strain*n_bio*3 - n_strain*3, n_strain*3-1)
                  )

mixed_all_temps_unc <- mixed_all_temps %>%
    summary() %>%
    .$varcor %>%
    as.data.frame() %>%
    mutate(temp = NA) %>%
    mutate(model = "mixed_allT") %>%
    mutate(source = ifelse(grp == "bio_rep", "bio",
                           ifelse(grp == "Residual", "exp",
                                  "strain"))) %>%
    select(temp, model, source, variance = vcov) %>%
    bind_rows(., my_vars) %>%
    # filter(!is.na(temp)) %>%  
    filter(source != "total") %>%
    full_join(dof_map) %>%
    mutate(chi0p05 = qchisq(0.05, dof),
           chi0p95 = qchisq(0.95, dof)) %>%
    mutate(var_left = variance/chi0p95*dof,
           var_right = variance/chi0p05*dof) %>%
    mutate(source = factor(source, levels = c("exp", "bio", "strain")))


mixed_all_temps_unc[mixed_all_temps_unc["model"]=="mcmc","var_left"] <- my_vars_mcmc_quant[,"2.5%"]
mixed_all_temps_unc[mixed_all_temps_unc["model"]=="mcmc","var_right"] <- my_vars_mcmc_quant[,"97.5%"]

# pdf(file=paste0("output/varplotunc_Dval_prior",eps_e,"_",eps_b,"_",eps_s,".pdf"))
mixed_all_temps_unc %>%
    ggplot(aes(x = factor(temp), y = variance, colour = model)) +
        geom_point() +
        geom_errorbar(aes(ymin = var_left, ymax = var_right)) +
        facet_wrap("source") +
        scale_y_sqrt(name = "Variance", 
                     sec.axis = sec_axis(~sqrt(.), name = "Standard deviation"))
# dev.off()













