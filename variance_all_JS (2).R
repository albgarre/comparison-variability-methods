eps_e=0.001
eps_b=0.001
eps_s=0.001

## Load libraries

library(tidyverse)
library(readxl)
library(lme4)
library(rjags)
library(ggplot2)
library(coda)

source("funcs_var_analysis_JS.R")

## Load the data

data_LA <- c("0mM", "3mM", "4mM") %>%
    set_names(., .) %>%
    map(., ~ read_excel("Result_combined_LA_var_fitting.xlsx", sheet = .)) %>%
    imap_dfr(., ~ mutate(.x, cond = .y)) %>%
    mutate(bio_rep = paste(strain, cond, day, sep = "-"),
           tec_rep = paste(strain, cond, day, rep, sep = "-"),
           bio_rep_c = bio_rep,
           tec_rep_c = tec_rep,
           log_D = mu
    )

data_aw <- c("aw_10", "aw_9", "aw_7p5", "aw_5", "aw_2p5", "aw_0") %>%
    set_names(., .) %>%
    map(., ~ read_excel("result combined_Aw_variability.xlsx", sheet = .)) %>%
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
    map(., ~ read_excel("result combined_pHvariability.xlsx", sheet = .)) %>%
    imap_dfr(., ~ mutate(.x, cond = .y)) %>%
    mutate(bio_rep = paste(strain, cond, day, sep = "-"),
           tec_rep = paste(strain, cond, day, rep, sep = "-"),
           bio_rep_c = bio_rep,
           tec_rep_c = tec_rep,
           log_D = mu
    )

data_T <- c("T5", "T10", "T20", "T30pH", "T30aw") %>%
    set_names(., .) %>%
    map(., ~ read_excel("result combined_T_var.xlsx", sheet = .)) %>%
    imap_dfr(., ~ mutate(.x, cond = .y)) %>%
    mutate(bio_rep = paste(strain, cond, day, sep = "-"),
           tec_rep = paste(strain, cond, day, rep, sep = "-"),
           bio_rep_c = bio_rep,
           tec_rep_c = tec_rep,
           log_D = mu
    )

all_data <- bind_rows(data_LA, data_aw, data_pH, data_T)

## Analyze the variances

all_data %>%
    group_by(cond) %>%
    summarize(max(day))

## Variances when bio == 3 (als er op max 3 dagen gemeten is)

n_strain <- 20
n_bio <- 3
n_exp <- 2

dof_strain <- n_strain - 1
dof_bio <- n_strain*n_bio - n_strain
dof_exp <- n_strain*n_bio*n_exp - n_strain*n_bio

fixed_bio3 <- all_data %>%
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
    mutate(model = "fixed")

mixed_bio_3 <- all_data %>%
    group_by(cond) %>%
    mutate(n_bio = max(day)) %>%
    filter(n_bio == 3) %>%
    split(.$cond) %>%
    map(., ~ lmer(mu ~ (1|strain)+ (1|bio_rep), data = ., REML = TRUE))

mcmc_bio_3 <- all_data %>% 
  group_by(cond) %>%
  mutate(n_bio = max(day)) %>%
  filter(n_bio == 3) %>%
  split(.$cond) %>%
  map(group_variance_mcmc, nstrain=n_strain, nbio=n_bio, nexp=n_exp, eps_e=eps_e, eps_b=eps_b, eps_s=eps_s)

## Variances when bio == 4

n_strain <- 20
n_bio <- 4
n_exp <- 2

dof_strain <- n_strain - 1
dof_bio <- n_strain*n_bio - n_strain
dof_exp <- n_strain*n_bio*n_exp - n_strain*n_bio

fixed_bio4 <- all_data %>%
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
    mutate(model = "fixed")

mixed_bio_4 <- all_data %>%
    group_by(cond) %>%
    mutate(n_bio = max(day)) %>%
    filter(n_bio == 4) %>%
    split(.$cond) %>%
    map(., ~ lmer(mu ~ (1|strain)+ (1|bio_rep), data = ., REML = TRUE))

mcmc_bio_4 <- all_data %>% 
  group_by(cond) %>%
  mutate(n_bio = max(day)) %>%
  filter(n_bio == 4) %>%
  split(.$cond) %>%
  map(group_variance_mcmc, nstrain=n_strain, nbio=n_bio, nexp=n_exp, eps_e=eps_e, eps_b=eps_b, eps_s=eps_s)

## Variances when bio == 6

n_strain <- 20
n_bio <- 6
n_exp <- 2

dof_strain <- n_strain - 1
dof_bio <- n_strain*n_bio - n_strain
dof_exp <- n_strain*n_bio*n_exp - n_strain*n_bio

fixed_bio6 <- all_data %>%
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
    mutate(model = "fixed")

mixed_bio_6 <- all_data %>%
    group_by(cond) %>%
    mutate(n_bio = max(day)) %>%
    filter(n_bio == 6) %>%
    split(.$cond) %>%
    map(., ~ lmer(mu ~ (1|strain)+ (1|bio_rep), data = ., REML = TRUE))

mcmc_bio_6 <- all_data %>% 
  group_by(cond) %>%
  mutate(n_bio = max(day)) %>%
  filter(n_bio == 6) %>%
  split(.$cond) %>%
  map(group_variance_mcmc, nstrain=n_strain, nbio=n_bio, nexp=n_exp, eps_e=eps_e, eps_b=eps_b, eps_s=eps_s)


## Compare them
my_vars_mcmc_bio_3 <- lapply(mcmc_bio_3, function(l) l[[1]])
my_vars_mcmc_bio_3_quant <- lapply(mcmc_bio_3, function(l) l[[2]])
my_vars_mcmc_bio_4 <- lapply(mcmc_bio_4, function(l) l[[1]])
my_vars_mcmc_bio_4_quant <- lapply(mcmc_bio_4, function(l) l[[2]])
my_vars_mcmc_bio_6 <- lapply(mcmc_bio_6, function(l) l[[1]])
my_vars_mcmc_bio_6_quant <- lapply(mcmc_bio_6, function(l) l[[2]])

my_vars_mcmc_bio_3 <- my_vars_mcmc_bio_3 %>% 
  map(as.data.frame) %>%
  imap_dfr(., ~ slice(.x, 2:n())) %>%
  mutate(std = sqrt(Mean)) %>%
  add_column(source=rep(c("bio", "exp", "strain"),length(my_vars_mcmc_bio_3)), .before = 1) %>%
  add_column(model=rep("mcmc",3*length(my_vars_mcmc_bio_3)), .before = 1) %>%
  add_column(cond=rep(names(my_vars_mcmc_bio_3),each=3), .before = 1) %>%
  select(source, model, std, cond)

my_vars_mcmc_bio_3_quant <- my_vars_mcmc_bio_3_quant %>% 
  map(as.data.frame) %>%
  imap_dfr(., ~ slice(.x, 2:n())) %>%
  add_column(source=rep(c("bio", "exp", "strain"),length(my_vars_mcmc_bio_3_quant)), .before = 1) %>%
  add_column(model=rep("mcmc",3*length(my_vars_mcmc_bio_3_quant)), .before = 1) %>%
  add_column(cond=rep(names(my_vars_mcmc_bio_3_quant),each=3), .before = 1)

my_vars_mcmc_bio_4 <- my_vars_mcmc_bio_4 %>% 
  map(as.data.frame) %>%
  imap_dfr(., ~ slice(.x, 2:n())) %>%
  mutate(std = sqrt(Mean)) %>%
  add_column(source=rep(c("bio", "exp", "strain"),length(my_vars_mcmc_bio_4)), .before = 1) %>%
  add_column(model=rep("mcmc",3*length(my_vars_mcmc_bio_4)), .before = 1) %>%
  add_column(cond=rep(names(my_vars_mcmc_bio_4),each=3), .before = 1) %>%
  select(source, model, std, cond)

my_vars_mcmc_bio_4_quant <- my_vars_mcmc_bio_4_quant %>% 
  map(as.data.frame) %>%
  imap_dfr(., ~ slice(.x, 2:n())) %>%
  add_column(source=rep(c("bio", "exp", "strain"),length(my_vars_mcmc_bio_4_quant)), .before = 1) %>%
  add_column(model=rep("mcmc",3*length(my_vars_mcmc_bio_4_quant)), .before = 1) %>%
  add_column(cond=rep(names(my_vars_mcmc_bio_4_quant),each=3), .before = 1)

my_vars_mcmc_bio_6 <- my_vars_mcmc_bio_6 %>% 
  map(as.data.frame) %>%
  imap_dfr(., ~ slice(.x, 2:n())) %>%
  mutate(std = sqrt(Mean)) %>%
  add_column(source=rep(c("bio", "exp", "strain"),length(my_vars_mcmc_bio_6)), .before = 1) %>%
  add_column(model=rep("mcmc",3*length(my_vars_mcmc_bio_6)), .before = 1) %>%
  add_column(cond=rep(names(my_vars_mcmc_bio_6),each=3), .before = 1) %>%
  select(source, model, std, cond)

my_vars_mcmc_bio_6_quant <- my_vars_mcmc_bio_6_quant %>% 
  map(as.data.frame) %>%
  imap_dfr(., ~ slice(.x, 2:n())) %>%
  add_column(source=rep(c("bio", "exp", "strain"),length(my_vars_mcmc_bio_6_quant)), .before = 1) %>%
  add_column(model=rep("mcmc",3*length(my_vars_mcmc_bio_6_quant)), .before = 1) %>%
  add_column(cond=rep(names(my_vars_mcmc_bio_6_quant),each=3), .before = 1)

my_vars_mcmc_quant <- my_vars_mcmc_quant %>% 
  map(as.data.frame) %>%
  imap_dfr(., ~ slice(.x, 2:n())) %>%
  add_column(source=rep(c("bio", "exp", "strain"),3), .before = 1) %>%
  add_column(model=rep("mcmc",9), .before = 1) %>%
  add_column(temp=rep(c(55, 60, 65),each=3), .before = 1)


## Compare the results
# by cond
# pdf(file=paste0("output/sdplot_all_cond_prior",eps_e,"_",eps_b,"_",eps_s,".pdf"))
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
#   bind_rows(., my_vars_mcmc_bio_3, my_vars_mcmc_bio_4, my_vars_mcmc_bio_6) %>%
#     bind_rows(., fixed_bio3, fixed_bio4, fixed_bio6) %>%
#     filter(source != "total") %>%
#     mutate(source = factor(source, levels = c("exp", "bio", "strain"))) %>%
#     ggplot() + 
#     geom_bar(aes(x = cond, y = std, fill = model), stat = "identity", position = "dodge") +
#     xlab("") + ylab("Estimate of the contribution to the standard deviation of mu") +
#     facet_wrap("source")
# dev.off()

# by source
# pdf(file=paste0("output/sdplot_all_source_prior",eps_e,"_",eps_b,"_",eps_s,".pdf"))
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
#   bind_rows(., my_vars_mcmc_bio_3, my_vars_mcmc_bio_4, my_vars_mcmc_bio_6) %>%
#     bind_rows(., fixed_bio3, fixed_bio4, fixed_bio6) %>%
#     filter(source != "total") %>%
#     mutate(source = factor(source, levels = c("exp", "bio", "strain"))) %>%
#     ggplot() + 
#     geom_bar(aes(x = source, y = std, fill = model), stat = "identity", position = "dodge") +
#     xlab("") + ylab("Estimate of the contribution to the standard deviation of mu") +
#     facet_wrap("cond")
# dev.off()

## Add the variances

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

# pdf(file=paste0("output/varplot_all_cont_prior",eps_e,"_",eps_b,"_",eps_s,".pdf"))
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
#     bind_rows(., fixed_bio3, fixed_bio4, fixed_bio6) %>%
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
# dev.off()

df_alldata_unc <- c(mixed_bio_3, mixed_bio_4, mixed_bio_6) %>%
    map(summary) %>%
    map(., ~ .$varcor) %>%
    map(as.data.frame) %>%
    imap_dfr(., ~ mutate(.x, cond = .y)) %>%
    mutate(model = "mixed") %>%
    mutate(source = ifelse(grp == "bio_rep", "bio",
                           ifelse(grp == "Residual", "exp",
                                  "strain"))) %>%
    select(source, model, std = sdcor, cond) %>%
    bind_rows(., fixed_bio3, fixed_bio4, fixed_bio6) %>%
  bind_rows(., my_vars_mcmc_bio_3, my_vars_mcmc_bio_4, my_vars_mcmc_bio_6) %>%
    filter(source != "total") %>%
    full_join(., dof_map) %>% 
    mutate(variance = std^2) %>%
    mutate(chi0p05 = qchisq(0.05, dof),
           chi0p95 = qchisq(0.95, dof)) %>%
    mutate(var_left = variance/chi0p95*dof,
           var_right = variance/chi0p05*dof) %>%
    mutate(source = factor(source, levels = c("exp", "bio", "strain")))

df_alldata_unc[df_alldata_unc["model"]=="mcmc","var_left"] <- c(my_vars_mcmc_bio_3_quant[,"2.5%"],
                                                                my_vars_mcmc_bio_4_quant[,"2.5%"],
                                                                my_vars_mcmc_bio_6_quant[,"2.5%"])
df_alldata_unc[df_alldata_unc["model"]=="mcmc","var_right"] <- c(my_vars_mcmc_bio_3_quant[,"97.5%"],
                                                                 my_vars_mcmc_bio_4_quant[,"97.5%"],
                                                                 my_vars_mcmc_bio_6_quant[,"97.5%"])

pdf(file=paste0("output/varplotunc_all_prior",eps_e,"_",eps_b,"_",eps_s,".pdf"))
df_alldata_unc %>%
    ggplot(aes(x = source, y = variance, colour = model)) +
    geom_point() +
    geom_errorbar(aes(ymin = var_left, ymax = var_right)) +
    facet_wrap("cond") +
    scale_y_sqrt(name = "Variance",
                 sec.axis = sec_axis(~sqrt(.), name = "Standard deviation")) +
    xlab("")
dev.off()

## Relative standard deviations

mean_mus <- all_data %>%
    group_by(cond) %>%
    summarize(m_mu = mean(mu, na.rm = TRUE))

c(mixed_bio_3, mixed_bio_4, mixed_bio_6) %>%
    map(summary) %>%
    map(., ~ .$varcor) %>%
    map(as.data.frame) %>%
    imap_dfr(., ~ mutate(.x, cond = .y)) %>%
    mutate(model = "mixed") %>%
    mutate(source = ifelse(grp == "bio_rep", "bio",
                           ifelse(grp == "Residual", "exp",
                                  "strain"))) %>%
    select(source, model, std = sdcor, cond) %>%
    left_join(., mean_mus) %>%
    mutate(rel_std = std/m_mu) %>%
    ggplot() +
        geom_col(aes(x = cond, y = rel_std, fill = source), position = "dodge") +
        ylab("Relative standard deviation") + ggtitle("Mixed model") 


bind_rows(my_vars_mcmc_bio_3, my_vars_mcmc_bio_4, my_vars_mcmc_bio_6) %>%
  mutate(model = "mixed") %>%
  left_join(., mean_mus) %>%
  mutate(rel_std = std/m_mu) %>%
  ggplot() +
  geom_col(aes(x = cond, y = rel_std, fill = source), position = "dodge") +
  ylab("Relative standard deviation") + xlab("") + ggtitle("MCMC model") 



c(mixed_bio_3, mixed_bio_4, mixed_bio_6) %>%
    map(summary) %>%
    map(., ~ .$varcor) %>%
    map(as.data.frame) %>%
    imap_dfr(., ~ mutate(.x, cond = .y)) %>%
    mutate(model = "mixed") %>%
    mutate(source = ifelse(grp == "bio_rep", "bio",
                           ifelse(grp == "Residual", "exp",
                                  "strain"))) %>%
    select(source, model, std = sdcor, cond) %>%
    left_join(., mean_mus) %>%
    mutate(rel_std = std/m_mu) %>%
    ggplot() +
    geom_col(aes(x = cond, y = std, fill = source), position = "dodge") +
    xlab("") + ylab("Standard deviation") + ggtitle("Mixed model") 

bind_rows(my_vars_mcmc_bio_3, my_vars_mcmc_bio_4, my_vars_mcmc_bio_6) %>%
  left_join(., mean_mus) %>%
  mutate(rel_std = std/m_mu) %>%
  ggplot() +
  geom_col(aes(x = cond, y = std, fill = source), position = "dodge") +
  xlab("") + ylab("Standard deviation") + ggtitle("MCMC model") 

## Look at T5 in detail

all_data %>%
    filter(cond == "T5") %>%
    ggplot() +
        geom_point(aes(x = strain, y = mu))