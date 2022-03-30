
## Load libraries

library(tidyverse)
library(GGally)
library(lme4)
library(cowplot)

source("funcs_var_analysis.R")

## Make the Monte Carlo simulations

n_strain <- 20
n_bio <- 3
n_exp <- 2

var_strain <- 1
var_bio <- 1
var_exp <- 1

mean_logD <- 10

dof_strain <- n_strain - 1
dof_bio <- n_strain*n_bio - n_strain
dof_exp <- n_strain*n_bio*n_exp - n_strain*n_bio

n_iter <- 1000

set.seed(19271)

my_simuls <- c(1:n_iter) %>%
  map(., ~ simulate_spread(n_strain, n_bio, n_exp,
                           var_strain, var_bio, var_exp,
                           mean_logD = mean_logD) 
      ) 

## Fit Diah's method

aryanis_pars <- my_simuls %>%
  map(., group_variance) %>%
  map_dfr(., ~ extract_estimates(., dof_strain, dof_bio, dof_exp)) %>%
  mutate(var_sum = strain + bio + exp,
         unexp_var = total - var_sum,
         method = "Aryani")

ggplot(aryanis_pars) +
  geom_density(aes(unexp_var/total)) 

aryanis_pars %>%
  select(strain, bio, exp, total) %>%
  gather(type, MSS) %>%
  ggplot() +
    geom_boxplot(aes(x = type, y = MSS)) +
    xlab("")
    # geom_density(aes(MSS, colour = type))

## Fit the mixed effects model

mixed_pars <- my_simuls %>%
  map(., ~ lmer(log_D ~ (1|strain)+ (1|strain:bio_rep), data = .,
                REML = TRUE)) %>%
  map(., summary) %>%
  map(., ~ .$varcor) %>%
  map(as_tibble) %>%
  map_dfr(., ~ mutate(., var_source = c("bio", "strain", "exp"))) %>%
  mutate(method = "mixed")


ggplot(mixed_pars) +
  geom_boxplot(aes(x = var_source, y = vcov)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "red") +
  ylab("Variance estimate") + xlab("")

## Compare the estimates

# p1 <- aryanis_pars %>%
#   select(strain, bio, exp, total, method) %>%
#   gather(type, MSS, -method) %>%
#   bind_rows(., select(mixed_pars, method, type = var_source, MSS = vcov)) %>%
#   filter(type != "total") %>%
#   mutate(type = ifelse(type == "bio", "Within-strain",
#                        ifelse(type == "exp", "Experimental",
#                               "Between-strain"))) %>%
#   mutate(type = factor(type, levels = c("Between-strain", "Within-strain", "Experimental"))) %>%
#   ggplot() +
#     geom_boxplot(aes(x = type, y = MSS, colour = method), notch = FALSE) +
#     geom_hline(yintercept = 1, linetype = 2, colour = "black") +
#   ylab(expression(Variance~estimate~(log~min^2))) + xlab("") +
#   xlab("") +
#     theme_bw() +
#     ylim(0, 4) +
#     theme(legend.position = "none",
#           axis.text = element_text(size = 14),
#           axis.title = element_text(size = 16)) 


p1 <- aryanis_pars %>%
  select(strain, bio, exp, total, method) %>%
  gather(type, MSS, -method) %>%
  bind_rows(., select(mixed_pars, method, type = var_source, MSS = vcov)) %>%
  filter(type != "total") %>%
  mutate(type = ifelse(type == "bio", "Within-strain",
                       ifelse(type == "exp", "Experimental",
                              "Between-strain"))) %>%
  mutate(type = factor(type, levels = c("Between-strain", "Within-strain", "Experimental"))) %>%
  ggplot() +
  geom_boxplot(aes(x = type, y = MSS, colour = method), notch = FALSE) +
  # geom_hline(yintercept = 1, linetype = 2, colour = "black") +
  geom_point(aes(x = x, y = y),
             inherit.aes = FALSE,
             shape = 10, 
             size = 8,
             data = tibble(
               x = c("Between-strain", "Within-strain", "Experimental"),
               y = c(1, 1, 1)
             )) +
  ylab(expression(Variance~estimate~(log~min^2))) + xlab("") +
  xlab("") +
  theme_bw() +
  ylim(0, 4) +
  theme(legend.position = "none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16)) 



## Make the Monte Carlo simulations (number 2)

n_strain <- 20
n_bio <- 3
n_exp <- 2

var_strain <- 1
var_bio <- .1
var_exp <- .1

mean_logD <- 10

dof_strain <- n_strain - 1
dof_bio <- n_strain*n_bio - n_strain
dof_exp <- n_strain*n_bio*n_exp - n_strain*n_bio

n_iter <- 1000

set.seed(19271)

my_simuls2 <- c(1:n_iter) %>%
  map(., ~ simulate_spread(n_strain, n_bio, n_exp,
                           var_strain, var_bio, var_exp,
                           mean_logD = mean_logD) 
  ) 

## Fit Diah's method

aryanis_pars2 <- my_simuls2 %>%
  map(., group_variance) %>%
  map_dfr(., ~ extract_estimates(., dof_strain, dof_bio, dof_exp)) %>%
  mutate(var_sum = strain + bio + exp,
         unexp_var = total - var_sum,
         method = "Aryani")

## Fit the mixed effects model

mixed_pars2 <- my_simuls2 %>%
  map(., ~ lmer(log_D ~ (1|strain)+ (1|strain:bio_rep), data = .,
                REML = TRUE)) %>%
  map(., summary) %>%
  map(., ~ .$varcor) %>%
  map(as_tibble) %>%
  map_dfr(., ~ mutate(., var_source = c("bio", "strain", "exp"))) %>%
  mutate(method = "mixed")

## Compare the estimates

p2 <- aryanis_pars2 %>%
  select(strain, bio, exp, total, method) %>%
  gather(type, MSS, -method) %>%
  bind_rows(., select(mixed_pars2, method, type = var_source, MSS = vcov)) %>%
  filter(type != "total") %>%
  mutate(type = ifelse(type == "bio", "Within-strain",
                       ifelse(type == "exp", "Experimental",
                              "Between-strain"))) %>%
  mutate(type = factor(type, levels = c("Between-strain", "Within-strain", "Experimental"))) %>%
  ggplot() +
  geom_boxplot(aes(x = type, y = MSS, colour = method)) +
  geom_point(aes(x = x, y = y),
             inherit.aes = FALSE,
             shape = 10, 
             size = 8,
             data = tibble(
               x = c("Between-strain", "Within-strain", "Experimental"),
               y = c(1, .1, .1)
             )) +
  ylab(expression(Variance~estimate~(log~min^2))) + xlab("") +
  theme_bw() +
  ylim(0, 4) + 
  theme(legend.position = "none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))

## Make the Monte Carlo simulations (number 3)

n_strain <- 20
n_bio <- 10
n_exp <- 10

var_strain <- 1
var_bio <- 1
var_exp <- 1

mean_logD <- 10

dof_strain <- n_strain - 1
dof_bio <- n_strain*n_bio - n_strain
dof_exp <- n_strain*n_bio*n_exp - n_strain*n_bio

n_iter <- 1000

set.seed(19271)

my_simuls3 <- c(1:n_iter) %>%
  map(., ~ simulate_spread(n_strain, n_bio, n_exp,
                           var_strain, var_bio, var_exp,
                           mean_logD = mean_logD) 
  ) 

## Fit Diah's method

aryanis_pars3 <- my_simuls3 %>%
  map(., group_variance) %>%
  map_dfr(., ~ extract_estimates(., dof_strain, dof_bio, dof_exp)) %>%
  mutate(var_sum = strain + bio + exp,
         unexp_var = total - var_sum,
         method = "Aryani")

## Fit the mixed effects model

mixed_pars3 <- my_simuls3 %>%
  map(., ~ lmer(log_D ~ (1|strain)+ (1|strain:bio_rep), data = .,
                REML = TRUE)) %>%
  map(., summary) %>%
  map(., ~ .$varcor) %>%
  map(as_tibble) %>%
  map_dfr(., ~ mutate(., var_source = c("bio", "strain", "exp"))) %>%
  mutate(method = "mixed")

## Compare the estimates

p3 <- aryanis_pars3 %>%
  select(strain, bio, exp, total, method) %>%
  gather(type, MSS, -method) %>%
  bind_rows(., select(mixed_pars3, method, type = var_source, MSS = vcov)) %>%
  filter(type != "total") %>%
  mutate(type = ifelse(type == "bio", "Within-strain",
                       ifelse(type == "exp", "Experimental",
                              "Between-strain"))) %>%
  mutate(type = factor(type, levels = c("Between-strain", "Within-strain", "Experimental"))) %>%
  ggplot() +
  geom_boxplot(aes(x = type, y = MSS, colour = method), notch = FALSE) +
  geom_point(aes(x = x, y = y),
             inherit.aes = FALSE,
             shape = 10, 
             size = 8,
             data = tibble(
               x = c("Between-strain", "Within-strain", "Experimental"),
               y = c(1, 1, 1)
             )) +
  ylab(expression(Variance~estimate~(log~min^2))) + xlab("") +
  theme_bw() +
  ylim(0, 4) + 
  theme(legend.position = "none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))

## Figure 1

plot_grid(p1, p2, p3, labels = "AUTO", nrow = 1)

## Means and sd

aryanis_pars %>%
  select(strain, bio, exp, total, method) %>%
  gather(type, MSS, -method) %>%
  bind_rows(., select(mixed_pars, method, type = var_source, MSS = vcov)) %>%
  filter(type != "total") %>%
  mutate(type = ifelse(type == "bio", "Within-strain",
                       ifelse(type == "exp", "Experimental",
                              "Between-strain"))) %>%
  mutate(type = factor(type, levels = c("Between-strain", "Within-strain", "Experimental"))) %>%
  group_by(method, type) %>%
  summarize(mean(MSS), sd(MSS))

aryanis_pars2 %>%
  select(strain, bio, exp, total, method) %>%
  gather(type, MSS, -method) %>%
  bind_rows(., select(mixed_pars2, method, type = var_source, MSS = vcov)) %>%
  filter(type != "total") %>%
  mutate(type = ifelse(type == "bio", "Within-strain",
                       ifelse(type == "exp", "Experimental",
                              "Between-strain"))) %>%
  mutate(type = factor(type, levels = c("Between-strain", "Within-strain", "Experimental"))) %>%
  group_by(method, type) %>%
  summarize(mean(MSS), sd(MSS))

aryanis_pars3 %>%
  select(strain, bio, exp, total, method) %>%
  gather(type, MSS, -method) %>%
  bind_rows(., select(mixed_pars3, method, type = var_source, MSS = vcov)) %>%
  filter(type != "total") %>%
  mutate(type = ifelse(type == "bio", "Within-strain",
                       ifelse(type == "exp", "Experimental",
                              "Between-strain"))) %>%
  mutate(type = factor(type, levels = c("Between-strain", "Within-strain", "Experimental"))) %>%
  group_by(method, type) %>%
  summarize(mean(MSS), sd(MSS))









