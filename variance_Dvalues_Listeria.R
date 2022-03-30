
## Load libraries

library(tidyverse)
library(readxl)
library(sofa)
library(lme4)

source("funcs_var_analysis.R")

## Data of variability (strain/rep) on L. monocytogenes (Figure 3.4)

aryani_data <- read_xlsx("./data/D-values.xlsx", skip = 2) %>%
    set_names(., paste(names(.), .[1,], .[2,], sep = "-")) %>%
    rename(., strain = "...1-NA-NA") %>%
    slice(., 3:n()) %>%
    gather(condition, D_val, -strain) %>%
    separate(condition, into = c("temp", "bio_rep", "tec_rep"), sep = "-") %>%
    separate(temp, into = c("temp", "foo"), sep = 2) %>%
    select(-foo) %>%
    mutate(temp = as.numeric(temp)) %>%
    mutate(bio_rep = LETTERS[as.numeric(bio_rep)]) %>%
    mutate(bio_rep = paste(strain, temp, bio_rep, sep = "-")) %>%
    mutate(tec_rep = paste(bio_rep, tec_rep, sep = "-")) %>%
    filter(!is.na(D_val))

## Some plot

aryani_data %>%
    ggplot(aes(x = temp, y = log10(D_val), colour = strain)) +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE, linetype = 2)

# ## Add to database
# 
# x <- Cushion$new(user = "test", pwd = "test")
# 
# db_name <- "test_inactivation3"
# # db_create(x, db_name)
# 
# aryani_data <- aryani_data %>%
#     mutate(Microorganism = "Listeria",
#            Species = "monocytogenes",
#            Ref = "Aryani",
#            Year = 2015,
#            response = "inactivation",
#            model_type = "primary",
#            temperature = "temp") %>%
#     rename(Strain = strain) 
# 
# aa <- apply(aryani_data, 1, function(y) {
#     doc_create(x, dbname = db_name, as.list(y))
# })
# 
# db_query(x, dbname = db_name,
#          selector = list(Ref = "Aryani"),
#          fields = c("Species", "Strain", "D_val"),
#          as = "json",
#          limit = 4000) %>%
#     # .$docs %>%
#     fromJSON() %>%
#     .$docs
# 
# db_query(x, dbname = db_name,
#          selector = list(`_id` = list(`$gt` = NULL)),
#          as = "json",
#          limit = 10000) %>%
#     fromJSON() %>%
#     .$docs %>%
#     group_by(Ref) %>%
#     summarize(n()) %>%
#     View()

## Analysis of variance

n_strain <- 20
n_bio <- 3
n_exp <- 2

dof_strain <- n_strain - 1
dof_bio <- n_strain*n_bio - n_strain
dof_exp <- n_strain*n_bio*n_exp - n_strain*n_bio

## Fit the models

all_fix <- aryani_data %>% 
    filter(temp %in% c(55, 60, 65)) %>%
    mutate(log_D = log10(D_val)) %>%
    rename(bio_rep_c = bio_rep) %>%
    split(.$temp) %>%
    map(group_variance) %>%
    map(., ~ extract_estimates(., dof_strain, dof_bio, dof_exp)) %>%
    # map(., ~ set_names(., paste0("var_", names(.)))) %>%
    imap_dfr(., ~ mutate(.x, temp = as.numeric(.y))) %>%
    mutate(model = "fixed") %>%
    gather(source, variance, -model, -temp)

all_mixed <- aryani_data %>% 
    filter(temp %in% c(55, 60, 65)) %>%
    mutate(log_D = log10(D_val)) %>%
    split(.$temp) %>%
    map(., ~ lmer(log_D ~ (1|strain)+ (1|bio_rep), data = ., REML = TRUE))

mixed_all_temps <- aryani_data %>%
    filter(temp %in% c(55, 60 ,65)) %>%
    mutate(log_D = log10(D_val)) %>%
    lmer(log_D ~ 1 + temp + (1|strain)+ (1|bio_rep), data = ., REML = TRUE)



## Compare them

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
    bind_rows(., all_fix)

mixed_all_temps %>%
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
    filter(source != "total") %>%
    mutate(source = factor(source, levels = c("exp", "bio", "strain"))) %>%
    ggplot() +
    geom_bar(aes(x = factor(temp), y = sqrt(variance), fill = model),
             position = "dodge", stat = "identity") + 
    facet_grid(~source) +
    xlab("Temperature (ºC)") + ylab("Attributed standard deviation")

## Add the uncertainties

dof_strain <- n_strain - 1
dof_bio <- n_strain*n_bio - n_strain
dof_exp <- n_strain*n_bio*n_exp - n_strain*n_bio

dof_map <- tibble(source = rep(c("exp", "bio", "strain"), 3),
                  model = rep(c("fixed", "mixed", "mixed_allT"), each = 3),
                  dof = c(rep(c(dof_exp, dof_bio, dof_strain), 2),
                          n_strain*n_bio*n_exp*3 - n_strain*n_bio*3, n_strain*n_bio*3 - n_strain*3, n_strain*3-1)
                  )

mixed_all_temps %>%
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
    mutate(source = factor(source, levels = c("exp", "bio", "strain"))) %>%
    ggplot(aes(x = factor(temp), y = variance, colour = model)) +
        geom_point() +
        geom_errorbar(aes(ymin = var_left, ymax = var_right)) +
        facet_wrap("source") +
        scale_y_sqrt(name = "Variance", 
                     sec.axis = sec_axis(~sqrt(.), name = "Standard deviation"))















