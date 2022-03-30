
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
