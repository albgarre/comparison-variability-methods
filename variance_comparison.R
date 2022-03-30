
library(tidyverse)

## Inactivation of L. plantarum

tibble(type_of = c("E", "R", "S"),
       RMSE = c(0.02, 0.05, 0.3),
       MSE = RMSE^2,
       dof = c(6*20 - 3*20, 3*20 - 1*20, 20 - 1),
       SE = MSE*dof
       ) %>%
    mutate(chi0p05 = qchisq(0.05, dof),
           chi0p95 = qchisq(0.95, dof)) %>%
    mutate(sigma2_left = SE/chi0p95,
           sigma2_right = SE/chi0p05) %>%
    ggplot(aes(x = type_of)) +
        geom_point(aes(y = MSE)) +
        geom_errorbar(aes(ymin = sigma2_left, ymax = sigma2_right)) +
        scale_y_sqrt() +
        ylab("variance") + xlab("")

## Growth of L. plantarum

# pH

tibble(type_of = c("E", "R", "S"),
       RMSE = c(0.04, 0.06, 0.07),
       MSE = RMSE^2,
       dof = c(6*20 - 3*20, 3*20 - 1*20, 20 - 1),
       SE = MSE*dof
       ) %>%
    mutate(chi0p05 = qchisq(0.05, dof),
           chi0p95 = qchisq(0.95, dof)) %>%
    mutate(sigma2_left = SE/chi0p95,
           sigma2_right = SE/chi0p05) %>%
    ggplot(aes(x = type_of)) +
    geom_point(aes(y = MSE)) +
    geom_errorbar(aes(ymin = sigma2_left, ymax = sigma2_right)) +
    scale_y_sqrt() +
    ylab("variance") + xlab("")

# tibble(dof = 10,
#        SE = seq(0, 5*dof, length = 1000),
#        var = 1,
#        x_chi = SE/var,
#        p = dchisq(x_chi, dof),
#        MSE = SE/dof
#        ) %>%
#     ggplot() +
#         geom_line(aes(MSE, p))

# aw

tibble(type_of = c("E", "R", "S"),
       RMSE = c(0.05, 0.08, 0.07),
       MSE = RMSE^2,
       dof = c(6*20 - 3*20, 3*20 - 1*20, 20 - 1),
       SE = MSE*dof
) %>%
    mutate(chi0p05 = qchisq(0.05, dof),
           chi0p95 = qchisq(0.95, dof)) %>%
    mutate(sigma2_left = SE/chi0p95,
           sigma2_right = SE/chi0p05) %>%
    ggplot(aes(x = type_of)) +
    geom_point(aes(y = MSE)) +
    geom_errorbar(aes(ymin = sigma2_left, ymax = sigma2_right)) +
    scale_y_sqrt() +
    ylab("variance") + xlab("")

# HLa

tibble(type_of = c("E", "R", "S"),
       RMSE = c(0.03, 0.05, 0.04),
       MSE = RMSE^2,
       dof = c(6*20 - 3*20, 3*20 - 1*20, 20 - 1),
       SE = MSE*dof
) %>%
    mutate(chi0p05 = qchisq(0.05, dof),
           chi0p95 = qchisq(0.95, dof)) %>%
    mutate(sigma2_left = SE/chi0p95,
           sigma2_right = SE/chi0p05) %>%
    ggplot(aes(x = type_of)) +
    geom_point(aes(y = MSE)) +
    geom_errorbar(aes(ymin = sigma2_left, ymax = sigma2_right)) +
    scale_y_sqrt() +
    ylab("variance") + xlab("")

# Temperature

tibble(type_of = c("E", "R", "S"),
       RMSE = c(0.03, 0.055, 0.05),
       MSE = RMSE^2,
       dof = c(6*20 - 3*20, 3*20 - 1*20, 20 - 1),
       SE = MSE*dof
) %>%
    mutate(chi0p05 = qchisq(0.05, dof),
           chi0p95 = qchisq(0.95, dof)) %>%
    mutate(sigma2_left = SE/chi0p95,
           sigma2_right = SE/chi0p05) %>%
    ggplot(aes(x = type_of)) +
    geom_point(aes(y = MSE)) +
    geom_errorbar(aes(ymin = sigma2_left, ymax = sigma2_right)) +
    scale_y_sqrt() +
    ylab("variance") + xlab("")

## Inactivation Listeria

# 55ºC

tibble(type_of = c("E", "R", "S"),
       RMSE = c(0.01, 0.04, 0.12),
       MSE = RMSE^2,
       dof = c(6*20 - 3*20, 3*20 - 1*20, 20 - 1),
       SE = MSE*dof
) %>%
    mutate(chi0p05 = qchisq(0.05, dof),
           chi0p95 = qchisq(0.95, dof)) %>%
    mutate(sigma2_left = SE/chi0p95,
           sigma2_right = SE/chi0p05) %>%
    ggplot(aes(x = type_of)) +
    geom_point(aes(y = MSE)) +
    geom_errorbar(aes(ymin = sigma2_left, ymax = sigma2_right)) +
    scale_y_sqrt() +
    ylab("variance") + xlab("")

# 60ºC

tibble(type_of = c("E", "R", "S"),
       RMSE = c(0.04, 0.055, 0.33),
       MSE = RMSE^2,
       dof = c(6*20 - 3*20, 3*20 - 1*20, 20 - 1),
       SE = MSE*dof
) %>%
    mutate(chi0p05 = qchisq(0.05, dof),
           chi0p95 = qchisq(0.95, dof)) %>%
    mutate(sigma2_left = SE/chi0p95,
           sigma2_right = SE/chi0p05) %>%
    ggplot(aes(x = type_of)) +
    geom_point(aes(y = MSE)) +
    geom_errorbar(aes(ymin = sigma2_left, ymax = sigma2_right)) +
    scale_y_sqrt() +
    ylab("variance") + xlab("")

## Growth Listeria

# A

tibble(type_of = c("E", "R", "S"),
       RMSE = c(0.035, 0.05, 0.045),
       MSE = RMSE^2,
       dof = c(6*20 - 3*20, 3*20 - 1*20, 20 - 1),
       SE = MSE*dof
) %>%
    mutate(chi0p05 = qchisq(0.05, dof),
           chi0p95 = qchisq(0.95, dof)) %>%
    mutate(sigma2_left = SE/chi0p95,
           sigma2_right = SE/chi0p05) %>%
    ggplot(aes(x = type_of)) +
    geom_point(aes(y = MSE)) +
    geom_errorbar(aes(ymin = sigma2_left, ymax = sigma2_right)) +
    scale_y_sqrt() +
    ylab("variance") + xlab("")

# B

tibble(type_of = c("E", "R", "S"),
       RMSE = c(0.03, 0.04, 0.042),
       MSE = RMSE^2,
       dof = c(6*20 - 3*20, 3*20 - 1*20, 20 - 1),
       SE = MSE*dof
) %>%
    mutate(chi0p05 = qchisq(0.05, dof),
           chi0p95 = qchisq(0.95, dof)) %>%
    mutate(sigma2_left = SE/chi0p95,
           sigma2_right = SE/chi0p05) %>%
    ggplot(aes(x = type_of)) +
    geom_point(aes(y = MSE)) +
    geom_errorbar(aes(ymin = sigma2_left, ymax = sigma2_right)) +
    scale_y_sqrt() +
    ylab("variance") + xlab("")

# C

tibble(type_of = c("E", "R", "S"),
       RMSE = c(0.025, 0.05, 0.055),
       MSE = RMSE^2,
       dof = c(6*20 - 3*20, 3*20 - 1*20, 20 - 1),
       SE = MSE*dof
) %>%
    mutate(chi0p05 = qchisq(0.05, dof),
           chi0p95 = qchisq(0.95, dof)) %>%
    mutate(sigma2_left = SE/chi0p95,
           sigma2_right = SE/chi0p05) %>%
    ggplot(aes(x = type_of)) +
    geom_point(aes(y = MSE)) +
    geom_errorbar(aes(ymin = sigma2_left, ymax = sigma2_right)) +
    scale_y_sqrt() +
    ylab("variance") + xlab("")

# D

tibble(type_of = c("E", "R", "S"),
       RMSE = c(0.015, 0.02, 0.018),
       MSE = RMSE^2,
       dof = c(6*20 - 3*20, 3*20 - 1*20, 20 - 1),
       SE = MSE*dof
) %>%
    mutate(chi0p05 = qchisq(0.05, dof),
           chi0p95 = qchisq(0.95, dof)) %>%
    mutate(sigma2_left = SE/chi0p95,
           sigma2_right = SE/chi0p05) %>%
    ggplot(aes(x = type_of)) +
    geom_point(aes(y = MSE)) +
    geom_errorbar(aes(ymin = sigma2_left, ymax = sigma2_right)) +
    scale_y_sqrt() +
    ylab("variance") + xlab("")




















