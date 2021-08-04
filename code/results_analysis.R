# Project:   compareHD_MVN
# Objective: Script to plot the results of the simulation
# Author:    Edoardo Costantini
# Created:   2021-07-30
# Modified:  2021-08-04

# barplots ----------------------------------------------------------------

# Define plot input
col_cond <- 1:6 # fixed columns that describe conditions
estimand <- unique(gg_shape$par)[1] # what parameter to plot
y_cond <- unique(gg_shape$pm)
x_cond <- unique(gg_shape$p)

# Bias
met_cond <- paste0(unique(gg_shape$method)[-9], # which methods to plot
                   collapse = "|")
col_plot <- which(colnames(out_ggready)
                    %in%
                    c("est", "ref"))

out_ggready %>%
  # Subset
  select(colnames(out_ggready)[c(col_cond, col_plot)]) %>%
  filter(grepl(estimand, par)) %>%
  filter(grepl(met_cond, method)) %>%
  # Compute Bias
  group_by(p, pm, par, method) %>%
  summarize(EQ_bar = mean(est),
            ref = mean(ref)) %>%
  mutate(RB = EQ_bar - ref,
         PB = 100 * abs((EQ_bar - ref)/ref)) %>%
  # Plot
  ggplot(aes(x = method,
             y = PB)) +
  geom_bar(stat = "identity") +
  facet_grid(rows = vars(factor(p,
                                labels = paste0("pm = ", y_cond))),
             cols = vars(factor(pm,
                                labels = paste0("p = ", x_cond))))

# Confidence Interval Coverage
met_cond <- paste0(unique(gg_shape$method), # which methods to plot
                   collapse = "|")
col_plot <- which(colnames(out_ggready)
                    %in%
                    "CIC")

out_ggready %>%
  # Subset
  select(colnames(out_ggready)[c(col_cond, col_plot)]) %>%
  filter(grepl(estimand, par)) %>%
  filter(grepl(met_cond, method)) %>%
  # Compute Bias
  group_by(p, pm, par, method) %>%
  summarize(CIC = mean(CIC)) %>%
  # Plot
  ggplot(aes(x = method,
             y = CIC)) +
  geom_bar(stat = "identity") +
  facet_grid(rows = vars(factor(p,
                                labels = paste0("pm = ", y_cond))),
             cols = vars(factor(pm,
                                labels = paste0("p = ", x_cond))))

# boxplots ----------------------------------------------------------------

# Define plot input
col_cond <- 1:6 # fixed columns that describe conditions
estimand <- unique(gg_shape$par)[1] # what parameter to plot
y_cond <- unique(gg_shape$pm)
x_cond <- unique(gg_shape$p)

# Bias
met_cond <- paste0(unique(gg_shape$method), # which methods to plot
                   collapse = "|")
col_plot <- which(colnames(out_ggready)
                    %in%
                    c("est", "ref"))

out_ggready %>%
  # Subset
  select(colnames(out_ggready)[c(col_cond, col_plot)]) %>%
  filter(grepl(estimand, par)) %>%
  filter(grepl(met_cond, method)) %>%
  # Plot
  ggplot(aes(x = method,
             y = est)) +
  geom_boxplot() +
  facet_grid(rows = vars(factor(p,
                                labels = paste0("pm = ", y_cond))),
             cols = vars(factor(pm,
                                labels = paste0("p = ", x_cond))))