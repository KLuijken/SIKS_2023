# ---------------------------------------------------------------------------- #
# Causal research questions - Exercise 1
# SIKS course Causal Inference
#
# Author: Kim Luijken, k.luijken@umcutrecht.nl
# ---------------------------------------------------------------------------- #

# Load required libraries
library( ggplot2) # for plotting

# load data
vaccine_data <- readRDS( "vaccine_data.rds")

# describe data ( simplisticly)
summary( vaccine_data)

# ---------------------------------------------------------------------------- #
# Example analysis 1
# ---------------------------------------------------------------------------- #

# estimate inverse probability of treatment weights via a logistic model
# fit propensity score model
ps_1 <- glm( vacc ~ sex +
                     age + I(age ^ 2) +
                     cvd +
                     pulm +
                     dm,
  family = binomial(),
  data = vaccine_data
)

vaccine_data$propensity_score_1 <- predict( ps_1, type = "response")

# Create density plot
unweighted_plot <- ggplot( vaccine_data,
        aes( x = propensity_score_1,
             fill = factor( vacc))) +
  geom_density( alpha = 0.5) +
  labs(title = "Unweighted Sample",
       x = "Propensity Score",
       y = "Density",
       fill = "Vaccination") +
  scale_fill_manual( values = c( "#F8766D", "#00BFC4")) +
  theme_classic() +
  ylim( 0, 5.5)
print( unweighted_plot)

# estimate weights
vaccine_data$weight_1 <- ifelse( vaccine_data$vacc == 1,
                                 1 / predict( ps_1, type = "response"),
                                 1 / (1 - predict( ps_1, type = "response"))
                                 )
# Create density plot
weighted_plot_1 <- ggplot( vaccine_data,
        aes(x = propensity_score_1,
            fill = factor( vacc),
            weight = weight_1)) +
  geom_density(alpha = 0.5) +
  labs(title = "Weighted Sample",
       x = "Propensity Score",
       y = "Density",
       fill = "Vaccination") +
  scale_fill_manual(values = c("#F8766D", "#00BFC4")) +
  theme_classic() +
  ylim( 0, 5.5)
print( weighted_plot_1)

# Display the plots
gridExtra::grid.arrange( unweighted_plot, weighted_plot_1, nrow = 1)

# inspect weights
summary( vaccine_data$weight_1)

# estimate effect of vaccination on mortality
mod_1 <- glm( mort ~ vacc,
              weights = weight_1,
              family = binomial(),
              data = vaccine_data)
# you may ignore the warning message showing up when running this command

summary( mod_1)

# compute marginal risk difference
# data under vacc == 1
vacc_yes_data <- vaccine_data
vacc_yes_data$vacc <- 1

# mortality under vacc == 1
mort_vacc_yes_1 <- predict( mod_1,
                        newdata = vacc_yes_data,
                        type = "response")

# data under vacc == 0
vacc_no_data <- vaccine_data
vacc_no_data$vacc <- 0

# mortality under vacc == 0
mort_vacc_no_1 <- predict( mod_1,
                        newdata = vacc_no_data,
                        type = "response")

# marginal risk difference
rd_1 <- mean( mort_vacc_yes_1 - mort_vacc_no_1)
rd_1

# compute standard error using bootstrap
## note: running this bootstrap takes around 3 minutes, but is not necessary
## to complete the practical. You may choose a number of bootstrap repetitions
## lower than 1000 to decrease computation time.
b_rep <- 1000
effect_estimates_1 <- matrix( NA, nrow = b_rep, ncol = 1)

for( i in 1:b_rep){
  # take bootstrap sample
  bootstrap_sample <- vaccine_data[ sample( 1:nrow( vaccine_data),
                                            size = nrow( vaccine_data),
                                            replace = T), ]
  # fit propensity score model
  ps_bs_1 <- glm( vacc ~ sex +
                 age + I(age ^ 2) +
                 cvd +
                 pulm +
                 dm,
               family = binomial(),
               data = bootstrap_sample
  )

  # estimate weights
  bootstrap_sample$weight_1 <- ifelse( bootstrap_sample$vacc == 1,
                                   1 / predict( ps_bs_1, type = "response"),
                                   1 / (1 - predict( ps_bs_1, type = "response"))
  )

  # estimate effect of vaccination on mortality
  mod_bs_1 <- glm( mort ~ vacc,
                weights = weight_1,
                family = binomial(),
                data = bootstrap_sample)

  # compute marginal risk difference
  # data under vacc == 1
  vacc_yes_bs_data <- bootstrap_sample
  vacc_yes_bs_data$vacc <- 1

  # mortality under vacc == 1
  mort_vacc_yes_1 <- predict( mod_bs_1,
                              newdata = vacc_yes_bs_data,
                              type = "response")

  # data under vacc == 0
  vacc_no_bs_data <- bootstrap_sample
  vacc_no_bs_data$vacc <- 0

  # mortality under vacc == 0
  mort_vacc_no_1 <- predict( mod_bs_1,
                             newdata = vacc_no_bs_data,
                             type = "response")

  # marginal risk difference
  effect_estimates_1[ i,] <- mean( mort_vacc_yes_1 - mort_vacc_no_1)

}

# estimate SE
se_1 <- sqrt( sum( (effect_estimates_1 - mean( effect_estimates_1))^2) / (b_rep - 1))


# effect estimate vaccination on mortality
cat( rd_1,
        "95% CI,",
        (rd_1 - 1.96 * se_1),
        "to",
        (rd_1 + 1.96 * se_1))


# ---------------------------------------------------------------------------- #
# Example analysis 2
# ---------------------------------------------------------------------------- #

# estimate inverse probability of treatment weights via a logistic model
# fit propensity score model
ps_2 <- glm( vacc ~ sex +
                age + I(age ^ 2) +
                cvd +
                pulm +
                dm,
              family = binomial(),
              data = vaccine_data
)

vaccine_data$propensity_score_2 <- predict( ps_2, type = "response")

# Create density plot
unweighted_plot <- ggplot( vaccine_data,
                           aes( x = propensity_score_2,
                                fill = factor( vacc))) +
  geom_density( alpha = 0.5) +
  labs(title = "Unweighted Sample",
       x = "Propensity Score",
       y = "Density",
       fill = "Vaccination") +
  scale_fill_manual( values = c( "#F8766D", "#00BFC4")) +
  theme_classic() +
  ylim( 0, 5.5)
print( unweighted_plot)

# estimate weights
vaccine_data$weight_2 <- ifelse( vaccine_data$vacc == 1,
                                 1,
                                 predict( ps_2, type = "response") / (1 - predict( ps_2, type = "response"))
                                 )

# Create density plot
weighted_plot_2 <- ggplot( vaccine_data,
                         aes(x = propensity_score_2,
                             fill = factor( vacc),
                             weight = weight_2)) +
  geom_density(alpha = 0.5) +
  labs(title = "Weighted Sample",
       x = "Propensity Score",
       y = "Density",
       fill = "Vaccination") +
  scale_fill_manual(values = c("#F8766D", "#00BFC4")) +
  theme_classic() +
  ylim( 0, 5.5)
print( weighted_plot_2)

# Display the plots
gridExtra::grid.arrange( unweighted_plot, weighted_plot_2, nrow = 1)

# estimate effect of vaccination on mortality
mod_2 <- glm( mort ~ vacc,
              weights = weight_2,
              family = binomial(),
              data = vaccine_data)

summary( mod_2)


# compute marginal risk difference
# data under vacc == 1
vacc_yes_data <- vaccine_data
vacc_yes_data$vacc <- 1

# mortality under vacc == 1
mort_vacc_yes_2 <- predict( mod_2,
                          newdata = vacc_yes_data,
                          type = "response")

# data under vacc == 0
vacc_no_data <- vaccine_data
vacc_no_data$vacc <- 0

# mortality under vacc == 0
mort_vacc_no_2 <- predict( mod_2,
                         newdata = vacc_no_data,
                         type = "response")

# marginal risk difference
rd_2 <- mean( mort_vacc_yes_2 - mort_vacc_no_2)
rd_2

# compute standard error using bootstrap
## note: running this bootstrap takes around 3 minutes, but is not necessary
## to complete the practical. You may choose a number of bootstrap repetitions
## lower than 1000 to decrease computation time.
b_rep <- 1000
effect_estimates_2 <- matrix( NA, nrow = b_rep, ncol = 1)

for( i in 1:b_rep){
  # take bootstrap sample
  bootstrap_sample <- vaccine_data[ sample( 1:nrow( vaccine_data),
                                            size = nrow( vaccine_data),
                                            replace = T), ]
  # fit propensity score model
  ps_bs_2 <- glm( vacc ~ sex +
                    age + I(age ^ 2) +
                    cvd +
                    pulm +
                    dm,
                  family = binomial(),
                  data = bootstrap_sample
  )

  # estimate weights
  bootstrap_sample$weight_2 <- ifelse( bootstrap_sample$vacc == 1,
                                       1,
                                       predict( ps_bs_2, type = "response") /
                                         (1 - predict( ps_bs_2, type = "response"))
                                       )

  # estimate effect of vaccination on mortality
  mod_bs_2 <- glm( mort ~ vacc,
                   weights = weight_2,
                   family = binomial(),
                   data = bootstrap_sample)

  # compute marginal risk difference
  # data under vacc == 1
  vacc_yes_bs_data <- bootstrap_sample
  vacc_yes_bs_data$vacc <- 1

  # mortality under vacc == 1
  mort_vacc_yes_2 <- predict( mod_bs_2,
                              newdata = vacc_yes_bs_data,
                              type = "response")

  # data under vacc == 0
  vacc_no_bs_data <- bootstrap_sample
  vacc_no_bs_data$vacc <- 0

  # mortality under vacc == 0
  mort_vacc_no_2 <- predict( mod_bs_2,
                             newdata = vacc_no_bs_data,
                             type = "response")

  # marginal risk difference
  effect_estimates_2[ i,] <- mean( mort_vacc_yes_2 - mort_vacc_no_2)
}

# estimate SE
se_2 <- sqrt( sum( (effect_estimates_2 - mean( effect_estimates_2))^2) / (b_rep - 1))


# effect estimate vaccination on mortality
cat( rd_2,
     "95% CI,",
     (rd_2 - 1.96 * se_2),
     "to",
     (rd_2 + 1.96 * se_2))



