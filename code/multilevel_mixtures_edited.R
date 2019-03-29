library(readr)
library(dplyr)
library(ggplot2)
library(brms)

# ===================================================
faithful <- read_csv('data/faithful.csv')
set.seed(101010)

# generate synthetic data
mixreg_df <- tibble(
  x = rnorm(150),
  y = c(x[1:100], 5 + 2*x[101:150]) + rnorm(150)
)

# ===================================================

# Plot the histogram of eruptions of old faithful
ggplot(faithful,
       aes(x = eruptions)
) + geom_histogram(col='white', bins=30) + theme_classic()

# Plot the scatterplot of the mixture of linear regressions
ggplot(mixreg_df,
       aes(x = x, y = y)
) + geom_point()

# ========  mixture of normals ===============================
# regression models with no predictors (intercept only)
M <- brm(eruptions ~ 1,
         data = faithful,
         family = mixture(gaussian, nmix = 2), # equal to mixture(gaussian, gaussian) 
         chains = 4,
         cores = 2)

M
# pop level effects: inferred means of the 2 normal distributions
# family specific parameters: sigmas = SD of the two distributions
# eff sample size = estimate based on autocorrelation -> can be higher than # of samples

# normal linear mixtures --------------------------------------------------

mixreg_prior <- c(
  prior(normal(0, 5), class=Intercept, dpar = mu1),
  prior(normal(0, 5), class=Intercept, dpar = mu2),
  prior(dirichlet(2, 2), class=theta)
)
M <- brm(y ~ x, 
         data=Df, 
         family = mixture(gaussian, nmix=2),
         prior = mixreg_prior, 
         chains = 2, 
         inits = 0)

summary(M)
# pop-level: intercepts and slopes of two components
# see mix-reg data generation above

# -------- zip -------------------------------------------
# zero-inflated poisson
owls <- read_csv('data/owls.csv')

ggplot(owls,
       aes(x = SiblingNegotiation, fill = SexParent)
) + geom_histogram(position = 'dodge', binwidth = 1) + facet_wrap(~ FoodTreatment, nrow = 2)
# SiblingNegotiation = begging for food

Mzip <- brm(SiblingNegotiation ~ FoodTreatment + 
              SexParent + offset(log(BroodSize)) + (1|Nest),
            data = owls, 
            cores = 2, 
            prior = set_prior('normal(0, 100)'), 
            family = zero_inflated_poisson())

Mzip

plot(Mzip, pars='^b_')
# family specific parameters: zi = probability of value coming from zero component

# non-zero-inflated Poisson model
M_p <- brm(SiblingNegotiation ~ FoodTreatment + 
              SexParent + offset(log(BroodSize)) + (1|Nest),
            data = owls, 
            cores = 2, 
            prior = set_prior('normal(0, 100)'), 
            family = poisson())

M_p

# begging bevarior varies by hungry or not
# offset only relevant to count, not to which distribution data comes from
# two formulas in brm command -> wrap in bf function
Mzip2 <- brm(bf(SiblingNegotiation ~ FoodTreatment + 
                  SexParent + offset(log(BroodSize)) + (1|Nest),
                zi ~ FoodTreatment + SexParent + (1|Nest)),
             data = owls, 
             cores = 2, 
             prior = set_prior('normal(0, 100)'), 
             family = zero_inflated_poisson())

Mzip2

plot(Mzip2, pars='^b_')

# negative binomial
Mzinb2 <- brm(bf(SiblingNegotiation ~ FoodTreatment + 
                  SexParent + offset(log(BroodSize)) + (1|Nest),
                zi ~ FoodTreatment + SexParent + (1|Nest)),
             data = owls, 
             cores = 2, 
             prior = set_prior('normal(0, 100)'), 
             family = zero_inflated_negbinomial())

Mzinb2
# shape = variance of negative binomial distribution

plot(Mzinb2, pars='^b_')

Mzip2_waic <- waic(Mzip2)
Mzinb2_waic <- waic(Mzinb2)
compare_ic(Mzinb2_waic, Mzip2_waic)
# neg binomial model = better fit than Poisson model