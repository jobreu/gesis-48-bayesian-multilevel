library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(brms)

mathach <- read_csv('data/mathach.csv')
sleepstudy <- read_csv('data/sleepstudy.csv') %>% 
  mutate(Subject = as.factor(Subject))
fake_rt <- read_csv('data/fake_rt.csv')
science_df <- read_csv('data/science.csv')
penicillin_df <- read_csv('data/penicillin.csv')
insulation_df <- read_csv('data/insulation.csv')

## ------------------------------------------------------------------------
ggplot(fake_rt,
       aes(x = reaction_time, fill = subject)
) + geom_histogram(col='white', binwidth = 10) +
  facet_wrap(~ subject) + guides(fill=F)


## ------------------------------------------------------------------------
ggplot(sleepstudy,
       aes(x=Days, y=Reaction, col=Subject)
) + geom_point() +
  stat_smooth(method='lm', se=F) +
  facet_wrap(~Subject) +
  guides(col=F)


## ------------------------------------------------------------------------
schools <- c('1308', '1946', '7276', '9198', 
             '2917', '9104', '4042', '7635', 
             '8188', '6464', '5838', '9158', 
             '5815', '6808', '6578', '4350')

mathach %>% 
  filter(school %in% schools) %>% 
  mutate(school = as.character(school)) %>% 
  ggplot(aes(x=ses, y=mathach, col=school)
) + geom_point() +
  stat_smooth(method='lm', se=F) +
  facet_wrap(~school) +
  guides(col=F)



# ----- fake_rt -----------------------------

M <- brm(reaction_time ~ 1 + (1|subject),
         cores = 2, 
         #prior = set_prior('normal(0, 100)'), 
         save_all_pars = T, # only necessary if we use the BayesFactor package for model comparison
         data = fake_rt)

prior_summary(M)
# 1st parameter of student_t = df (t-distribution with df of 1 = Cauchy)
# other parts of default student_t priors come from the data
fake_rt %>% pull(reaction_time) %>% median() # robust measure of central tendency
fake_rt %>% pull(reaction_time) %>% mad() # robust measure of spread

summary(M)
# Rhat should be close to 1 and not >= 1.1
# Pop-level Intercept = mean of means
# Group level sd = mean variation between subjects
# family specific sd: within-subject

# ------- mathach --------------------
M <- brm(mathach ~ ses + (ses|school),
         cores = 2, 
         prior = set_prior('normal(0, 100)'), 
         save_all_pars = T,
         data = mathach)

prior_summary(M)

# -------- sleep ---------------------
# random intercepts and slopes: Reaction ~ 1 + Days + (1 + Days|Subject) = # Reaction ~ Days + (Days|Subject)
# random intercepts: Reaction ~ Days + (1|Subject)

M <- brm(Reaction ~ Days + (Days|Subject),
         cores = 2, 
         #prior = set_prior('normal(0, 100)'), 
         save_all_pars = T,
         data = sleepstudy)

prior_summary(M)
# blank field in prior column (first for category!) = uniform distribution
# uniform distribution = improper prior for b
# family specific parameters: sigma = sd of residuals across all individuals
# group-level effects: sd(Days) -> inter-individual variability in the slopes
# group-level effects: sd(Intercept) -> inter-individual variability in the intercepts
# cor(Intercept, Days): relationship between intercept and slope

summary(M)
# pop level effects (fixed effects) = avg. individual

plot(M)

M <- brm(Reaction ~ Days + (Days|Subject),
         cores = 2, 
         prior = set_prior('normal(0, 100)'), # specify prior for slope term
         save_all_pars = T,
         data = sleepstudy)

prior_summary(M)
summary(M)
plot(M)

# Random intercepts model
# M_lmer_ri <- lmer(Reaction ~ Days + (1|Subject),
#                   data = Df)

M_ri <- brm(Reaction ~ Days + (1|Subject),
            cores = 2,               
            prior = set_prior('normal(0, 100)'), # flat prior on coefs
            save_all_pars = T,
            data = sleepstudy)

# Random intercepts and random slopes model
# M_lmer <- lmer(Reaction ~ Days + (Days|Subject),
#                data = Df)

M <- brm(Reaction ~ Days + (Days|Subject),
         cores = 2,               
         #prior = set_prior('normal(0, 100)'), # flat prior on coefs
         save_all_pars = T,
         data = sleepstudy)


# Model comparison
waic(M_ri, M)
loo(M_ri, M)
bayes_factor(M_ri, M)




# Nested models -----------------------------------------------------------

glimpse(science_df)

# Each class is within a school. 
group_by(science_df, school, class) %>% 
  summarise(n = n()) %>% 
  xtabs(n ~ school + class, data=.)


# To model this nesting, we'd do the following:
M <- brm(like ~ sex + PrivPub + (1|school/class), # random intercepts for schools & classes
           cores = 2,               
           prior = set_prior('normal(0, 100)'),  
           save_all_pars = T,
           data = science_df)

summary(M)
# ~school: inter-school variability (in this case after controlling for school type & student gender)
# ~school:class inter-class variability

# which is identical to
M <- brm(like ~ sex + PrivPub + (1|school) + (1|school:class), 
           cores = 2,               
           prior = set_prior('normal(0, 100)'),  
           save_all_pars = T,
           data = science_df)

# however, "Class" is unambiguous 

M <- brm(like ~ sex + PrivPub + (1|school) + (1|Class), 
         cores = 2,               
         prior = set_prior('normal(0, 100)'),  
         save_all_pars = T,
         data = science_df)

# here's a hierarchical model with random slopes

M <- brm(like ~ sex + PrivPub + (sex|school) + (sex|Class), 
         cores = 2,               
         prior = set_prior('normal(0, 100)'),  
         save_all_pars = T,
         data = science_df)


# Crossed designs ---------------------------------------------------------

# Crossed designs
group_by(penicillin_df, plate, sample) %>% 
  summarize(n=n()) %>% 
  xtabs(n ~ plate + sample, data=.)

M <- brm(diameter ~ 1 + (1|plate) + (1|sample), # random intercepts for plates & samples (no predictor besides grouping variables)
         cores = 2,               
         prior = set_prior('normal(0, 100)', class='Intercept'), # note: no predictor var (b) 
         save_all_pars = T,
         data = penicillin_df)

summary(M)

# varying intercept and slopes, varying by two groups
ggplot(insulation_df,
       mapping = aes(x = Temp, y = Gas, col = Insul)
) + geom_point() +
  stat_smooth(method = 'lm', se = F)


# Varying slopes by insul (non-multilevel)
# interaction between continuous (Temp) & binary variable (Insul: after = 0, before = 1)
M_vsvi <- brm(Gas ~ Temp*Insul, # equivalent to Gas ~ Temp + Insul + Temp:Insul 
              data = insulation_df,
              cores = 2, 
              prior = set_prior('normal(0, 100)'), 
              save_all_pars = T 
)

summary(M_vsvi)
# after: Temp 0 -> gas 4.73, slope = -0.28
# before: Intercept = 4.73 + 2.12 = 6.85; slope = -0.28 - 0.11 = -0.39

# Multilevel slopes by insul
M_vsvi_ml <- brm(Gas ~ Temp + (Temp|Insul), 
                 data = insulation_df,
                 cores = 2, 
                 prior = set_prior('normal(0, 100)'), 
                 save_all_pars = T 
)

summary(M_vsvi_ml)


