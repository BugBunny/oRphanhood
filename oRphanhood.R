
#oRphanhood

# This script calculates the proportions of people with living mothers or
# fathers by age in a set of stable population models (SPMs) defined by
# 2-parameter relational logit model life tables, 2-parameter relational
# Gompertz fertility models, and their growth rates, r.

# These populations are then used to model and predict:

#  i.    the relationship between life table survivorship (lm+n/lm, where m is
#        25 for women and 35 for men) and the proportion of the population by
#        five-year age group with a living father or mother (Sn) controlling
#        for M, the unstandardised mean age at childbearing. These regression
#        models are used for estimating adult mortality.

#   ii.  the relationship in single-year age groups between parental survival
#        and summary indices of adult survivorship (lm+n/lm) given the timing
#        of childbearing (M). These models are useful when estimating
#        fertility using the Own-Children Method as the estimate of orphans
#        can be removed from the count of children living apart from their
#        mother, bringing the numerators and denominators of the fertility
#        rates into correspondence.

# Options are available to replicate the Timaeus (1992) coefficients and to
# calculate multipliers for converting orphanhood to survivorship and vice
# versa in a specific population.

# Ian M. Timaeus (LSHTM)
# 1 June, 2021

# Parameters ------------------------------------------------------------------
sex <- 2 # Mothers = 1; Fathers = 2
# If one2one is TRUE, the scripts will calculate a single set of 'tailor-made'
# multipliers for converting proportions surviving to life table survivorship
# in a population with the specified characteristics
one2one <- FALSE
# If replicating is TRUE, the analysis uses the set of stable populations that
# Timaeus (1992) adopted to produce coefficients for estimating adult mortality
# by the orphanhood method. (The results duplicate those that were first
# calculated using a spreadsheet, Supercalc, in 1989)
replicating <- FALSE
if (!replicating) {
   # Initialise the input parameters for the simulated populations
   alpha_mortality <- c(-1.4, -1.0, -0.6, -0.2, 0.2)
   beta_mortality <- c(0.7, 1.0, 1.1, 1.3)
   alpha_fertility_f <- c(-0.5, -0.2, 0.1, 0.4)
   alpha_fertility_m <- c(-0.4, -0.1, 0.2, 0.6)
   beta_fertility <- c(0.75, 0.95, 1.15)
   growth_rate <- c(0.025, 0.005)
} else {
   # Do not alter these parameters by mistake
   alpha_mortality <- c(-1.0, -0.6, -0.2, 0.2)
   beta_mortality <- c(0.8, 1.1)
   alpha_fertility_f <- c(-0.5, -0.2, 0.1, 0.4)
   alpha_fertility_m <- c(-0.4, -0.1, 0.2, 0.6)
   beta_fertility <- c(0.7, 1.0, 1.3)
   growth_rate <- c(0.02)
}

# Calculate orphanhood by age in a specified SPM ==============================
SimulateOrphans <- function(a_m, b_m, a_f, b_f, growth_rate, sex) {

   # Logits of Brass's General Standard model life table
   Ysx <- c(NA, -0.8670, -0.7152, -0.6552, -0.6219, -0.6015, -0.5879, -0.5766,
            -0.5666, -0.5578, -0.5498, -0.5431, -0.5365, -0.5296, -0.5220,
            -0.5131, -0.5043, -0.4941, -0.4824, -0.4694, -0.4551, -0.4401,
            -0.4248, -0.4103, -0.3963, -0.3829, -0.3686, -0.3549, -0.3413,
            -0.3280, -0.3150, -0.3020, -0.2889, -0.2759, -0.2627, -0.2496,
            -0.2364, -0.2230, -0.2094, -0.1956, -0.1816, -0.1674, -0.1530,
            -0.1381, -0.1229, -0.1073, -0.0911, -0.0745, -0.0574, -0.0396,
            -0.0212, -0.0021, 0.0177, 0.0383, 0.0598, 0.0821, 0.1055, 0.1299,
            0.1554, 0.1821, 0.2100, 0.2394, 0.2701, 0.3024, 0.3364, 0.3721,
            0.4097, 0.4494, 0.4912, 0.5353, 0.5818, 0.6311, 0.6832, 0.7385,
            0.7971, 0.8591, 0.9253, 0.9962, 1.0714, 1.1513, 1.2377, 1.3298,
            1.4287, 1.5346, 1.6496, 1.7717, 1.9043, 2.0501, 2.2054, 2.3737,
            2.5550, 2.7587, 2.9748, 3.2181, 3.4534, 3.7090, 4.0557, 4.2585,
            4.6051, 5.1270, 5.5550)
   # Gompits of Booth's standard relational Gompertz fertility model for women
   Zsx_f <- c(NA, -3.17091, -2.74255, -2.36854, -2.04079, -1.75210, -1.49286,
              -1.25061, -1.04479, -0.85927, -0.69130, -0.53325, -0.38524,
              -0.24423, -0.10783, 0.02564, 0.15853, 0.29147, 0.42515, 0.56101,
              0.70000, 0.84272, 0.99014, 1.14407, 1.30627, 1.47872, 1.66426,
              1.86597, 2.08894, 2.33992, 2.62602, 2.95500, 3.32873, 3.75984,
              4.25499, 4.80970, 5.41311, 6.12864, 7.07022, 8.64839, NA)
   # Gompits of Paget's standard relational Gompertz fertility model for men
   Zsx_m <- c(NA, -3.94388, -3.60472, -3.29098, -3.00073, -2.73207, -2.48300,
              -2.25180, -2.03701, -1.83621, -1.64711, -1.46475, -1.28997,
              -1.13178, -0.98688, -0.84932, -0.71966, -0.59719, -0.47931,
              -0.36513, -0.25533, -0.14925, -0.04603, 0.05488, 0.15404,
              0.25174, 0.34846, 0.44494, 0.54133, 0.63775, 0.73432, 0.83123,
              0.92890, 1.02762, 1.12760, 1.22900, 1.33191, 1.43658, 1.54324,
              1.65239, 1.76445, 1.87982, 1.99896, 2.12229, 2.25053, 2.38436,
              2.52442, 2.67168, 2.82734, 2.99295, 3.17016, 3.36064, 3.56611,
              3.78876, 4.02974, 4.28837, 4.56749, 4.87088, 5.19998, 5.55620,
              5.93800, 6.34357, 6.76448, 7.21465, 7.70720, 8.26447, 8.91782,
              9.78653, 10.96787, 12.55913, NA)

   # Calculate life table survivorship
   lx <- 1 / (1 + exp(2*(a_m + b_m * Ysx)))
   lx[1] <- 1
   Lx = (Ysx[1:(length(Ysx)-1)] + Ysx[2:length(Ysx)]) / 2
   Lx <- 1 / (1 + exp(2*(a_m + b_m * Lx)))
   Lx[1] = 0.25 + 0.75 * lx[2]

   # Calculate the age-specific fertility rates for either men or women
   if (sex == 2) {
      Zsx <- Zsx_m
      base_age = 35
      # Survivorship of fathers from conception to the birth of their children
      last <- 11 + (length(Zsx) - 1) - 1
      prepartum_survivors <- 2 * Lx[11:last] / (lx[11:last] + Lx[10:(last-1)])
   } else {
      Zsx <- Zsx_f
      base_age = 25
      prepartum_survivors <- replicate(length(Zsx) - 1, 1)
   }
   Fx <- exp(-exp(-a_f - b_f * Zsx))
   Fx[1] <- 0
   last <- length(Fx)
   Fx[last] <- 1
   fx <- Fx[2:last] - Fx[1:(last-1)]
   last <- last - 1

   # Future survivorship of women/men having children at each age 10->49/79
   Lxy <- t(Lx[11:length(Lx)])
   for (i in 2:last) Lxy <- rbind(Lxy, c(Lxy[nrow(Lxy), 2:ncol(Lxy)], 0))
   # Survivorship (in cols) of women/men who had children at age x (in rows)
   age <- seq(10, last+9) + 0.5
   Sxy <- exp(-growth_rate * age) * Lxy * fx
   # Women's mean age of childbearing
   births <- sum(Sxy[, 1])
   mbar <- sum(Sxy[, 1] * age) / births
   # Calculate what proportion of fathers survive the 9 months before a birth
   Sxy <- Sxy * prepartum_survivors
   # Proportions of all mothers/fathers surviving by exact age of the child
   Sy <- colSums(Sxy)
   # Age-weighted parental survival by single- and 5-year age group of child
   age_distrib <- exp(-growth_rate * seq(0, length(Sy)-1)) * lx[1:length(Sy)]
   Sy <- age_distrib * Sy / births
   S5y <- replicate(9, NA)
   for (i in 1:9) S5y[i] <- (0.5 * Sy[i*5+1]
                              + sum(Sy[(i*5+2):(i*5+5)])
                              + 0.5 * Sy[i*5+6]) /
                            (0.5 * age_distrib[i*5+1]
                              + sum(age_distrib[(i*5+2):(i*5+5)])
                              + 0.5 * age_distrib[i*5+6])
   Sy <- (Sy + c(Sy[2:length(Sy)], 0)) /
         (age_distrib + c(age_distrib[2:length(Sy)], 0))

   # Life table survivorship from exact age m to age group m+y
   Py <- replicate(length(Sy), 0)
   P5y <- replicate(9, NA)
   Py[1:(length(Lx)-base_age)] <- Lx[base_age+1:length(Lx)] / lx[base_age+1]
   for (i in 1:9) P5y[i] <- lx[base_age + i*5 + 6] / lx[base_age+1]
   # Return everything to the main program as columns of a single matrix
   c(mbar, Sy, Py, S5y, P5y)
}

# Main program ================================================================
# Set up a data frame with every combination of the input parameters ----------
if (sex == 2) {
   alpha_fertility <- alpha_fertility_m
} else {
   alpha_fertility <-  alpha_fertility_f
}
in_puts <- expand.grid(r = growth_rate,
                       b_f = beta_fertility,
                       a_f = alpha_fertility,
                       b_m = beta_mortality,
                       a_m = alpha_mortality)
# Drop early, yet wide, and late, yet narrow, male fertility distributions
if (sex==2 & replicating) {
   in_puts <- subset(in_puts, (a_f<0.6 | b_f!=0.7) & (a_f>-0.4 | b_f!=1.3))
}
# Calculate orphanhood for each combination of the input parameters -----------
out_puts <- mapply(SimulateOrphans,
                    in_puts$a_m,
                    in_puts$b_m,
                    in_puts$a_f,
                    in_puts$b_f,
                    in_puts$r, sex)
out_puts <- data.frame(t(out_puts))
last <- (length(out_puts) - 1 - 2 * 9) / 2
names(out_puts)[1] <- 'mbar'
names(out_puts)[2:(last+1)] <- paste0('S', seq(0, last-1))
names(out_puts)[(last+2):(2*last+1)] <- paste0(seq(0, last-1), 'Pm')
names(out_puts)[(2*last+2):(2*last+10)] <- paste0('5S', seq(5, 45, 5))
names(out_puts)[(2*last+11):(2*last+19)] <- paste0(seq(10, 50, 5), 'Pm-exact')

# Model relationship between living mothers or fathers and survivorship -------
o_coeffs <- matrix(nrow = last, ncol = ifelse(one2one, 1, 3))
for (i in 1:last) {
   if (one2one) {
      models <- lm(out_puts[, i+1] ~ 0 + out_puts[, i + 1 + last])
   } else {
      models <- lm(out_puts[, i+1] ~ out_puts[, 1] + out_puts[, i + 1 + last])
   }
   # Strip out and save the coefficients of the fitted models
   o_coeffs[i, ] <- models$coefficients
}
o_coeffs <- data.frame(o_coeffs) # o_coeffs are for estimating orphanhood
names(o_coeffs)[1] <- ifelse(one2one, 'Multiplier', 'Intercept')
names(o_coeffs)[2] <- 'M_bar'
names(o_coeffs)[3] <- 'Pn'
# Model relationship between survivorship and living mothers or fathers -------
first <- last * 2 + 2
# Men's suriviorship is modelled using data from two adjoining age groups
m_coeffs <- matrix(nrow = 9, ncol =  ifelse(one2one, 1, 2 + sex))
for (i in first:(first+9-sex)) {
   if (one2one) {
      models <- lm(out_puts[, i+9] ~ 0 + out_puts[, i])
   } else {
      if (sex == 2) {
         models <- lm(out_puts[, i+9] ~ out_puts[, 1]
                                       + out_puts[, i]
                                       + out_puts[, i+1])
      } else {
         models <- lm(out_puts[, i+9] ~ out_puts[, 1] + out_puts[, i])
      }
   }
   # Strip out and save the coefficients of the fitted models
   m_coeffs[i-first+1, ] <- models$coefficients
}
m_coeffs <- data.frame(m_coeffs) # m_coeffs are for estimating adult mortality
names(m_coeffs)[1] <- ifelse(one2one, 'Multiplier', 'Intercept')
names(m_coeffs)[2] <- 'M_bar'
names(m_coeffs)[3] <- '5Sn-5'
if (sex == 2) names(m_coeffs)[4] <- '5Sn'