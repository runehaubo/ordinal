#################################
## test profile and confint methods:
library(ordinal)
data(wine)
fm1 <- clm(rating ~ contact + temp, data = wine)
summary(fm1)

## profile.clm and confint.clm:
pr1 <- profile(fm1)
confint(pr1)
pr1 <- profile(fm1, which.beta = 1:2)
confint(pr1)
pr1 <- profile(fm1, which.beta = 2:1)
confint(pr1)
pr1 <- profile(fm1, which.beta = 1)
confint(pr1)
pr1 <- profile(fm1, which.beta = 2)
confint(pr1)
pr1 <- try(profile(fm1, which.beta = 0), silent = TRUE) ## error
pr1 <- try(profile(fm1, which.beta = "no.par"), silent = TRUE) ## error
pr1 <- try(profile(fm1, which.beta = -1), silent = TRUE) ## error
pr1 <- profile(fm1, which.beta = "tempwarm") 
confint(pr1)
pr1 <- profile(fm1, alpha = 0.1)
confint(pr1) ## should give NA in this case?
pr1 <- profile(fm1, max.steps = 9)
pr1 <- profile(fm1, step.warn = 7)
pr1 <- profile(fm1, nsteps = 6)
pr1 <- profile(fm1, trace = 1)
pr1 <- profile(fm1, control = list(gradTol = .1))
confint(pr1) ## not at all unreliable...

## single regression coef setting:
fm2 <- clm(rating ~ contact, data = wine)
summary(fm2)
pr2 <- profile(fm2)
confint(pr2)

## confint.clm:
confint(fm1)
confint(fm1, 2)
confint(fm1, 1)
confint(fm1, "tempwarm")
confint(fm1, type = "profile")
confint(fm1, type = "Wald")
confint(fm1, 2, type = "Wald")
confint(fm1, level = 0.5)
confint(fm1, level = 1 - 1e-6)
confint(fm1, level = 1 - 1e-10) ## extreme, but it works
confint(fm1, trace = 1)

## plot.profile:
pr1 <- profile(fm1, which.beta=1:2, alpha = 1e-3)
par(mfrow = c(1,2))
plot(pr1)
plot(pr1, 1)
plot(pr1, "contactyes")
plot(pr1, level = .97)
plot(pr1, Log = TRUE)
plot(pr1, relative = FALSE)
plot(pr1, root = TRUE)
plot(pr1, approx = TRUE)
plot(pr1, n=10)
plot(pr1, ylim = c(0,2))
plot(pr1, las = 1)
plot(pr2)

