library(ordinal)

## Testing that the profile remains the same - that the model object
## is not 'distorted' by update(object/fitted, doFit=FALSE)
set.seed(1234)
wts <- runif(nrow(wine), 0, 2)
fm3 <- clm(rating ~ temp + contact, data=wine,
           weights=wts)
pr <- profile(fm3)

set.seed(1234)
fm3 <- clm(rating ~ temp + contact, data=wine,
           weights=runif(nrow(wine), 0, 2))
pr3 <- profile(fm3)
## > set.seed(1234)
## > fm3 <- clm(rating ~ temp + contact, data=wine,
## +            weights=runif(nrow(wine), 0, 2))
## > pr3 <- profile(fm3)
## Warning messages:
## 1: In profile.clm.beta(fitted, which.beta, alpha, max.steps, nsteps,  :
##   profile may be unreliable for tempwarm because only 1
##   steps were taken down
## 2: In profile.clm.beta(fitted, which.beta, alpha, max.steps, nsteps,  :
##   profile may be unreliable for tempwarm because only 1
##   steps were taken up
## 3: In profile.clm.beta(fitted, which.beta, alpha, max.steps, nsteps,  :
##   profile may be unreliable for contactyes because only 1
##   steps were taken down
## 4: In profile.clm.beta(fitted, which.beta, alpha, max.steps, nsteps,  :
##   profile may be unreliable for contactyes because only 1
##   steps were taken up
##
stopifnot(isTRUE(all.equal(pr, pr3, check.attributes=FALSE)))
stopifnot(
    isTRUE(all.equal(pr$tempwarm[, "lroot"], pr3$tempwarm[, "lroot"])),
    isTRUE(all.equal(pr$contactyes[, "lroot"], pr3$contactyes[, "lroot"])))
