# test.makeThresholds.R

library(ordinal)

# Prvious bug which is now fixed:
res <- ordinal:::makeThresholds(letters[1:3], "symmetric")
stopifnot(length(res$alpha.names) == res$nalpha)
# length(res$alpha.names) used to be 4

# Real data example:
wine <- within(wine, {
  rating_comb3b <- rating
  levels(rating_comb3b) <- c("1-2", "1-2", "3", "4-5", "4-5")
}) 
wine$rating_comb3b[1] <- "4-5" # Need to remove the zero here to avoid inf MLE
ftable(rating_comb3b ~ temp + contact, data=wine)

fm.comb3_c <- clm(rating_comb3b ~ contact, #scale=~contact, 
                  threshold = "symmetric", data=wine) # no error
