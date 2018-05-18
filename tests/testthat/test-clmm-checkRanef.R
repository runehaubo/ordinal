context("Testing error-warning-message from clmm via checkRanef")

## Make example with more random effects than observations:
wine$fake <- factor(c(1:65, 1:65)[1:nrow(wine)])
wine$fakeToo <- factor(1:nrow(wine))

## Check warning, error and 'message' messages:
expect_warning(
  fmm2 <- clmm(rating ~ temp + contact + (1|judge) + (1|fake), data=wine)
  , "no. random effects")

expect_warning(
  fmm2 <- clmm(rating ~ temp + contact + (1|judge) + (1|fake), data=wine,
               checkRanef="warn")
  , "no. random effects")

expect_error(
  fmm2 <- clmm(rating ~ temp + contact + (1|judge) + (1|fake), data=wine,
               checkRanef="error")
  , "no. random effects")

expect_message(
  fmm2 <- clmm(rating ~ temp + contact + (1|judge) + (1|fake), data=wine,
               checkRanef="message")
  , "no. random effects")

expect_error(
  fmm2 <- clmm(rating ~ temp + contact + (1|fakeToo), data=wine,
               checkRanef="error")
  , "no. random effects")

expect_error(
  fmm2 <- clmm(rating ~ temp + contact + (1|judge) + (1|fakeToo), data=wine,
               checkRanef="error")
  , "no. random effects")
