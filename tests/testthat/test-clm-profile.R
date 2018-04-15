context("Testing error message from profile.clm")

expect_warning(
    fm2 <- clm(rating ~ contact, scale=~contact, nominal=~contact,
               data=wine)
    , "\\(1\\) Hessian is numerically singular")

expect_error(profile(fm2)
             , "Cannot get profile when vcov\\(fitted\\) contains NAs")

