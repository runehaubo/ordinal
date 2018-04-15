## This file contains:
## Functions that set control parameters for clm() and clmm().

clm.control <-
  function(method = c("Newton", "model.frame", "design", "ucminf", "nlminb",
             "optim"), ...,  trace = 0L, maxIter = 100L, gradTol = 1e-6,
           maxLineIter = 15L, relTol = 1e-6, tol = sqrt(.Machine$double.eps),
           maxModIter = 5L,
           convergence=c("warn", "silent", "stop", "message"))
{
  method <- match.arg(method)
  convergence <- match.arg(convergence)

  if(!all(is.numeric(c(maxIter, gradTol, maxLineIter, relTol, tol,
                       maxModIter))))
      stop("maxIter, gradTol, relTol, tol, maxModIter and maxLineIter should all be numeric")

  ctrl <- list(method = method,
               convergence = convergence,
               trace = as.integer(trace),
               maxIter = as.integer(maxIter),
               gradTol = as.numeric(gradTol),
               relTol = as.numeric(relTol),
               tol = as.numeric(tol),
               maxLineIter = as.integer(maxLineIter),
               maxModIter = as.integer(maxModIter))
if(method %in% c("ucminf", "nlminb", "optim"))
    ctrl$ctrl <- list(trace = as.integer(abs(trace)), ...)

  return(ctrl)
}

clmm.control <-
  function(method = c("nlminb", "ucminf", "model.frame"),
           ..., trace = 0, maxIter = 50, gradTol = 1e-4,
           maxLineIter = 50, useMatrix = FALSE,
           innerCtrl = c("warnOnly", "noWarn", "giveError"))
{
  method <- match.arg(method)
  innerCtrl <- match.arg(innerCtrl)
  useMatrix <- as.logical(useMatrix)
  stopifnot(is.logical(useMatrix))
  ctrl <- list(trace=if(trace < 0) 1 else 0,
               maxIter=maxIter,
               gradTol=gradTol,
               maxLineIter=maxLineIter,
               innerCtrl=innerCtrl)
  optCtrl <- list(trace = abs(trace), ...)

  if(!is.numeric(unlist(ctrl[-5])))
    stop("maxIter, gradTol, maxLineIter and trace should all be numeric")
  if(any(ctrl[-c(1, 5)] <= 0))
    stop("maxIter, gradTol and maxLineIter have to be > 0")
  if(method == "ucminf" && !"grtol" %in% names(optCtrl))
    optCtrl$grtol <- 1e-5
  if(method == "ucminf" && !"grad" %in% names(optCtrl))
    optCtrl$grad <- "central"

  list(method = method, useMatrix = useMatrix, ctrl = ctrl,
       optCtrl = optCtrl)
}

