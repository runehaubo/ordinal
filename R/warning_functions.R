givesWarnings <- function(expr) countWarnings(expr) > 0L

countWarnings <- function(expr)
{
    .number_of_warnings <- 0L
    frame_number <- sys.nframe()
    ans <- withCallingHandlers(expr, warning = function(w)  {
        assign(".number_of_warnings", .number_of_warnings + 1L,
               envir = sys.frame(frame_number))
        invokeRestart("muffleWarning")
    })
    .number_of_warnings
}

