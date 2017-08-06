.normArgRanges <- function(ranges_, errorOnEmpty=FALSE, warnOnEmpty=TRUE) {
    var_ <- deparse(substitute(ranges_))
    if (is(ranges_, "GRanges") || is(ranges_, "IRanges")) {
        if (length(ranges_) == 0)
            if (!errorOnEmpty) {
                if (warnOnEmpty)
                    warning(paste0("'", var_, "' is empty"))
            } else {
                stop(paste0("'", var_, "' cannot be empty"))
            }
        return(ranges_)
    }
    stop(paste0("'", var_, "' must be a GRanges or an IRanges object"))
}
