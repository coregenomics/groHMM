## Workaround for system.file() failing with devtools::check() and Travis CI:
## https://github.com/hadley/testthat/issues/86#issuecomment-143878211
Sys.setenv("R_TESTS" = "")

library(testthat)
library(groHMM)

test_check("groHMM")
