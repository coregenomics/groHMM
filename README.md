## groHMM

[![Build Status](https://api.travis-ci.org/coregenomics/groHMM.svg)](https://travis-ci.org/coregenomics/groHMM)
[![codecov.io](https://codecov.io/gh/coregenomics/groHMM/coverage.svg)](https://codecov.io/gh/coregenomics/groHMM)

## Usage

Install the groHMM latest version using:

``` R
devtools::install_github("coregenomics/groHMM",
                         ref = "1.99.x",
                         repos = BiocInstaller::biocinstallRepos())
```

If the above command fails, install Bioconductor and the `devtools` package.

``` R
source("https://bioconductor.org/biocLite.R")
install.packages(devtools)
```

## Hacking

### This repo

Fork or clone the git repository,
enter the git directory,
then install the dependencies:

``` R
source("https://bioconductor.org/biocLite.R")
install.packages(devtools)
devtools::install(repos = BiocInstaller::biocinstallRepos(),
                  dependencies = c("Imports", "Suggests"))
```

### Unit tests and coverage

Run the unit tests with:

``` R
devtools::test()
```

Check your test coverage with [`covr`](https://github.com/jimhester/covr):

``` R
library(IRanges)
library(covr)

## Collapse integer ranges in string form.
## 
## Show problem lines in a collapsed range form. For example, if
## indices c(1:3, 9:11) have no coverage, show as 1-3,9-11 instead of
## 1 2 3 9 10 11.
str_range <- function(x) {
	ir <- IRanges::IRanges(start = x, width = rep(1, length(x)))
	ir <- IRanges::reduce(ir)
	ranges <- paste(IRanges::start(ir), IRanges::end(ir), sep = "-")
	single_lines <- IRanges::start(ir) == IRanges::end(ir)
	ranges[single_lines] <- as.character(IRanges::start(ir)[single_lines])
	paste(ranges, collapse = ",")
}

## Collapse coverage output to filenames and lines.
lines <- function(coverage, pattern = NULL) {
	df <- zero_coverage(coverage)
	df_grouped <- dplyr::group_by(df, filename)
	df <- dplyr::summarise(df_grouped,
                           lines = str_range(line))
	if (!is.null(pattern))
	    df <- dplyr::filter(df, grepl(filename, pattern = pattern))
    message(paste(df$filename, df$lines, collapse = "\n\n"))
	invisible(df)
}
```

Then run it with:

``` R
cov <- package_coverage()
cov
lines(cov)
```

### Package quality checks

Run R's standard check followed by BiocCheck:

``` R
devtools::check()
BiocCheck::BiocCheck(".")
```
