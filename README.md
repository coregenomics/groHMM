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
install.packages("devtools")
```

## Hacking

The place where groHMM needs the most attention
is increasing the test coverage
by adding unit tests to `./inst/tests/`.
One can see where work is needed from the web reports of
[codecov](https://codecov.io/gh/coregenomics/groHMM).

Otherwise, one can check
[Travis CI](https://travis-ci.org/coregenomics/groHMM)
for notes and warnings to fix.

To tackle this and more,
it's best not to rely on the web CI services
(which are very thorough but take an hour to run)
and to instead
run the tests, coverage and checks locally (which take tens of seconds)
as explained below.

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

The web integration reports take about ~1 hour to generate;
the pacing item being
[Travis CI](https://travis-ci.org/coregenomics/groHMM).

One can therefore run
[`covr`](https://github.com/jimhester/covr) and the checks locally:

``` R
library(covr)
cov <- package_coverage()
cov
zero_coverage(cov)
```

### Package quality checks

Besides testing and coverage,
one can catch a broader range of quality issues using
R's standard check followed by BiocCheck:

``` R
devtools::check(build_args="--no-build-vignettes")
BiocCheck::BiocCheck(".")
```

### R daily build

Strictly speaking,
Bioconductor's development process requires using a recent
[R daily build](http://bioconductor.org/developers/how-to/useDevel/).
There are a few different approaches for compiling from source:
for example one can use the package manager to install build-time dependencies;
on Debian one could run `apt-get build-dep r-base`.
However I suggest using `spack` and environmental modules instead
to allows one to easily switch between R versions
and keep a separate R library for each installation.

``` bash
cd
git clone https://github.com/llnl/spack.git
spack/bin/spack install r@$(date -d -I yesterday)
```

Install environmental-modules to load our new r module:

``` bash
aptitude install environment-modules
```

Note the version of r that you have loaded.
We will add the version to our `R-devel` launcher script.

``` bash
spack/bin/spack find r
```

Create a file `/usr/local/bin/R-devel`:

``` bash
#!/bin/bash
. /etc/profile.d/modules.sh
. ~/spack/bin/spack/setup-env.sh
spack load r@2017-08-01  # Change version to match `spack find r`
exec R "$@"
```

In emacs `ess-mode`
as long as `R-devel` is in your PATH,
one can launch it with `M-x R-devel`
as your associated `R` shell.
If you use RStudio or RStudio server,
to use `R-devel` as your interpreter
with environmental variables or Rprofile hacks
you're on your own `(^_~)`
