# rewind
[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)

An R Package with Biased Resampling Strategies for Time Series (and other Dependency-Oriented Data in the future)

For Imbalanced Domain Learning, resampling strategies are some of the most known tools. This package focuses on the problem of predicting extreme numerical values in dependency-oriented data such as time series. The package rewind presents a set of biased resampling strategies which are focused on applying such methods while taking into account the dependencies between consecutive data points in time.

**To install from github use the following command lines in R:**

    library(devtools)  # You need to install this package!
    install_github("nunompmoniz/rewind",ref="master")

After installation the package can be used loaded by doing:

     library(rewind)
