# GreenSimulation
Nested Green Simulation Algorithm.

## Description
This package is an implementation of Nested Green Simulation Algorithm.

## Installation
Install `GreenSim` from GitHub:
```r
install.packages("devtools")
library(devtools)
devtools::install_github("chenqi57/GreenSim")
library(GreenSim)
```
or directly from the `GreenSim_0.1.0.tar.gz`:
```r
install.packages("GreenSim_0.1.0.tar.gz", repos = NULL, type = "source")
library(GreenSim)
```
## Function
Three functions are contained in this package, where two functions are the implementation of the MMDS algorithm.

* `eigen_centered`: Center the sample matrix *X* to *Y* and return *Y*'s corresponding distance matrix.
* `MMDS`: Modified Multi-Dimensional Scaling Algorithm using eigen-decompision methods provided by R.
* `MMDS`: Modified Multi-Dimensional Scaling Algorithm using eigen-decompision methods provided by C++.

## Documentation
Documentations of the above functions can be accesed by typing `?` before each function's name at the R command. 
For instance, the user can read the function `MMDS`'s argument, output and examples in detail by typing `?MMDS`.

## Example
Suppose you have a 1000×1500 sample matrix `sample`, where the sample consists of 1500 data points from five guassian distributions with same covariance matrix 0.45I(1000) but different means. The information about which group each data point belongs to is stored in `Labels`.

`MMDS`
```r
library(mmds)

sample_MMDS = MMDS(X = t(sample), MM = 2, sigma = sqrt(0.45), centered = FALSE)
data = data.frame(sample_MMDS)
data$label = Labels

ggplot(data,aes(x = data[, 1], y = data[, 2], colour = Labels)) + 
geom_point(size = 1) + xlab("") + ylab("") + theme(legend.position = "none")
```
![image](https://imgur.com/md5H41j.png#center)

`MMDS.cpp`
```r
sample_MMDS_cpp = MMDS.cpp(X = t(sample), MM = 2, sigma = sqrt(0.45), centered = FALSE)
data2 = data.frame(sample_MMDS_cpp)
data2$label = Labels

ggplot(data2,aes(x = data2[, 1], y = data2[, 2], colour = Labels)) +
geom_point(size = 1) + xlab("") + ylab("") + theme(legend.position = "none")
```
![image](https://imgur.com/ZPkpacl.png#center)
