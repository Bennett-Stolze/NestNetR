# NestNetR

NestNetR is an R-package that provides an end-to-end workflow to **preprocess light-level geolocator data** and **classify breeding behaviour** (e.g. incubation and brooding activity) using a **pre-trained deep-learning model** based on **light intensity** and **temperature** recordings.

## Installation

Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("Bennett-Stolze/NestNetR")
```

## Deep-learning backend
NestNetR uses keras3 with a TensorFlow backend for behaviour classification. Please follow the instructions to install keras3 and TensorFlow in R: https://cran.r-project.org/web/packages/keras3/vignettes/getting_started.html
```r
install.packages("keras3")
keras3::install_keras(backend = "tensorflow")
```

## Documentation & Usage
A comprehensive package vignette is included and provides a step-by-step introduction to the complete workflow, including data preprocessing, breeding-period detection, segmentation, and behaviour classification.
```r
library(calibrar)
vignette("NestNetR")
```