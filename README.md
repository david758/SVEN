# SVEN (Support Vector Elastic Net)

Trains an elastic net model using the method described in Zhou et al. 2014

## How to Use

Put the files in the same directory. Load the function "SVEN" into R.

You can load test data with:

```
load("heart.rda")
```

### Parameters

[B] A nxd matrix of training examples
[y] A n length vector of training labels
[t] L1 norm constraint
[lambda] L2 regularization coefficient

