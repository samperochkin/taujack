# taujack
Code computing (Generalized) Kendall's tau and its jackknife variance under serial dependence. It is tied to the paper:

Perreault, S. (2024) Simultaneous computation of Kendall’s tau and its jackknife variance. Statistics & Probability Letters, 213, 110181.
doi.org/10.1016/j.spl.2024.110181 

### Usage
The folder src/ contains the C++ implementations of the algorithms and the functions.R file contains wrappers.
It should be clear in the R file which C++ functions are used.

```
library(Rcpp)

sourceCpp("src/ms.cpp")       # for Knight's extended alg.
sourceCpp("src/dac_seq.cpp")  # divide-and-conquer alg.
sourceCpp("src/bf.cpp")       # brute force alg.
source("functions.R")         # wrappers (performs re-ordering if necessary)
```
