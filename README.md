# hestonmcmc
An R package `hestonmcmc` for Estimating Heston Stochastic Volatility Models using MCMC method

Heston Stochastic Volatility Model

Heston, S. L. (1993). A closed-form solution for options with stochastic volatility with applications
to bond and currency options. *The review of financial studies*, 6(2), 327-343.
[https://doi.org/10.1093/rfs/6.2.327](https://doi.org/10.1093/rfs/6.2.327)

## Sample Data

### Option 1: using existing sample data

`Y_series` is a sample trajectory with 1,000 observations which are generated
by Euler approximation. The true parameters are
|mu | k | theta | sigma_v | rho|
|:---|:---|:---|:---|:----|:---|
|0.125|0.1|0.25|0.1|-0.7|
and the time unit h = 1.

```r
library(hestonmcmc)
y = Y_series
```

### Option 2: generate a sample

Just generate a random sample trajectory using crHeston()

```r
S0 = c(0.125,0.1,0.25,0.1,-0.7)
y_series_0 = crHeston(v_0=S0[3], n_segment=10, par=S0, N=1000, h=1)
```


### Option 3: Generate and Write to Files Large Number of Sample Trajectories

```r
par0 = c(0.125,0.1,0.25,0.1,-0.7)
gen_data('par0', par0, N=100000, N_rep=200, h=1, n_segment=10)
```

## Estimate the parameters

### Example Using `mcmc()`

```{r}
parameters_record_v = mcmc(y,echo = TRUE)
# record_v is the last 100 trajectories
```

The priors and other setting refers to mcmc(), by typing `?mcmc` + Enter in the
Console.

### Example Using `cmcmc()`

```r
parameters = cmcmc(Y_series,g=90000,G=100000)
parameters
```

Run records:

|     |mu| k | theta | sigma_v | rho|
|:---:|:----|:---|:----------|:-----------|:------|
|true values|0.125|0.1 |0.25       |0.1         |-0.7   |
|run1(g=2K,G=10K) |0.05727125 |0.01633974 |0.10742531 |0.01051467 |-0.06455093|
|run2(g=2K,G=10K) |0.05790192 |0.01503496 |0.10796416 |0.01009222 |-0.05384306|
|run3(g=2K,G=10K) |0.05849966 |0.01498290 |0.10895028 |0.01038532 |-0.04910205|
|run4(g=5K,G=10K) |0.06004212 |0.01577060 |0.11097616 |0.01056999 |-0.08140302|
|run5(g=15K,G=20K)|0.06744194 |0.02566914 |0.12225802 |0.01139778 |-0.21733187|
|run6(g=90K,G=100K)|0.11382411|0.05106944 |0.21238676 |0.02905102 |-0.71997816|

### Example Using `mcmc2()`

```r
parameters = cmcmc2(Y_series, g=2000, G =10000, G_sub=10)
parameters
```

Run records:

|     |mu| k | theta | sigma_v | rho|
|:---:|:----|:---|:----------|:-----------|:------|
|true values|0.125|0.1 |0.25       |0.1         |-0.7   |
|run1(g=2K,G=10K,G_sub=10)|0.09760627 |0.04726157 |0.17978683 |0.02281419 |-0.58624698|


## Notes

The cpp version: all R functions have a equivalent cpp implemented version, 
such as following:

* `mcmc()` v.s. `cmcmc()`
* `rmu()` v.s. `crmu()`
* `rk()`  v.s. `crk()`
* and so on.
