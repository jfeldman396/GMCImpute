The Gaussian Mixture Copula for Inference with Missing Data
================
Joe Feldman
7/10/2022

# Installation

``` r
# library(devtools)
# remotes::install_github('jfeldman396/GMCImpute')

library(GMCImpute)
```

# Background: Missing Data

Missing data is commonplace in survey and fused data sets, necessitating
sophisticated methods for dealing with missingness when deriving
inference. First, these data sets may be comprised of mixed data types,
such as continuous, count, and nominal variables. Second, the missing
data may bias certain properties of the data set. As a result, a
complete case analysis, where the analyst drops any observations with
missing values, would lead to misleading insights on the data.

Consider the following simulated example of a continuous variable
![X\_{1}](https://latex.codecogs.com/png.latex?X_%7B1%7D "X_{1}"), a
count variable
![X\_{2}](https://latex.codecogs.com/png.latex?X_%7B2%7D "X_{2}"), and a
binary variable
![X\_{3}](https://latex.codecogs.com/png.latex?X_%7B3%7D "X_{3}")
simulated under the following dependence structure:

![ X\_{1} \\sim N(0,1) \\\\
X\_{2} \\mid X\_{1} = x\_{1} \\sim Poisson(5 \* \\lvert x\_{1} \\rvert) \\\\
X\_{3} \\mid X\_{1} = x\_{1}, X\_{2} = x\_{2} \\sim Bernoulli(scale(x\_{2}))](https://latex.codecogs.com/png.latex?%20X_%7B1%7D%20%5Csim%20N%280%2C1%29%20%5C%5C%0AX_%7B2%7D%20%5Cmid%20X_%7B1%7D%20%3D%20x_%7B1%7D%20%5Csim%20Poisson%285%20%2A%20%5Clvert%20x_%7B1%7D%20%5Crvert%29%20%5C%5C%0AX_%7B3%7D%20%5Cmid%20X_%7B1%7D%20%3D%20x_%7B1%7D%2C%20X_%7B2%7D%20%3D%20x_%7B2%7D%20%5Csim%20Bernoulli%28scale%28x_%7B2%7D%29%29 " X_{1} \sim N(0,1) \\
X_{2} \mid X_{1} = x_{1} \sim Poisson(5 * \lvert x_{1} \rvert) \\
X_{3} \mid X_{1} = x_{1}, X_{2} = x_{2} \sim Bernoulli(scale(x_{2}))")

``` r
num= 500
X1<-rnorm(num)
X2<- rpois(num,5*abs(X1))
X3<- as.factor(rbernoulli(num,pnorm(-.5+scale(X2))))
X<- data.frame(X1,X2,X3)
```

Next, we introduce a missingness mechanism
![\\boldsymbol{R}](https://latex.codecogs.com/png.latex?%5Cboldsymbol%7BR%7D "\boldsymbol{R}")
that creates bias across the data set by correlating missing values in
![X\_{2}](https://latex.codecogs.com/png.latex?X_%7B2%7D "X_{2}") and
![X\_{3}](https://latex.codecogs.com/png.latex?X_%7B3%7D "X_{3}") with
![X\_{1}](https://latex.codecogs.com/png.latex?X_%7B1%7D "X_{1}").
Specifically, if
![R\_{ij} = 1](https://latex.codecogs.com/png.latex?R_%7Bij%7D%20%3D%201 "R_{ij} = 1"),
then
![X\_{ij} = missing](https://latex.codecogs.com/png.latex?X_%7Bij%7D%20%3D%20missing "X_{ij} = missing").
We use the following missingness mechanism in this example:

![ P(R\_{ij} = 1 \\mid X\_{1} = x\_{1}) = Bernoulli(-0.5 + 0.5\*\\Phi(x\_{1}))](https://latex.codecogs.com/png.latex?%20P%28R_%7Bij%7D%20%3D%201%20%5Cmid%20X_%7B1%7D%20%3D%20x_%7B1%7D%29%20%3D%20Bernoulli%28-0.5%20%2B%200.5%2A%5CPhi%28x_%7B1%7D%29%29 " P(R_{ij} = 1 \mid X_{1} = x_{1}) = Bernoulli(-0.5 + 0.5*\Phi(x_{1}))")

Where ![\\Phi](https://latex.codecogs.com/png.latex?%5CPhi "\Phi") is
the standard normal cdf.

``` r
R = t(sapply(1:num, function(x)rbernoulli(2, p = pnorm(-.5 + .5*abs(X1[x]))))) 
X_noMis = X
X[which(R[,1] == T),2] = NA
X[which(R[,2] == T),3] = NA
```

We can visualize bias that the mechanism creates in
![X\_{2}, X\_{3}](https://latex.codecogs.com/png.latex?X_%7B2%7D%2C%20X_%7B3%7D "X_{2}, X_{3}")
with the following plots.
![\\boldsymbol{R}](https://latex.codecogs.com/png.latex?%5Cboldsymbol%7BR%7D "\boldsymbol{R}")
removes values of
![X\_{2}](https://latex.codecogs.com/png.latex?X_%7B2%7D "X_{2}") and
![X\_{3}](https://latex.codecogs.com/png.latex?X_%7B3%7D "X_{3}") for
large values of
![X\_{1}](https://latex.codecogs.com/png.latex?X_%7B1%7D "X_{1}") in
absolute value.

    ## Warning: Removed 216 rows containing missing values (geom_point).

![](README_files/figure-gfm/unnamed-chunk-2-1.png)<!-- --> As a result,
both margins are affected. We show a comparison of the empirical cdfs of
![X\_{2}](https://latex.codecogs.com/png.latex?X_%7B2%7D "X_{2}") before
and after inputting missing values, while the incidence of positive
indicators is greatly reduced for
![X\_{3}](https://latex.codecogs.com/png.latex?X_%7B3%7D "X_{3}")

![](README_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

# The Gaussian Mixture Copula with Margin Adjustment

The function `GMC_Impute` allows users to fit a Gaussian mixture copula
to data comprised of unordered categorical, binary, count and continuous
data types with missing values. This is done through utilization of the
, which enables copula estimation on the aforementioned data types. The
function then produces a user specified number of multiple imputations.

The model is given by

![z\_{i} \\sim N(\\Lambda \\eta\_{i}, \\Sigma),\\ \\eta\_{i} \\sim \\sum\_{h=1}^{H^{\*}\\leq H} \\pi\_{h}N(\\mu\_{h},\\Delta\_{h})\\\\
\\Sigma \\sim Gamma(a\_{\\sigma},b\_{\\sigma}),\\lambda\_{jh} \\sim N(0, \\phi\_{jh}^{-1}\\tau\_{h}^{-1})\\\\
\\phi\_{jh} \\sim Gamma(\\nu/2,\\nu/2), \\tau\_{h} = \\prod\_{l=1}^{h}\\delta\_{h}^{\\tau} \\\\
\\text{with} \\ \\delta^{\\tau}\_{1} \\sim Gamma(a\_{1},1) \\ \\text{and} \\ \\delta^{\\tau}\_{l} \\sim Gamma(a\_{2},1), \\ l\\geq 2, \\ a\_{2} \\geq 1](https://latex.codecogs.com/png.latex?z_%7Bi%7D%20%5Csim%20N%28%5CLambda%20%5Ceta_%7Bi%7D%2C%20%5CSigma%29%2C%5C%20%5Ceta_%7Bi%7D%20%5Csim%20%5Csum_%7Bh%3D1%7D%5E%7BH%5E%7B%2A%7D%5Cleq%20H%7D%20%5Cpi_%7Bh%7DN%28%5Cmu_%7Bh%7D%2C%5CDelta_%7Bh%7D%29%5C%5C%0A%5CSigma%20%5Csim%20Gamma%28a_%7B%5Csigma%7D%2Cb_%7B%5Csigma%7D%29%2C%5Clambda_%7Bjh%7D%20%5Csim%20N%280%2C%20%5Cphi_%7Bjh%7D%5E%7B-1%7D%5Ctau_%7Bh%7D%5E%7B-1%7D%29%5C%5C%0A%5Cphi_%7Bjh%7D%20%5Csim%20Gamma%28%5Cnu%2F2%2C%5Cnu%2F2%29%2C%20%5Ctau_%7Bh%7D%20%3D%20%5Cprod_%7Bl%3D1%7D%5E%7Bh%7D%5Cdelta_%7Bh%7D%5E%7B%5Ctau%7D%20%5C%5C%0A%5Ctext%7Bwith%7D%20%5C%20%5Cdelta%5E%7B%5Ctau%7D_%7B1%7D%20%5Csim%20Gamma%28a_%7B1%7D%2C1%29%20%5C%20%5Ctext%7Band%7D%20%5C%20%5Cdelta%5E%7B%5Ctau%7D_%7Bl%7D%20%5Csim%20Gamma%28a_%7B2%7D%2C1%29%2C%20%5C%20l%5Cgeq%202%2C%20%5C%20a_%7B2%7D%20%5Cgeq%201 "z_{i} \sim N(\Lambda \eta_{i}, \Sigma),\ \eta_{i} \sim \sum_{h=1}^{H^{*}\leq H} \pi_{h}N(\mu_{h},\Delta_{h})\\
\Sigma \sim Gamma(a_{\sigma},b_{\sigma}),\lambda_{jh} \sim N(0, \phi_{jh}^{-1}\tau_{h}^{-1})\\
\phi_{jh} \sim Gamma(\nu/2,\nu/2), \tau_{h} = \prod_{l=1}^{h}\delta_{h}^{\tau} \\
\text{with} \ \delta^{\tau}_{1} \sim Gamma(a_{1},1) \ \text{and} \ \delta^{\tau}_{l} \sim Gamma(a_{2},1), \ l\geq 2, \ a_{2} \geq 1")

For the mixing weights, we use the truncated dirichlet process prior,
which specifies an upper bound for the number of clusters and allows the
data to inform which
![H^{\*}](https://latex.codecogs.com/png.latex?H%5E%7B%2A%7D "H^{*}")
get populated

![\\pi\_{h} = V\_{h}\\prod\_{l&lt;h} (1-V\_{l}),\\  V\_{h}\\sim Beta(1, \\alpha\_{\\pi})\\\\
        \\alpha\_{\\pi} \\sim Gamma(a\_{\\alpha}, b\_{\\alpha})](https://latex.codecogs.com/png.latex?%5Cpi_%7Bh%7D%20%3D%20V_%7Bh%7D%5Cprod_%7Bl%3Ch%7D%20%281-V_%7Bl%7D%29%2C%5C%20%20V_%7Bh%7D%5Csim%20Beta%281%2C%20%5Calpha_%7B%5Cpi%7D%29%5C%5C%0A%20%20%20%20%20%20%20%20%5Calpha_%7B%5Cpi%7D%20%5Csim%20Gamma%28a_%7B%5Calpha%7D%2C%20b_%7B%5Calpha%7D%29 "\pi_{h} = V_{h}\prod_{l<h} (1-V_{l}),\  V_{h}\sim Beta(1, \alpha_{\pi})\\
        \alpha_{\pi} \sim Gamma(a_{\alpha}, b_{\alpha})")

For the prior on
![\\eta](https://latex.codecogs.com/png.latex?%5Ceta "\eta"), we use a
Normal-Inverse Wishart:
![(\\boldsymbol{\\mu\_{h}, \\Delta\_{h}}) \\sim Normal-InverseWishart(\\boldsymbol{\\mu\_{0}}, \\delta^{2}\*\\boldsymbol{I\_{k}},\\kappa\_{0}, \\nu\_{mix})](https://latex.codecogs.com/png.latex?%28%5Cboldsymbol%7B%5Cmu_%7Bh%7D%2C%20%5CDelta_%7Bh%7D%7D%29%20%5Csim%20Normal-InverseWishart%28%5Cboldsymbol%7B%5Cmu_%7B0%7D%7D%2C%20%5Cdelta%5E%7B2%7D%2A%5Cboldsymbol%7BI_%7Bk%7D%7D%2C%5Ckappa_%7B0%7D%2C%20%5Cnu_%7Bmix%7D%29 "(\boldsymbol{\mu_{h}, \Delta_{h}}) \sim Normal-InverseWishart(\boldsymbol{\mu_{0}}, \delta^{2}*\boldsymbol{I_{k}},\kappa_{0}, \nu_{mix})")

Key to gaussian mixture copula are the marginal distributions of each
variable in the data, as latent variables, modeled with finite mixture
are linked to the observed scale using the inverse marginal distribution
function. Previous work estimates these margins empirically, which is
problematic given that the missing data clearly biases these estimates.
The margin adjustment corrects these biases, yielding proper inference
with missing data.

# Fitting the model

To use the function, the user can specify a number of properties of the
model:

-   `nImp`: The number of imputations to create
-   `H`: The upper bound for the number of clusters in the truncated DP
    mixture
-   `k.star`: The dimension of the latent factors, defult is
    `ceiling(0.7*p)`
-   `nsamp`: number of interations in the MCMC
-   `burn`: burn-in before posterior samples are saved
-   `hyperparams`:
    -   `delta`: to the precision of the prior covariance. This
        parameter has been the most influential in the discovery of new
        clusters. Lower to find more clusters. Default value is 10
    -   `k.star`: change to increase or decrease dimension of latent
        factors
    -   `a_alpha`
    -   `b_alpha`
    -   `nu_mix`
    -   `kappa_0`
    -   `nu`
    -   `a1`
    -   `a2`
    -   `a.sigma`
    -   `b.sigma`
    -   `D_0`: k.star dimensional identity

Default values are included in the function documentation, but we
recommend altering
![\\delta](https://latex.codecogs.com/png.latex?%5Cdelta "\delta") to
improve model fit. The function is called below:

``` r
hyperparams = list(delta = 10,
                   k.star = 2,
                   plugin.threshold = 100,
                   a_alpha = 1,
                   b_alpha = 3,
                   nu_mix = 4,
                   kappa_0 = .001,
                   nu = 3,
                   a1 = 2,
                   a2 = 3,
                   a.sigma = 1,
                   b.sigma = .3,
                   D_0 = diag(1,2))
imps<-GMC_Impute(Data = X, nImp = 5,hyperparams = hyperparams, burn = 5000,nsamp = 10000, seed = 47)
```

`GMC_Impute` returns `nImp` imputations, as well as posterior samples of
Copula parameters which may be used for simulation of posterior
predictive data sets or posterior inference. See documentation for
format.

# Plotting Results

## Visualizing Imputations:

    ## Warning: Removed 2 rows containing missing values (geom_point).

![](README_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

## Plotting posterior samples from the margin adjustment

We can plot posterior samples of the marginal distribution of
![X\_{2}](https://latex.codecogs.com/png.latex?X_%7B2%7D "X_{2}"), as
well as point-wise posterior means:

``` r
  par(mar = c(5,6,4,2))
range = range(imps$Support[[3]][1,2:18]) # get support
quantiles = imps$Quantiles[[3]] #get quantiles
 plot(range[1]:range[2],quantiles[1,2:(range[2]+2)],
      col = "gray",
      type = 'b',
       xlab = expression(X[2]),
       ylab = expression(P(X[2] <= x)),
       main = "Posterior Samples of F_X2",
       cex.lab = 1.5,
       cex.main = 2)
  sapply(2:5000,function(x)points(range[1]:range[2],quantiles[x,2:(range[2]+2)], col ='gray', type = 'b'))
  lines(ecdf(X_noMis$X2), cex = 2)
  points(range[1]:range[2],colMeans(quantiles[2:5000,2:(range[2]+2)]),pch =2, cex = 2)
  lines(ecdf(X$X2), col = 2, cex = 2)
    legend("bottomright",c("ECDF w/o Mis","ECDF w/ Mis","Posterior Mean"), pch = c(16,16,2),col = c(1,2,1),bty = 'n', cex = 1.3, text.font = 2)
```

![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

# Simulating Predictive Data Sets

Finally, we can use the posterior samples to generative posterior
predictive data sets for checks and inference. This is done by using the
returned samples from GMC\_Impute

``` r
#use the 2000th sample as an example

#GMC parameters
mu = lapply(1:25, function(x) return(imps$mus[[x]][2000,]))
alpha = rep(0,3)
Delta = lapply(1:25, function(x) return(imps$Deltas[[x]][2000,,]))
Lambda = imps$Lambdas[2000,,]
pi_h = imps$pis[2000,]
Sigma.diag= imps$Sigmas[2000,]
#column names and memberships
col_mem = imps$col_mem
cat_col_names = imps$cat_col_names
bin_col_names = imps$bin_col_names
count_col_names = imps$count_col_names
cont_col_names = imps$con_col_names
H = 25
Y = imps$dat
z = imps$zs[2000,]
Fns = vector('list', 3) #format marginal distributions
for(i in 1:3){
  support = imps$Support[[i]][2000,]
  qs = imps$Quantiles[[i]][2000,]
  Fj = cbind(support,rep(0,length(support)),qs)
  Fns[[i]] = Fj
  
}

Y = imps$dat #for 
H =25

pred<- get_predictive_Y(n = dim(Y)[1],
                          Lambda = Lambda,
                          mu = mu,
                          alpha = alpha,
                          Delta = Delta,
                          pi_h = pi_h,
                          Sigma.diag = Sigma.diag,
                          col_mem = col_mem,
                          cat_col_names = cat_col_names,
                          bin_col_names = bin_col_names, 
                          count_col_names = count_col_names,
                          cont_col_names = cont_col_names,
                          H = H,
                          Y = Y,
                          z = z,
                          Fns = Fns,
                          seed = 10)
```

![](README_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->
