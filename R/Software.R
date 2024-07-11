#' Fit Gaussian mixture copula with the margin adjustment and impute missing data
#'
#' Fit the Gaussian mixture copula to a data set with or without missing data. The user does not have to deal with one-hot-encoding of categorical variables,
#' which is accomplished by the function. The MCMC non-parametrically estimates marginal distributions through the margin adjustment.
#' The user is given control over whether or not to impute missing data through the logical \code{Impute} and then how many completed data sets to produce with \code{nImp}.
#' In addition, the user can specify several hyperparameters under the latent mixture, including the upper bound for the number of clusters in the truncated DP,
#' the number of unique values of a variable for RPL re-sampling to take place, and prior precision of cluster specific covariance matrices.
#'
#' The function returns posterior samples of model parameters for posterior predictive checks, as well multiple completed data sets.
#' @param Data Data frame that is the input into the Gaussian Mixture Copula. Categorical and  binary variables should be factors
#' @param nImp Number of imputations to Create
#' @param Impute logical; Whether or not to create multiple imputations
#' @param H Upper bound for the number of components of the mixtures
#' @param k.star Dimension of latent factor in mixture
#' @param nsamp Number of samples for the MCMC
#' @param burn Number of burn-in samples before imputations
#' @param hyperparams Optional; Hyperparameter specifications for finite mixuture and RPL
#'  with default values given by:
#' \itemize{
#' \item \code{k.star}: dimension of latent factors, default is \code{floor(.7*p))}
#' \item \code{plugin.threshold}: threshold for number of unique values in each column for the re-sampling step to be conducted,
#' default is 350
#' \item \code{a_alpha}: prior for the stick breaking precision parameter, where the prior is given by \code{gamma(a_alpha,b_alpha)},
#' default value is 1
#' \item \code{b_alpha}: prior for the stick breaking precision parameter, where the prior is given by \code{gamma(a_alpha,b_alpha)},
#' default value is 3
#' \item \code{nu_mix}: hyper parameter for the component specific prior \eqn{N-IW(0,I_{k.star},nu_mix,kappa_0)} on \eqn{\eta},
#' default is k.star +2
#' \item \code{kappa_0}: hyper parameter for the component specific prior \eqn{N-IW(0,I_{k.star},nu_mix,kappa_0)} on \eqn{\eta},
#' default is .001
#' \item \code{nu}: hyper parameter for local shrinkage of MGP prior on components of Lambda, default is 3
#' \item \code{a1}: hyper parameter for global shrinkage of MGP prior on components of Lambda, default is 2
#' \item \code{a2}: hyper parameter for global shrinkage of MGP prior on components of Lambda, default is 3
#' \item \code{a.sigma,b.sigma}: hyper parameter for gamma prior on error variance, default value is 1 and 0.3
#' \item \code{delta}: prior precision parameter for covariance matrix in latent mixture. Increase to discover fewer clusters.
#' default value is 10
#' \item \code{D_0}: prior covariance for N-IW. Default is a k.star dimenional identity
#' }

#'
#' @return A list with the following elements:
#' \itemize{
#' \item \code{Imputations} nImp complete data sets
#' \item \code{Deltas}: List of length H, each containing an \code{(nsamp - burn x k.star x k.star)} array of
#'  posterior samples of cluster specific covariance matrices
#' \item \code{mus}: List of length H, each containing an \code{(nsamp-burn x k.star)} array of posterior samples of cluster
#' specific covariance matrices
#' \item \code{zs}: Indicies of occupied clusters of dimension \code{(nsamp-burn x H)}. Used for posterior predictive sampling to identify which components of the mixture to use
#' \item \code{clustering}: Posterior samples of size \code{(nsamp-burn x n)} with discovered clusterings
#' \item \code{Sigmas}: array of dimension \code{(nsamp-burn x p)} containing posterior samples of the error covariance
#' \item \code{Quantiles}:List of length p, each containing an array of dimension \code{(nsamp-burn x length(support))} with posterior samples of point-wise cumulative probabilities.
#' \item \code{Support}: List of length p, each containing an array of dimension \code{(nsamp-burn x length(support))} with support of each variable
#' \item \code{col_mem}: array of length p containing variable names
#' \item \code{cat_col_names}: array containing names of categorical variables, with the number of repititions indicating eachs' number of levels
#' \item \code{bin_col_names}: array containing names of binary variables
#' \item \code{count_col_names}: array containing names of count-valued variables
#' \item \code{cont_col_names}: array containing names of continous variables
#' \item \code{dat}: Original data frame, useful for simulating predictive data sets
#' \item \code{seed}: Seed for replication of results
#' \item \code{Y_aug}: augmented data set with binary representation of categorical variables.
#' }
#'
#'
#' @import tmvtnorm
#' @import mvtnorm
#' @import LaplacesDemon
#' @import truncnorm
#' @import purrr
#' @import MASS
#' @import dplyr
#' @import tidyr
#' @import data.table
#' @export
#'
#' @examples
#' num = 500 #sample size
#' X1<-rnorm(num) # simulate data set
#' X2<- rpois(num,5*abs(X1))
#' X3<- rbernoulli(num,1,pnorm(-.5+scale(X2)))
#' #remove values with missingness mechanism R
#' R = t(sapply(1:num, function(x)rbernoulli(2, p = pnorm(-.5 + .5*abs(X1[x]))))) #change to -.25 for more missing
#' X<- data.frame(X1,X2,X3)
#' X_noMis = X
#' X[which(R[,1] == T),2] = NA
#' X[which(R[,2] == T),3] = NA
#' mcmc<-GMC.mcmc(Data = X, nsamp = 1000) #fit model


GMC.mcmc<- function(Data, nImp = 10, Impute = T,H = 25, k.star = NULL, nsamp = 10000, burn = floor(.5*nsamp),
                      hyperparams = NULL, seed = NULL){

  if(any(sapply(colnames(Data), function(x) is.character(Data[,x])))){
    stop('Nominal variables must be factors')
  }

  if(is.null(seed)){

    seed = sample(1,100000,1)
    set.seed(seed)
  }else{
    set.seed(seed)
  }



  ########### format data
  Y = Data
  #Pull out columns that are factors
  factor.cols<- which(sapply(Y,is.factor))

  #Arrange Factors in order of number of Factors
  if(length(factor.cols)>1){
    inc_f = Y[,names(sort(sapply(Y[,sapply(Y, is.factor)], nlevels), decreasing = TRUE))]
    Y<- cbind(inc_f, Y[,-factor.cols])
    factor.cols <- which(sapply(Y,is.factor) == TRUE)
  }


  # Generate augmented data data set for RPL
  options(na.action='na.pass')
  cat<- NULL
  cat_col = 1

  for(i in factor.cols){
    if(nlevels(Y[,i]) >2){
      one.hot.cat<- model.matrix(~Y[,i] - 1)
      colnames(one.hot.cat)<- levels(Y[,i])
      cat<- cbind(cat,one.hot.cat)
      cat_col = cat_col + nlevels(Y[,i])



    }
    else{
      one.hot.bin<- model.matrix(~Y[,i])[,2]
      cat<-cbind(cat,one.hot.bin)
      colnames(cat)[cat_col] = colnames(Y)[i]
      cat_col = cat_col + 1
    }

  }



  one.hot.cols<- ncol(cat)
  Y_mod<- cbind(cat,Y[,-factor.cols])
  # dimension of Y_mod
  n<-dim(Y_mod)[1]
  p<-dim(Y_mod)[2]
  is_cat_bin<- c(rep(1,one.hot.cols), rep(0,p - one.hot.cols))


  ########## data
  vnames<-colnames(Y_mod)
  Y_mod<-as.matrix(Y_mod)
  colnames(Y_mod)<-vnames

  #initialize latent Z
  Z<-NULL
  R<- NULL
  for(j in 1:p) { Z<-cbind(Z, scale(Y_mod[,j])) }
  for(j in 1:p) { R<-cbind(R, match(Y_mod[,j],sort(unique(Y_mod[,j])))) }
  Rlevels<-apply(R,2,max,na.rm=TRUE)
  Z[which(is.na(Y_mod),arr.ind = T)] = 0

  # Hyperparameters for DP mixture

  if(is.null(hyperparams)){
    ### dimension of latent factors
    k.star = floor(.7*p)

    ### Marginal threshooldd for re-sampling under RPL
    #marginal threshold for RL
    plugin.threshold = 350
    plugin.marginal = (apply(Y_mod, 2, function(x) {
      length(unique(x))
    }) > plugin.threshold)

    ### Cluster specific parameters
    #hyper parameteres for stick breaking process
    a_alpha = 1; b_alpha = 3;

    #hyperparameters for cluster parameters
    nu_mix = k.star + 2; kappa_0=.001

    ### MGP
    #local shrinkage
    nu = 3
    #global shrinkage for MGP
    a1 = 2; a2 = 3;

    ### Error Variance
    a.sigma = 1; b.sigma = 0.3

    ### Mixture
    #cluster covariance
    D_0<- diag(1,k.star)
    #precision parameter
    delta = diag(rep(10,k.star))



  }else{
    k.star = hyperparams$k.star

    plugin.threshold = hyperparams$plugin.threshold
    plugin.marginal = (apply(Y_mod, 2, function(x) {
      length(unique(x))
    }) > plugin.threshold)


    a_alpha = hyperparams$a_alpha; b_alpha =hyperparams$b_alpha;

    nu_mix = hyperparams$nu_mix; kappa_0=hyperparams$kappa_0

    nu = hyperparams$nu

    a1 = hyperparams$a1; a2 = hyperparams$a2;

    a.sigma = hyperparams$a.sigma; b.sigma = hyperparams$b.sigma

    D_0<- hyperparams$D_0

    delta = diag(rep(hyperparams$delta,k.star))


  }

  ###Initialize factors, factor loadings, error covariance, alpha_star under stick-breaking process
  # Use SVD for simplicity
  svd0 = svd(Z);

  # Factor loadings (p x k.star):
  Lambda = Lambda.s = as.matrix(svd0$v[,1:k.star])
  ##MGP initialization
  #Local and global precision parameters:
  phi.jh = 1/Lambda^2; tau.h = delta.h = rep(1, k.star)
  # Factors (n x k.star):
  eta = svd0$u[,1:k.star]%*%diag(svd0$d[1:k.star], k.star) #y%*%Lambda

  # Residuals (n x p):
  eps = Z - (tcrossprod(eta, Lambda))

  # Diagonal error variance (p-dimensional vector):
  Sigma.diag = apply(eps, 2, var)

  alpha = rep(0,p)

  # For stick-breaking
  alpha_star = 1

  #### mixture initialization
  # Storage for probability weights:
  prob_h =  array(0, c(n, H))

  # Naive initialization for labels:
  z = sample(1:H, n, replace = TRUE)

  # Initialize the other parameters, given labels:
  n_h = sapply(1:H, function(h) sum(z==h))
  V_h = rep(1, H)
  V_h[1:(H-1)] = sapply(1:(H-1), function(h) (1 + n_h[h])/(1 + n_h[h] + alpha_star + sum(n_h[-(1:h)])))
  pi_h = rep(1, H); pi_h[1] = V_h[1]; pi_h[-1] = sapply(2:H, function(h) V_h[h]*prod(1 - V_h[1:(h-1)]))

  #initialize eta cluster paramters
  mu = lapply(1:H, function(h) colMeans(eta[z==h,]))
  Delta = lapply(1:H, function(h) diag(1,k.star))

  #for arranging categorical variables
  factor_order = format_data(Y, factor.cols)
  col_mem = factor_order$col_mem
  cat_col_names = factor_order$cat_col_names
  bin_col_names = factor_order$bin_col_names
  count_col_names = unlist(sapply(colnames(Y[,-factor.cols]),function(x) if(all(Y[,x]%%1 == 0,na.rm = T)){return(x)}))
  cont_col_names = setdiff(colnames(Y[,-factor.cols]),count_col_names)

  #storage

  Deltas<- vector('list', H)
  mus<-vector('list',H)
  Sigmas<- array(0, c(nsamp- burn, p))
  Lambdas<- array(0, c(nsamp-burn, p, k.star))
  pis<- array(0, c(nsamp-burn,H))
  zs<-array(0, c(nsamp-burn,H))
  clustering<-array(0,c(nsamp-burn,n))
  quantiles= support = vector('list', p)
  Imputations<- vector('list', nImp)
  indexer = 1
  for(i in 1:H){
    Deltas[[i]] = array(0, c(nsamp-burn,k.star, k.star))
    mus[[i]] = array(0, c(nsamp-burn,k.star))
  }

  #for indexing imputations
  impevery = floor((nsamp -burn)/nImp)

  for(ns in  1:nsamp){



    # 1) Sample the cluster indicators
    for(h in 1:H) prob_h[,h] = pi_h[h]*dmvnorm(eta, mean = mu[[h]], sigma = Delta[[h]])


    z = apply(prob_h, 1, function(p_i) sample(1:H, 1, prob = p_i))


    # Update the cluster counts:
    n_h = sapply(1:H, function(h) sum(z==h))


    # 2) Sample the probability weights:
    V_h[1:(H-1)] = sapply(1:(H-1), function(h)
      rbeta(n = 1,
            shape1 = 1 + n_h[h],
            shape2 = alpha_star + sum(n_h[-(1:h)]))
    )


    # And update the probabilities of each cluster:
    pi_h[1] = V_h[1];

    pi_h[-1] = sapply(2:H, function(h)
      V_h[h]*prod(1 - V_h[1:(h-1)])
    )

    # 3) Update Cluster Specific Components
    for(h in 1:H){
      #update u_j
      ind<- which(z == h)
      n_i = n_h[h]
      eta_h = eta[ind,]
      # mean eta, S
      if(n_i ==0){ # if there is no membership, sample from prior
        Delta[[h]] = rinvwishart(nu_mix,delta^2%*%D_0)
        mu[[h]] = mvrnorm(1, mu = rep(0,k.star), Sigma = Delta[[h]]/kappa_0)
      }else{
        if(length(eta_h)>dim(eta)[2]){
          eta_bar = colMeans(eta_h)
          S = crossprod(t(apply(eta_h,1,function(x) x- eta_bar)))
        }else{
          eta_bar = eta_h
          S = array(0, c(k.star,k.star))
        }

        cross_eta_bar = tcrossprod(eta_bar)

        #IW updates
        nu_h = nu_mix + n_i
        psi_h =(delta^2%*%D_0 + S + (kappa_0*n_i)/(kappa_0 + n_i)*cross_eta_bar)

        Delta[[h]] = rinvwishart(nu_h, psi_h)

        #mean updates
        alpha_h = (n_i/(kappa_0 +n_i))*eta_bar
        kappa_h = kappa_0 + n_i

        mu[[h]] = mvrnorm(1, mu = alpha_h, Sigma = Delta[[h]]/kappa_h)

      }

    }

    # 4) sample alpha_star
    alpha_star = rgamma(n = 1,
                        shape = a_alpha + H - 1,
                        rate = b_alpha - sum(log(1-V_h[-H])))

    # 5)  sample the factors

    for(h in 1:H){
      n_i = n_h[h]
      mu_h = mu[[h]]
      Delta_h = Delta[[h]]
      inds = which(z == h)
      Z_sub = Z[inds,]

      if(n_i >1){
        chQeta = chol(solve(Delta_h)+ crossprod(Lambda, diag(1/Sigma.diag))%*%Lambda)
        leta = apply(t(tcrossprod(crossprod(Lambda, diag(1/Sigma.diag)), Z_sub-alpha)),1,function(x) x + (crossprod(mu_h, solve(Delta_h))))
        eta[inds,] = t(backsolve(chQeta,forwardsolve(t(chQeta), leta) + rnorm(n_i*k.star))) #for(i in 1:n) eta[i,]= backsolve(chQeta,forwardsolve(t(chQeta), leta[,i]) + rnorm(k.star))
      }else if(n_i == 1){
        chQeta = chol(solve(Delta_h)+ crossprod(Lambda, diag(1/Sigma.diag))%*%Lambda)
        leta = apply(tcrossprod(crossprod(Lambda, diag(1/Sigma.diag)), t(Z_sub-alpha)),2,function(x) x + t(crossprod(mu_h, solve(Delta_h))))
        eta[inds,] = t(backsolve(chQeta,forwardsolve(t(chQeta), leta) + rnorm(n_i*k.star))) #for(i in 1:n) eta[i,]= backsolve(chQeta,forwardsolve(t(chQeta), leta[,i]) + rnorm(k.star))
      }
    }

    # 6) Sample lambda_jh
    cp.eta = crossprod(eta)
    for(j in 1:p){
      chQj = chol(diag(phi.jh[j,]*tau.h, k.star) + cp.eta/Sigma.diag[j])
      lj = crossprod(eta, Z[,j] - alpha[j])/Sigma.diag[j]
      Lambda[j,] = backsolve(chQj,forwardsolve(t(chQj), lj) + rnorm(k.star))
    }

    # 7) Sample the error variances
    eps =  Z - (tcrossprod(eta, Lambda)+alpha)

    Sigma.diag = apply(eps, 2, function(x) 1/rgamma(n = 1, shape = a.sigma + n/2,
                                                    rate = b.sigma + 1/2*sum(x^2)))


    # 8) sample phi.jh
    phi.jh = matrix(rgamma(n = p*k.star, shape = (nu + 1)/2,
                           rate = (nu + Lambda^2*matrix(rep(tau.h, each = p), nr = p))/2), nr = p) #for(h in 1:k.star){for(j in 1:p) phi.jh[j,h] = rgamma(n = 1, shape = (nu + 1)/2, rate = (nu + Lambda[j,h]^2*tau.h[h])/2)

    # 9) sample tau.h via delta.h
    delta.h = sampleMGP(theta.jh = sqrt(phi.jh)*Lambda, delta.h = delta.h, a1 = a1, a2 = a2)
    tau.h = cumprod(delta.h)


    # 10) impute categorical data for RPL

    Y_mod_dummy = array(0, c(n,p))
    colnames(Y_mod_dummy) = colnames(Y_mod)
    for(col in unique(c(cat_col_names,bin_col_names))){
      inds = which(col_mem == col) # for deciding whether the variable has more than one lev
      ir_NA = (1:n)[is.na(Y[,col])] # for finding null values

      if(length(ir_NA)>0  ){

        cop_mean=  alpha+ t(Lambda%*%t(eta))
        sds<-sqrt(Sigma.diag[inds])
        #
        if(length(inds)>1){ #categorical variable

          DO_center =(cop_mean[ir_NA,inds])

          #compute categorical probs
          probs = t(apply(DO_center,1,function(x)compute_probs(x, sds))) +.000001
          probs = probs/rowSums(probs)

          #impute the category as the highest probability
          impute_cat<- apply(probs, 1, which.max)

          replace<- t(sapply(impute_cat,function(x) DO_sample(x,inds)))
          Y_mod_dummy[ir_NA, inds] = replace



        }

      }

    }


    # 11) Re-sample Z
    for(j in 1:p){

      muj = alpha[j] + Lambda[j,]%*%t(eta) # univariate moments
      sdj = sqrt(Sigma.diag[j])
      ir_NA <- (1:n)[is.na(Y_mod[,j])]

      if(is_cat_bin[j] == 0){
        # sample missing Z missing values
        Z[ir_NA, j] <- rnorm(length(ir_NA), muj[ir_NA], sdj)

        if(!plugin.marginal[j]){ # re-sample latent Z according to RPL
          for(r in 1:Rlevels[j]){

            ir<-(1:n)[R[,j] == r & !is.na(R[,j])]
            lb<- suppressWarnings(max(Z[R[,j] == r-1,j], na.rm = T))
            ub <- suppressWarnings(min(Z[R[,j] == r+1,j], na.rm = T))
            Z[ir,j] = qnorm(runif(length(ir), pnorm(lb, muj[ir],sdj),pnorm(ub,muj[ir],sdj)),muj[ir],sdj)
          }
        }
      }
      else{ #categorical imputation/re-sampling
        ir = (1:n)[!is.na(Y_mod[,j])]
        Z[ir[which(Y_mod[ir,j] == 1)],j] =   rtruncnorm(n = 1,a =0 ,mean = muj[ir[which(Y_mod[ir,j] == 1)]],sd = sdj)
        Z[ir[which(Y_mod[ir,j] == 0)], j] = rtruncnorm(n = 1,b = 0,mean = muj[ir[which(Y_mod[ir,j] == 0)]],sd = sdj)

        if(!(colnames(Y_mod)[j]%in% bin_col_names)){
          if(length(muj[ir_NA[which(Y_mod_dummy[ir_NA,j] == 1)]])>0){
            Z[ir_NA[which(Y_mod_dummy[ir_NA,j] == 1)],j] =   rtruncnorm(n = 1,a =0 ,mean = muj[ir_NA[which(Y_mod_dummy[ir_NA,j] == 1)]],sd = sdj)

          }
          if(length( Z[ir_NA[which(Y_mod_dummy[ir_NA,j] == 0)], j]) >0){
            Z[ir_NA[which(Y_mod_dummy[ir_NA,j] == 0)],j] =   rtruncnorm(n = 1,b =0 ,mean = muj[ir_NA[which(Y_mod_dummy[ir_NA,j] == 0)]],sd = sdj)

          }
        }else{
          Z[ir_NA, j] <- rnorm(length(ir_NA), muj[ir_NA], sdj)
        }

      }
    }

    #progress
    if(ns%%100 == 0 ){
      cat(round(100 * ns/nsamp), "percent done ", date(),
          "\n")
      print(paste("Current number of clusters/occupancy:"))
      print(table(z))


    }


    if(ns> burn){
      #compute margins
      Fns = vector("list",p)
      for( x in 1:p){

        if(length(unique(Y_mod[,x]))>4){
          bounds<-get_bound(R[,x], Rlevels[x], Z[,x])
          wts<-compute_mixprobs(bounds,Y_mod[,x],pi_h,z,Delta, mu, alpha = alpha,Lambda, Sigma.diag, x)
          if(!plugin.marginal[x] ){

            if(all(Y_mod[!is.na(Y_mod[,x]),x]%%1 == 0)){
              if(min(Y_mod[,x], na.rm = T)>=0)
                Fns[[x]] = fill_in_gaps_cdf(max(min(Y_mod[,x], na.rm = T),0),max(Y_mod[,x], na.rm = T),wts)$cdf
            }else{
              Fns[[x]] = fill_in_gaps_cdf(min(Y_mod[,x], na.rm = T)-5,max(Y_mod[,x], na.rm = T),wts)$cdf
            }
          }else{
            Fns[[x]] = wts
          }
          quantiles[[x]] = rbind(quantiles[[x]], Fns[[x]][,3]) #save margin estimates
          support[[x]] = rbind(support[[x]], Fns[[x]][,1])
        }

      }


      for(i in 1:H){
        Deltas[[i]][ns-burn,,] = Delta[[i]]
        mus[[i]][ns-burn,] = mu[[i]]
      }
      Sigmas[ns-burn,] = Sigma.diag
      Lambdas[ns-burn,,] = Lambda

      pis[ns-burn,] = pi_h
      occupied = as.numeric(names(table(z)))
      zs[ns-burn,occupied] = 1
      clustering[ns-burn,] = z


      if(ns%%impevery == 0 & Impute == T){
        Y_impute <- Y

        for(col in unique(col_mem)){
          inds = which(col_mem == col)
          ir_NA = (1:n)[is.na(Y[,col])]

          if(length(ir_NA)>0  ){

            if(length(inds)>1){ #categorical variable

              Y_impute[ir_NA,col] = levels(Y[,col])[apply(Z[ir_NA,inds],1,which.max)]

            }

            else if(col %in% bin_col_names){

              Y_impute[ir_NA,col] = ifelse(Z[ir_NA,inds]<0,levels(Y[,col])[1], levels(Y[,col])[2])

            }
            else if (col %in% count_col_names){

              if(min(Y_mod[which(!is.na(Y_mod[,inds])),inds], na.rm = T)>= 0){
                #if values are strictly positive

                Y_impute[ir_NA,col] = inverse_cdf(Fns[[inds]][,3],probs = compute_mixprobs1(Z[ir_NA,inds],
                                                                                            unique(z),
                                                                                            pi_h = pi_h,
                                                                                            Delta = Delta,
                                                                                            mu = mu,
                                                                                            alpha,
                                                                                            Lambda,
                                                                                            Sigma.diag,
                                                                                            inds),

                                                  min = Fns[[inds]][2,1], max =max(Fns[[inds]][,1]) )


              }
              else{

                Y_impute[ir_NA,col] = inverse_cdf(Fns[[inds]][,3],probs =compute_mixprobs1(Z[ir_NA,inds],
                                                                                           unique(z),
                                                                                           pi_h = pi_h,
                                                                                           Delta = Delta,
                                                                                           mu = mu,
                                                                                           alpha,
                                                                                           Lambda,
                                                                                           Sigma.diag,
                                                                                           inds),
                                                  min = min(Fns[[inds]][,1])+1, max =max(Fns[[inds]][,1]) )

              }

            }
            else{

              Fn_inv = splinefun(Fns[[inds]][,3],
                                 Fns[[inds]][,1],
                                 method = c("monoH.FC"))


              Y_impute[ir_NA,col] = Fn_inv(compute_mixprobs1(Z[ir_NA,inds],
                                                             unique(z),
                                                             pi_h = pi_h,
                                                             Delta = Delta,
                                                             mu = mu,
                                                             alpha = alpha,
                                                             Lambda,
                                                             Sigma.diag, inds))

            }

          }
        }
        Imputations[[indexer]] = Y_impute
        indexer = indexer +1

      }

    }
  }

  return(list(Deltas = Deltas,
              Imputations = Imputations,
              mus = mus,
              pis = pis,
              zs = zs,
              clustering = clustering,
              Sigmas = Sigmas,
              Lambdas = Lambdas,
              Quantiles = quantiles,
              Support = support,
              col_mem = col_mem,
              cat_col_names = cat_col_names,
              bin_col_names = bin_col_names,
              count_col_names = count_col_names,
              cont_col_names = cont_col_names,
              dat = Y,
              Y_aug = Y_mod,
              seed = seed))


}



#' Generate a posterior predictive data set
#'
#' Utilize posterior samples of Gaussian mixture copula parameters to generate a posterior
#' predictive data set for missing data inference. Each posterior predictive data set is
#' constructed with one radomly selected sample from the posterior of GMC parameters.
#'
#' @param \code{mcmcobj}: fitted GMC mcmc object from the function \code{GMC.mcmc}
#' @param \code{nobs}: number of observations in the posterior predictive data set, default is \code{nrow(data)}
#' @param \code{nsets}: number of posterior predictive data sets to create
#' @param \code{seed}: seed for replication
#'
#'
#' @return
#' \itemize{
#' \item \code{Y_pred}: A list of length \code{nsets} containing posterior predictive data sets with \code{nobs} observations.
#' These data sets are constructed \emph{with} the margin adjustment
#' \item \code{Y_pred_noMA}:  A list of length \code{nsets} containing posterior predictive data sets with \code{nobs} observations.
#' These data sets are constructed \emph{without} the margin adjustment. In this case empirical distribution functions are used
#' for count variables, and smoothed empirical distributions are used for continuous variables.
#' \item \code{seed}: seed for replication
#' }
#' @import tmvtnorm
#' @import mvtnorm
#' @import LaplacesDemon
#' @import truncnorm
#' @import purrr
#' @import MASS
#' @import dplyr
#' @import tidyr
#' @import data.table
#' @export
#'
#' @examples

#' mcmc<-GMC.mcmc(Data = X, nsamp = 1000)
#'pred<- get_predictive_Y(mcmc,
#'                        nsets = 1,
#'                        nobs = dim(X)[1],
#'                        seed = 2)
#' plot(pred$Y_pred[[1]])
#'
get_predictive_Y<-function(mcmc,
                           nobs,
                           nsets,
                           seed = NULL){

  if(is.null(seed)){
    seed =sample(1:10000000,1)
    set.seed(seed)
  }else{
    set.seed(seed)
  }

  Y = mcmc$dat #for data formatting
  Y_aug = mcmc$Y_aug
  p = dim(Y_aug)[2]
  H = dim(mcmc$zs)[2] # upper bound for number of clusters in mixtures
  n = nobs
  nsamps = dim(mcmc$Sigmas)[1] #number of posterior samples from mcmc fit
  if(nsamps < nsets){
    stop('Number of posterior predictive data sets exceeds number of posterior samples')
  }
  postsamps = sample(1:nsamps,size = nsets, replace = F)

  pred_datasets = pred_datasets.noMA = vector('list', nsets)

  for(i in 1:nsets){

    # get posterior sample of mixture model parameters
    mu = lapply(1:H, function(x) return(mcmc$mus[[x]][postsamps[i],]))
    alpha = rep(0,p)
    Delta = lapply(1:H, function(x) return(mcmc$Deltas[[x]][postsamps[i],,]))
    Lambda = mcmc$Lambdas[postsamps[i],,]
    pi_h = mcmc$pis[postsamps[i],]
    Sigma.diag= mcmc$Sigmas[postsamps[i],]
    z = mcmc$zs[postsamps[i],]
    #column names and memberships
    col_mem = mcmc$col_mem
    cat_col_names = mcmc$cat_col_names
    bin_col_names = mcmc$bin_col_names
    count_col_names = mcmc$count_col_names
    cont_col_names = mcmc$con_col_names

    Fns = vector('list', p) #format marginal distributions
    for(j in 1:p){
      support = mcmc$Support[[j]][postsamps[i],]
      qs = mcmc$Quantiles[[j]][postsamps[i],]
      Fj = cbind(support,rep(0,length(support)),qs)
      Fns[[j]] = Fj

  }



  unique_z = which(z == 1)

  if(length(unique(z))>1){
    cs<- sample(unique_z, n, replace = T, prob =  pi_h[unique_z])
    etas = NULL
    #need to fix this

    etas = do.call(rbind, t(sapply(unique(cs),function(x) return(rbind(etas,
                                                                       mvrnorm(sum(cs == x),mu[[x]],
                                                                               Delta[[x]]))))))


  }else{
    etas = mvrnorm(n,mu[[unique_z]], Delta[[unique_z]])
  }

  means = t(sapply(1:n,function(x) alpha + Lambda%*%(etas[x,])))

  Y_pred = Y_pred_noMA = data.frame(array(0,c(n,dim(Y)[2])))

  colnames(Y_pred) = colnames(Y_pred_noMA) = colnames(Y)

  for(col in unique(col_mem)){

    inds = which(col_mem == col)
    ir_NA = (1:n)[is.na(Y[,col])]

    if(col %in% cat_col_names){ #get predictive category by sampling from diagonal orthant probabilities

      sds = sqrt(Sigma.diag[inds])

      DO_center =means[,inds]

      probs = t(apply(DO_center,1,function(x)compute_probs(x, sds)))+.0000001
      probs = probs/rowSums(probs)

      impute_cat<- apply(probs, 1,function(x) sample(1:length(inds),1,replace = T,x))
      preds = factor(levels(Y[,col])[impute_cat], levels = levels(Y[,col]))
      Y_pred[,col] = preds
      Y_pred_noMA[,col]= preds

    }

    else if(col %in% bin_col_names){ #binary predictive sampling; whether or not the latent simulation is above zero
      latent_z<- rnorm(n,means[,inds],sqrt(Sigma.diag[inds]))
      pred = factor(ifelse(latent_z<0,levels(Y[,col])[1], levels(Y[,col])[2]), levels = levels(Y[,col]))
      Y_pred[,col] = pred
      Y_pred_noMA[,col] = pred
    }
    else if (col %in% count_col_names){

      latent_z<- rnorm(n,means[,inds],sqrt(Sigma.diag[inds]))

      Y_pred[,col] = inverse_cdf(Fns[[inds]][,3],probs = compute_mixprobs1(latent_z,
                                                                           unique_z,
                                                                           pi_h = pi_h,
                                                                           Delta = Delta,
                                                                           mu = mu,
                                                                           alpha = alpha,
                                                                           Lambda,
                                                                           Sigma.diag, inds),
                                 min = Fns[[inds]][2,1], max =max(Fns[[inds]][,1]))
      Y_pred_noMA[,col] = quantile(Y[,col],
                                   probs = (n/(n+1))*compute_mixprobs1(latent_z,
                                                                       unique_z,
                                                                       pi_h = pi_h,
                                                                       Delta = Delta,
                                                                       mu = mu,
                                                                       alpha = alpha,
                                                                       Lambda,
                                                                       Sigma.diag, inds),
                                   type = 1,
                                   na.rm= T)



    }


    else{ #continuous
      latent_z<- rnorm(n,means[,inds],sqrt(Sigma.diag[inds]))

      Fn_inv = splinefun(Fns[[inds]][,3],
                         Fns[[inds]][,1],
                         method = c("monoH.FC"))
      #smoothed empirical estimate for noMA
      range = sort(unique(Y[,col]))
      Fn_inv_noMA = splinefun(ecdf(Y[,col])(range),range, method = c("monoH.FC"))


      Y_pred[,col]= Fn_inv(compute_mixprobs1(latent_z,
                                             unique_z,
                                             pi_h = pi_h,
                                             Delta = Delta,
                                             mu = mu,
                                             alpha = alpha,
                                             Lambda,
                                             Sigma.diag,
                                             inds))
      Y_pred_noMA[,col] = Fn_inv_noMA(compute_mixprobs1(latent_z,
                                                        unique_z,
                                                        pi_h = pi_h,
                                                        Delta = Delta,
                                                        mu = mu,
                                                        alpha = alpha,
                                                        Lambda,
                                                        Sigma.diag,
                                                        inds))


      }
    }

  pred_datasets[[i]] = Y_pred
  pred_datasets.noMA[[i]] = Y_pred_noMA

  }



  return(list(Y_pred = pred_datasets, Y_pred_noMA = pred_datasets.noMA, seed = seed))

}

