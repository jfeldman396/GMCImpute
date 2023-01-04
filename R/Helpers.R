## Description: formats data and returns the column names. Important for multiple imputation and predictive sampling

## Inputs:
# Y: data frame
# factor.cols: which columns are factors

## Outputs:
# col_mem: column names and which categorical variable it belongs to
# bin_col_names: binary column names
# cat_col_names: categorical (nlevels>2) column names

#' Title
#'
#' @param Y
#' @param factor.cols
#'
#' @return
#'
#' @examples
format_data<-function(Y,factor.cols){


  cat<- NULL
  cat_col = 1
  cat_col_names = NULL
  bin_col_names = NULL
  for(i in factor.cols){
    if(nlevels(Y[,i]) >2){
      one.hot.cat<- model.matrix(~. -1,data =model.frame(~.,data = Y[,i], na.action = na.pass))
      cat<- cbind(cat,one.hot.cat)
      cat_col_names = c(cat_col_names, rep(names(Y)[i], nlevels(Y[,i])))
    }
    else{
      one.hot.bin<- model.matrix(~. , data = model.frame(~.,data = Y[,i], na.action = na.pass))[,2]
      cat<-cbind(cat,one.hot.bin)
      bin_col_names = c(bin_col_names, rep(names(Y)[i],1))
    }

  }

  one.hot.cols<- ncol(cat) # number of one hot encoded categories

  Y_mod<- cbind(cat,Y[,-factor.cols])
  col_mem = c(cat_col_names,bin_col_names,colnames(Y_mod)[(one.hot.cols+1):ncol(Y_mod)])# data frame used for sampler
  #Y_mod_comp = Y_mod %>% tidyr::drop_na() #complete observations
  return(list( cat_col_names = cat_col_names,
               bin_col_names = bin_col_names, col_mem = col_mem))
}

## Description: Samples MGP parameters for factor model

## Inputs:
# theta.jh, delta.h: local and global shrikage parameters

## Outputs
# new sample of delta.h

#' Title
#'
#' @param theta.jh
#' @param delta.h
#' @param a1
#' @param a2
#'
#' @return

#'
#' @examples
sampleMGP<-function(theta.jh, delta.h, a1 = 3, a2 = 5){
  # Store the dimensions locally
  p = nrow(theta.jh); K = ncol(theta.jh)

  # Sum over the (squared) replicates:
  sum.theta.l = colSums(theta.jh^2)

  # h = 1 case is separate:
  tau.not.1 = cumprod(delta.h)/delta.h[1]
  delta.h[1] = rgamma(n = 1, shape = a1 + p*K/2,
                      rate = 1 + 1/2*sum(tau.not.1*sum.theta.l))
  # h > 1:
  if(K > 1){for(h in 2:K){
    tau.not.h = cumprod(delta.h)/delta.h[h]
    delta.h[h] = rgamma(n = 1, shape = a2 + p/2*(K - h + 1),
                        rate = 1 + 1/2*sum(tau.not.h[h:K]*sum.theta.l[h:K]))
  }}
  delta.h #list(tau.h = cumprod(delta.h), delta.h = delta.h)
}

## Description: Returns a predictive categorical observation with an augmented one-hot encoding

##Inputs:
# inds: the column numbers in the data frame corresponding to levels of categorical vbl
# max: index of 'hot' level

##Outputs:
# one row of predictive one-hot categorical vbl

#' Title
#'
#' @param max
#' @param inds
#'
#' @return

#'
#' @examples
DO_sample<-function(max,inds){
  #create one-hot-encoding for imputed data
  sample<- rep(0,length(inds))

  #replace predictive category with 1
  sample[max] = 1
  return(sample)

}


## Description: Compute diagonal orthant probabilities for categorical imputation

##Inputs:
# means: observation specific mean, should be n_mis x nlevels(j)
# sds: observation specific sd should be n_mis X 1

##Outputs:
# probabilities for each level of the categorical variable in question
#' Title
#'
#' @param mean
#' @param sds
#'
#' @return

#'
#' @examples
compute_probs = function(mean,sds){

  return(sapply(1:length(mean),function(x) (1-pnorm(0,mean[x], sds[x]))*prod(pnorm(0,mean[-x], sds[-x]))))
}

## Description: get margin adjustment bounds for margin estimation

##Inputs:
# lvls: the level of each observation, should be (n_obs x p)
# num_lvls: number of unique levels for each column
# Zs: Z_j

##Outputs:
# bounds on Z_j for margin adjustment calculation

#' Title
#'
#' @param lvls
#' @param num_lvls
#' @param Zs
#'
#' @return

#'
#' @examples
get_bound<- function(lvls,num_lvls, Zs){
  bounds = c(-Inf,max(Zs[which(lvls == 1)]))
  bounds = rbind(bounds,t(sapply(2:(num_lvls-1), function(x) c(max(Zs[which(lvls == x-1)]), max(Zs[which(lvls == x)])))))
  bounds = rbind(bounds, c(max(Zs[which(lvls == max(num_lvls) - 1)]),Inf))
  return(bounds)
}

## Description: Computes cdf estimate for each margin, to be used within the MCMC

##Inputs:
# bounds: bounds for computation of cumulative probabilities
# Y: Y_j, for inclusion of support
# pi_h: mixture weights
# z: mixture membership
# Delta: list of length H with cluster specific parameters
# mu: list of length H with cluster specific parameters
# alpha: array of length 0, keep this way
# Lambda: factor loadings
# Sigma: error covariance
# inds: indexer for column

##Outputs:
# cdf: point-wise estimate for cdf of Y_j

#' Title
#'
#' @param bounds
#' @param Y
#' @param pi_h
#' @param z
#' @param Delta
#' @param mu
#' @param alpha
#' @param Lambda
#' @param Sigma
#' @param inds
#'
#' @return

#'
#' @examples
compute_mixprobs<- function(bounds,Y,pi_h,z,Delta, mu,alpha,Lambda, Sigma, inds){

  #cluster specific marginal means
  means<-t(sapply(unique(z),function(x) alpha +Lambda %*% mu[[x]]))

  #cluster specific marginal sds
  sds<- sqrt(t(sapply(unique(z) ,function(x) diag(tcrossprod(Lambda%*%Delta[[x]],Lambda) + diag(Sigma)))))

  pis = pi_h[unique(z)]
  pis = pis/sum(pis)
  #compute cumulative probabilities using predictive bounds
  cdf<- cbind(sort(unique(Y[!is.na(Y)])),
              apply(bounds,1,function(x) sum(pis*(pnorm(x[2], means[,inds],sds[,inds]) - pnorm(x[1], means[,inds],sds[,inds])))),
              cumsum(apply(bounds,1,function(x) sum(pis*(pnorm(x[2], means[,inds],sds[,inds]) - pnorm(x[1], means[,inds],sds[,inds]))))))


  return(cdf)
}

## Description: Computes cumulative probability for transforming latent variables via the copula link, to be used in immputation or predictive sampling

##Inputs:
# bounds: bounds for computation of cumulative probabilities
# Zs: Z_j,
# unique_z: mixture memberships, which of the H clusters are populated
# pi_h: mixture weights
# Delta: list of length H with cluster specific parameters
# mu: list of length H with cluster specific parameters
# alpha: array of length 0, keep this way
# Lambda: factor loadings
# Sigma: error covariance
# inds: indexer for column

##Outputs:
# \psi_j(Z_{j})
#' Title
#'
#' @param Zs
#' @param unique_z
#' @param pi_h
#' @param Delta
#' @param mu
#' @param alpha
#' @param Lambda
#' @param Sigma
#' @param inds
#'
#' @return

#'
#' @examples
compute_mixprobs1<- function(Zs,unique_z,pi_h,Delta, mu, alpha,Lambda, Sigma, inds){
  means<-t(sapply(unique_z ,function(x) alpha + Lambda %*% mu[[x]]))
  #means = matrix(0,nrow = 2311, ncol = p)
  sds<- sqrt(t(sapply(unique_z, function(x) diag(tcrossprod(Lambda%*%Delta[[x]],Lambda) + diag(Sigma)))))


  return(sapply(Zs,function(x)sum(pi_h[unique_z]*(pnorm(x, means[,inds],sds[,inds]) ))))



}

## Description: fills in gaps in the support of observed Y_{j} if Y_{j} is discrete

##Inputs:
# min: smallest value in the support of Y_{j}
# max: largest value in the support of Y_{j}
# cdf: cdf estimate from compute_mixprobs

##Outputs:
# Adjusted with mass approximated at missing points in support for discrete variables cdf (full_cdf)
#' Title
#'
#' @param min
#' @param max
#' @param cdf
#'
#' @return
#'
#' @examples
fill_in_gaps_cdf<-function(min,max, cdf){
  cdf = data.frame(cdf)

  colnames(cdf) = c('val','prob', 'cum_prob')
  obs_val = cdf$val
  support = data.frame(val = seq(min,max,by = 1))

  full_cdf = full_join(support,cdf, by = c('val' = 'val'))
  #find the gaps, assign 1/n probability
  num_gaps = sum(is.na(full_cdf$prob))
  loc_gaps = which(is.na(full_cdf$prob))

  full_cdf$cum_prob[which(full_cdf$cum_prob == max(full_cdf$cum_prob, na.rm = T))] =NA

  full_cdf$cum_prob[dim(full_cdf)[1]] = 1

  full_cdf = rbind(c(-1,0,0), full_cdf)

  spline = splinefun(full_cdf$val, full_cdf$cum_prob, method = 'monoH.FC')
  full_cdf$cum_prob[which(is.na(full_cdf$cum_prob))] = spline(full_cdf$val[which(is.na(full_cdf$cum_prob))])

  full_cdf$prob[2:length(full_cdf$prob)] = diff(full_cdf$cum_prob)
  return(list(cdf = full_cdf, num_gaps = num_gaps))
}

## Description: Inverse link function from latent to observed, given margin estimates for multiple imputation or predictive sampling


## Inputs:
# cdf: cdf estimate from compute_mixprobs or fill_in_gaps_cdf, depending on whether Y_{j} is discrete or continuous
# probs: cumulative probabilities under the mixture from compute_mixprobs1
# min: lower boundary of support for Y_{j}
# max: upper boundary of support for Y_{j}
#' Title
#'
#' @param cdf
#' @param probs
#' @param min
#' @param max
#'
#' @return

#'
#' @examples
inverse_cdf<- function(cdf, probs, min, max){
  range = min:max
  return(range[cut(probs, breaks = cdf, labels = min:max)])
}




