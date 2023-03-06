#' @param n sample size
#' @param m number of traits
#' @param p number of variants
#' @param L levels
#' @param varY variance of Y
#' @param scaled_prior_variance prior variance relative to V[Y]
#'
#' @return optimization_state
#'
susie_initialize <- function(n, m, p, L, varY,
                             scaled_prior_variance = .1){

    .zero <- function(d1, d2) { matrix(0, nrow=d1, ncol=d2) }

    .var0 <- function(l) { scaled_prior_variance }

    list(alpha        = matrix(1/p, nrow=L, ncol=p),
         mu           = lapply(1:L, function(l) .zero(p,m)),
         mu2          = lapply(1:L, function(l) .zero(p,m)),
         lbf_list     = lapply(1:L, function(l) .zero(p,m)),
         KL           = rep(as.numeric(NA), L),
         lbf          = rep(as.numeric(NA), L),
         lbf_variable = matrix(as.numeric(NA), L, p),
         V            = sapply(1:L, .var0),
         fitted       = matrix(0, nrow=n, ncol=m),
         L            = L,
         m            = m,
         p            = p)
}


#' @param xx
#' @param yy
#' @param .rescale
#' @return t(xx) %*% yy
safe.xty <- function(xx, yy, .rescale = FALSE) {

    x <- xx; x[is.na(xx)] <- 0
    y <- yy; y[is.na(yy)] <- 0

    if(.rescale){
        n.obs <- crossprod(!is.na(xx), !is.na(yy))
        return(crossprod(x, y) / n.obs * nrow(xx))
    } else {
        return(crossprod(x, y))
    }
}


#' @param xx
#' @param .rescale
#' @return apply(xx, 2, sum)
safe.colsum <- function(xx, .rescale = FALSE) {
    if(.rescale){
        .sum <- function(x) {
            nn <- sum(!is.na(x))
            sum(x, na.rm=TRUE) / nn * length(x)
        }
        return(apply(xx, 2, .sum))
    } else {
        return(apply(xx, 2, sum, na.rm=TRUE))
    }
}


#' @rdname single_effect_shared
#'
#' @title Bayesian single-effect linear regression
#'
#' @param Y An n by m output matrix
#' @param X An n by p design matrix
#' @param V A scalar giving the (initial) prior variance
#' @param residual_variance The residual variance.
#'
#' @return A list with the following elements:
#'
#' \item{alpha}{Posterior probability of variant across k traits; \code{alpha[j]}}
#' \item{mu}{Matrix of posterior means; E[b[j,k] | z[j] = 1]}
#' \item{mu2}{Matrix of posterior second moments; E[b[j,k]^2 | z[j] = 1]}
#' \item{lbf}{Vector of log-Bayes factors for each variable.}
#' \item{lbf_model}{Log-Bayes factor for the single effect regression.}
#' \item{V}{Prior variance by the EM method.}
#' \item{loglik}{The log-likelihood, \eqn{\log p(y | X, V)}.}
#'
#' @importFrom stats dnorm
#'
single_effect_shared <- function(Y, X, V, residual_variance = NULL){

    if(is.null(residual_variance)){
        residual_variance <- apply(Y, 2, var, na.rm=TRUE)
    }

    ## 1. Compute single variant effects
    XtY <- safe.xty(X, Y)
    x2 <- safe.colsum(X * X)

    betahat <- sweep(XtY, 1, x2, `/`)
    shat2 <- matrix(1/x2, ncol=1) %*% residual_variance

    ## 2. log(bf) for each SNP and trait
    lbf.mat <- (
        dnorm(betahat, 0, sqrt(shat2 + V), log=TRUE) -
        dnorm(betahat, 0, sqrt(shat2), log=TRUE)
    )

    ## 3. Combine LBF across the columns (outputs)
    lbf.mat[is.infinite(shat2)] <- -Inf

    .mat <- lbf.mat
    .mat[!is.finite(lbf.mat)] <- NA
    lbf <- apply(.mat, 1, sum, na.rm=TRUE)

    maxlbf <- max(lbf)
    w <- exp(lbf - maxlbf)
    w_sum <- sum(w)
    alpha <- w/w_sum

    post_var <- (1/shat2 + 1/V)^(-1)
    post_mean <- XtY * sweep(post_var, 2, residual_variance, `/`)
    post_mean2 <- post_var + post_mean^2

    lbf_model <- w_sum + maxlbf
    stuff <- sum(dnorm(Y, 0, sqrt(residual_variance), log=TRUE), na.rm=TRUE)
    loglik <- lbf_model + stuff

    ## prior variance by MLE
    V <- sum(sweep(post_mean2, 1, alpha, `*`), na.rm=TRUE)

    list(alpha = alpha,
         mu = post_mean,
         mu2 = post_mean2,
         lbf = lbf,
         lbf.mat = lbf.mat,
         lbf_model = lbf_model,
         V = V,
         loglik = loglik)
}

#' @param X an n by p design matrix
#' @param Y an n by m output matrix
#' @param vR the residual variance
#' @param alpha the posterior inclusion (p x 1)
#' @param mu the posterior mean of b (p x m)
#' @param mu2 the posterior second moment of b (p x m)
posterior_loglik <- function(X, Y, vR, alpha, mu, mu2){
    n <- nrow(X)
    hat <- safe.xty(t(X), sweep(mu, 1, alpha, `*`))
    hat2 <- safe.xty(t(X * X), sweep(mu2, 1, alpha, `*`))
    stuff <- apply(Y*Y - 2*Y*hat + hat2, 2, sum, na.rm=TRUE)
    -0.5 * n * sum(log(2*pi*vR), na.rm=TRUE) - 0.5 * sum(stuff / vR, na.rm=TRUE)
}

#' @title Estimate single effect regression with shared variable selection
#' 
#' @param X an n by p design matrix
#' @param Y an n by m output matrix
#' @param state variational inference state
update_shared_effect <- function(X, Y, state){

    ## update residual variance
    llik <- 0
    for(l in 1:state$L){

        ## Remove l-th effect from fitted values
        theta.l <- sweep(state$mu[[l]], 1, state$alpha[l, ], `*`)
        state$fitted <- state$fitted - safe.xty(t(X), theta.l)

        ## Take partial residuals
        ## for the l-th model
        R <- Y - state$fitted
        v0 <- state$V[l]

        ## Update by shared single effect regression
        vR <- apply(R, 2, var, na.rm=TRUE)
        res <- single_effect_shared(R, X, v0, vR)

        state$alpha[l, ] <- res$alpha
        state$mu[[l]] <- res$mu
        state$mu2[[l]] <- res$mu2
        state$lbf_list[[l]] <- res$lbf.mat
        state$V[l] <- res$V
        state$lbf[l] <- res$lbf_model
        state$lbf_variable[l, ] <- res$lbf

        theta.l <- sweep(state$mu[[l]], 1, state$alpha[l, ], `*`)
        state$fitted <- state$fitted + safe.xty(t(X), theta.l)

        post.llik <- posterior_loglik(X,R,vR,
                                      res$alpha,
                                      res$mu,
                                      res$mu2)

        state$KL[l] <- (post.llik -res$loglik)

        llik <- llik + post.llik
    }
    state$llik <- llik
    return(state)
}

#' @title Calibrate shared effect loading across multiple levels
#'
#' @param state
#' @param EPS
#'
shared_effect_loading <- function(state, EPS = 1e-8){
    loading <- matrix(0, state$L, state$m)
    for(l in 1:state$L){
        alpha.l <- state$alpha[l, ]
        Z <- state$lbf_list[[l]]
        w.l <- apply(Z, 2, function(z) exp(z - max(z)))
        w.l <- apply(w.l, 2, function(z) z/pmax(sum(z), EPS))
        lambda <- apply(sweep(w.l, 1, alpha.l, `*`), 2, sum, na.rm=TRUE)
        loading[l, ] <- lambda
    }
    return(loading)
}

#' @title Compute posterior statistics
#' @param state
#'
#' @return A list with the following elements:
#'
#' \item{mean}{Posterior Mean combined across all the levels}
#' \item{sd}{Posterior Standard Deviation combined across all the levels}
#' \item{lfsr}{local False Sign Rate}
#' 
#' @importFrom states pnorm
#'
posterior_statistics <- function(state) {

    m <- state$m
    L <-  state$L
    lfsr <- matrix(NA, L, m)
    mean.list <- list()
    sd.list <- list()
    for(l in 1:L){
        .mu <- state$mu[[l]]
        .mu2 <- state$mu2[[l]]
        .alpha <- state$alpha[l, ]
        ## calibrate local false sign rates
        pos <- pnorm(0, .mu, sqrt(.mu2 - .mu^2))
        neg <- 1 - pos
        lfsr[l, ] <- 1 - apply(sweep(pmax(pos, neg), 1, .alpha, `*`), 2, sum)

        ## calibrate the mean and std dev
        .mean <- sweep(.mu, 1, .alpha, `*`)
        .var <- sweep(.mu2, 1, .alpha, `*`) - .mean^2

        mean.list[[l]] <- .mean
        sd.list[[l]] <- sqrt(.var)
    }

    list(mean = mean.list,
         sd = sd.list,
         lfsr = lfsr)
}

#' @param state
#' @param coverage
take_credible_sets <- function(state, coverage){
    p <- state$p
    .fun <- function(l){
        oo <- order(state$alpha[l, ], decreasing=TRUE)
        nn <- sum(cumsum(state$alpha[l, oo]) <= coverage) + 1
        k <- oo[1:min(nn, p)]
        ret <- data.frame()
        .alpha <- state$alpha[l, k]
        .lbf <- state$lbf_list[[l]]
        for(tt in 1:state$m){
            .beta <- state$mu[[l]][k, tt]
            .sigma <- sqrt(state$mu2[[l]][k, tt] - .beta^2)
            ret <- rbind(ret, data.frame(l, x.col = k,
                                         alpha = .alpha,
                                         lbf = .lbf[k, tt],
                                         beta = .beta,
                                         se = .sigma,
                                         y.col = tt))
        }
        rownames(ret) <- NULL
        return(ret)
    }
    ret <- do.call(rbind, lapply(1:state$L, .fun))
}

#' @rdname fit_mt_susie
#'
#' @title Fit multi-trait SuSiE across several regression models
#'
#' @param X n x p design matrix
#' @param Y n x m output matrix
#' @param L expected # of independent effects
#' @param max.iter maximum iterations
#' @param tol tolerance level check convergence
#'
#' @param clamp handle outliers by winsorization
#'
fit_mt_susie <- function(X, Y, L = 15, max.iter=100, tol = 1e-8, clamp = 6, coverage = .9) {

    xx <- apply(X, 2, scale)
    yy <- apply(Y, 2, scale)

    if(!is.null(clamp)){
        xx[xx > clamp] <- clamp
        xx[xx < -clamp] <- -clamp
        yy[yy > clamp] <- clamp
        yy[yy < -clamp] <- -clamp

        xx <- apply(xx, 2, scale)
        yy <- apply(yy, 2, scale)

        message("winsorization with the clamping value =", 6)
    }

    ret <- susie_initialize(n = nrow(yy), m = ncol(yy), p = ncol(xx), L)

    message("Initialized the Mt SuSiE model")

    loglik <- c()

    for(iter in 1:max.iter){

        ret <- update_shared_effect(xx, yy, ret)

        if(iter > 5){
            .prev <- tail(loglik, 1)
            .curr <- ret$llik
            diff <- abs(.prev - .curr)/abs(.curr)

            if(diff < tol){
                loglik <- c(loglik, .curr)
                message("converged at ", iter, ", ", .curr)
                break
            }
        }
        message("Iteration [", iter, "] ", ret$llik)
        loglik <- c(loglik, ret$llik)
    }

    ret$loglik_trace <- loglik

    ret$stats <- posterior_statistics(ret)
    ret$cs <- take_credible_sets(ret, coverage)

    return(ret)
}
