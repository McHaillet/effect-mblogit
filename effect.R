# Effects package implementation adapted for mblogit models from mclogit package
# by Marten Chaillet (@McHaillet)

eff.mul <- function(x0, B, se, m, p, r, V){
    mu <- exp(x0 %*% B)
    mu <- mu/(1 + sum(mu))
    mu[m] <- 1 - sum(mu)
    logits <- log(mu/(1 - mu))
    # if (!se) return(list(p=mu, logits=logits))
    d <- array(0, c(m, m - 1, p))
    exp.x0.B <- as.vector(exp(x0 %*% B))
    sum.exp.x0.B <- sum(exp.x0.B)
    for (j in 1:(m-1)){
        d[m, j,] <- - exp.x0.B[j]*x0
        for (jj in 1:(m-1)){
            d[j, jj,] <- if (jj != j) {
                - exp(as.vector(x0 %*% (B[,jj] + B[,j])))*x0
            } else {
                exp.x0.B[j]*(1 + sum.exp.x0.B - exp.x0.B[j])*x0
            }
        }
    }
    d <- d/(1 + sum.exp.x0.B)^2
    V.mu <- rep(0, m)
    for (j in 1:m){
        dd <- as.vector(t(d[j,,]))
        for (s in 1:r){
            for (t in 1:r){
                V.mu[j] <- V.mu[j] + V[s,t]*dd[s]*dd[t]
            }
        }
    }
    V.logits <- V.mu/(mu^2 * (1 - mu)^2)
    list(p=mu, std.err.p=sqrt(V.mu), logits=logits,
       std.error.logits=sqrt(V.logits))
}
                                         
logit2p <- function(logit) 1/(1 + exp(-logit))

p2logit <- function(p) log(p/(1 - p))

Effect.mblogit <- function(mod, confidence.level=0.95){
    
    # get predictors/response names from model object and put the baseline last
    resp.names <- rownames(mod$D)
    resp.names <- c(resp.names[-1], resp.names[1])
    pred.names <- mod$xlevels$state
    
    # get size of predictor and response
    n_pred <- length(pred.names)
    n_resp <- length(resp.names) - 1
    
    # prepare predictor grid
    predictors <- expand.grid(state = pred.names)
    X0 <- model.matrix(~ state, predictors)
    
    B <- t(matrix((coef(mod)), n_resp, n_pred))
    # test
    V <- vcov(mod)
    id <- rep(-1, length = n_pred * n_resp)
    for (i in 1:n_resp) {
        for (j in 1:n_pred) {
            id[(i - 1) * n_pred + j] <- i + (j - 1) * n_resp
        }
    }
    V <- V[id, id]
    
    m <- ncol(B) + 1
    p <- nrow(B)
    r <- p*(m - 1)
    n <- nrow(X0)
    P <- Logit <- matrix(0, n, m)
    colnames(P) <-  paste("prob.", resp.names, sep="")
    colnames(Logit) <-  paste("logit.", resp.names, sep="")


    z <- qnorm(1 - (1 - confidence.level)/2)
    Lower.P <- Upper.P <- Lower.logit <- Upper.logit <- SE.P <- SE.logit <- matrix(0, n, m)
    colnames(Lower.logit) <-  paste("L.logit.", resp.names, sep="")
    colnames(Upper.logit) <-  paste("U.logit.", resp.names, sep="")
    colnames(Lower.P) <-  paste("L.prob.", resp.names, sep="")
    colnames(Upper.P) <-  paste("U.prob.", resp.names, sep="")
    colnames(SE.P) <-  paste("se.prob.", resp.names, sep="")
    colnames(SE.logit) <-  paste("se.logit.", resp.names, sep="")

    for (i in 1:n){
        res <- eff.mul(X0[i,], B, se, m, p, r, V) # compute effects
        #        P[i,] <- prob <- res$p # fitted probabilities
        P[i,] <- res$p # fitted probabilities
        Logit[i,] <- logit <- res$logits # fitted logits

        #            SE.P[i,] <- se.p <- res$std.err.p # std. errors of fitted probs
        SE.P[i,] <- res$std.err.p # std. errors of fitted probs
        SE.logit[i,] <- se.logit <- res$std.error.logits # std. errors of logits
        Lower.P[i,] <- logit2p(logit - z*se.logit)
        Upper.P[i,] <- logit2p(logit + z*se.logit)
        Lower.logit[i,] <- logit - z*se.logit
        Upper.logit[i,] <- logit + z*se.logit
    }
    resp.levs <- c(m, 1:(m-1)) # restore the order of the levels
    P <- P[, resp.levs]
    Logit <- Logit[, resp.levs]

    Lower.P <- Lower.P[, resp.levs]
    Upper.P <- Upper.P[, resp.levs]
    Lower.logit <- Lower.logit[, resp.levs]
    Upper.logit <- Upper.logit[, resp.levs]
    SE.P <- SE.P[, resp.levs]
    SE.logit <- SE.logit[, resp.levs]

    result <- list(formula=formula(mod),
                 model.matrix=X0, model="mblogit",
                 prob=P, logit=Logit)
    result <- c(result, list(se.prob=SE.P, se.logit=SE.logit,
                                   lower.logit=Lower.logit, upper.logit=Upper.logit,
                                   lower.prob=Lower.P, upper.prob=Upper.P,
                                   confidence.level=confidence.level))    
    result
}

