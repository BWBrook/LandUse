AICc.lmer <- function(...) {
 models <- list(...)
 num.mod <- length(models)
 AICcs <- numeric(num.mod)
 ns <- numeric(num.mod)
 ks <- numeric(num.mod)
 AICc.vec <- rep(0,num.mod)
 for (i in 1:num.mod) {
  n <- models[[i]]@XtX@x[1]
  k <- length(models[[i]]@fixef)+sum(models[[i]]@nc)+1
  AICcs[i] <- (-2*logLik(models[[i]])[1]) + ((2*k*n)/(n-k-1))
  ns[i] <- n
  ks[i] <- k
  AICc.vec[i] <- AICcs[i]}
 return(AICc.vec)}



BIC.lmer <- function(...) {
	models <- list(...)
	num.mod <- length(models)
	BICs <- numeric(num.mod)
	ns <- numeric(num.mod)
	ks <- numeric(num.mod)
	BIC.vec <- rep(0,num.mod)
	for (i in 1:num.mod) {
    n <- models[[i]]@XtX@x[1]
    k <- length(models[[i]]@fixef)+sum(models[[i]]@nc)+1
		BICs[i] <- (-2*logLik(models[[i]])) + k*log(n)
		ns[i] <- n
		ks[i] <- k
		BIC.vec[i] <- BICs[i]}
	return(BIC.vec)}




AICc.glm <- function(...) {
	models <- list(...)
	num.mod <- length(models)
	AICcs <- numeric(num.mod)
	ns <- numeric(num.mod)
	ks <- numeric(num.mod)
	AICc.vec <- rep(0,num.mod)
	for (i in 1:num.mod) {
		if (length(models[[i]]$df.residual) == 0) n <- models[[i]]$dims$N else n <- length(models[[i]]$residuals)
		if (length(models[[i]]$df.residual) == 0) k <- sum(models[[i]]$dims$ncol) else k <- (length(models[[i]]$coeff))+1
		AICcs[i] <- (-2*logLik(models[[i]])) + ((2*k*n)/(n-k-1))
		ns[i] <- n
		ks[i] <- k
		AICc.vec[i] <- AICcs[i]}
	return(AICc.vec)}



BIC.glm <- function(...) {
	models <- list(...)
	num.mod <- length(models)
	BICs <- numeric(num.mod)
	ns <- numeric(num.mod)
	ks <- numeric(num.mod)
	BIC.vec <- rep(0,num.mod)
	for (i in 1:num.mod) {
		if (length(models[[i]]$df.residual) == 0) n <- models[[i]]$dims$N else n <- length(models[[i]]$residuals)
		if (length(models[[i]]$df.residual) == 0) k <- sum(models[[i]]$dims$ncol) else k <- (length(models[[i]]$coeff))+1
		BICs[i] <- (-2*logLik(models[[i]])) + k*log(n)
		ns[i] <- n
		ks[i] <- k
		BIC.vec[i] <- BICs[i]}
	return(BIC.vec)}
                 


k.lmer <- function(x) if(length(x$df.residual) == 0) length(x@fixef)+sum(x@nc)+1


k.glm <- function(x) if(length(x$df.residual) == 0) k <- sum(x$dims$ncol) else k <- (length(x$coeff)+1)
 

delta.IC <- function(x) x - min(x) ## where x is a vector of an IC


weight.IC <- function(x) (exp(-0.5*x))/sum(exp(-0.5*x)) ## Where x is a vector of dIC


chdev.lmer <- function(x) ((( as.numeric(dev.null) - as.numeric(deviance(x)))/ as.numeric(dev.null)*100))


chdev.glm <- function(x) ((( as.numeric(x[12]) - as.numeric(x[10]))/ as.numeric(x[12]))*100) ## % change in deviance, where x is glm object


linreg.ER <- function(x,y) { # where x and y are vectors of the same length; calls AICc, delta.AIC, weight.AIC functions
  fit.full <- lm(y ~ x); fit.null <- lm(y ~ 1)
  AIC.vec <- c(AICc.glm(fit.full),AICc.glm(fit.null))
  dAIC.vec <- delta.IC(AIC.vec); wAIC.vec <- weight.IC(dAIC.vec)
  ER <- wAIC.vec[1]/wAIC.vec[2]
  r.sq.adj <- as.numeric(summary(fit.full)[9])
  return(ER,r.sq.adj)
}


tl.full.func <- function(x)
   {
   	r.fit <- x[1]*(1-(N.fix/x[2])^x[3])
	  r.dev <- (r.fix-r.fit)^2
	  r.sum <- sum(r.dev)
    return(r.sum)
   }


tl.drop.func <- function(x)
   {
   	r.fit <- x[1]*(1-(N.fix[-i]/x[2])^x[3])
	  r.dev <- (r.fix[-i]-r.fit)^2
	  r.sum <- sum(r.dev)
    return(r.sum)
   }

 aicW <- function(model.list, finite = TRUE, null.model = NULL, N, order = TRUE){
	if(hasArg(N)) { N <- N }
        else N <- NULL
	if(is.null(null.model)) {
		null.model <- model.list[[length(model.list)]]
	}
	LLlist <- vector()
	for(index in 1:length(model.list)){
		if(class(model.list[[index]])[1] == "mer") {
			LLlist <- c(LLlist, as.numeric(logLik(model.list[[index]])))
		}
		else {LLlist <- c(LLlist, as.numeric(logLik(model.list[[index]])))}
	}
	aiclist <- vector()
	for(index in 1:length(model.list)){
		aicname <- "AIC"
		if(finite==T) {
			aiclist <- c(aiclist, AIC.finite.lm(model.list[[index]], N))
			aicname <- "AICc"}
		else {aiclist <- c(aiclist, AIC.lm(model.list[[index]]))}
	}
	biclist <- vector()
	for(index in 1:length(model.list)){
		if(class(model.list[[index]])[1] == "mer") {
			biclist <- c(biclist, BIC.lm(model.list[[index]], N))
		}
		else {biclist <- c(biclist, BIC.lm(model.list[[index]], N))}
	}
	dflist <- vector()
	for(index in 1:length(model.list)){
		if(class(model.list[[index]])[1] == "mer") {
			dflist <- c(dflist, as.numeric(attributes(logLik(model.list[[index]]))$df))
		}
		else {dflist <- c(dflist, length(coef(model.list[[index]])))}
	}
	delist <- vector()
	for(index in 1:length(model.list)){
		if(class(model.list[[index]])[1] == "mer") {
			delist <- c(delist, lmer.de(model.list[[index]], null.model, N))
		}
		else {delist <- c(delist, lmer.de(model.list[[index]], null.model, N))}
	}
	modform <- vector()	
	if(is.null(names(model.list))){
		for(index in 1:length(model.list)){
			modform <- c(modform, formula(model.list[[index]]))
		}
	}	
	if(!is.null(names(model.list))){
		models <- names(model.list)
		delta.i <- aiclist - min(aiclist)
		aicw <- exp(delta.i * -0.5)/sum(exp(delta.i * -0.5))
		delta.BIC <- biclist - min(biclist)
    bicw <- exp(delta.BIC * -0.5)/sum(exp(delta.BIC * -0.5))
		return.matrix <- matrix(c(round(LLlist, 3), dflist, round(aiclist, 3), round(delta.i, 3), round(aicw, 4), round(biclist,3), round(delta.BIC, 3), round(bicw,4), round(delist, 2)), ncol = 9)
		colnames(return.matrix) <- c("logLik", "df", aicname, paste("d", aicname, sep = ""), "wAICc", "BIC", "dBIC", "wBIC", "%DE")
		# rownames(return.matrix) <- sapply(strsplit(models, NULL), function(x) paste(x, collapse = " + "))
		rownames(return.matrix) <- models
	}
	else { 
		delta.i <- aiclist - min(aiclist)
		aicw <- exp(delta.i * -0.5)/sum(exp(delta.i * -0.5))
		delta.BIC <- biclist - min(biclist)
    bicw <- exp(delta.BIC * -0.5)/sum(exp(delta.BIC * -0.5))
		return.matrix <- matrix(c(round(LLlist, 3), dflist, round(aiclist, 3), round(delta.i, 3), round(aicw, 4), round(biclist,3), round(delta.BIC, 3), round(bicw,4), round(delist, 2)), ncol = 9)
#		return.matrix <- matrix(c(round(LLlist, 3), dflist, round(aiclist, 3), round(delta.i, 3), round(aicw, 4), round(delta.BIC, 3), round(delist, 2)), ncol = 7)
		colnames(return.matrix) <- c("logLik", "df", aicname, paste("d", aicname, sep = ""), "wAICc", "BIC", "dBIC", "wBIC", "%DE")
#		colnames(return.matrix) <- c("logLik", "df", aicname, paste("d", aicname, sep = ""), "weight", "dBIC", "%DE")
		rownames(return.matrix) <- modform
	}
	if(order) {
	return.matrix <- return.matrix[order(return.matrix[, 4]), ]
	}
	return(return.matrix)
}

AIC.finite.lm <- function(model, N){
	LL = as.numeric(logLik(model))
	p = as.numeric(attributes(logLik(model))$df)
	if(class(model)[1] == "mer"){
	    if(is.null(N)){
	        N <- model@A@Dim[2]
	    }
		else N <- N
	}
	else {
	    if(is.null(N)){
	        N <- length(model$data[, 1])
	    }
		else N <- N
	}
	AICc <- AIC(logLik(model)) + (2 * p * (p + 1))/(N - p - 1)
	return(AICc)
}

AIC.lm <- function(object) {
	ret <- AIC(logLik(object))
	return(ret)
}

BIC.lm <- function(model, N) {
	LL = as.numeric(logLik(model))
	p = as.numeric(attributes(logLik(model))$df)
	if(class(model)[1] == "mer"){
	    if(is.null(N)){
	        N <- model@A@Dim[2]
	    }
		else N <- N
	}
	else {
	    if(is.null(N)){
	        N <- length(model$data[, 1])
	    }
		else N <- N
	}
	BIC <- (-2 * LL) + (p * log(N))
	return(BIC)
}

lmer.de <- function(model, null.model, N) {
	if(class(model)[1] == "mer"){
	    if(is.null(N)){
	        N <- model@A@Dim[2]
	    }
		else N <- N
	}
	else {
	    if(is.null(N)){
	        N <- length(model$data[, 1])
	    }
		else N <- N
	}
    if(class(null.model)[1] == "mer"){
        null.dev <- as.vector(null.model@deviance["ML"])
	    model.lr <- as.vector(null.model@deviance["ML"]) - as.vector(model@deviance["ML"])
    }
    else {
		if(class(model)[1] == "mer"){
			null.dev <- deviance(null.model)
		    model.lr <- (deviance(null.model)) - as.vector(model@deviance["ML"])
		}
		else {
			null.dev <- deviance(null.model)
		    model.lr <- deviance(null.model) - deviance(model)
		}
	}
	de <- 100 * (model.lr/null.dev)
    return(de)
}
