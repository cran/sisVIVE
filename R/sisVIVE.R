sisVIVE <- function(Y,D,Z,intercept=TRUE,normalize=TRUE) {
  # Error checking
  if( (!is.vector(Y) && !is.matrix(Y)) | !is.numeric(Y)) stop("Y is not a numeric vector.")
  if( (!is.vector(D) && !is.matrix(D)) | !is.numeric(D)) stop("D is not a numeric vector.")
  if( !is.matrix(Z) | !is.numeric(Z)) stop("Z is not a numerical matrix.")
  if(nrow(Z) != length(Y)) stop("The dimension of Y differs from the row dimension of Z")
  if(nrow(Z) != length(D)) stop("The dimension of D differs from the row dimension of Z")
  if(length(D) != length(Y)) stop("The dimension of Y differs from the dimension of D")
  if(intercept) {
    if( (ncol(Z) + 1) > nrow(Z)) stop("The number of instruments is greater than the sample size.")
  } else {
    if( ncol(Z) > nrow(Z)) stop("The number of instruments is greater than the sample size.")
  }
  
  # Define constants 
  n = nrow(Z); L = ncol(Z);
  
  if(intercept) {
    meanY = mean(Y); meanD = mean(D); meanZ = colMeans(Z);
	Y = Y - mean(Y); D = D - meanD; Z = scale(Z,center=TRUE,scale=FALSE);
  } 
  if(normalize) {
  	normZ = sqrt(colSums(Z^2))
    Z = scale(Z,center=FALSE,scale=TRUE) / sqrt(n-1) 
  } else {
  	normZ = rep(1,L)
  }
  
  QR = qr(Z); Yhat = qr.fitted(QR,Y); Dhat = qr.fitted(QR,D);
  DhatZ = t(Dhat) %*% Z; Z.transformed = Z - Dhat %*% DhatZ / sum( Dhat^2)
  Y.transformed = Yhat - Dhat * sum(Dhat * Yhat) / sum( Dhat^2)
  fit = lars(Z.transformed,Y.transformed,intercept = FALSE,normalize=TRUE)
  estEffect = drop(t(Dhat) %*% (as.numeric(Y) - Z %*% t(coef(fit)) )) / sum(Dhat^2)
  
  # Packaging Object
  alpha = coef(fit)[-nrow(coef(fit)),,drop=FALSE]; alpha = scale(alpha,FALSE,normZ); 
  beta = estEffect[-length(estEffect)];
  alphaSuppSize = apply(alpha,1,function(x){sum(x != 0)}); indexProperSupp = which(alphaSuppSize < L);
  alpha = alpha[indexProperSupp,,drop=FALSE]; beta = estEffect[indexProperSupp]; lambda = fit$lambda[indexProperSupp];
  
  whichInvalid = apply(alpha,1,function(x){paste(which(abs(x) > 0),sep="",collapse=",")})
  dimnames(alpha) = list(paste(0:(nrow(alpha)-1)),dimnames(Z)[[2]]);
  names(beta) = paste(0:(nrow(alpha)-1));
  names(whichInvalid) = paste(0:(nrow(alpha)-1))
  
  object <- list(call = match.call(), alpha = alpha, beta = beta, whichInvalid = whichInvalid,lambda = lambda,larsObject = fit,Y=Y,D=D,Z=Z,Dhat=Dhat, normZ = normZ)
  class(object) <- "sisVIVE"
  return(object)
}

print.sisVIVE <- function(x,...) {
  object = x
  cat("\nCall:\n")
  dput(object$call)
  printOut = data.frame(object$beta,object$whichInvalid,object$lambda)
  dimnames(printOut) = list(paste(0:(nrow(object$alpha)-1)),c("Estimates of Beta","Invalid Instruments","Lambda Value"))
  print(printOut)
  invisible(object)
}

plot.sisVIVE <- function(x,...) {
  object = x
  alpha = object$alpha; beta = object$beta; lambda = object$lambda
  xplot = lambda; yplot = beta
  
  plot(xplot,yplot,ylab="Estimate of Beta",xlab="Lambda",type="b",lty=1,main="Estimate of beta for different lambda",...)
  invisible()
}

predict.sisVIVE <- function(object,lambda,type=c("coefficients","instruments"),...) {
  type = match.arg(type)
  if(missing(lambda)) lambda = object$lambda
  alpha = predict(object$larsObject,s = lambda,type="coefficients",mode="lambda",...)$coefficient
  if(is.vector(alpha)) alpha = matrix(alpha,1,length(alpha))
  beta = drop(t(object$Dhat) %*% (as.numeric(object$Y) - object$Z %*% t(alpha)))/sum(object$Dhat^2) 
  
  alpha = scale(alpha,FALSE,object$normZ); attr(alpha,"scaled:scale") = NULL #rescaling back to original
  alphaSuppSize = apply(alpha,1,function(x){sum(x != 0)}); indexImproperSupp = which(alphaSuppSize == ncol(alpha));
  alpha[indexImproperSupp,] = NA; beta[indexImproperSupp] = NA;
  return(switch(type, coefficients = list(lambda = lambda, alpha = alpha, beta = beta),
					  instruments = list(lambda = lambda, instruments = apply(alpha,1,function(x){paste(which(abs(x) > 0),sep="",collapse=",")}))))			
}

summary.sisVIVE <- function(object,...) {
  print.sisVIVE(object,...)
}

cv.sisVIVE <- function(Y,D,Z,lambdaSeq,K = 10,intercept=TRUE,normalize=TRUE) {
  # Error checking
  if( (!is.vector(Y) && !is.matrix(Y)) | !is.numeric(Y)) stop("Y is not a numeric vector.")
  if( (!is.vector(D) && !is.matrix(D)) | !is.numeric(D)) stop("D is not a numeric vector.")
  if( !is.matrix(Z) | !is.numeric(Z)) stop("Z is not a numerical matrix.")
  if(nrow(Z) != length(Y)) stop("The dimension of Y differs from the row dimension of Z")
  if(nrow(Z) != length(D)) stop("The dimension of D differs from the row dimension of Z")
  if(length(D) != length(Y)) stop("The dimension of Y differs from the dimension of D")
  if(intercept) {
    if( (ncol(Z) + 1) > nrow(Z)) stop("The number of instruments is greater than the sample size.")
  } else {
    if( ncol(Z) > nrow(Z)) stop("The number of instruments is greater than the sample size.")
  }
  if(!is.numeric(K) | K > length(Y)) stop("K is not a proper numeric number")
  fitall = sisVIVE(Y,D,Z,intercept=intercept,normalize=normalize)
  if(missing(lambdaSeq) || all(is.na(lambdaSeq))) {
    #warning("Lambda sequence not provided; defaulting to using lambdas provided by sisVIVE")
    lambdaSeq = c(fitall$lambda,seq(from=min(fitall$lambda,na.rm=TRUE),to=2*max(fitall$lambda,na.rm=TRUE),length.out = 100))
    lambdaSeq = sort(unique(lambdaSeq))
  }
  if(length(lambdaSeq) < 2) stop("Only one lambda provided. Please provide multiple lambdas")
  if(any(is.na(lambdaSeq))) {
    warning("Some lambda values are missing. Ignoring these lambda values for cross-validation")
	lambdaSeq = lambdaSeq[!is.na(lambdaSeq)]
  }
  lambdaSeq = sort(lambdaSeq)
  
  # Define constants
  n = nrow(Z); L = ncol(Z);
 
  # Cross validation
  Kfolds = split(sample(1:n),rep(1:K, length = n))
  errormat = matrix(0,K,length(lambdaSeq))
  for(i in seq(K)) {
    testSet = Kfolds[[i]]
	fit = sisVIVE(Y[-testSet], D[-testSet], Z[-testSet, , drop=FALSE], intercept = intercept,normalize = normalize)
	coefs = predict(fit, lambda = lambdaSeq, type = "coefficients")
	Y.test = Y[testSet]; D.test = D[testSet]; Z.test = Z[testSet,]
	if(intercept) {
      meanY.test = mean(Y.test); meanD.test = mean(D.test); meanZ.test = colMeans(Z.test)
	  Y.test = Y.test - meanY.test; D.test = D.test - meanD.test; Z.test = scale(Z.test,center=TRUE,scale=FALSE)
    } 
	
	QR = qr(Z.test)
	residTest = (as.numeric(Y.test) - Z.test %*% t(coefs$alpha) - D.test %*% t(coefs$beta))
	badLambda = which(is.na(residTest[1,]) == TRUE); goodLambda = which(is.na(residTest[1,]) == FALSE)
	errormat[i,badLambda] = NA
	if(length(goodLambda) == 1) {
	  errormat[i,goodLambda] = sum(qr.fitted(QR,residTest[,goodLambda])^2)
	} else {
	  errormat[i,goodLambda] = colSums(qr.fitted(QR, residTest[,goodLambda])^2) #does (PZ %*% residual)^2 summed across observations
	}
  }
  cv = colMeans(errormat)
   
  if(all(is.nan(cv))) {
    warning("All lambdas were invalid. Please try different values of lambda for cross-validation to work")
	return(list(lambda = rep(NA,length(lambdaSeq)),estCVError = NA, alpha = rep(NA,ncol(Z)),beta = NA, whichInvalid = NA))
  } else {
    stderror = apply(errormat,2,function(x){sd(x)/sqrt(K)})
    mincv.index = which.min(cv)
    onestderr.rule.index = max(which(cv <= (cv[mincv.index] + stderror[mincv.index]) & cv[mincv.index] <= cv)) #this is never empty vector 
    returnOut = predict(fitall,lambdaSeq[onestderr.rule.index],type="coefficients")
    return(list(lambda = returnOut$lambda, estCVError = cv[onestderr.rule.index], alpha = drop(returnOut$alpha), beta = drop(returnOut$beta),whichInvalid = paste(which(abs(returnOut$alpha) > 0),sep="",collapse=",")))
  }
}
