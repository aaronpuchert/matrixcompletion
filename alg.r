# svd algorithm with constant k
alg.svd <- function(df, init=mean(df$stars), k=1, eps=0.01)
{
	mat <- matrix(init, max(df$user), max(df$movie));
	given <- mask <- matrix(0, max(df$user), max(df$movie));
	given[(df$movie-1)*nrow(mat)+df$user] <- df$stars;
	mask[(df$movie-1)*nrow(mat)+df$user] <- 1;

	errvec <- max(df$stars) - min(df$stars);
	errvec <- c(errvec, errvec/(1-eps)^2);
	while (errvec[1]/errvec[2] < 1-eps) {
		errvec <- c(sqrt(sum((mat*mask-given)^2) / nrow(df)), errvec);
		mat <- mat*(1-mask)+given;
		sing <- La.svd(mat, k, k);
		mat <- sing$u[,1:k] %*% diag(sing$d[1:k],k) %*% sing$v[1:k,]
	};

	print(errvec);
	return(mat)
}

# Hazan's algorithm
alg.hazan <- function(df, t=1, eps=0.01)
{
	n <- max(df$user); m <- max(df$movie); l <- nrow(df);
	Y <- matrix(NA, n+m, n+m);
	apply(df, 1, function(x) {Y[x[[1]]+m, x[[2]]] <<- x[[3]]} );

	errvec <- max(df$stars) - min(df$stars);	# how to initialize?
	errvec <- c(errvec, errvec/(1-eps)^2);
	i <- 0; v <- runif(n+m);
	X <- t * v %*% t(v);
	while (errvec[1]/errvec[2] < 1-eps) {
		errvec <- c(sum(ifelse(Y!=NA, (X-Y)^2, 0))/l, errvec);
		alpha <- 2/(i+2);	i <- i+1;
		Nabla <- ifelse(Y!=NA, 2*(X-Y), 0);
		eig <- eigen(Nabla, symmetric=TRUE);
		X <- (1-alpha)*X + alpha*t*eig$vectors[1,] %*% t(eig$vectors[1,]);
	};

	print(errvec);
	return(X[(m+1):(n+m),1:m])
}

# some post processing
post.clamp <- function(data, rg) { ifelse(data<rg[1], rg[1], ifelse(data>rg[2], rg[2], data)) }
post.linear <- function(data, rg) { old<-range(data); (data-old[1])/(old[2]-old[1])*(rg[2]-rg[1])+rg[1] }
post.logistic <- function(data, rg)
{
	data <- (2*data-(rg[1]+rg[2]))/(rg[2]-rg[1]);	# scale to [-1,1]
	data <- 1/(1+exp(-2*data));						# apply logistic function with derivative 1/2 in 0.
	data*(rg[2]-rg[1]) + rg[1]
}
