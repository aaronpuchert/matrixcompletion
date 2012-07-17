# svd algorithm with constant k
alg.svd <- function(df, init=mean(df$stars), k=1, eps=0.01)
{
	mat <- matrix(init, max(df$user), max(df$movie));
	given <- mask <- matrix(0, max(df$user), max(df$movie));
	given[(df$movie-1)*nrow(mat)+df$user] <- df$stars;
	mask[(df$movie-1)*nrow(mat)+df$user] <- 1;

	errvec <- (max(df$stars) - min(df$stars)) * c(1, 1/(1-eps)^2);
	while (errvec[1]/errvec[2] < 1-eps) {
		errvec <- c(sqrt(sum((mat*mask-given)^2) / nrow(df)), errvec);
		mat <- mat*(1-mask)+given;
		sing <- La.svd(mat, k, k);
		mat <- sing$u[,1:k] %*% diag(sing$d[1:k],k) %*% sing$v[1:k,]
	};

	print(errvec);
	return(mat)
}

# Hazan's algorithm with assumed contraction speed c of von Mises iteration
alg.hazan <- function(df, tr=1, alpha=0.05, eps=0.01, c=0.5)
{
	n <- max(df$user); m <- max(df$movie); len <- nrow(df);
	Y <- matrix(NA, n+m, n+m);
	Y[(df$movie-1)*nrow(mat)+df$user+m] <- df$stars;

	i <- 0; v <- runif(n+m);
	X <- tr * (v %*% t(v)) / (t(v) %*% v);
	errvec <- (sum(ifelse(Y!=NA, (X-Y)^2, 0))/len) * c(1/(1-eps)^2, 1/(1-eps)^4);

	w <- runif(n+m); Z <- tr * (w %*% t(w)) / (t(w) %*% w)
	Cf <- sum((X-Z)^2 * (Y!=NA))/len;	# this is just a lower bound (if necessary, enhance)
	while (errvec[1]/errvec[2] < 1-eps) {
		errvec <- c(sum(ifelse(Y!=NA, (X-Y)^2, 0))/len, errvec);
		alpha <- 2/(i+2);	i <- i+1;
		Nabla <- ifelse(Y!=NA, 2*(X-Y), 0);

		# von Mises iteration
		v <- runif(n+m); l<-2; oldl<-1;	# better initialization
		while (l/oldl > 1 + c*alpha*Cf/t) {
			v <- Nabla %*% v;
			l <- sqrt(sum(v*v));
			v <- v/l;
		}

		X <- (1-alpha)*X + alpha*tr * v %*% t(v);
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
