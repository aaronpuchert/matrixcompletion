### ALGORITHMS
# SVD algorithm with rank k approximation, stops when error decrease is below eps.
alg.svd <- function(df, pl=alg.svd.pl(df), debug=FALSE)
{
	mat <- matrix(pl$init, max(df$user), max(df$movie))

	# Matrix of given values
	given <- mask <- matrix(0, max(df$user), max(df$movie))
	given[(df$movie-1)*nrow(mat)+df$user] <- df$stars

	# Mask containing a '1' for each given value
	mask[(df$movie-1)*nrow(mat)+df$user] <- 1

	# Error history of length 2, initialize with dummy values making no trouble.
	errvec <- (max(df$stars) - min(df$stars)) * c(1, 1/(1-pl$eps)^2)
	while (errvec[1]/errvec[2] < 1-pl$eps) {
		# Compute error
		err <- sqrt(sum((mat*mask-given)^2) / nrow(df))
		if (debug) print(err)
		errvec <- c(err, errvec[1])
		# Then overwrite with given values ...
		mat <- mat*(1-mask)+given
		sing <- La.svd(mat, pl$k, pl$k)
		# ... and compute rank k approximation
		mat <- sing$u %*% diag(sing$d[1:pl$k],pl$k) %*% sing$v
	}

	return(mat)
}

alg.svd.pl <- function(df, init=mean(df$stars), k=1, digits=2)
	list(init=init, k=k, eps=10^-digits)

# Hazan's algorithm with target trace tr, curvature constant Cf and eps as above,
#	averaged on error history of maximum length maxhist.
alg.hazan <- function(df, pl=alg.hazan.pl(df, mean(df$stars)*(max(df$user) + max(df$movie))), debug=FALSE)
{
	n <- max(df$user); m <- max(df$movie); len <- nrow(df)

	# "Given" matrix Y we want to approximate
	Y <- matrix(0, n+m, n+m)
	Y[(df$movie-1)*(n+m)+df$user+m] <- df$stars
	Y[(df$user-1+m)*(n+m)+df$movie] <- df$stars

	# Initialize X = v*v^T with an eigenvector belonging to the greatest eigenvalue
	i <- 0; v <- power.method(Y, pl$Cf/pl$tr)
	X <- pl$tr * (v %*% t(v)) / sum(v*v)
	errvec <- (sum(ifelse(Y!=0, (X-Y)^2, 0))/len) * c(1/(1-pl$eps)^2, 1/(1-pl$eps)^4)

	decr <- 1;	# Average error decrease
	while (decr > pl$eps | i < 100) {
		# Compute error and average error decrease
		err <- sum((Y!=0) * (X-Y)^2)/len
		if (debug) print(err)
		errvec <- c(err, errvec)
		hist <- min(pl$maxhist, length(errvec))
		decr <- lm(errvec[1:hist] ~ as.numeric(1:hist))$coefficients[[2]]

		alpha <- 2/(i+2);	i <- i+1
		# Compute "symmetricized" gradient matrix of
		# f(X) = \sum_(\Omega+(0,m)) (X_ij-Y_ij)^2
		Nabla <- 2 * (Y != 0) * (X-Y)

		# Compute an eigenvector corresponding to the greatest eigenvalue
		v <- power.method(Nabla, alpha*pl$Cf/pl$tr)

		# Blend old X with tr*v*v^T
		X <- (1-alpha)*X + alpha*pl$tr * v %*% t(v)
	}

	return(X[(m+1):(n+m),1:m])
}

alg.hazan.pl <- function(df, tr, digits=2, Cf=curvature(df[c(1,2)]), maxhist=50)
	list(tr=tr, eps=10^-digits, Cf=tr^2 * Cf, maxhist=maxhist)

# "Power method" to compute an eigenvector corresponding to the greatest eigenvalue
power.method <- function(A, eps)
{
	# Start with random normalized vector
	v <- runif(dim(A)[1])
	v <- v / sqrt(sum(v*v))

	l<-2; oldl<-1;	# not "correct", but works in high dimensions
	while (l/oldl > 1 + eps) {
		v <- A %*% v
		oldl <- l; l <- sqrt(sum(v*v))
		v <- v/l
	}
	return(v)
}

# Compute lower bound for the "curvature constant" C_f
curvature <- function(df, samples=10)
{
	n <- max(df$user); m <- max(df$movie); len <- nrow(df)
	Y <- matrix(0, n+m, n+m)
	Y[(df$movie-1)*nrow(Y)+df$user+m] <- 1

	samp <- vector(mode="list", length=samples)
	for (i in 1:samples) {
		# Generate n independent vectors
		A <- matrix(rnorm((n+m)^2), nrow=n+m)
		# Normalize them
		A <- apply(A, 2, function(x) {x/sqrt(sum(x*x))})
		# Generate weights
		w <- rexp(n+m, n+m); w <- w / sum(w)
		# Generate matrix
		samp[[i]] <- A %*% (t(A)*w)
	};

	maximum <- 0;
	for (i in 1:samples)
		for (j in 1:samples) {
			diff <- sum((samp[[i]]-samp[[j]])^2 * Y)/len
			maximum <- max(maximum, diff)
		}
	return(maximum)
}

### POST PROCESSING ###
post.clamp <- function(data, rg)
{
	ifelse(data<rg[1],
		rg[1],
		ifelse(data>rg[2],
			rg[2],
			data
		)
	)
}

post.linear <- function(data, rg)
{
	old<-range(data);
	(data-old[1])/(old[2]-old[1])*(rg[2]-rg[1])+rg[1]
}

post.logistic <- function(data, rg)
{
	# Scale to [-1,1]
	data <- (2*data-(rg[1]+rg[2]))/(rg[2]-rg[1])
	# Apply logistic function with derivative 1/2 in 0.
	data <- 1/(1+exp(-2*data))
	data*(rg[2]-rg[1]) + rg[1]
}
