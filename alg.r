# SVD algorithm with rank k approximation, stops when error decrease is below eps.
alg.svd <- function(df, pl=alg.svd.pl(df))
{
	mat <- matrix(pl$init, max(df$user), max(df$movie))
	given <- mask <- matrix(0, max(df$user), max(df$movie))
	given[(df$movie-1)*nrow(mat)+df$user] <- df$stars		# matrix of given values
	mask[(df$movie-1)*nrow(mat)+df$user] <- 1				# contains a '1' for each given value

	# error history, initialize with dummy values making no trouble.
	errvec <- (max(df$stars) - min(df$stars)) * c(1, 1/(1-pl$eps)^2)
	while (errvec[1]/errvec[2] < 1-pl$eps) {
		errvec <- c(sqrt(sum((mat*mask-given)^2) / nrow(df)), errvec) # compute error
		mat <- mat*(1-mask)+given										# then overwrite with given values 
		sing <- La.svd(mat, pl$k, pl$k)
		mat <- sing$u %*% diag(sing$d[1:pl$k],pl$k) %*% sing$v			# and compute rank k approximation
	}

	print(errvec)
	return(mat)
}

alg.svd.pl <- function(df, init=mean(df$stars), k=1, digits=2) list(init=init, k=k, eps=10^-digits)

# Hazan's algorithm with target trace tr, curvature constant Cf and eps as above,
#	averaged on error history of maximum length maxhist.
alg.hazan <- function(df, pl=alg.hazan.pl(df, mean(df$stars)*(max(df$user) + max(df$movie))))
{
	n <- max(df$user); m <- max(df$movie); len <- nrow(df)
	Y <- matrix(0, n+m, n+m)
	Y[(df$movie-1)*(n+m)+df$user+m] <- df$stars	# "given" matrix Y we want to approximate 
	Y[(df$user-1+m)*(n+m)+df$movie] <- df$stars

	# initialize X = v*v^T with an eigenvector belonging to the greatest eigenvalue 
	i <- 0; v <- power.method(Y, pl$Cf/pl$tr)
	X <- pl$tr * (v %*% t(v)) / sum(v*v)
	errvec <- (sum(ifelse(Y!=0, (X-Y)^2, 0))/len) * c(1/(1-pl$eps)^2, 1/(1-pl$eps)^4)

	decr <- 1;	# average error decrease
	while (decr > pl$eps | i < 10) {
		# compute error and average error decrease
		errvec <- c(sum(ifelse(Y!=0, (X-Y)^2, 0))/len, errvec)
		hist <- min(pl$maxhist, length(errvec))
		decr <- lm(errvec[1:hist] ~ as.numeric(1:hist))$coefficients[[2]]

		alpha <- 2/(i+2);	i <- i+1
		# compute "symmetricized" gradient matrix of f(X) = \sum_(\Omega+(0,m)) (X_ij-Y_ij)^2
		Nabla <- ifelse(Y!=0, 2*(X-Y), 0)

		# compute an eigenvector corresponding to the greatest eigenvalue
		v <- power.method(Nabla, alpha*pl$Cf/pl$tr)

		# blend old X with tr*v*v^T
		X <- (1-alpha)*X + alpha*pl$tr * v %*% t(v)
	}

	print(errvec)
	return(X[(m+1):(n+m),1:m])
}

alg.hazan.pl <- function(df, tr, digits=2, Cf=curvature(df[c(1,2)]), maxhist=5) list(tr=tr, eps=10^-digits, Cf=tr^2 * Cf, maxhist=maxhist);

# "power method" to compute an eigenvector corresponding to the greatest eigenvalue
power.method <- function(A, eps)
{
	v <- runif(dim(A)[1]); v <- v/sqrt(sum(v*v))
	l<-2; oldl<-1;	# not "correct", but works in high dimensions
	while (l/oldl > 1 + eps) {
		v <- A %*% v
		oldl <- l; l <- sqrt(sum(v*v))
		v <- v/l
	}
	return(v)
}

# compute lower bound for the "curvature constant" C_f
curvature <- function(df, samples=10)
{
	n <- max(df$user); m <- max(df$movie); len <- nrow(df)
	Y <- matrix(0, n+m, n+m)
	Y[(df$movie-1)*nrow(Y)+df$user+m] <- 1

	samp <- vector(mode="list", length=samples)
	for (i in 1:samples) {
		A <- matrix(rnorm((n+m)^2), nrow=n+m)		# generate n independent vectors
		A <- apply(A, 2,
			function(x) {x/sqrt(sum(x*x))})		# normalize them
		w <- rexp(n+m, n+m); w <- w / sum(w)		# generate weights
		samp[[i]] <- A %*% (t(A)*w)					# generate matrix
	};

	maximum <- 0;
	for (i in 1:samples) for (j in 1:samples) { diff<-sum((samp[[i]]-samp[[j]])^2 * Y)/len; maximum <- max(maximum, diff) }
	return(maximum)
}

# some post processing
post.clamp <- function(data, rg) { ifelse(data<rg[1], rg[1], ifelse(data>rg[2], rg[2], data)) }
post.linear <- function(data, rg) { old<-range(data); (data-old[1])/(old[2]-old[1])*(rg[2]-rg[1])+rg[1] }
post.logistic <- function(data, rg)
{
	data <- (2*data-(rg[1]+rg[2]))/(rg[2]-rg[1])	# scale to [-1,1]
	data <- 1/(1+exp(-2*data))						# apply logistic function with derivative 1/2 in 0.
	data*(rg[2]-rg[1]) + rg[1]
}
