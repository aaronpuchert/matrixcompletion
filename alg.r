# svd algorithm with constant k
svd.alg <- function(df, init=1, k=1, eps=0.05)
{
	mat <- matrix(init, max(df$user), max(df$movie));
	given <- mask <- matrix(0, max(df$user), max(df$movie));
	apply(df, 1, function(x) {given[x[[1]],x[[2]]] <<- x[[3]]; mask[x[[1]],x[[2]]] <<- 1} );

	errvec <- max(df$stars) - min(df$stars);
	errvec <- c(errvec, errvec/(1-eps)^2);
	while (errvec[1]/errvec[2] < 1-eps) {
		errvec <- c(sqrt(sum((mat*mask-given)^2) / sum(mask)), errvec);
		mat <- mat*(1-mask)+given;
		sing <- La.svd(mat, k, k);
		mat <- sing$u[,1:k] %*% diag(sing$d[1:k],k) %*% sing$v[1:k,]
	};

	print(errvec);
	return(mat)
}