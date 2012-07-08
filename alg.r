# svd algorithm with constant k
svd.alg <- function(df, init=1, k=1, eps=0.05)
{
	mat <- matrix(init, max(df[[1]]), max(df[[2]]));
	err <- max(df$stars) - min(df$stars); olderr <- err/(1-eps)^2;
	while (err/olderr < 1-eps) {
		olderr <- err; err <- 0;
		apply(df, 1, function(x)
			{ err <<- err + (mat[x[[1]],x[[2]]] - x[[3]])^2;
			mat[x[[1]],x[[2]]] <<- x[[3]];}
		);
		err <- sqrt(err / length(df[[1]]));
		print(err);
		sing <- La.svd(mat, k, k);
		mat <- sing$u[,1:k] %*% diag(sing$d[1:k],k) %*% sing$v[1:k,]
	};
	return(mat)
}