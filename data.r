### INPUT PROCESSING ###
# read csv data
read.data <- function (filename, ...)
{
	tab <- read.table(filename, sep=",", header=FALSE, ...);
	if (dim(tab)[2] == 2) col.names(tab) = c("user", "movie") else col.names(tab) = c("user", "movie", "stars");
	class(tab) <- c("data", "data.frame");
	tab
}

# reduce user and movie domain
reduce.data <- function (df, feat=-length(df))
{
	ui <- sort(union(df$user, NULL));
	mi <- sort(union(df$movie, NULL));
	df$user <- sapply(df$user, function(x) {which(ui==x)});
	df$movie <- sapply(df$movie, function(x) {which(mi==x)});
	list(data=df,ui=ui,mi=mi)
}

# plot in rainbow colors
plot.data <- function(df, sr=5, ...)
{plot.default(df[-3], pch=".", col = rainbow(sr)[df[[3]]], ...)}

# create data matrix
matrix.data <- function (df, init=0)
{
	mat <- matrix(init, max(df[[1]]), max(df[[2]]));
	apply(df, 1, function(x) {mat[x[[1]],x[[2]]] <<- x[[3]]});
	mat
}

### EVALUATION & OUTPUT ###
# evaluate matrix on test data set
eval <- function(df.test, mat) {
	cbind(df.test, est=0);
	apply(df.test, 1, function(x) {x$est <- mat[x$user,x$movie];})
}

# write to result file
write.data <- function(df.test, filename, ...) {write.table(df.test, filename, sep=",", col.names=FALSE, row.names=FALSE, ...}

# compute RMSE
error.data <- function(df.test) {sqrt(mean((df.test$stars - df.test$est)^2))}

# cross-validation
crossval <- function(df, k, alg, ...)
{
	err <- 0; len <- length(df[[1]]);
	for (i in 1:k) {
		// split data set
		df.test <- df[((i-1)*len/k):(i*len/k),];
		df.learn <- df[-(((i-1)*len/k):(i*len/k)),];
		// learn
		mat <- alg(df.learn, ...);
		// estimate on test data
		df.test <- eval(df.test, mat);
		// compute RMSE
		err <- err + error.data(df.test)^2
	}
	sqrt(err/k)
}