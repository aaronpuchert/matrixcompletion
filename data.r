### INPUT PROCESSING ###
# read csv data
read.data <- function (filename, ...)
{
	tab <- read.table(filename, sep=",", header=FALSE, ...);
	if (dim(tab)[2] == 2) colnames(tab) <- c("user", "movie") else colnames(tab) <- c("user", "movie", "stars");
	class(tab) <- c("data", "data.frame");
	tab
}

# reduce user and movie domain
reduce.data <- function (df, ui = sort(union(df$user, NULL)), mi = sort(union(df$movie, NULL)))
{
	df$user <- sapply(df$user, function(x) {which(ui==x)});
	df$movie <- sapply(df$movie, function(x) {which(mi==x)});
	list(data=df,ui=ui,mi=mi)
}

expand.data <- function (reduced, ui, mi)
{
	data.frame(user=ui[reduced$user], movie=mi[reduced$movie], reduced[-(1:2)]);
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
eval <- function(df.test, mat)
{
	new <- apply(df.test, 1, function(x) {c(x, est=mat[x[[1]],x[[2]]]);});
	as.data.frame(t(new))
}

# write to result file
write.data <- function(df.test, filename, ...) {write.table(df.test, filename, sep=",", col.names=FALSE, row.names=FALSE, ...)}

# compute RMSE
error.data <- function(df.test) {sqrt(mean((df.test$stars - df.test$est)^2))}

# cross-validation
crossval <- function(df, m, alg, ...)
{
	err <- 0; len <- length(df[[1]]);
	for (i in 1:m) {
		# split data set
		df.test <- df[((i-1)*len/m):(i*len/m),];
		df.learn <- df[-(((i-1)*len/m):(i*len/m)),];
		# learn
		mat <- alg(df.learn, ...);
		# estimate on test data
		df.test <- eval(df.test, mat);
		# compute RMSE
		err <- err + error.data(df.test)^2
	}
	sqrt(err/m)
}