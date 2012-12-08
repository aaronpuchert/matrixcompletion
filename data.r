# read csv data
read.data <- function (filename, ...)
{
	tab <- read.table(filename, sep=",", header=FALSE, col.names = c("user", "movie", "stars"), ...);
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