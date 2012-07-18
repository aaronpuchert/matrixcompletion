### INPUT PROCESSING ###
# read csv data
read.data <- function (filename, ...)
{
	tab <- read.table(filename, sep=",", header=FALSE, ...);
	if (ncol(tab) == 2) colnames(tab) <- c("user", "movie") else colnames(tab) <- c("user", "movie", "stars");
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
plot.data <- function(df, rg=c(1,5), col=3, ...)
{plot.default(df[c(1,2)], pch=".", col = hsv(h=(df[[col]]-rg[1])/(rg[2]-rg[1])), ...)}

# create data matrix
matrix.data <- function (df, init=NA)
{
	mat <- matrix(init, max(df[[1]]), max(df[[2]]));
	mat[(df$movie-1)*nrow(mat)+df$user] <- df$stars;
	mat
}

### EVALUATION & OUTPUT ###
# evaluate matrix on test data set
eval <- function(df, mat)
{
	df <- cbind(df, est=NA);
	df$est <- mat[(df$movie-1)*nrow(mat)+df$user]; df
}

# write to result file
write.data <- function(df.test, filename, ...) {write.table(df.test, filename, sep=",", col.names=FALSE, row.names=FALSE, ...)}

# compute error (& analyze it)
error.data <- function(df.test)
{
	list(error = sqrt(mean((df.test$stars - df.test$est)^2)),
		lm = lm(est ~ stars, data=df.test)$coefficients)
}

# cross-validation
crossval <- function(df, m, alg, ...)
{
	len <- nrow(df); df <- cbind(df, est=NA);
	perm <- sample(len);		# permute randomly to gain fair results
	for (i in 1:m) {
		# select test data
		ind <- perm[ceiling((i-1)*len/m):floor(i*len/m)];
		# learn
		mat <- alg(df[-ind,], ...);
		# estimate on test data
		df$est[ind] <- mat[(df$movie[ind]-1)*nrow(mat) + df$user[ind]];
	}

	# return cross-validated data set
	return(df)
}

### META-ROUTINES ###
read.all <- function(trainfile, qualfile)
{
	train <- read.data(trainfile);
	newtrain <- reduce.data(train);
	qual <- read.data(qualfile);
	newqual <- reduce.data(qual, newtrain$ui, newtrain$mi);
	list(train=newtrain$data, qual=newqual$data, ui=newtrain$ui, mi=newtrain$mi)
}

write.res <- function(all, mat, resfile)
{
	qual <- eval(all$qual, mat);
	qual <- expand.data(qual, all$ui, all$mi);
	write.data(qual, resfile);
}
