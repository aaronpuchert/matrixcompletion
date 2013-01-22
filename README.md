MATRIX COMPLETION
=================

Purpose
-------
This is an implementation of different algorithms for matrix completion. Given a partially known
matrix, we want to guess the remaining entries.

This can be used for recommendation systems, by predicting user ratings. It was partially inspired
by the [Netflix Prize][Netflix]. In this case, we have a matrix of ratings *(a_ij)* of user *i* for
movie *j*, of which only some are given. Computing the others gives good hints: just recommend user
*i* the best-rated movies *j* he has not yet watched.

Theory
------
There is no theory for the SVD algorithm, but the idea is rather straightforward: we want the best
rank *k* approximation. We fill the matrix with some arbitrary, but constant value (like the mean
value of the given entries) and then compute the singular value decomposition. By taking only the
singular vectors belonging to the *k* largest singular values, we compute the best rank *k*
approximation. Now we iterate these steps to fill the matrix with appropriate values. Seems to work
though.

The theory for Hazan's algorithm is explained in [Giesen et al][Hazan]. The actual algorithm is
much simpler than the math behind, surprisingly.


Usage
-----
The data should be provided in CSV format, each line encoding one matrix entry:

    i,j,a_ij

To evaluate, just provide another data set with just `i,j` lines, you'll then get the computed `a_ij`
as output. Read in the data with `read.data`. You can also use `read.all` for reading in both files
at once. If your data is in another format, or already read into R, don't bother. Just do what
`read.data` does, then there will be no problems. If many columns or rows in your matrix are unused,
use `reduce.data` to 'compress' it (and later `expand.data` to decompress). `read.all` also does this.

To use the algorithms, first set their parameters. There are functions `alg.*.pl` supporting you here.
They produce 'parameter lists', which then can be passed on to the algorithm. The algorithms return the
completed matrix, which you might save. Or evaluate it on the test set by `eval`.

For parameter tuning, validation etc., there is a built-in cross validation routine. (It is in fact a
'cross evaluation'.) The results can be analyzed later. (by `error.data`, for example) There is also a
neat little tool, which I called bisection. This can be used to tweak parameters. It takes a long time
though, and some parameters are just easily 'guessed'.

Results
-------
On real world data, we got a root mean squared error of about 0.4 to 0.5 on a range of ratings from 1
to 5. The best results came from combining the output of both algorithms. Especially the unsophisticated
SVD algorithm showed surprisingly good results. (about 24% of the entries were given after compression)

A minimal example
-----------------

    > source("data.r")
    > source("alg.r")
    > data<-read.all("train.csv", "qualifiying.csv")
    > params<-alg.svd.pl(data$train, k=20)
    > mat<-alg.svd(data$train, params)
    [1] 0.4468885 0.4513283 0.4576820 0.4675882 0.4852974 0.5279047 0.8372528 4.0000000 4.0812162
    > write.res(data, mat, "result.csv")

The output is the (reversed) error history on the training data set. Thus we needed 7 iterations in this
case.

[Netflix]: http://www.netflixprize.com//index "The Netflix Prize"
[Hazan]:   http://theinf2.informatik.uni-jena.de/theinf2_multimedia/Publications/grow_specM.pdf "Joachim Giesen, Martin Jaggi, and Sören Laue: Optimizing over the Growing Spectrahedron"