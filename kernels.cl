#pragma OPENCL EXTENSION cl_khr_fp64 : enable

/**
 * Compute ifelse(target != 0, (matrix-target)^2, 0) for two square matrices and
 * compute the sum over the rows/columns.
 */
__kernel void error(
	__global double* output,
	const unsigned int outputsize,
	__global const double* matrix,
	__global const double* target)
{
	size_t i = get_global_id(0);

	double sum = 0.0;
	__global const double *matrixCol = matrix + i*outputsize;
	__global const double *targetCol = target + i*outputsize;
	for (int j = 0; j < outputsize; ++j)
		if (targetCol[j] != 0) {
			double diff = matrixCol[j] - targetCol[j];
			sum += diff * diff;
		}

	output[i] = sum;
}

/**
 * Compute the matrix 2 * ifelse(target != 0, matrix - target, 0).
 */
__kernel void nabla(
	__global double* output,
	const unsigned int outputsize,
	__global const double* matrix,
	__global const double* target)
{
	size_t i = get_global_id(0);

	if (target[i] != 0)
		output[i] = 2 * (matrix[i] - target[i]);
	else
		output[i] = 0.0;
}

/**
 * Blend the given matrix with an outer product of a vector.
 */
__kernel void blend(
	__global double* output,
	const unsigned int outputsize,
	__global const double* matrix,
	double matrixFactor,
	__global const double* vector,
	double vectorFactor,
	unsigned int length)
{
	size_t i = get_global_id(0);
	size_t row = i / length, col = i % length;

	output[i] = matrixFactor * matrix[i] + vectorFactor * vector[row] * vector[col];
}

__kernel void matvecmul(
	__global double* output,
	const unsigned int outputsize,
	__global const double* matrix,
	__global const double* vector,
	const unsigned int inputsize)
{
	size_t i = get_global_id(0);

	double sum = 0.0;
	__global const double *column = matrix + i*inputsize;
	for (int j = 0; j < inputsize; ++j)
		sum += column[j] * vector[j];

	output[i] = sum;
}

__kernel void vecdiv(
	__global double* output,
	const unsigned int outputsize,
	__global const double* input,
	double factor)
{
	size_t i = get_global_id(0);

	output[i] = input[i] / factor;
}
