#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel void matvecmul(
	__global double* output,
	const unsigned int outputsize,
	__global const double* matrix,
	__global const double* vector,
	const unsigned int inputsize)
{
	size_t i = get_global_id(0);

	double sum = 0.0f;
	__global double *column = matrix + i*inputsize;
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
