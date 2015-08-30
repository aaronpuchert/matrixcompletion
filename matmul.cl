__kernel void matvecmul(
	__global float* output,
	const unsigned int outputsize,
	__global float* matrix,
	__global float* vector,
	const unsigned int inputsize)
{
	int i = get_global_id(0);

	float sum = 0.0f;
	__global float *column = matrix + i*inputsize;
	for (int j = 0; j < inputsize; ++j)
		sum += column[j] * vector[j];

	output[i] = sum;
}
