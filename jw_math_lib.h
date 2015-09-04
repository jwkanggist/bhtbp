#ifdef __JW_MATH_LIB__
#define __JW_MATH_LIB__

// the basic vector operation 
//-----------------------------------------------------------------------------
// element by element vector addition
extern inline void add_vec(const double *in1, const double *in2, double *out, unsigned L);

// element by element vector multiplication 
extern inline void mul_vec(const double *in1, const double *in2, double *out, unsigned L);

//=========================================================================
// prod_vec(): This function for product of element by element vector multiplication
// **in1 : each of in1[1],...,in[n] indicate the input vectors (in1[i] is L x 1 vector)
// *out: output vector ( out[i] is L x 1 vector)
// L   : the number of vector in the product
// N   : the length of each vector in the product
//=========================================================================
extern inline void prod_vec(const double **in1, double *out, unsigned L, unsigned N);

extern inline void find_less_eq(const double *in, const double th_value, const unsigned in_length, unsigned* out_index, unsigned out_length);// <=
extern inline void find_eq     (const double *in, const double th_value, const unsigned in_length, unsigned* out_index, unsigned out_length);// ==
extern inline unsigned find_max_vec(const double *in, const unsigned in_length);
extern inline unsigned find_min_vec(const double *in, const unsigned in_length);

//=========================================================================
// function name find_eq()
/* inputs :1) object vector X
 *         2) size of the vector X: N
 *         3) the target value: value
 * outputs: 1) index of X corresponding to "value"
 * objective : the function try to find the index of X corresponding to "value"
========================================================================== */
extern inline int find_eq(const int N, const unsigned int *in, const unsigned int value);

//=========================================================================
// function name find_nonzero()
/* inputs :1) object vector X
 *         2) size of the vector X: N
 * outputs: 1) output array of nonzero index 
 *         - if output_index[0]=0 then X does not include nonzero elements
 *           2) actual size of output_index[]
 * objective: the function try to find the index array of nonzero element in X
 =============================================================================*/
extern inline int find_nonzero(const double *in, const int N, unsigned int *output_index);

extern inline void linear_rightshift(double* in_array, const unsigned in_length, const unsigned shift, double blank_value);
extern inline void linear_leftshift(double* in_array, const unsigned in_length, const unsigned shift, double blank_value);

/*==========================================================================
* FFT code from the book Numerical Recipes in C *
* FFT/IFFT routine. (see pages 507-508 of Numerical Recipes in C)
 Inputs:
	data[] : array of complex* data points of size 2*NFFT+1.
		data[0] is unused,
		* the n'th complex number x(n), for 0 <= n <= length(x)-1, is stored as:
			data[2*n+1] = real(x(n))
			data[2*n+2] = imag(x(n))
		if length(Nx) < NFFT, the remainder of the array must be padded with zeros

	nn : FFT order NFFT. This MUST be a power of 2 and >= length(x).
	isign:  if set to 1, 
				computes the forward FFT
			if set to -1, 
				computes Inverse FFT - in this case the output values have
				to be manually normalized by multiplying with 1/NFFT.
 Outputs:
	data[] : The FFT or IFFT results are stored in data, overwriting the input.
================================================================================*/
extern inline void fft_verNRC(double data[], int nn, int isign);

extern void gen_Qmatix_inC_Ver2(const double *matrixA, const int N, const int M, const int MaxdegV, const int MaxdegC, unsigned int *Q1, unsigned int *Q2);


     
# endif //__JW_MATH_LIB__
