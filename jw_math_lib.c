/* -------------------------------------------------------------------------- *
	FILE: jw_math_lib.c
	AUTH:JwKang
	2011.10.30 

	Description: This library includes some math. function to implement BHT-BP in C


	Revision: 2013.June
	
* -------------------------------------------------------------------------- */
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include "mex.h"
#include "jw_math_lib.h"

#define _USE_MATH_DEFINES
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr // for FFT
// The following line must be defined before including math.h to correctly define M_PI
#define PI	3.14159265358979323	/* pi to machine precision, defined in math.h */
#define TWOPI	(2.0*PI)

////////// Basic vector operation ///////////////////////////////
// element by element vector addition
 void add_vec(const double *in1, const double *in2, double *out, unsigned L)
{
  unsigned i=0;
  for ( i=0;i< L;i++)
    out[i]=(in1[i]+in2[i]);
}
// element by element vector multiplication 
 void mul_vec(const double *in1, const double *in2, double *out, unsigned L)
{ 
  unsigned i=0;
  for (i=0;i<L;i++)
    out[i]=in1[i]*in2[i]; 
}

/********************************************************************
/ prod_vec(): This function for product of element by element vector multiplication
/ **in1 : each of in1[1][],...,in[N][] indicate the input vectors (in1 is N x L vector)
/ *out: output vector ( out is L x 1 matrix)
/ L   : the number of vector in the product
/ N   : the length of each vector in the product
  in[0][0] ,.................,in[0][L]
  in[1][0],..................,in[1][L]
       .                          .
       .                          .
       .                          .
X in[N][0],..................,in[N][L]
------------------------------------------
  out[0],.....................,out[L]

*********************************************************************/
void prod_vec(const double **in1, double *out, unsigned L, unsigned N)
{
    unsigned i=0, j=0;  
    for (i=0; i < L ; i++)
    {
        out[i]=1.0;
        for ( j=1; j < N ; j++)
            out[i]=out[i]*in1[j][i];
    }
}

 void find_less_eq(const double *in, const double th_value, const unsigned in_length, unsigned* out_index, unsigned out_length)
{
    unsigned i=0, j=0;
    
    for (i=0;i < in_length ; i++)
    {
        if (in[i] <= th_value)
            out_index[j++]=i;
    }
    out_length=j;
}
  
 void find_less(const double *in, const double th_value, const unsigned in_length, unsigned* out_index, unsigned out_length)
{
    unsigned i=0, j=0;
    
    for (i=0;i < in_length ; i++)
    {
        if (in[i] == th_value)
            out_index[j++]=i;
    }
    out_length=j;
}

unsigned find_max_vec(const double *in, const unsigned in_length)
{
    unsigned i=0, max_index=1;
    double max_value=0.0;
    
    max_value=in[0];


    for (i=1 ; i < in_length ; i++)
    {
        if (in[i] > max_value)
        {
            max_value=in[i];
            max_index=i+1;
        }
    }
    return (max_index);
}


 unsigned find_min_vec(const double *in, const unsigned in_length)
{
    unsigned i=0, min_index=1;
    double min_value=0.0;
    
    min_value=in[0];


    for (i=1 ; i < in_length ; i++)
    {
        if (in[i] < min_value)
        {
            min_value=in[i];
            min_index=i+1;
        }
    }
    return (min_index);
}

/* ========================================================================*/
 void linear_rightshif(double* in_array, const unsigned in_length, const unsigned shift, double blank_value)
{
    unsigned i=0;
    
    for (i=in_length ; i > in_length - shift -1 ; i--)
        in_array[i-1]=in_array[i-shift-1];
   
    for (i=0; i < shift ; i++)
        in_array[i]=blank_value;    
}

    
 void linear_leftshif(double* in_array, const unsigned in_length, const unsigned shift, double blank_value)
{
    unsigned i=0;
    
    for (i=0; i < in_length - shift ; i++)
        in_array[i]=in_array[i+shift];
   
    for (i=shift+1; i < in_length ;i++)
        in_array[i]=blank_value;    
}
//=========================================================================
// function name find_nonzero()
/* inputs :1) object vector X
 *         2) size of the vector X: N
 * outputs: 1) output array of nonzero index 
 *         - if output_index[0]=0 then X does not include nonzero elements
 *           2) actual size of output_index[]
 * objective: the function try to find the index array of nonzero element in X
 =============================================================================*/
int find_nonzero(const double *in, const unsigned N, unsigned int *output_index)
{
    unsigned i=0,k=0;
    for (i=0;i<N;i++)
    {
        if (in[i]!=0)
            output_index[k++]=i;
    }
    return (k);
}
//=========================================================================
// function name find_eq()
/* inputs :1) object vector X
 *         2) size of the vector X: N
 *         3) the target value: value
 * outputs: 1) index of X corresponding to "value"
 * objective : the function try to find the index of X corresponding to "value"
========================================================================== */
int find_eq(const unsigned N, const unsigned int *in, const unsigned int value)
{
    unsigned i=0;
    int re=-1;
    
    for (i=0;i<N;i++)
    {
        if (in[i]==(int)value)
        {
            re=i;
            break;
        }
    }
    return(re);// there are no index mathcing with value" return -1
}
//=========================================================================        
// function name gen_Qmatrix()
/* inputs : 1) M by N matrix A
 *          2) size of matrix A: M,N
 *          3) maximun degree of A: MaxdegV, MaxdegC
 * outputs : 1) Q1 matrix to represent graphical relation of A (MaxdegV by N matrix )
 *           2) Q2 matrix to represetn graphical relatio of A (MaxdegC by M matrix )
 * objective : the function try to generate Q1 and Q2 function given the matrix A
 */     

void gen_Qmatix_inC_Ver2(const double *matrixA, const int N, const int M, const int MaxdegV, const int MaxdegC, unsigned int *Q1, unsigned int *Q2)
{
    double *temp_colA;
    double *temp_rowA;
    
    unsigned int *indexofcolA;
    unsigned int *indexofrowA;
    int i=0,j=0;
    int nonzero_num_colA=0,nonzero_num_rowA=0;

    
    /* temporal memory allocation */
    // the mxCalloc function should be replaced by ANSI C when you transplant this code to the other applications
    temp_colA=(double*)mxCalloc(M,sizeof(double)); 
    temp_rowA=(double*)mxCalloc(N,sizeof(double)); 
    
    indexofcolA=(unsigned int*)mxCalloc(MaxdegV,sizeof(unsigned int));
    indexofrowA=(unsigned int*)mxCalloc(MaxdegC,sizeof(unsigned int));

       //  printf("%f",matrixA[24]);
    for(i=0;i<N;i++)
    {
        /* ----------- generation of Q1 matrix ---------------------*/
        for (j=0;j<M;j++){
            temp_colA[j]=matrixA[j+M*i]  ; // creat the i-th colunm vector of A
        }
        
        nonzero_num_colA=find_nonzero(temp_colA, M, indexofcolA);// find nonzero index of the i-th col of A
        
        
        for (j=0;j<nonzero_num_colA;j++){
            Q1[j+MaxdegV*i]=indexofcolA[j]+1;}
        /* ----------- generation of Q2 matrix ---------------------*/    
        if (i < M)
        {
            for (j=0;j<N;j++)
                temp_rowA[j]=matrixA[i+M*j]  ; // creat the j-th row vector of A
            
            nonzero_num_rowA=find_nonzero(temp_rowA, N, indexofrowA);// find nonzero index of the i-th row of A
            for (j=0;j<nonzero_num_rowA;j++)
               Q2[j+MaxdegC*i]=indexofrowA[j]+1;
        }
    }
    
    mxFree(temp_colA);
    mxFree(temp_rowA);
    mxFree(indexofcolA);
    mxFree(indexofrowA);
        
} 
 
 
//----------------------------------------------------------------------------------------------------------------
/************************************************
* FFT code from the book Numerical Recipes in C *
* Visit www.nr.com for the licence.             *
************************************************/
/*
 FFT/IFFT routine. (see pages 507-508 of Numerical Recipes in C)

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
*/

void fft_verNRC(double data[], int nn, int isign)
{
    int n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;
    
    n = nn << 1;
    j = 1;
    for (i = 1; i < n; i += 2) {
	if (j > i) {
	    tempr = data[j];     data[j] = data[i];     data[i] = tempr;
	    tempr = data[j+1]; data[j+1] = data[i+1]; data[i+1] = tempr;
	}
	m = n >> 1;
	while (m >= 2 && j > m) {
	    j -= m;
	    m >>= 1;
	}
	j += m;
    }
    mmax = 2;
    while (n > mmax) {
	istep = 2*mmax;
	theta = TWOPI/(isign*mmax);
	wtemp = sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi = sin(theta);
	wr = 1.0;
	wi = 0.0;
	for (m = 1; m < mmax; m += 2) {
	    for (i = m; i <= n; i += istep) {
		j =i + mmax;
		tempr = wr*data[j]   - wi*data[j+1];
		tempi = wr*data[j+1] + wi*data[j];
		data[j]   = data[i]   - tempr;
		data[j+1] = data[i+1] - tempi;
		data[i] += tempr;
		data[i+1] += tempi;
	    }
	    wr = (wtemp = wr)*wpr - wi*wpi + wr;
	    wi = wi*wpr + wtemp*wpi + wi;
	}
	mmax = istep;
    }
}

