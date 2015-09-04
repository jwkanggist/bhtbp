/*======================================================================
 * BPprocess_verC5.c - codes for BPprocess in C language
 *
 * Input parameters (N,M,Niter,MaxdegV,MaxdegC,Nd,zeropos,H, Zindex, degofC, prior_pdf,noise_pdf)
 *      1) N: The signal length
 *      2) M: The number of measurements
 *      3) maxiter: The maximum number of the iterations
 *      4) MaxdegV: the number of the maximum variable node degree with H
 *      5) MaxdegC: the number of the maximum factor node degree with H
 *      6) Nd     : The number of samples for the BP message sampling
 *      7) zeropos: The index of zero in message PDF
 *      8) H:  the sensing matrix H
 *      9) quant_Z: The quantized noisy measurements with Nd
 *      10) degofC: The degree of factor nodes wrt. the sensing matrix H
 *      11) prior_pdf: the sampled density of prior PDF 
 *      12) noise_pdf: the sampled density of zeromean Gaussiance PDF with a certain variance
 * Output: LVp: Marginal posterior PDFs obtained from the BP process
 *
 * copyright@ Jaewook Kang with Gwangju Institute of Science and Technology 2013
 * feedback: jwkkang@gist.ac.kr
 *=======================================================================*/

#include "mex.h"
#include <math.h>
#include <stdio.h>
#include "jw_math_lib.h"
#include "jw_math_lib.c"
#define VTH_CONST 0.0000000001
#define EPSILON 0.00001



/*-----------------main function for BPprocess_verC-------------------------*/    
void BPprocess_verC5(const unsigned N, const unsigned M, const unsigned maxiter, const unsigned MaxdegV, const unsigned MaxdegC, const unsigned Nd, const unsigned zeropos, 
                    const double *H, const double *quant_Z, const double *degofC, const double *prior_pdf, const double *noise_pdf, double *LVp, double *terminate_cond)                    
{
    const unsigned SIZE_OF_DOUBLE=sizeof(double);
    const unsigned SIZE_OF_UNSIGNED=sizeof(unsigned);
    
    double *LVr, *LVq;//  message PDFs
    double *tempPDF;
    unsigned int  *Q1, *Q2;
    unsigned int *a_Q1col,*a_Q2col;

    /*-----------------------------------------------------------------------------------*/
    double *FFTbuff, *FFTacc, *IFFTout;
    
    unsigned targetindex=0, f_index=0, v_index=0;
    unsigned i=0,j=0,k=0,l=0, t=0,iter_index=0;
    unsigned FFTshift=0,maxindex=0, minindex=0, Fdeg=0;
    double Z=0.0,min_of_LVr=0.0,  LVp_diff_norm=0.0, LVp_norm=0.0;
    int  measurement_index=0; 
    
// dynamic memory allocation
    LVq=(double *)mxCalloc(N*MaxdegV*Nd, SIZE_OF_DOUBLE);//  N frames of MaxdegV by Nd matrix
    LVr=(double *)mxCalloc(M*MaxdegC*Nd, SIZE_OF_DOUBLE);//  M frames of MaxdegC by Nd matrix

    
    tempPDF=(double *)mxCalloc(Nd, SIZE_OF_DOUBLE);
    IFFTout=(double *)mxCalloc(Nd, SIZE_OF_DOUBLE);
    
    FFTacc=(double *)mxCalloc(2*Nd+1, SIZE_OF_DOUBLE); // for FFT operations
    FFTbuff=(double *)mxCalloc(2*Nd+1, SIZE_OF_DOUBLE);// for FFT operations

    Q1=(unsigned int*)mxCalloc(MaxdegV*N, SIZE_OF_UNSIGNED);//  N by MaxdegV  matrix 
    Q2=(unsigned int*)mxCalloc(MaxdegC*M, SIZE_OF_UNSIGNED);// M by MaxdegC matrix 
    a_Q2col=(unsigned int*)mxCalloc(MaxdegC,SIZE_OF_UNSIGNED);// 1 by MaxdegC matrix
    a_Q1col=(unsigned int*)mxCalloc(MaxdegV,SIZE_OF_UNSIGNED);// 1 by MaxdegV matrix
    
   /*--------------------- algorithm description------------------------ */

    // initialization of LVr
    for(t=0;t<Nd;t++){for(i=0;i<MaxdegC;i++){for(j=0;j<M;j++){LVr[MaxdegC*Nd*j + Nd*i + t]=1.0;}}}
                
    //1) Generate Q1 Q2 matrices
    gen_Qmatix_inC_Ver2(H, N,  M,  MaxdegV,  MaxdegC, Q1, Q2);

    //2) message passing iterations
    for (iter_index=0;iter_index<maxiter;iter_index++)
    { 
        *terminate_cond=0.0; 
     //   variable nodes to factor nodes message paasing
        for (i=0;i<N;i++)
        {   
            for(j=0;j<MaxdegV;j++)
            {    // matrixA[rowindex+M*colindex] "/*remember MATLAB store matrices in their transposed form*/"                 
                if (Q1[j+MaxdegV*i]==0){break;}// the case of variable nodes which do not have neighoring factor nodes
                else{targetindex=Q1[j+MaxdegV*i];}
               
                for(t=0;t<Nd;t++){tempPDF[t]=prior_pdf[t];}// copy the knowledge of prior PDF first to tempPDF, " tempPDF = prior_pdf"
                
                for(k=0;k<MaxdegV;k++)
                {
                    f_index=Q1[k+MaxdegV*i];// index of the  factor nodes sending message
                                   
                    if (f_index ==0){break;}
                    else
                    {
                        if (f_index != targetindex)
                        {// for summation of message from factor nodes except the target node//
                            //------------------ the part of "v_index = find(Q2(:,f_index)==i);" in MATLAB -----------------------------------//
                            for(l=0;l<MaxdegC;l++){a_Q2col[l]=Q2[l+MaxdegC*(f_index-1)];} // extract "f_index"-th col vector from Q2 matrix 
                            v_index=find_eq(MaxdegC,a_Q2col,i+1); // find index of the variablen node  from the col vector
                            //----------------------------------------------------------------------------------------------------------------//
                            for(t=0;t<Nd;t++){tempPDF[t]=tempPDF[t]*LVr[MaxdegC*Nd*(f_index-1)  + Nd*v_index + t];}// 나중에 mul_vec()함수를 써서 고치기
                            //------------------------------------------------------------------------//
                        }
                    }
                      
                }
                /// normalization of LVq ////
                Z=0.0;
                for(t=0; t<Nd;t++) { Z = Z + tempPDF[t];}// 나중에 sum_vec() 함수 만들어서 고치기
                if(Z!=1.0){for(t=0; t<Nd;t++) {LVq[MaxdegV*Nd*i + Nd*j + t]= tempPDF[t]/Z;}}// LVq update //
            }
            
            // here marginal posterior computation //
            f_index=Q1[1+MaxdegV*i];// index of the  factor nodes sending message
            //------------------ the part of "v_index = find(Q2(:,f_index)==i);" in MATLAB -----------------------------------//
            for(l=0;l<MaxdegC;l++){a_Q2col[l]=Q2[l+MaxdegC*(f_index-1)];} // extract "f_index"-th col vector from Q2 matrix 
            v_index=find_eq(MaxdegC,a_Q2col,i+1); // find index of the variablen node  from the col vector
            //----------------------------------------------------------------------------------------------------------------//
            for(t=0;t<Nd;t++){ tempPDF[t]=LVq[MaxdegV*Nd*i + Nd*1 + t]*LVr[MaxdegC*Nd*(f_index-1)+ Nd*v_index+t]; }
            /// normalization of LVp ////
            Z=0.0;
            for(t=0; t<Nd;t++) { Z = Z + tempPDF[t];}// 나중에 sum_vec() 함수 만들어서 고치기 
            if (Z!=1.0) {for(t=0; t<Nd;t++) { tempPDF[t]= tempPDF[t]/Z;}   }
            
           // termination condition calculation //
             LVp_diff_norm=0.0;
             LVp_norm=0.0;
            for (t=0; t<Nd;t++)
            {  // LVp_diff_norm  = norm(prev_LVp -next_LVp,2)^2
                LVp_diff_norm =  LVp_diff_norm+pow(LVp[Nd*i + t] -tempPDF[t],2);
                LVp_norm      =  LVp_norm + pow(tempPDF[t],2);
                LVp[Nd*i + t]= tempPDF[t];
            }  
            *terminate_cond = *terminate_cond + LVp_diff_norm/LVp_norm;
        }
        // ---------- termination condition check ------------//
        *terminate_cond=  *terminate_cond /N;
//         printf("\n terminate condition %f\n",terminate_cond); 
        if ((*terminate_cond < EPSILON )&&(iter_index > 5)){break;}
// ---------factor nodes to variable nodes message passing ---------------//    
        for(i=0;i<M;i++)
        {
            measurement_index=(int)quant_Z[i];  
            Fdeg=(unsigned)degofC[i];
            
            for(j=0;j<MaxdegC;j++)
            {
                if (Q2[j+MaxdegC*i]==0){break;}
                else{targetindex=Q2[j+MaxdegC*i];}// index of the variable nodes receving message
                

               for(t=0;t<Nd;t++)
               {
                    FFTacc[2*t+1]=noise_pdf[t];
                    FFTacc[2*t+2]=0.0;// the imaginary part should be set to zero
               }//noise_pdf --> FFTacc

               fft_verNRC(FFTacc, Nd, 1);// FFTacc=fft (noise_pdf)
             
                for(k=0;k<MaxdegC;k++)
                { // ----------------------------------- FFT domain calculation --------------------------------------------------------//
                    v_index=Q2[k+MaxdegC*i]; // index of the variable nodes sending message

                    if ( v_index == 0){break;}// except the target variable node
                    else
                    {
                        if(v_index != targetindex)
                        {
                            //------------------ the part of "f_index = find(Q1(:,v_index)==i);" in MATLAB -----------------------------------//
                            for(l=0;l<MaxdegV;l++){ a_Q1col[l]=Q1[l+MaxdegV*(v_index-1)];} // extract "v_index"-th col vector from Q1 matrix
                            f_index=find_eq(MaxdegV,a_Q1col,i+1); // find index of the factor node  from the col vector
                            //------------------ the part of "f_index = find(Q1(:,v_index)==i);" in MATLAB -----------------------------------//
                           
                            for (t=0;t<Nd;t++) 
                            {
                                FFTbuff[2*t+1]=LVq[MaxdegV*Nd*(v_index-1) + Nd*f_index+t];
                                FFTbuff[2*t+2]=0.0;
                            }// LVq --> FFTbuff

                            fft_verNRC(FFTbuff, Nd, 1);// FFTbuff=fft (LVq) by from the Numerical recipes in C

                            for(t=0;t<Nd;t++)
                            {// elementwise multiplication  of only real part
                                FFTacc[2*t+1]= FFTacc[2*t+1]*FFTbuff[2*t+1]; // FFTacc = FFTacc * FFTbuff
                            }
                        }
                    }
                }
//--------------------------------------- The end of FFT domain calculation  --------------------------------------------------------//  
              
                 fft_verNRC(FFTacc, Nd, -1);// ifft function from the Numerical recipes in C
                 for(t=0;t<Nd;t++)
                 {
                     IFFTout[t]=FFTacc[2*t+1]/Nd;//normalization for IFFT
                     tempPDF[t]=0.0;// initialization to zeros
                 }

                FFTshift = (unsigned)!(Fdeg%2)*(zeropos-1);
                if (FFTshift != 0)
                {//------- shift adjustment due to FFT convolution --------//
                    for (t=0;t<FFTshift;t++)
                    {
                        tempPDF[Nd-FFTshift+t]=IFFTout[t];
                        tempPDF[t]=IFFTout[FFTshift+t];
                    }
                }
                else{for(t=0;t<Nd;t++){tempPDF[t]=IFFTout[t];}}

               //  ---------------------------------------------------------//   
 
                // confind dynamic range of message PDF
                for (t=0;t<Nd;t++){ if (tempPDF[t] <= VTH_CONST ){tempPDF[t] = VTH_CONST;}}
//                 minindex=find_min_vec(tempPDF, Nd);
//                 min_of_LVr=tempPDF[minindex];
                if (min_of_LVr < VTH_CONST ){min_of_LVr=VTH_CONST;}
                           
//                 // shift the peak of message PDF given noisy measurements Z
                if (measurement_index  > 0)
                {
                     for (t=Nd-1;t>measurement_index-1;t--){tempPDF[t]=tempPDF[t-measurement_index];}
//                      for (t=0;t<measurement_index;t++){tempPDF[t]=min_of_LVr;}
                     for (t=0;t<measurement_index;t++){tempPDF[t]=tempPDF[0];}
                }
                else if (measurement_index < 0 )
                {
                    for (t=0; t< Nd+measurement_index;t++){tempPDF[t]= tempPDF[t-measurement_index];}
//                     for (t=Nd+measurement_index; t<Nd ;t++) {tempPDF[t]=min_of_LVr;}
                    for (t=Nd+measurement_index; t<Nd ;t++) {tempPDF[t]=tempPDF[Nd-1];}
                }

//                 // LVr update -- 0.5 times old value + 0.5 times updated value //
                for (t=0;t<Nd;t++) {LVr[MaxdegC*Nd*i+ Nd*j+t]=tempPDF[t]*0.5 + LVr[MaxdegC*Nd*i+ Nd*j+t]*0.5;}
            }

        }

   }
 //   printf("\n The number of performing BP iterations : %d \n",iter_index); 
       
    // Dynamic memory release    
    mxFree(LVq);
    mxFree(LVr);
    
    mxFree(tempPDF);
    mxFree(IFFTout);
    
    mxFree(FFTacc);
    mxFree(FFTbuff);
  
           
    mxFree(Q1);
    mxFree(Q2);
    mxFree(a_Q2col);      
    mxFree(a_Q1col);
   
}
 

         
/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    /* Remember !!! MATLAB stores matrices in their transposed form
     * i.e., columnwise. 
     */
    
/* nlhs : the number of output parameters
 * plhs : outputs
 * nrhs : the number of input parameters
 * prhs : input
 **/
    // input variable declear 
    unsigned MaxdegV, MaxdegC;             /* input scalar */
    unsigned maxiter,zeropos;             /* input scalar */

    double *matrixH;                   /* MxN input measurement matrix */
    double *prior_pdf;                 /* 1xNd quantized prior density */ 
    double    *quant_Z;            /* 1xM quantized measurement vector */
    double    *degofC;            /* 1xM the array of degree of factor nodes */
    double *noise_pdf;     /* 1xNd quantized noise density */ 
    mwSize M,N,Nd;                         /*Size of  matrix H */
    
    // output variable declear
    double *LVp;
    double *terminate_cond;  

    /* check for proper number of arguments */
    if(nrhs!=9) {//inputs
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","Nine inputs required.");
    }
    if(nlhs!=2) {// outputs
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","two outputs required.");
    }
     /*(maxiter,MaxdegV,MaxdegC,Nd,zeropos, quant_Z,degofC,f_x,noise_pdf_FFT)*/
    //---------------- For inputs -------------------------------------//
    maxiter = (unsigned)mxGetScalar(prhs[0]);/* create a pointer to the real data in the input matrix  */
    MaxdegV =(unsigned)mxGetScalar(prhs[1]);
    MaxdegC =(unsigned)mxGetScalar(prhs[2]);
    zeropos =(unsigned)mxGetScalar(prhs[3]);
    
    matrixH  =  mxGetPr(prhs[4]);
    quant_Z=   mxGetPr(prhs[5]);
    degofC=    mxGetPr(prhs[6]);
    prior_pdf=  mxGetPr(prhs[7]);
    noise_pdf = mxGetPr(prhs[8]);//??????
   
    N=(mwSize)mxGetN(prhs[4]);/* get row dimension of the inputX */ 
    M=(mwSize)mxGetM(prhs[4]);/* get col dimension of the inputX */
    Nd=(mwSize)mxGetM(prhs[7]);
    

    //---------------- For outputs ----------------------------------//
    plhs[0] = mxCreateDoubleMatrix(Nd,N,mxREAL);/* create the output matrix */
    plhs[1] = mxCreateDoubleScalar(mxREAL);/* create the output matrix */ 

    /* get a pointer to the real data in the output matrix */
    LVp = mxGetPr(plhs[0]); 
    terminate_cond = mxGetPr(plhs[1]);
    //-----------------call the computational routine--------------------//
    BPprocess_verC5((unsigned)N, (unsigned)M, maxiter, MaxdegV, MaxdegC, (unsigned)Nd,  zeropos, matrixH, quant_Z, degofC, prior_pdf, noise_pdf, LVp, terminate_cond);
}
