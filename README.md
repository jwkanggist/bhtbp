# bhtbp
 This project is for performing compressed sensing reconstruction of sparse signal via the BHT-BP solver.
 
 
 README for BHT_BP solver by Jaewook Kang (jwkkang@gist.ac.kr) 
=============================================================

Description: 
=========================================================================
 This package contains methods for performing compressed sensing reconstruction of sparse signal via the BHT-BP solver.

-Jaewook Kang, Heung-No Lee, Kiseon Kim, "Bayesian Hypothesis Test using Nonparametric Belief 
Propagation for Noisy Sparse Recovery,"  IEEE Trans. on Signal process., vol. 63, no. 4,  pp. 935-948, Feb. 2015

 This package shows a numerical comparison of solvers for compressed sensing recovery. 
 The solvers included in this comparison are as given below:

	-CS-BP : D. Baron, S. Sarvotham, and R. Baraniuk, "Bayesian compressive sensing via belief propagation,
          IEEE Trans. Signal Process., vol. 58, no. 1, pp. 269-280, Jan. 2010.
          
	-BCS: Shihao Ji, Ya Xue, and Lawrence Carin, "Bayesian compressive sensing," 
       IEEE Trans. Signal Process., vol. 56, no. 6, pp. 2346-2356, June. 2008.
       
 -SuPrEM: M. Akcakaya, J. Park, and V. Tarokh, “A coding theory approach to noisy compressive 
       sensing using low density frame,” IEEE Trans. Signal Process., vol. 59, no. 12, pp. 5369-5379, Nov. 2011.
       
 -L1-DS: E. Candes and T. Tao, “The Dantzig selector: Statistical estimation when p is much larger than n,” 
            Ann. Statist., vol. 35, no. 6, pp. 2313?2351, 2007
            
==========================================================================

The BHT-BP solver includes many input and output parameter, which is described below. 
For the other solvers, please see each author's website or the corresponding papers.


Matlab function for the BHT-BP solver :
=======================================
[estX_BHT,estX_CSBP,Merror_CSBP,Merror_BHT,state_err_BHT,state_err_CSBP,telapsed]
=BHT_BP_ver3(trueX,H,Z,q,maxiter,sigma_N,sigmaX,Xstate,Nd,Xmin)
   
 Input parameters
====================
* trueX - This is the true value of the sparse signal X, which is used to calculate MSE      performance of the reconstruction
* H - This is a measurement matrix (M by N)
* Z- This is a noisy measurment vector (M by 1).  
* q - The sparsity rate which is a parameter of the prior PDF
* maxiter - The maxmum number of belief propagation iterations in BHT-BP
* sigma_N - The standard deviation for additive measurement noise
* sigma_X - The standard deviation for signal, which is a parameter of the prior PDF
* Xstate - The binary state vector indicating the location of nonzero value in the signal. 
  
-Output parameters
====================
* estX_BHT- The BHT-BP signal estimate
* estX_CSBP -The CS-BP signal estimate 
* Merror_CSBP -  The MSE performance of the BHT-BP recovery  
* Merror_BHT - The MSE performance of the CS-BP recovery
* state_err_BHT-  The number of state errors in the detected signal support by the BHT-BP recovery
* state_err_CSBP - The number of state errors in the detected signal support by the CS-BP recovery
* telapsed - the running time expending for the belief propagation iterations
=========================================================================

Size: 
The total size of the file is 1.12 MB
Player Information:
"	MATLAB Version: 8.2.0.701 (R2013b)
"	Operating System: Microsoft Windows 7 Version 6.1 (Build 7601: Service Pack 1)
"	Java Version: Java 1.7.0_11-b21 with Oracle Corporation Java HotSpot(TM) 64-Bit Server VM mixed mode

Packing List: 
 	demo_BHT_BP.m: provides an exemplary sparse signal   recovery among BHT-BP, CS-BP, and BCS solvers.
 	demo_ImagRev.m: providing image (cameraman128x128) recovery demonstration using BHT-BP, CS-BP, and BCS solvers.
    BHT_BP_Fig10_MSE_wrt_Nd_ver1.m: The experment setup for Fig 10 in the corresponding paper.
    BHT_BP_Fig11a_MSE_wrt_L_M512.m : The experment setup for Fig 11-(a) in the corresponding paper.  
    BHT_BP_Fig11b_MSE_wrt_L_M768.m : The experment setup for Fig 11-(a) in the corresponding paper. 
    BHT_BP_Fig8n9_MSE_wrt_SNR_ver2.m: The experment setup for Fig 8 and Fig 9 in the corresponding paper. 
   
Contact Information:

"	Lab. phone: +82-62-715-2264
"	E-mail        :  jwkkang@gist.ac.kr, jwkang10@gmail.com
"   Final Update @ May 2015
