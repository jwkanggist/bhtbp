%-------------------------------------------------------------------------
% filename :BHT_BP_ver3.m
% objective :This is a compressive sensing recovery algorithm named
% Bayesian hypothesis test and belief propagation (BHT-BP)BHT_BP
%
% - input parameters
%   trueX - This is the true value of the sparse signal X, which is used to calculate MSE performance of the reconstruction 
%   H - This is a measurement matrix (M by N)
%   Z- This is a noisy measurment vector (M by 1).  
%   q - The sparsity rate which is a parameter of the prior PDF
%   maxiter - The maxmum number of belief propagation iterations in BHT-BP
%   sigma_N - The standard deviation for additive measurement noise 
%   sigma_X - The standard deviation for signal, which is a parameter of the prior PDF  
%   Xstate - The binary state vector indicating the location of nonzero value in the signal. 
%       
% - output parameters
%   estX_BHT- The BHT-BP signal estimate 
%   estX_CSBP -The CS-BP signal estimate 
%   Merror_CSBP -  The MSE performance of the BHT-BP recovery    
%   Merror_BHT - The MSE performance of the CS-BP recovery  
%   state_err_BHT-  The number of state errors in the detected signal support by the BHT-BP recovery  
%   state_err_CSBP - The number of state errors in the detected signal support by the CS-BP recovery
%   telapsed - the running time expending for the belief propagation iterations
% --------------------------------------------------------------------------------
% Written by: Jaewook Kang in GIST Korea
% Email: jwkkang@gist.ac.kr
% Created: July 2011, updated using C code at 2013 Apr. updated at 2014 July
%--------------------------------------------------------------------------
function  [estX_BHT,estX_CSBP,Merror_CSBP,Merror_BHT,state_err_BHT,state_err_CSBP,telapsed]...
    =BHT_BP_ver3(trueX,H,Z,q,maxiter,sigma_N,sigmaX,Xstate,Nd,Xmin)

% ----------- preallocation  ----------------------%
[M,N]=size(H);
Psupp=zeros(N,1);
Pnonsupp=zeros(N,1);
Xstate_hatBHT=zeros(N,1);
Xstate_hatBP=zeros(N,1);
estX_BHT=zeros(N,1);
estX_CSBP=zeros(N,1);
f_nonsupp=zeros(Nd,1);
f_supp=zeros(Nd,1);
f_x=zeros(Nd,1);
ref_function_supp=zeros(Nd,1);
ref_function_nonsupp=zeros(Nd,1);
Xindex=zeros(Nd,1);
%------------------------------------------------%


% parameters for message sampling and truncation 
sigmagate=sigmaX*3;% use the 3 sigma rule;
t=2*sigmagate/Nd;% sampling stepsize
Xindex=[-sigmagate:t:t*(Nd/2-1)]';
zeropos=Nd/2+1;% the location of zero in the pdf
absmax_Z=(Nd-1)*t;
BHT_th=log(q/(1-q));% The threshold for the likelitest of Bayesian support detection



% ------priori distirbution on x------
f_supp=pdf2('Normal',Xindex,0,sigmaX);

pXmin_index=round(Xmin/t)+zeropos-1;
nXmin_index=round(-Xmin/t)+zeropos+1;

while (pXmin_index-zeropos)*t > Xmin
    pXmin_index=pXmin_index-1;
    nXmin_index=nXmin_index+1;
end
f_supp(nXmin_index:pXmin_index)=1e-3;

f_nonsupp(zeropos)=1;
f_x=f_supp*q+f_nonsupp*(1-q);
f_x=f_x/sum(f_x);
% f_x(zeropos)=1-q;
% 
% figure(1)
% semilogy(Xindex,f_x)
%---- reference function for BHT----%
if Xmin==0
    epsilon=eps;
else
    epsilon=Xmin;
end
    
f_nonsuppBHT=pdf2('Normal',Xindex,0,epsilon/12);% use 6-sigma gate for prior for nonsupport set
f_xBHT=f_supp*q+f_nonsuppBHT*(1-q);
f_xBHT=f_xBHT/sum(f_xBHT);



ref_function_supp=f_supp./f_xBHT;
ref_function_nonsupp=f_nonsuppBHT./f_xBHT;

%----------------------------------------------------------------------

% pdf for noise distribution
noise_pdf=pdf2('Normal',Xindex,0,sigma_N+sqrt(t^2/12/sigmaX^2));
noise_pdf=noise_pdf/sum(noise_pdf);

% get the information on the sensing matrix
degofV=sum(H)';degofC=sum(H')';
MaxdegV=max(degofV);MaxdegC=max(degofC);

% saturate the value of the noisy measurements
tempindex=find(abs(Z)> absmax_Z);
Z(tempindex)=sign(Z(tempindex))*absmax_Z;

% generate the index of the noisy measurement Z
quant_Z=round(Z/t);% This process quantization of the measurements

V=diag(ones(M,1)*sigma_N^2);
%%-----------------------------------------------------------------------%%
% Belief propagation in C code implementation
tstart=tic;
 [LVp,BPterminate_cond]=BPprocess_verC5(maxiter,MaxdegV,MaxdegC,zeropos,H, quant_Z,degofC,f_x,noise_pdf);
 telapsed=toc(tstart);
%%-----------------------------------------------------------------------%%


% MAP estimation for CS-BP
Xmin_BP=Xmin-mod(Xmin,t);
[v, est_Xindex]=max(LVp);
estX_CSBP=((est_Xindex-zeropos)').*t;
estX_CSBP(find(abs(estX_CSBP) < Xmin_BP))=0;
if Xmin_BP ==0
    est_suppBP=find(abs(estX_CSBP)> Xmin_BP)';
else
    est_suppBP=find(abs(estX_CSBP)>= Xmin_BP)';
end
Xstate_hatBP(est_suppBP)=1;

% the number of state error
state_err_CSBP=sum(bitxor(Xstate,Xstate_hatBP));

%--------------------- Bayesian support detection ------------------------%
for z=1:N
    Psupp(z)=sum(ref_function_supp.*LVp(:,z));
    Pnonsupp(z)=sum(ref_function_nonsupp.*LVp(:,z));
end
% likelifood function for the support detect
Lx=log(Pnonsupp./Psupp);
% likehood test for the support set estimation
est_supp=find(Lx<BHT_th);% the index set of the support 
Xstate_hatBHT(est_supp)=1;

% the number of state error
state_err_BHT=sum(bitxor(Xstate,Xstate_hatBHT));

% MMSE estimation to find value of the support
if length(est_supp)~=0
    estX_BHT(est_supp)=lscov2(H(:,est_supp),Z,V);
else

end
% Calculate normalized Mean square error 
energy_x=norm(trueX)^2;
Merror_CSBP=norm(estX_CSBP - trueX)^2/energy_x;
Merror_BHT=norm(estX_BHT - trueX)^2/energy_x;


[index,value]=find(estX_BHT);
Xstate_BHT=zeros(N,1);
Xstate_BHT(index)=1;
state_err_BHT=sum(bitxor(Xstate_BHT,Xstate));

end


