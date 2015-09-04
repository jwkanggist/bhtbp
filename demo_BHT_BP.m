%-------------------------------------------------------------------------
% filename :demo_BHT_BP.m
% objective :Testbench to  test out the BHT_BP (sparse recovery algorithm via
% Bayesian support detection)
%
% This MATLAB corresponds to the paper:
% Jaewook Kang, Heung-No Lee, Kiseon Kim, 
% "Bayesian Hypothesis Test using Nonparametric Belief Propagation 
% for Noisy Sparse Recovery,"  
% to appear in IEEE Transaction on Signal processing.
%
% Written by: Jaewook Kang in GIST Korea
% Email: jwkkang@gist.ac.kr
% Created: July 2011, modified  9th Feb. 2012
%--------------------------------------------------------------------------
clc
clear all
close all
% put key subdirectories in path if not already there
% path(path, './Sensing_matrix_data');
path(path, './solver');
path(path, './etc_tool');
path(path, './solver/C_codes');

% setup system parameters
N=200; % Signal length of X
M=round(N*0.5);% the number of measurements Y
q=0.1;
K=round(q*N);
Nd=256;

sigmaX=5;% standard deviation for support of signal
Ts=2*3*sigmaX/Nd;

SNRdB=20;% SNR over measurement Y

Xmin=sigmaX/4;
maxiter=30;
%-----0------------------------------------------------------
% % load a good M by N LDPC-like measurement matrix 
L=4;
H=gen_uLDPC_H(N,M,L);

%----------------- Intro ------------------------------
disp('%-----------------------------------------------------------------------------------------%');
disp('% Sparse signal recovery via Bayesian support detection (BHT_BP) via belief propagation');
disp('% Copyright:INFONET-CSNL in GIST Korea.');
disp('%');
disp('% Written by Jaewook Kang, Phd student in GIST-DIC, jwkkang@gist.ac.kr');
disp('% July, 2011');
disp('%');
disp('%------------------------------------------------------------------------------------------%');
disp('<Experiment setup>');
disp(sprintf('System scale: N=%d, M=%d,',N,M));
disp(sprintf('Sparsity : q=%f,',q));
disp(sprintf('SNR=%d <dB>',SNRdB));
disp(sprintf('The number of iteration for BP process =%d ',maxiter));
disp('%------------------------------------------------------------------------------------------%');
% signal generation
disp('<The sparse signal generation>');
    
%[Xstate,X,supp]=Ksparse_signal_gen(N,K,sigmaX,Xmin);
[Xstate,X,supp]=sparse_signal_gen(N,q,sigmaX,Xmin);
 %[Xstate,X,supp]=signed_sparse_signal_gen(N,q,sigmaX);

% random projection (sensing) with binary sparse matrix H 
disp('<Measurement generation>');
Y=H*X;
energy_Y=norm(Y)^2; 

sigma_N=sqrt(energy_Y*10^(-SNRdB/10)/M);% calculate noise power accoding to the given SNR
noise=sigma_N*randn(M,1);
% add measurement noise
Z=Y+noise;

noiseE=norm(noise)^2;
measuredSNR=10*log10(energy_Y/noiseE);

% signal recoveyr
[estX_BHT,estX_CSBP,Merror_CSBP,Merror_BHT,state_err_BHT,state_err_CSBP,telapsed]=BHT_BP_ver3(X,H,Z,q,maxiter,sigma_N,sigmaX,Xstate,Nd,Xmin);



% BCS
initsigma2 = std(Z)^2/1e2;
[XsuppBCS,supportBCS,sigma2,errbars] = BCS_fast_rvm(H,Z,initsigma2,1e-8);
estX_BCS = zeros(N,1); 
estX_BCS(supportBCS) = XsuppBCS;
state_hat_BCS=zeros(N,1);
state_hat_BCS(supportBCS) = 1;
state_err_BCS=sum(bitxor(Xstate,state_hat_BCS));

Merror_BCS=norm(estX_BCS - X)^2/norm(X)^2;




disp('%------------------------------------------------------------------------------------------%');
disp('<Resulting output>')
disp(sprintf('CS-BP: Total Mean square error  = %8.7f',Merror_CSBP));
disp(sprintf('BHT_BP: Total Mean square error  = %8.7f',Merror_BHT));
disp(sprintf('BCS: Total Mean square error  = %8.7f',Merror_BCS));
disp(sprintf('      : The number of the state mismatch using BP = %d',state_err_CSBP));
disp(sprintf('      : The number of the state mismatch using BHT = %d',state_err_BHT));
disp(sprintf('      : The number of the state mismatch using BCS = %d',state_err_BCS));
disp(sprintf('Telapsed Time  = %8.7f sec',telapsed));
disp('%------------------------------------------------------------------------------------------%');


figure(10)
subplot(3,1,1); plot(1:N,X,1:N,estX_CSBP,'rO');axis([1 N -max(abs(X))-0.5 max(abs(X))+0.5]); title(['(a) Reconstruction with CS-BP']);legend('Ground Truth','CS-BP');
subplot(3,1,2);plot(1:N,X,1:N,estX_BHT,'rO'); axis([1 N -max(abs(X))-0.5 max(abs(X))+0.5]);title(['(c) Reconstruction with BHT-BP']); legend('Ground Truth','BHT-BP');
subplot(3,1,3);plot(1:N,X,1:N,estX_BCS,'rO'); axis([1 N -max(abs(X))-0.5 max(abs(X))+0.5]);title(['(c) Reconstruction with BCS']); legend('Ground Truth','BCS');box on;


