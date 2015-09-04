%-------------------------------------------------------------------------
% filename :demo_ImageRev.m
% objective :Testbench to compare the performance of CS recovery via BCS, BHT-BP, and
% CS-BP
%
% This MATLAB corresponds to the paper:
% Jaewook Kang, Heung-No Lee, Kiseon Kim, 
% "Bayesian Hypothesis Test using Nonparametric Belief Propagation 
% for Noisy Sparse Recovery,"  
% to appear in IEEE Transaction on Signal processing.
%
% Written by: Jaewook Kang in GIST Korea
% Email: jwkkang@gist.ac.kr
% Created: July 2012
%--------------------------------------------------------------------------
clc
clear all
close all
% put key subdirectories in path if not already there
path(path, './solver');
path(path, './etc_tool');
path(path, './image');

%-----------------------------------------------
% load an image
% [Xmat,map] = imread('cameraman256.tif');
[Xmat,map] = imread('cameraman128.tif');

% X_temp=double(X_temp(:,:,1));
Xmat=double(Xmat(:,:,1));
%-----Discrete wavelet transform (Haar)------------%
[cA, cH, cV, cD] =dwt2(Xmat);
[cAA, cAH, cAV, cAD]=dwt2(cA);
[cAAA, cAAH, cAAV, cAAD]=dwt2(cAA);
cAA2=[cAAA cAAV;
      cAAH,cAAD];
cA2=[cAA2 cAV;
     cAH,cAD];
Xmat_wavelet=[cA2 cV;cH cD];
% figure(2);imagesc(Xmat_wavelet);
% rev_image=cast(Xmat_wavelet,'uint8');
% imwrite(rev_image,'./image\lena128_wavelet.tif','tif')
clear cA cH cV cD cAA cAH cAD cA2 cAA2
X=Xmat_wavelet(:);
%-----------------------------------------------------------%


%-----------setup system parameters----------------------------%
th_value=100; % threshold to belong to the signal support
q=length(find(abs(X)>th_value))/(128*128);
N=length(X); % Signal length of X
M=N*0.5;% the number of measurements Y

Nd=256;
supp=find(abs(X)>th_value);
nonsupp=find(abs(X)<=th_value);
X(nonsupp)=0;% K-term approximation
Xstate=zeros(N,1);
Xstate(supp)=1;
%[Xstate,X,supp]=quantize_signal(X,N,sigmaX,q_level);
%-----------------------------
%load sensing matrix pre-generated 
load H_8192_N16384.mat

SNRdB=20;% SNR over measurement Y
sigmaX=max(abs(X))/2;% standard deviation for support of signal
Xmin=min(abs(X(supp)));
maxiter=30;

disp('%-----------------------------------------------------------------------------------------%');
disp('% Sparse signal recovery via Bayesian support detection (BHT_BP) via belief propagation');
disp('% The cameramen image recovery testing');
disp('% Copyright:INFONET-CSNL in GIST Korea.');
disp('%');
disp('% Written by Jaewook Kang, Phd student in GIST-DIC, jwkkang@gist.ac.kr');
disp('% July, 2013');
disp('%');
disp('%------------------------------------------------------------------------------------------%');
disp('<Experiment setup>');
disp(sprintf('System scale: N=%d, M=%d,',N,M));
disp(sprintf('Sparsity : q=%f,',q));
disp(sprintf('SNR=%d <dB>',SNRdB));
disp(sprintf('The number of iteration for BP process =%d ',maxiter));
disp('%------------------------------------------------------------------------------------------%');

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
%--------------- CS-BP and BHT-BP recovery -------------------------------%
[estX_BHT,estX_CSBP,Merror_CSBP,Merror_BHT,state_err_BHT,state_err_CSBP,telapsed]=BHT_BP_ver3(X,H,Z,q,maxiter,sigma_N,sigmaX,Xstate,Nd,Xmin);
%-----------------BCS recovery ------------------------------------------%
initsigma2 = std(Z)^2/1e2;
[XsuppBCS,supportBCS,sigma2,errbars] = BCS_fast_rvm(H,Z,initsigma2,1e-8);
estX_BCS = zeros(N,1); 
estX_BCS(supportBCS) = XsuppBCS;
state_hat_BCS=zeros(N,1);
state_hat_BCS(supportBCS) = 1;
NumOfState_err_BCS=sum(bitxor(Xstate,state_hat_BCS));

Merror_BCS=norm(estX_BCS - X)^2/norm(X)^2;
clear H
disp('%------------------------------------------------------------------------------------------%');
disp('<Resulting output>')
disp(sprintf('CS-BP: Total Mean square error  = %8.7f',Merror_CSBP));
disp(sprintf('BHT_BP: Total Mean square error  = %8.7f',Merror_BHT));
disp(sprintf('BCS: Total Mean square error  = %8.7f',Merror_BCS));
disp(sprintf('      : The number of the state mismatch using BP = %d',state_err_CSBP));
disp(sprintf('      : The number of the state mismatch using BHT = %d',state_err_BHT));
disp(sprintf('      : The number of the state mismatch using BCS = %d',NumOfState_err_BCS));
disp(sprintf('Telapsed Time  = %8.7f sec',telapsed));
disp('%------------------------------------------------------------------------------------------%');

figure(10)
subplot(3,1,1); plot(1:N,X,1:N,estX_CSBP,'rO');axis([1 N -max(abs(X))-0.5 max(abs(X))+0.5]); title(['(a) Reconstruction with CS-BP']);legend('Ground Truth','CS-BP');
subplot(3,1,2);plot(1:N,X,1:N,estX_BHT,'rO'); axis([1 N -max(abs(X))-0.5 max(abs(X))+0.5]);title(['(b) Reconstruction with BHT_BP']); legend('Ground Truth','BHT_BP');
subplot(3,1,3);plot(1:N,X,1:N,estX_BCS,'rO'); axis([1 N -max(abs(X))-0.5 max(abs(X))+0.5]);title(['(c) Reconstruction with BCS']); legend('Ground Truth','BCS');box on;


%---------------- reconstruction to image format--------------------------%
rev_image_BP=zeros(sqrt(N));
rev_image_BHT=zeros(sqrt(N));
rev_image_BCS=zeros(sqrt(N));

cAAA_BP=zeros(sqrt(N)/8);cAAH_BP=zeros(sqrt(N)/8);cAAV_BP=zeros(sqrt(N)/8);cAAD_BP=zeros(sqrt(N)/8);
cAAA_BHT=zeros(sqrt(N)/8);cAAH_BHT=zeros(sqrt(N)/8);cAAV_BHT=zeros(sqrt(N)/8);cAAD_BHT=zeros(sqrt(N)/8);
cAAA_BCS=zeros(sqrt(N)/8);cAAH_BCS=zeros(sqrt(N)/8);cAAV_BCS=zeros(sqrt(N)/8);cAAD_BCS=zeros(sqrt(N)/8);

cAA_BP=zeros(sqrt(N)/4);cAH_BP=zeros(sqrt(N)/4);cAV_BP=zeros(sqrt(N)/4);cAD_BP=zeros(sqrt(N)/4);
cAA_BHT=zeros(sqrt(N)/4);cAH_BHT=zeros(sqrt(N)/4);cAV_BHT=zeros(sqrt(N)/4);cAD_BHT=zeros(sqrt(N)/4);
cAA_BCS=zeros(sqrt(N)/4);cAH_BCS=zeros(sqrt(N)/4);cAV_BCS=zeros(sqrt(N)/4);cAD_BCS=zeros(sqrt(N)/4);

cA_BP=zeros(sqrt(N)/2);cH_BP=zeros(sqrt(N)/2);cV_BP=zeros(sqrt(N)/2);cD_BP=zeros(sqrt(N)/2);
cA_BHT=zeros(sqrt(N)/2);cH_BHT=zeros(sqrt(N)/2);cV_BHT=zeros(sqrt(N)/2);cD_BHT=zeros(sqrt(N)/2);
cA_BCS=zeros(sqrt(N)/2);cH_BCS=zeros(sqrt(N)/2);cV_BCS=zeros(sqrt(N)/2);cD_BCS=zeros(sqrt(N)/2);

for i=1:sqrt(N)
    rev_image_BP(:,i)=round(estX_CSBP(sqrt(N)*(i-1)+1:sqrt(N)*i));
    rev_image_BHT(:,i)=round(estX_BHT(sqrt(N)*(i-1)+1:sqrt(N)*i));
    rev_image_BCS(:,i)=round(estX_BCS(sqrt(N)*(i-1)+1:sqrt(N)*i));
end

% inverse wavelet transform
cAAA_BP=rev_image_BP(1:sqrt(N)/8,1:sqrt(N)/8);           cAAV_BP=rev_image_BP(1:sqrt(N)/8,sqrt(N)/8+1:sqrt(N)/4);
cAAH_BP=rev_image_BP(sqrt(N)/8+1:sqrt(N)/4,1:sqrt(N)/8); cAAD_BP=rev_image_BP(sqrt(N)/8+1:sqrt(N)/4,sqrt(N)/8+1:sqrt(N)/4); 

cAAA_BHT=rev_image_BHT(1:sqrt(N)/8,1:sqrt(N)/8);           cAAV_BHT=rev_image_BHT(1:sqrt(N)/8,sqrt(N)/8+1:sqrt(N)/4);
cAAH_BHT=rev_image_BHT(sqrt(N)/8+1:sqrt(N)/4,1:sqrt(N)/8); cAAD_BHT=rev_image_BHT(sqrt(N)/8+1:sqrt(N)/4,sqrt(N)/8+1:sqrt(N)/4); 

cAAA_BCS=rev_image_BCS(1:sqrt(N)/8,1:sqrt(N)/8);           cAAV_BCS=rev_image_BCS(1:sqrt(N)/8,sqrt(N)/8+1:sqrt(N)/4);
cAAH_BCS=rev_image_BCS(sqrt(N)/8+1:sqrt(N)/4,1:sqrt(N)/8); cAAD_BCS=rev_image_BCS(sqrt(N)/8+1:sqrt(N)/4,sqrt(N)/8+1:sqrt(N)/4); 

cAA_BP=idwt2(cAAA_BP, cAAH_BP, cAAV_BP, cAAD_BP);
cAA_BHT=idwt2(cAAA_BHT, cAAH_BHT, cAAV_BHT, cAAD_BHT);
cAA_BCS=idwt2(cAAA_BCS, cAAH_BCS, cAAV_BCS, cAAD_BCS);

cAV_BP=rev_image_BP(1:sqrt(N)/4,sqrt(N)/4+1:sqrt(N)/2);
cAH_BP=rev_image_BP(sqrt(N)/4+1:sqrt(N)/2,1:sqrt(N)/4); 
cAD_BP=rev_image_BP(sqrt(N)/4+1:sqrt(N)/2,sqrt(N)/4+1:sqrt(N)/2); 

cAV_BHT=rev_image_BHT(1:sqrt(N)/4,sqrt(N)/4+1:sqrt(N)/2);
cAH_BHT=rev_image_BHT(sqrt(N)/4+1:sqrt(N)/2,1:sqrt(N)/4); 
cAD_BHT=rev_image_BHT(sqrt(N)/4+1:sqrt(N)/2,sqrt(N)/4+1:sqrt(N)/2); 

cAV_BCS=rev_image_BCS(1:sqrt(N)/4,sqrt(N)/4+1:sqrt(N)/2);
cAH_BCS=rev_image_BCS(sqrt(N)/4+1:sqrt(N)/2,1:sqrt(N)/4); 
cAD_BCS=rev_image_BCS(sqrt(N)/4+1:sqrt(N)/2,sqrt(N)/4+1:sqrt(N)/2); 

cA_BP=idwt2(cAA_BP, cAH_BP, cAV_BP, cAD_BP);
cA_BHT=idwt2(cAA_BHT, cAH_BHT, cAV_BHT, cAD_BHT);
cA_BCS=idwt2(cAA_BCS, cAH_BCS, cAV_BCS, cAD_BCS);

cV_BP=rev_image_BP(1:sqrt(N)/2,sqrt(N)/2+1:sqrt(N));
cH_BP=rev_image_BP(sqrt(N)/2+1:sqrt(N),1:sqrt(N)/2); 
cD_BP=rev_image_BP(sqrt(N)/2+1:sqrt(N),sqrt(N)/2+1:sqrt(N)); 

cV_BHT=rev_image_BHT(1:sqrt(N)/2,sqrt(N)/2+1:sqrt(N));
cH_BHT=rev_image_BHT(sqrt(N)/2+1:sqrt(N),1:sqrt(N)/2); 
cD_BHT=rev_image_BHT(sqrt(N)/2+1:sqrt(N),sqrt(N)/2+1:sqrt(N)); 

cV_BCS=rev_image_BCS(1:sqrt(N)/2,sqrt(N)/2+1:sqrt(N));
cH_BCS=rev_image_BCS(sqrt(N)/2+1:sqrt(N),1:sqrt(N)/2); 
cD_BCS=rev_image_BCS(sqrt(N)/2+1:sqrt(N),sqrt(N)/2+1:sqrt(N)); 

Xmat_hat_BP=idwt2(cA_BP, cH_BP, cV_BP, cD_BP);
Xmat_hat_BHT=idwt2(cA_BHT, cH_BHT, cV_BHT, cD_BHT);
Xmat_hat_BCS=idwt2(cA_BCS, cH_BCS, cV_BCS, cD_BCS);

Xmat_hat_BP=cast(Xmat_hat_BP,'uint8');
Xmat_hat_BHT=cast(Xmat_hat_BHT,'uint8');
Xmat_hat_BCS=cast(Xmat_hat_BCS,'uint8');

imwrite(Xmat_hat_BP,'./image\rev_cameraman_BP.tif','tif');
imwrite(Xmat_hat_BHT,'./image\rev_cameraman_BHT.tif','tif');
imwrite(Xmat_hat_BCS,'./image\rev_cameraman_BCS.tif','tif');
 
figure(1)
subplot(2,2,1);imagesc(Xmat); title(['(a) Origianl']);
subplot(2,2,2);imagesc(Xmat_hat_BHT); title(['(b) Reconstruction via BHTBP']);
subplot(2,2,3);imagesc(Xmat_hat_BP); title(['(c)  Reconstruction via CSBP']);
subplot(2,2,4);imagesc(Xmat_hat_BCS); title(['(d)  Reconstruction via BCS']);
box on;


 
