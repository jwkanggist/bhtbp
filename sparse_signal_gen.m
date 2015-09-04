%-------------------------------------------------------------------------
% filename :sparse_signal_gen.m
% objective :To  generate stochastically sparse signals with N, q, and
% sigmaX
%
% Input : 
%         N - signal length
%         q - the sparsity rate
%         sigmaX - standard deviation for supporting values in the signal
%          Xmin  - the minimum value of X
% Output : Xstate : information on the sparse pattern on the signal
%          X      : the generated signals
%
% Written by: Jaewook Kang in GIST Korea
% Email: jwkkang@gist.ac.kr
% Created: July 2011
%--------------------------------------------------------------------------
function [Xstate,X,supp]=sparse_signal_gen(N,q,sigmaX,Xmin)

% Numofsupp=binornd(N,q,1,NumofX);
Numofsupp=length(find(rand(N,1)<q));

if Numofsupp==0
    Numofsupp=round(q*N);
end

Xstate=randerr2(1,N,Numofsupp)';
supp=find(Xstate==1);
Xsupp=sigmaX*randn(Numofsupp,1);


% saturate the value of the signal
Maxvalue_X=sigmaX*2.9;
tempindex=find(abs(Xsupp)> Maxvalue_X);
Xsupp(tempindex)=sign(Xsupp(tempindex))*Maxvalue_X;

% remove too small value
tempindex=find(abs(Xsupp) < Xmin); 
Xsupp(tempindex)=sign(Xsupp(tempindex))*Xmin;

X=zeros(N,1);
X(supp)=Xsupp;
