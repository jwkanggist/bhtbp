%-------------------------------------------------------------------------
% filename :gen_Qmatrix.m
% objective :To generate Q matrix for bipartite graph composition with H
%
%
% Written by: Jaewook Kang in GIST Korea
% Email: jwkkang@gist.ac.kr
% Created: July 2011
%--------------------------------------------------------------------------
function [Q1,Q2]=gen_Qmatrix(N,M,H,MaxdegV,MaxdegC)

Q1 = zeros(MaxdegV,N);
Q2 = zeros(MaxdegC,M);

for i=1:N
    temp=find(H(:,i));
    temp_length=length(temp);
    Q1(1:temp_length,i)=temp;
    
    if (i<=M)
        temp=find(H(i,:));
        temp_length=length(temp);
        Q2(1:temp_length,i)=temp';
    end
end
