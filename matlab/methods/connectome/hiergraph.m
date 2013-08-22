function  [CIJ,HM,K] = hiergraph(N, E, mx_lvl, nc_lvl, model)
% function  [CIJ,K] = hiergraph(N, E, mx_lvl, nc_lvl, model)
% inputs:
%           N        network size
%           E        desired number of edges
%           mx_lvl   number of hierarchical levels, N = 2^mx_lvl
%           nc_lvl   number of clusters at each level
%           model    edge density model 
% outputs:
%           CIJ      connection matrix
%           K        number of connections present in the output CIJ
%
%
% Author: Marcus Kaiser     Date: 30 July 2009

% initialization
k = nc_lvl;
CIJ = zeros(N);
HM = zeros(N);
E = round(E);

if k^(mx_lvl-1)>N, error('Not possible!'); end

for i=1:mx_lvl
    switch model
        case 1,
            E_i = E/mx_lvl;
        case 2,
            E_i = E * (k-1)/k^i;
        case 3,
            E_i = E / (i*(i+1));
        case 4,
            E_i = E * 2 / ((2*i-1) * (2*i+1));
        case 5,
            E_i = E * 4 / (i * (i+1) * (i+2));        
        case 6, % control I vary E_i (increase with level)
            ff = 3/2;
            fk = (ff^(mx_lvl+1) - ff) / (ff - 1) / mx_lvl; 
            f = (ff^i) / fk;
            E_i = f * E/mx_lvl;
        case 7, % control I vary E_i (decrease with level)
            ff = 2/3;
            fk = (ff^(mx_lvl+1) - ff) / (ff - 1) / mx_lvl; 
            f = (ff^i) / fk;
            E_i = f * E/mx_lvl;           
        case 8, % control II vary k (increase with level)
            ff = 1.1; % was 3/2
            %fk = (ff^(mx_lvl+1) - ff) / (ff - 1) / mx_lvl; 
            k = ceil( nc_lvl * (ff^i) / 1 );
            E_i = E/mx_lvl;            
        case 9, % control II vary k (decrease with level)
            ff = 0.9; % was 2/3
            %fk = (ff^(mx_lvl+1) - ff) / (ff - 1) / mx_lvl; 
            k = ceil( nc_lvl * (ff^i) / 1 );
            E_i = E/mx_lvl;            
    end; % switch
    A_i = (k-1)/k^i;
    p_i = E_i / (N^2 * A_i); % edge density within cluster
    N_i = floor(N * (1/k)^(i-1));   % size of current cluster
    N_c = floor(N/N_i);

    for j=1:N_c
        if j==N_c, n_i = N -(N_c-1)*N_i;
        else n_i = N_i;
        end
        c_i = (rand(n_i)<p_i);
        r0 = 1+(j-1)*N_i;
        r1 = r0 + n_i - 1;
        while length(unique(unique(HM(r0:r1,r0:r1))))~=1 && r1<N, r0=r0+1;r1=r1+1; end
        CIJ(r0:r1, r0:r1) = c_i;
        HM(r0:r1, r0:r1) = i;
    end;
end; % for i

for i=1:N, CIJ(i,i)=0; end

missE = E-nnz(CIJ);

while missE ~= 0
    r = 1+floor(rand(1,2)*N);
    if (CIJ(r(1),r(2))==1) && (missE<0)
        CIJ(r(1),r(2)) = 0;
        missE = missE + 1;
    elseif (CIJ(r(1),r(2))==0) && (missE>0)
        CIJ(r(1),r(2)) = 1;
        missE = missE - 1;
    end;
end;

K = nnz(CIJ);
