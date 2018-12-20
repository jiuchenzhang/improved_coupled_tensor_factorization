clear;
clc;
addpath('../ctf-ac','cmtf_toolbox_v1_1','poblano_toolbox_v1_1','tensor_toolbox_2.5')
%% Parameters
%Number of tensors
NumOfT              = 2;
%Dimension
dim                 = 3;
%Number of shared factors
rCf                 = 1;
noise               = 0.1;
%Size of each dimension
Num                 = 50;
size(1:dim)              = Num;
flag_sparse(1:NumOfT)         = 0;

rRank               = 5;

%The sampling ratio of missing values
% M                   = [0.5 0.5];   %Average sampling
M                   = [0 0.95];      %Differential sampling

for n = 1:NumOfT
    lambdas{n}(1:rRank)   = 1;
end

%% initialization
options.Display ='final';
options.MaxFuncEvals = 100000;
options.MaxIters     = 10000;
options.RelFuncTol   = 1e-6;

%% generate synthetic data  and add noise
Z = createSysData('size', size, 'lambdas', lambdas, 'noise', noise, 'rRank', rRank, 'rCf', rCf);

%%  Initialization  of  Algorithm
%rank of the factorization
R                   = 10;

pCf                 = 1;

csize               = zeros ( NumOfT * dim , 1);
csize(1:end)        = repmat (size,1,NumOfT);

Z.csize             = csize;
Z.pCf               = pCf;
Z.R                 = R;
Z.N                 = NumOfT;

%% add missing value
Zt = addmv(Z, M, flag_sparse);
Zt.m = M;

Zt.size(1:(dim-1)*NumOfT+1) = Num; 
Zt.modes = {[1,2,3],[1,4,5]};

runTime = 20;
%% fit model
% initialization
[mtf,ctfac,cp] = initResult(NumOfT,runTime);

    path =strcat( 'result/synthetic-data/','size',num2str(size),'-rR',num2str(rRank),'-R',num2str(R),'-Miss',num2str(M));

    for t = 1: runTime
        ctfac_out = ctf_opt(Zt,R,'alg_options',options);
        ctfac.Iter(t) = ctfac_out.Iter;
        ctfac.ExitFlag(t) = ctfac_out.ExitFlag;
        ctfac.F{t} = ctfac_out.F;
        ctfac.TCS(t,:) = ctfac_out.TCS;
        ctfac.RMSE(t,:) = ctfac_out.RMSE;
        fprintf('%s\n','finish ctf-ac');

        for i = 1:length(M)
            cp_out = c_cp_wopt(Zt.object{i},Zt.miss{i},R,Zt,i,'alg_options',options);
            cp.Iter(t,i) = cp_out.Iter;
            cp.ExitFlag(t,i) = cp_out.ExitFlag;
            cp.F{t,i} = cp_out.F;
            if(M(i)~=0)
                cp.TCS(t,i) = cp_out.TCS;
            end
            cp.RMSE(t,i) = cp_out.RMSE;
        end
        fprintf('%s\n','finish cp');

        mtf_out = c_cmtf_opt(Zt,R,options);
        mtf.Iter(t) = mtf_out.Iter;
        mtf.ExitFlag(t) = mtf_out.ExitFlag;
        mtf.F{t} = mtf_out.F;
        mtf.TCS(t,:) = mtf_out.TCS;
        mtf.RMSE(t,:) = mtf_out.RMSE;
        fprintf('%s\n','finish mtf');
    end
    fprintf('%s\n','save the result')

