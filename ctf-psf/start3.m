%%
% An illustration of the joint decomposition methods with a third-order
% tensor and a matrix with R components. X and Y share one dimension. 
%%
clear;
clc;
addpath('../ctf-psf','cmtf_toolbox_v1_1','poblano_toolbox_v1_1','tensor_toolbox_2.5')
%% Parameters
size        = [50,30,20,100];
lambdas     = {[1,1,1,1],[1,1,1,1]};
modes       = {[1,2,3],[1,4]};
R           = 4;
flag_sparse = [0 0];                                                                               
noise       = 0.1;
M           = [0.9  0.4];
c =3;

%% Parameters for CTF-PSF
eachSize = [50,30,20,50,100];
structSize = {[1,2,3],[4,5]};
sameDim = {[1,4]};

%% initialization
options = ncg('defaults');
options.Display ='final';
options.MaxFuncEvals = 100000;
options.MaxIters     = 10000;
options.StopTol      = 1e-7;
options.RelFuncTol   = 1e-6;

%% set parameter of acmtf
beta_cp  = 0.001;
beta_pca = 0.001;

%% initializations
cp_result =[];
acmtf_result = [];
cmtf_result = [];
psfCTF_result = [];

for R  = 4:4
    path =strcat( '../result/synthetic-data-two/','size',num2str(size),'-Rr',num2str(length(lambdas{1})),'-R',num2str(R),'-Rc',num2str(c),...
        '-Noise',num2str(noise),'-Miss',num2str(M));
    for i = 1:10
        [Z,A] = createData('eachSize',eachSize,'structSize',structSize,'c',c,'sameDim',sameDim,'lambdas',lambdas,'noise',noise);
        Z.modes = modes;
        Z.size = size;
        Z.r = c;
        Z.m  = M;

        Zt = ssingVal(Z,M,flag_sparse);
        
         psfCTF_res = psfCTF_opt(Zt,R,0,'alg_options',options);
         acmtf_res = c_acmtf_opt(Zt,R,'alg_options',options,'beta_cp',beta_cp,'beta_pca',beta_pca);
         cmtf_res = c_cmtf_opt(Zt,R,'alg_options',options);
         
        if(~isempty(acmtf_res))
            acmtf_result.ExitFlags(i) = acmtf_res.ExitFlag;
            acmtf_result.Fits(:,i) = acmtf_res.Fit;
            acmtf_result.FxValues(i) = acmtf_res.FunV;
            acmtf_result.Iters(i) = acmtf_res.Iter;
        end

        if(~isempty(cmtf_res))
            cmtf_result.ExitFlags(i) = cmtf_res.ExitFlag;
            cmtf_result.Fits(:,i) = cmtf_res.Fit;
            cmtf_result.FxValues(i) = cmtf_res.FunV;
            cmtf_result.Iters(i) = cmtf_res.Iter;
        end

        if(~isempty(psfCTF_res))
            psfCTF_result.ExitFlags(i) = psfCTF_res.ExitFlag;
            psfCTF_result.Fits(:,i) = psfCTF_res.Fit;
            psfCTF_result.FxValues(i) = psfCTF_res.FunV;
            psfCTF_result.Iters(i) = psfCTF_res.Iter;
        end
    end

    %% save the results 
    fprintf('%s\n','save the result');
    
    psfCTF_result_Path = strcat(path,'/psfCTF_result.mat');
    acmtf_result_Path = strcat(path,'/acmtf_result.mat');
    cmtf_result_Path = strcat(path,'/cmtf_result.mat');

    mkdir(path);
    save(acmtf_result_Path,'acmtf_result');
    save(cmtf_result_Path,'cmtf_result');
    save(psfCTF_result_Path,'psfCTF_result');
end