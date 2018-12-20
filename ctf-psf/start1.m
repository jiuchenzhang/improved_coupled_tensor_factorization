%%
% The TCS for each baseline with different amounts of missing data (estimate rank = 4 when Rc=3,Rs=1)
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
M           = [0  0];
c = 3;

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

mp = [0.3,0.4,0.5,0.6,0.7,0.8,0.85,0.9,0.93,0.95];
[a,b] = ndgrid([4],mp);
parameter  =[a(:),b(:)];
parfor p = 1:length(parameter)
        M = [0 0];
        M(1)  = parameter(p,2);
        R = parameter(p,1);
        
        if R==c
            shared = 1;
        else
            shared = 0;
        end 
        
        path =strcat( '../result/synthetic-data-1-9/','size',num2str(size),'-Rr',num2str(length(lambdas{1})),'-R',num2str(R),'-Rc',num2str(c),...
            '-Noise',num2str(noise),'-Miss',num2str(M));
        %% initializations
        cp_result =[];
        acmtf_result = [];
        cmtf_result = [];
        psfCTF_result = [];
        for i = 1:10
            [Z,A] = createData('eachSize',eachSize,'structSize',structSize,'c',c,'sameDim',sameDim,'lambdas',lambdas,'noise',noise);
            Z.modes = modes;
            Z.size = size;
            Z.r = c;
            Z.m  = M;

            Zt = addMissingVal(Z,M,flag_sparse);
             % factorization 
            [psfCTF_res,acmtf_res,cmtf_res,cp_res]=jointGet8(Zt,R,beta_cp,beta_pca,shared,options);

            if(~isempty(acmtf_res))
                acmtf_result.ExitFlags(i) = acmtf_res.ExitFlag;
                acmtf_result.Fits(:,i) = acmtf_res.Fit;
                acmtf_result.FxValues(i) = acmtf_res.FunV;
                acmtf_result.Iters(i) = acmtf_res.Iter;
                acmtf_result.Times(i) = acmtf_res.Time;
            end
 
            if(~isempty(cmtf_res))
                cmtf_result.ExitFlags(i) = cmtf_res.ExitFlag;
                cmtf_result.Fits(:,i) = cmtf_res.Fit;
                cmtf_result.FxValues(i) = cmtf_res.FunV;
                cmtf_result.Iters(i) = cmtf_res.Iter;
                cmtf_result.Times(i) = cmtf_res.Time;
            end

            if(~isempty(cp_res))
                cp_result.ExitFlags(:,i) = cp_res.ExitFlag;
                cp_result.Fits(:,i) = cp_res.Fit;
                cp_result.FxValues(:,i) = cp_res.FunV;
                cp_result.Iters(:,i) = cp_res.Iter;
                cp_result.Times(:,i) = cp_res.Time;
            end

            if(~isempty(psfCTF_res))
                psfCTF_result.ExitFlags(i) = psfCTF_res.ExitFlag;
                psfCTF_result.Fits(:,i) = psfCTF_res.Fit;
                psfCTF_result.FxValues(i) = psfCTF_res.FunV;
                psfCTF_result.Iters(i) = psfCTF_res.Iter;
                psfCTF_result.Times(i) = psfCTF_res.Time;
            end
        end

        %% save the results 
        fprintf('%s\n','save the result');
        saveRes(path,cp_result,acmtf_result,cmtf_result,psfCTF_result);
end