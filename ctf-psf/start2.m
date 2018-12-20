%%
% The TCS for CTF-PSF at different missing value rations (SR% where Rc={1,3,5},\hat{R}=R=6, noise=0.1)
%%
clear;
clc;
addpath('../ctf-psf','cmtf_toolbox_v1_1','poblano_toolbox_v1_1','tensor_toolbox_2.5')
%% Parameters

size        = [50,30,20,100];
lambdas     = {[1,1,1,1,1,1],[1,1,1,1,1,1]};
modes       = {[1,2,3],[1,4]};
R           = 6;
flag_sparse = [0 0];
noise       = 0.1;
M           = [0 0];

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

mp = [0.1,0.3,0.5,0.7,0.8,0.85];
[a,b] = ndgrid([1,3,5],mp);
parameter  =[a(:),b(:)];
parfor p = 1:length(parameter)
        M = [0 0];
        M(1)  = parameter(p,2);
        c = parameter(p,1);
        
        if R==c
            shared = 1;
        else
            shared = 0;
        end 
        
        path =strcat( '../result/synthetic-data-1-9-dif_Rc/','size',num2str(size),'-Rr',num2str(length(lambdas{1})),'-R',num2str(R),'-Rc',num2str(c),...
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

            Zt = ssingVal(Z,M,flag_sparse);

            psfCTF_res=jointGet_dif_Rc(Zt,R,shared,options);
 
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
        saveRes_Rc(path,psfCTF_result);
end