%%
% The TCSs for coupled factorization methods at noise = 0.1 and 0.3,
% estimate rank = 4, Rc={1,2,3,4}, SR(%)={20,40,60,80}
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

%% Parameters for CTF-PSF
eachSize = [50,30,20,50,100];
structSize = {[1,2,3],[4,5]};
sameDim = {[1,4]};

%% CMTF-initialization
options = ncg('defaults');
options.Display ='final';
options.MaxFuncEvals = 100000;
options.MaxIters     = 10000;
options.StopTol      = 1e-7;
options.RelFuncTol   = 1e-6;

%% set parameter of acmtf
beta_cp  = 0.001;
beta_pca = 0.001;

runTimes = 10;
%% main loop
[a,b,c]= ndgrid(1:2:3,1:R,1:4);
parameters = [a(:),b(:),c(:)];

parfor i = 1:length(parameters)
    noise = parameters(i,1)*0.1;

    %shared components
    c = parameters(i,2);
    
    m = parameters(i,3)*0.2;
    
   % create coupled data
    [Z,A] = createData('eachSize',eachSize,'structSize',structSize,'c',c,'sameDim',sameDim,'lambdas',lambdas,'noise',noise);
    Z.modes = modes;
    Z.size = size;
    Z.r = c;
    
    M = zeros(1,2);
    M(1) = m;
    
    if c==R
        shared = 1;
    else
        shared = 0;
    end
    
    msize = length(M);
    mi = 0;
    for k = 1:msize
        if M(k)~=0
            mi = k;
        end
    end
    
    % save path
    path =strcat( '../result/synthetic-data/','size',num2str(size),'-Rr',num2str(length(lambdas{1})),'-R',num2str(R),'-Rc',num2str(c),...
        '-Noise',num2str(noise),'-Miss',num2str(M(1)));
    
    % fit models
    fitModels(Z,M,R,mi,beta_cp,beta_pca,shared,path,runTimes,flag_sparse,options);
end