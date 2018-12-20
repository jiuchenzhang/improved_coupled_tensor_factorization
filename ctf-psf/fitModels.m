function fitModels(Z,M,R,mi,beta_cp,beta_pca,shared,path,runTimes,flag_sparse,options)
%% initializations
acmtf_result = [];
cmtf_result = [];
psfCTF_result = [];

%% main loop, fit models and choose the best point
for i = 1:runTimes
    % add missing value and missing pattern
    Zt = addMissingVal(Z,M,flag_sparse);
    
    % choose the best starting poiont
    psfCTF = getMinFofpsfCTF(Zt,R,shared,mi,options);
    acmtf = getMinFofAcmtf(Zt,R,beta_cp,beta_pca,mi,options);
    cmtf = getMinFofCmtf(Zt,R,mi,options);
    
    cmtf_result.ExitFlags(i) = cmtf.ExitFlag;
    cmtf_result.Fits(i) = cmtf.Fit;
    cmtf_result.FxValues(i) = cmtf.FunV;
    cmtf_result.Iters(i) = cmtf.Iter;
    cmtf_result.Times(i) = cmtf.Time;

    acmtf_result.ExitFlags(i) = acmtf.ExitFlag;
    acmtf_result.Fits(i) = acmtf.Fit;
    acmtf_result.FxValues(i) = acmtf.FunV;
    acmtf_result.Iters(i) = acmtf.Iter;
    acmtf_result.Times(i) = acmtf.Time;

    psfCTF_result.ExitFlags(i) = psfCTF.ExitFlag;
    psfCTF_result.Fits(i) = psfCTF.Fit;
    psfCTF_result.FxValues(i) = psfCTF.FunV;
    psfCTF_result.Iters(i) = psfCTF.Iter;
    psfCTF_result.Times(i) = psfCTF.Time;
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