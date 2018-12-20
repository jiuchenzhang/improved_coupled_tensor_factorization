function [psfCTF,acmtf,cmtf,cp]=jointGet8(Z,R,beta_cp,beta_pca,shared,options)
runTimes = 10;
% initialization
psfCTF =[];
psfCTF_flag = zeros(1,runTimes);
psfCTF_funv = zeros(1,runTimes);
psfCTF_iter =zeros(1,runTimes);
psfCTF_fit = zeros(1,runTimes);
psfCTF_time  = zeros(1,runTimes);

% initialization
acmtf =[];
acmtf_flag = zeros(1,runTimes);
acmtf_funv = zeros(1,runTimes);
acmtf_iter =zeros(1,runTimes);
acmtf_fit = zeros(1,runTimes);
acmtf_time  = zeros(1,runTimes);

% initialization
cmtf =[];
cmtf_flag = zeros(1,runTimes);
cmtf_funv = zeros(1,runTimes);
cmtf_iter =zeros(1,runTimes);
cmtf_fit = zeros(1,runTimes);
cmtf_time  = zeros(1,runTimes);

% initialization
cp =[];
cp_flag = zeros(1,runTimes);
cp_funv = zeros(1,runTimes);
cp_iter =zeros(1,runTimes);
cp_fit = zeros(1,runTimes);
cp_time  = zeros(1,runTimes);


%% main loop, fit models and choose the best point
parfor i = 1:runTimes

    tic;
    psfCTF_res = psfCTF_opt(Z,R,shared,'alg_options',options);
    psfCTF_time(i)=toc;
    psfCTF_flag(i) = psfCTF_res.ExitFlag;
    psfCTF_funv(i) = psfCTF_res.FunV;
    psfCTF_iter(i) = psfCTF_res.Iter;

   tic;
   acmtf_res = c_acmtf_opt(Z,R,'alg_options',options,'beta_cp',beta_cp,'beta_pca',beta_pca);
   acmtf_time(i) = toc;
   acmtf_flag(i) = acmtf_res.ExitFlag;
   acmtf_funv(i) = acmtf_res.FunV; 
   acmtf_iter(i)= acmtf_res.Iter;

   tic;
   cmtf_res = c_cmtf_opt(Z,R,'alg_options',options);
   cmtf_time(i) = toc;
   cmtf_flag(i) = cmtf_res.ExitFlag;
   cmtf_funv(i) = cmtf_res.FunV;
   cmtf_iter(i)= cmtf_res.Iter;

    acmtf_fit(i) = acmtf_res.Fit(1);
    cmtf_fit(i) = cmtf_res.Fit(1);
    psfCTF_fit(i) = psfCTF_res.Fit(1);

   cp_res = c_cp_wopt(Z.object{1},Z.miss{1},R,Z.initData{1},'alg_options',options);
   cp_time(1,i) = toc;
   cp_flag(1,i) = cp_res.ExitFlag;
   cp_funv(1,i) = cp_res.FunV;
   cp_fit(1,i) = cp_res.Fit;
   cp_iter(1,i)= cp_res.Iter;

end
% save the result
[f,index] = min(psfCTF_funv);
if((psfCTF_flag(index)==0 || psfCTF_flag(index)==3))
    psfCTF.ExitFlag = psfCTF_flag(index);
    psfCTF.FunV = psfCTF_funv(index);
    psfCTF.Iter = psfCTF_iter(index);
    psfCTF.Time = psfCTF_time(index);
    psfCTF.Fit = psfCTF_fit(:,index);
end

[f,index] = min(acmtf_funv);
if((acmtf_flag(index) ==0||acmtf_flag(index)==3))
    acmtf.ExitFlag = acmtf_flag(index);
    acmtf.FunV = acmtf_funv(index);
    acmtf.Iter = acmtf_iter(index);
    acmtf.Time = acmtf_time(index);
    acmtf.Fit = acmtf_fit(:,index);
end

[f,index] = min(cmtf_funv);
if(cmtf_flag(index)==3||cmtf_flag(index)==0)
    cmtf.ExitFlag = cmtf_flag(index);
    cmtf.FunV = cmtf_funv(index);
    cmtf.Iter = cmtf_iter(index);
    cmtf.Time = cmtf_time(index);
    cmtf.Fit = cmtf_fit(:,index);
end

[f,index] = min(cp_funv);
if(cp_flag(index)==0 || cp_flag(index)==3)
    cp.ExitFlag = cp_flag(:,index);
    cp.FunV = cp_funv(:,index);
    cp.Iter = cp_iter(:,index);
    cp.Time = cp_time(:,index);
    cp.Fit = cp_fit(:,index);
end