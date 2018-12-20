function psfCTF=jointGet_dif_Rc(Z,R,shared,options)
runTimes = 10;
% initialization
psfCTF =[];
psfCTF_flag = zeros(1,runTimes);
psfCTF_funv = zeros(1,runTimes);
psfCTF_iter =zeros(1,runTimes);
psfCTF_fit = zeros(1,runTimes);
psfCTF_time  = zeros(1,runTimes);

%% main loop, fit models and choose the best point
for i = 1:runTimes

    tic;
    psfCTF_res = psfCTF_opt(Z,R,shared,'alg_options',options);
    psfCTF_time(i)=toc;
    psfCTF_flag(i) = psfCTF_res.ExitFlag;
    psfCTF_funv(i) = psfCTF_res.FunV;
    psfCTF_iter(i) = psfCTF_res.Iter;
    psfCTF_fit(i) = psfCTF_res.Fit(1);
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
