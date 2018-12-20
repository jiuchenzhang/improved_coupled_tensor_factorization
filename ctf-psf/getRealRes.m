function getRealRes(Zt,R,M,beta_cp,beta_pca,shared,flag_sparse,path,options)
for i = 1:10
        % add missing values
        Zt = addMissingVal(Zt,M,flag_sparse);

        psfCTF_res = psfCTF_opt(Zt,R,shared,'alg_options',options);
        psfCTF_result.ExitFlags(i) = psfCTF_res.ExitFlag;
        psfCTF_result.Fits(i) = psfCTF_res.Fit(1);
        psfCTF_result.FxValues(i) = psfCTF_res.FunV;
        psfCTF_result.Iters(i) = psfCTF_res.Iter;

        cmtf_res = c_cmtf_opt(Zt,R,'alg_options',options);
        cmtf_result.ExitFlags(i) = cmtf_res.ExitFlag;
        cmtf_result.Fits(i) = cmtf_res.Fit(1);
        cmtf_result.FxValues(i) = cmtf_res.FunV;
        cmtf_result.Iters(i) = cmtf_res.Iter;

        acmtf_res = c_acmtf_opt(Zt,R,'alg_options',options,'beta_cp',beta_cp,'beta_pca',beta_pca);
        acmtf_result.ExitFlags(i) = acmtf_res.ExitFlag;
        acmtf_result.Fits(i) = acmtf_res.Fit(1);
        acmtf_result.FxValues(i) = acmtf_res.FunV;
        acmtf_result.Iters(i) = acmtf_res.Iter;

end
    

    psfCTF_result_Path = strcat(path,'/psfCTF_result.mat');
    acmtf_result_Path = strcat(path,'/acmtf_result.mat');
    cmtf_result_Path = strcat(path,'/cmtf_result.mat');
    mkdir(path);

    save(acmtf_result_Path,'acmtf_result');
    save(cmtf_result_Path,'cmtf_result');
    save(psfCTF_result_Path,'psfCTF_result');