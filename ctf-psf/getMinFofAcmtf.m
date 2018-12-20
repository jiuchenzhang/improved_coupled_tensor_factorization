function acmtf = getMinFofAcmtf(Z,R,beta_cp,beta_pca,mi,options)
fprintf('%s\n','ACMTF begins to chose');

m_acmtf= 1000;
times = 10;

% initialization
acmtf =[];
acmtf_flag = zeros(1,times);
acmtf_funv = zeros(1,times);
acmtf_iter =zeros(1,times);
acmtf_fit = zeros(1,times);
acmtf_time  = zeros(1,times);

for i = 1:times
   tic;
   acmtf_res = c_acmtf_opt(Z,R,'alg_options',options,'beta_cp',beta_cp,'beta_pca',beta_pca);
   
   acmtf_time(i) = toc;
   acmtf_flag(i) = acmtf_res.ExitFlag;
   acmtf_funv(i) = acmtf_res.FunV; 
   acmtf_fit(i) = acmtf_res.Fit(mi);
   acmtf_iter(i)= acmtf_res.Iter;
end

for i = 1:times
    if(m_acmtf > acmtf_funv(i) && (acmtf_flag(i) ==0||acmtf_flag(i)==3))
        m_acmtf = acmtf_funv(i);
        min = i;
    end
end

% save the result
acmtf.ExitFlag = acmtf_flag(min);
acmtf.FunV = acmtf_funv(min);
acmtf.Iter = acmtf_iter(min);
acmtf.Time = acmtf_time(min);
acmtf.Fit = acmtf_fit(min);
fprintf('%s\n','ACMTF finish');


