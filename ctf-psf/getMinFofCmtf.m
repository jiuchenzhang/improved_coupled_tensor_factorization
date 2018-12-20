function cmtf = getMinFofCmtf(Z,R,mi,options)
fprintf('%s\n','CMTF begins to chose');

m_cmtf = 1000;
times = 10;

% initialization
cmtf =[];
cmtf_flag = zeros(1,times);
cmtf_funv = zeros(1,times);
cmtf_iter =zeros(1,times);
cmtf_fit = zeros(1,times);
cmtf_time  = zeros(1,times);

for i = 1:times
   tic;
   cmtf_res = c_cmtf_opt(Z,R,'alg_options',options);
   
   cmtf_time(i) = toc;
   cmtf_flag(i) = cmtf_res.ExitFlag;
   cmtf_funv(i) = cmtf_res.FunV;
   cmtf_fit(i) = cmtf_res.Fit(mi);
   cmtf_iter(i)= cmtf_res.Iter;
end

for i = 1:times
    if(m_cmtf > cmtf_funv(i) && (cmtf_flag(i)==3||cmtf_flag(i)==0))
       m_cmtf = cmtf_funv(i);
       min =i;
    end
end
% save the result
cmtf.ExitFlag = cmtf_flag(min);
cmtf.FunV = cmtf_funv(min);
cmtf.Iter = cmtf_iter(min);
cmtf.Time = cmtf_time(min);
cmtf.Fit = cmtf_fit(min);
fprintf('%s\n','CMTF finish');