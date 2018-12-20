function cp = getMinFofCp(Z,R,mi,options)
fprintf('%s\n','CP begins to chose');

m_cp = 1000;
times = 10;

% initialization
cp =[];
cp_flag = zeros(1,times);
cp_funv = zeros(1,times);
cp_iter =zeros(1,times);
cp_fit = zeros(1,times);
cp_time  = zeros(1,times);

for i = 1:times
   tic;
   cp_res = c_cp_wopt(Z.object{mi},Z.miss{mi},R,Z.initData{mi},'alg_options',options);
   
   cp_time(i) = toc;
   cp_flag(i) = cp_res.ExitFlag;
   cp_funv(i) = cp_res.FunV;
   cp_fit(i) = cp_res.Fit;
   cp_iter(i)= cp_res.Iter;
end

for i = 1:times
    if(m_cp > cp_funv(i) && (cp_flag(i)==0 || cp_flag(i)==3))
        m_cp = cp_funv(i);
        min = i;
    end
end
% save the result 
cp.ExitFlag = cp_flag(min);
cp.FunV = cp_funv(min);
cp.Iter = cp_iter(min);
cp.Time = cp_time(min);
cp.Fit = cp_fit(min);
fprintf('%s\n','CP finish');


