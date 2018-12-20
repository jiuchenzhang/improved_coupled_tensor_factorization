function psfCTF = getMinFofpsfCTF(Z,R,shared,mi,options)
fprintf('%s\n','psfCTF begins to chose');

m_psfCTF = 1000;
times = 10;

% initialization
psfCTF =[];
psfCTF_flag = zeros(1,times);
psfCTF_funv = zeros(1,times);
psfCTF_iter =zeros(1,times);
psfCTF_fit = zeros(1,times);
psfCTF_time  = zeros(1,times);

for i = 1:times
   tic;
   psfCTF_res = psfCTF_opt(Z,R,shared,'alg_options',options);
   
   psfCTF_time(i)=toc;
   psfCTF_flag(i) = psfCTF_res.ExitFlag;
   psfCTF_funv(i) = psfCTF_res.FunV;
   psfCTF_fit(i) = psfCTF_res.Fit(mi);
   psfCTF_iter(i) = psfCTF_res.Iter;
end

% save the result
for i = 1:times
    if(m_psfCTF > psfCTF_funv(i)&&(psfCTF_flag(i)==0 || psfCTF_flag(i)==3))
       m_psfCTF = psfCTF_funv(i);
       min =i;
    end
end

psfCTF.ExitFlag = psfCTF_flag(min);
psfCTF.FunV = psfCTF_funv(min);
psfCTF.Iter = psfCTF_iter(min);
psfCTF.Time = psfCTF_time(min);
psfCTF.Fit = psfCTF_fit(min);
fprintf('%s\n','psfCTF finish');




