function fitMode(M,path,Zt,R,options)
%% fit model
AMTC = ctf_opt(Zt,R,'alg_options',options);

for i = 1:length(M)
    if(Zt.m(i)~=0)
        options.num = length(Zt.object{i});
        CP{i} = c_cp_wopt(Zt.object{i},Zt.miss{i},R,Zt,i,'alg_options',options);
    end
end

Zt.size = [378 ,378,378];
Zt.modes = {[1,2],[1,3]};
MTF = c_cmtf_opt(Zt,R,options);
fprintf('%s\n','save the result');
saveRes(path,MTF,AMTC,CP);