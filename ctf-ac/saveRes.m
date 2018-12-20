function saveRes(path,mtf,ctf,cp)
mtf_Path = strcat(path,'/mtf.mat');
ctf_Path = strcat(path,'/ctf.mat');
cp_Path = strcat(path,'/cp.mat');
mkdir(path);
save(mtf_Path,'mtf');
save(ctf_Path,'ctf');
save(cp_Path,'cp');