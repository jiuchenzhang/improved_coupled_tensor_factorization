function saveRes_Rc(path,psfCTF_result)

psfCTF_result_Path = strcat(path,'/psfCTF_result.mat');
mkdir(path);
save(psfCTF_result_Path,'psfCTF_result');