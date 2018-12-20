function saveRes(path,cp_result,acmtf_result,cmtf_result,psfCTF_result)
cp_result_Path = strcat(path,'/cp_result.mat');
psfCTF_result_Path = strcat(path,'/psfCTF_result.mat');
acmtf_result_Path = strcat(path,'/acmtf_result.mat');
cmtf_result_Path = strcat(path,'/cmtf_result.mat');
mkdir(path);
save(cp_result_Path,'cp_result');
save(acmtf_result_Path,'acmtf_result');
save(cmtf_result_Path,'cmtf_result');
save(psfCTF_result_Path,'psfCTF_result');