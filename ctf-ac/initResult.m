function [mtf,ctf,cp] = initResult(NumOfT,runTime)
mtf.Iter = zeros(runTime,1);
mtf.ExitFlag = zeros(runTime,1);
mtf.F = cell(runTime,1);
mtf.TCS = zeros(runTime,NumOfT);
mtf.RMSE = zeros(runTime,NumOfT);

ctf.Iter = zeros(runTime,1);
ctf.ExitFlag = zeros(runTime,1);
ctf.F = cell(runTime,1);
ctf.TCS = zeros(runTime,NumOfT);
ctf.RMSE = zeros(runTime,NumOfT);

cp.Iter = zeros(runTime,NumOfT);
cp.ExitFlag = zeros(runTime,NumOfT);
cp.F = cell(runTime,NumOfT);
cp.TCS = zeros(runTime,NumOfT);
cp.RMSE = zeros(runTime,NumOfT);





