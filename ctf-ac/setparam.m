function params = setparam(params)
params.addParamValue('Display','iter', @(x) ismember(x,{'iter','final','off'}));
params.addParamValue('DisplayIters', 1, @(x) x >= 1);
params.addParamValue('MaxIters', 1000, @(x) x >= 0);
params.addParamValue('RelFuncTol', 1e-6, @isnumeric);
params.addParamValue('num', 2, @(x) x >= 0);