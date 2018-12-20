function output = c_cmtf_opt(Z,R,varargin)
%% Error checking
cmtf_check(Z);

if (nargin < 2)
    error('Error: invalid input arguments');
end

%% Set parameters
params = inputParser;
params.addParamValue('alg' , 'ncg', @(x) ismember(x,{'ncg','tn','lbfgs'}));
params.addParamValue('init', 'random', @(x) (iscell(x) || ismember(x,{'random','nvecs'})));
params.addOptional('alg_options', '', @isstruct);
params.parse(varargin{:});

%% Set up optimization algorithm
switch (params.Results.alg)
    case 'ncg'
        fhandle = @ncg;
    case 'tn'
        fhandle = @tn;
    case 'lbfgs'
        fhandle = @lbfgs;
end

%% Set up optimization algorithm options
if isempty(params.Results.alg_options)
    options = feval(fhandle, 'defaults');
else
    options = params.Results.alg_options;
end       

%% Initialization
sz = Z.size;
N = length(sz);

if iscell(params.Results.init)
    G = params.Results.init;
elseif strcmpi(params.Results.init,'random')
    G = cell(N,1);
    for n=1:N
        G{n} = randn(sz(n),R);
        for j=1:R
            G{n}(:,j) = G{n}(:,j) / norm(G{n}(:,j));
        end
    end
elseif strcmpi(params.Results.init,'nvecs')
    G = cell(N,1);
    for n=1:N
        G{n} = cmtf_nvecs(Z,n,R);
    end
else
    error('Initialization type not supported')
end

%% Fit CMTF using Optimization
P = numel(Z.object);
Znormsqr = cell(P,1);
for p = 1:P
    if isa(Z.object{p},'tensor') || isa(Z.object{p},'sptensor')
        Znormsqr{p} = norm(Z.object{p})^2;
    else
        Znormsqr{p} = norm(Z.object{p},'fro')^2;
    end
end

out = feval(fhandle, @(x)cmtf_fun(x,Z,Znormsqr), tt_fac_to_vec(G), options);

%% Compute factors and model fit
P = ktensor(cmtf_vec_to_fac(out.X, Z));
F = P.U;
output.Init = G;
output.F  = F;
output.ExitFlag  = out.ExitFlag;
output.FunV = out.F;
output.Iter = out.Iters;
M = Z.m;

TCS =zeros(1,length(M));

for s = 1:length(M)
    if M(s)~=0
        D = full(ktensor(F(Z.modes{s})));

        trueval = Z.initData{s}(find(Z.miss{s}==0));
        estvalCMTF =  D(find(Z.miss{s}==0));
        TCS(1,s) = norm(trueval-estvalCMTF)/norm(trueval);
    end
end
output.Fit = TCS;

