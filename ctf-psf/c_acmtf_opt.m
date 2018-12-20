function output = c_acmtf_opt(Z,R,varargin)
%% Error checking
cmtf_check(Z);

if (nargin < 2)
    error('Error: invalid input arguments');
end

%% Set parameters
params = inputParser;
params.addParamValue('alg', 'ncg', @(x) ismember(x,{'ncg','tn','lbfgs'}));
params.addParamValue('beta_cp', 0.001, @(x) x >= 0);
params.addParamValue('beta_pca',0.001, @(x) x >= 0);
params.addParamValue('init', 'random', @(x) (isstruct(x) || ismember(x,{'random','nvecs'})));
params.addOptional('alg_options', '', @isstruct);
params.parse(varargin{:});
P = numel(Z.object);

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

if isstruct(params.Results.init)
    G.fac   = params.Results.init.fac;
    G.norms = params.Results.init.norms;
elseif strcmpi(params.Results.init,'random')
    G.fac = cell(N,1);
    for n=1:N
        G.fac{n} = randn(sz(n),R);
        for j=1:R
            G.fac{n}(:,j) = G.fac{n}(:,j) / norm(G.fac{n}(:,j));
        end
    end
    for p=1:P
        G.norms{p} =ones(R,1);
    end
elseif strcmpi(params.Results.init,'nvecs')
    G.fac = cell(N,1);
    for n=1:N
        G.fac{n} = cmtf_nvecs(Z,n,R);
    end
    for p=1:P
        G.norms{p} =ones(R,1);
    end             
else
    error('Initialization type not supported')
end

%% Fit ACMTF using Optimization
Znormsqr = cell(P,1);
for p = 1:P
    if isa(Z.object{p},'tensor') || isa(Z.object{p},'sptensor')
        Znormsqr{p} = norm(Z.object{p})^2;
    else
        Znormsqr{p} = norm(Z.object{p},'fro')^2;
    end
end
out = feval(fhandle, @(x)acmtf_fun(x,Z,R,Znormsqr, params.Results.beta_cp, params.Results.beta_pca),  acmtf_struct_to_vec(G), options);

%% Compute factors 
output.ExitFlag  = out.ExitFlag;
output.FunV = out.F;
output.Iter = out.Iters;

Temp = acmtf_vec_to_struct(out.X, Z, R);
Zhat = cell(P,1);
M = Z.m;
TCS = zeros(1,length(M));

for p=1:P
    Zhat{p} = ktensor(Temp.norms{p},Temp.fac(Z.modes{p}));
    if M(p)~=0
        D = full(Zhat{p});
        
        trueval = Z.initData{p}(find(Z.miss{p}==0));
        estvalACMTF =  D(find(Z.miss{p}==0));
        TCS(1,p) = norm(trueval-estvalACMTF)/norm(trueval);
        
    end
end
output.F = Temp;
output.Init = G;
output.Fit = TCS;


