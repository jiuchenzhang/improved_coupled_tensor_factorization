function output = ctf_opt(Z,R,varargin)
%% Error checking
%cmtf_check(Z);
% 
% if (nargin < 2)
%     error('Error: invalid input arguments');
% end

%% Set parameters
params = inputParser;
params.addParamValue('alg' , 'ctf_ncg', @(x) ismember(x,{'ctf_ncg','tn','lbfgs','sgd'}));
params.addParamValue('init', 'random', @(x) (iscell(x) || ismember(x,{'random','nvecs'})));
params.addOptional('alg_options', '', @isstruct);
params.parse(varargin{:});

%% Set up optimization algorithm
switch (params.Results.alg)
    case 'sgd'
        fhandle = @sgdUpdate;
    case 'ctf_ncg'
        fhandle = @ctf_ncg;
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
sz      = Z.csize;
N       = length(Z.csize)/Z.N;
P      = length(Z.object);

if iscell(params.Results.init)
    G = params.Results.init;
elseif strcmpi(params.Results.init,'random')
    G = cell(P,1);
    for t = 1:P
        G{t} = cell(N,1);
        for n = 1:N
            G{t}{n} = randn(sz((t-1)*N+n),R);
            for j=1:R                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
                G{t}{n}(:,j) = G{t}{n}(:,j) / norm(G{t}{n}(:,j));
            end
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
Znormsqr = cell(P,1);
for p = 1:P
    if isa(Z.object{p},'tensor') || isa(Z.object{p},'sptensor')
        Znormsqr{p} = norm(Z.object{p})^2;
    else
        Znormsqr{p} = norm(Z.object{p},'fro')^2;
    end
end

%% compute factors and model fit
out = feval(fhandle, @(x)ctf_fun(x,Z,Znormsqr), ctf_fac_to_vec(G), options);
Ai = ctf_vec_to_fac(out.X, Z);
for p  = 1:P
    F = ktensor(Ai{p});
    A{p} = F.U;
end
output.Init = G;
output.F  = A;
output.ExitFlag  = out.ExitFlag;
output.FunV = out.F;
output.Iter = out.Iters;


TCS = zeros(1,P);
RMSE = zeros(1,P);

for s = 1:P
    D = full(ktensor(A{s}));
   if(Z.m(s)~=0)
        trueval = Z.initData{s}(find(Z.miss{s}==0));
        estvalCMTF =  D(find(Z.miss{s}== 0));
        TCS(1,s) = norm(trueval-estvalCMTF)/norm(trueval);
        output.TCS = TCS;
        RMSE(1,s) = norm(Z.miss{s}.*(D-Z.initData{s}))/norm(Z.miss{s});
   else
       RMSE(1,s) = norm(Z.miss{s}.*(D-Z.initData{s}))/norm(Z.miss{s});
   end
end
output.RMSE = RMSE;
