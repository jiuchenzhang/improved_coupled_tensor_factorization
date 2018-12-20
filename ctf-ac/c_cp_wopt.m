function output = c_cp_wopt(Zcp,Wcp,R,Zt,i,varargin)
%% Check for POBLANO
if ~exist('poblano_params','file')
    error(['CP_WOPT requires Poblano Toolbox for Matlab. This can be ' ...
           'downloaded at http://software.sandia.gov/trac/poblano.']);
end

%% Set parameters
params = inputParser;
params.addParamValue('alg','ncg', @(x) ismember(x,{'ncg','tn','lbfgs'}));
params.addParamValue('init','random', @(x) (iscell(x) || ismember(x,{'random','nvecs'})));
params.addParamValue('fun','auto', @(x) ismember(x,{'auto','default','sparse','sparse_lowmem'}));
params.addParamValue('alg_options', '', @isstruct);
params.parse(varargin{:});

%% Set up optimization algorithm
switch (params.Results.alg)
    case 'ncg'
        opthandle = @ncg;
    case 'tn'
        opthandle = @tn;
    case 'lbfgs'
        opthandle = @lbfgs;
end

%% Set up optimization algorithm options
if isempty(params.Results.alg_options)
    options = feval(opthandle, 'defaults');
else
    options = params.Results.alg_options;
end

%% Set up function handle
normZsqr = norm(Zcp)^2;
funtype = params.Results.fun;

if (isequal(funtype,'auto') && isa(Zcp,'tensor')) || isequal(funtype,'default')
    funhandle = @(x) tt_cp_wfun(Zcp,Wcp,x,normZsqr);
else
    if ~isa(Zcp,'sptensor') || ~isa(Wcp,'sptensor')
        warning('Converting dense tensor to sparse');
        Zcp = sptensor(Zcp);
        Wcp = sptensor(Wcp);
    end
    Zvals = tt_cp_wfg_sparse_setup(Zcp,Wcp);
    fflag = ~isequal(funtype,'sparse_lowmem');
    funhandle = @(x) tt_cp_wfun(Zvals,Wcp,x,normZsqr,fflag);
end
    
%% Initial guess
sz = size(Zcp);
N = length(sz);

if iscell(params.Results.init)
    P0 = params.Results.init;
elseif strcmpi(params.Results.init,'random')
    P0 = cell(N,1);
    for n=1:N
        P0{n} = randn(sz(n),R);
        for j=1:R
            P0{n}(:,j) = P0{n}(:,j) / norm(P0{n}(:,j));
        end
    end
elseif strcmpi(params.Results.init,'nvecs')
    P0 = cell(N,1);
    for n=1:N
        P0{n} = nvecs(Zcp,n,R);
    end
else
    error('Initialization type not supported')
end

%% Fit CP using CP_WOPT by ignoring missing entries
out = feval(opthandle, funhandle, tt_fac_to_vec(P0), options);

P  = full(ktensor(tt_cp_vec_to_fac(out.X,Zcp)));

if Zt.m(i)~=0
    TCS = norm((1-Wcp).*(Zt.initData{i}-P))/norm((1-Wcp).*Zt.initData{i});
    output.TCS = TCS;
    RMSE = norm(Wcp.*(Zt.initData{i}-P))/norm(Wcp);
else
    RMSE = norm(Wcp.*(Zt.initData{i}-P))/norm(Wcp);
end
output.RMSE = RMSE;
output.Init = P0;
output.F = tt_cp_vec_to_fac(out.X,Zcp);
output.ExitFlag  = out.ExitFlag;
output.Iter = out.Iters;
output.FunV = out.F;
  
