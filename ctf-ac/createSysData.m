%% generate synthetic data
function Z = createSysData (varargin)

params = inputParser;
params.addParamValue('size',  @isnumeric);
params.addParamValue('lambdas',  @iscell);
params.addParamValue('rRank', @isnumeric);
params.addParamValue('rCf', @isnumeric);
params.addParamValue('noise', @(x) x > 0);

params.parse(varargin{:});
%%Parameters
size        = params.Results.size;
lambdas     = params.Results.lambdas;
nlevel      = params.Results.noise;
rRank       = params.Results.rRank;
rCf          = params.Results.rCf;

TenNum  = length(lambdas);
dim = length(size);
factor = cell(TenNum,1);

for n = 1:TenNum
    factor{n} = cell(dim,1);
    for d = 1:dim
        factor{n}{d} = rand(size(d),rRank);
    end
end

if rCf>0
    for cf = 1:rCf
         for n = 1:TenNum
            if n==1
                factor{n}{cf} = normrnd(0,1,size(cf),rRank);
            else
                factor{n}{cf} =rand(size(cf),rRank).*factor{1}{cf}+ rand(size(cf),rRank);
                %factor{n}{cf} = normrnd(n*0.5,1,size(cf),rRank);
            end
        end
    end
end

X  = cell(TenNum,1);
for n = 1:TenNum
    X{n} = full(ktensor(lambdas{n}',factor{n}));
    N    = tensor(randn(size));
    X{n} = X{n}/norm(X{n});
    T{n} = X{n} + nlevel*norm(X{n})*N/norm(N);
end
Z.noNoiseData   = X;
Z.initData      = T;
Z.size          = size;
Z.cf            = rCf;