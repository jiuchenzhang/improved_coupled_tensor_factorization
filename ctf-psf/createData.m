
%% Simulate dataset generation
function [Z,Ac] = createData(varargin)

params = inputParser;
params.addParamValue('eachSize',  @isnumeric);
params.addParamValue('structSize',  @iscell);
params.addParamValue('lambdas',  @iscell);
params.addParamValue('flag_sparse', @isnumeric);
params.addParamValue('noise', @(x) x > 0);
params.addParamValue('c', @(x) x > 0);
params.addParamValue('sameDim', @(x) x > 0);
params.parse(varargin{:});
%%Parameters
lambdas     = params.Results.lambdas;
%Dimension structure
structSize  = params.Results.structSize;
%Dimension size
eachSize          = params.Results.eachSize;
%Noise level
nlevel      = params.Results.noise;
%Number of shared components
c           = params.Results.c;
%Shared dimension index
sameDim     = params.Results.sameDim;

%% Check parameters
if length(lambdas)~=length(structSize)
    error('There should be weights for each data set');
end
P = length(structSize);
for p=1:P
    l(p) = length(lambdas{p});
end
if length(unique(l))>1
    error('There should be the same number of weights for each data set');
end
%According to lambda to generate the rank of data
Rtotal    = length(lambdas{1});

%% Generate factor matrices
nb_modes  = length(eachSize);
Ac        = cell(nb_modes,1);

for n = 1:nb_modes
    Ac{n} = randn(eachSize(n),c);
    for r=1:c
        Ac{n}(:,r)=Ac{n}(:,r)/norm(Ac{n}(:,r));
    end
end

%Make Shared parts consistent
mz = sameDim;
D = length(mz);
for d = 1:D
    Ac{mz{d}(1)} = Ac{mz{d}(2)};
end

%Generate tensors of Rs
As        = cell(nb_modes,1);
Rs = Rtotal-c;
if Rs > 0
    for n = 1:nb_modes
        As{n} = rand(eachSize(n),Rtotal-c);
        for r=1:Rs
            As{n}(:,r)=As{n}(:,r)/norm(As{n}(:,r));
        end
    end
end

%Form a complete factor matrix
A        = cell(nb_modes,1);
for n = 1:nb_modes
    A{n} = [Ac{n},As{n}];
end

%% Generate tensors and add noise
P  = length(structSize);
X  = cell(P,1);
for p = 1:P
    X{p} = full(ktensor(lambdas{p}',A(structSize{p})));
    N    = tensor(randn(size(X{p})));
    X{p} = X{p} + nlevel*norm(X{p})*N/norm(N);
end

%Assignment parameters
Z.initData = X;
Z.structSize = structSize;
Z.eachSize  = eachSize;
Z.sameDim = sameDim;
