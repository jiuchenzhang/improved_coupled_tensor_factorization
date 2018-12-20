function [f,G] = psfCTF_fg(Z,A,Znormsqr)

if isa(A,'ktensor')
    A = tocell(A);
end

if ~iscell(A)
    error('A must be a cell array');
end

P = numel(Z.object);

if ~exist('Znormsqr','var')
    Znormsqr = cell(P,1);
    for p = 1:P
        if isa(Z.object{p},'tensor') || isa(Z.object{p},'sptensor')
            Znormsqr{p} = norm(Z.object{p})^2;
        else
            Znormsqr{p} = norm(Z.object{p},'fro')^2;
        end
    end
end

fp = cell(P,1);
Gp = cell(P,1);
for p = 1:P   
    if length(size(Z.object{p}))>=3
        % Tensor
        if isfield(Z,'miss') && ~isempty(Z.miss{p})            
            [fp{p},Gp{p}] = tt_cp_wfg(Z.object{p}, Z.miss{p}, A(Z.structSize{p}), Znormsqr{p});                
        else            
            [fp{p},Gp{p}] = tt_cp_fg(Z.object{p}, A(Z.structSize{p}), Znormsqr{p});
        end       
    elseif length(size(Z.object{p}))==2
        % Matrix
        if isfield(Z,'miss') && ~isempty(Z.miss{p})
            [fp{p},Gp{p}] = pca_wfg(Z.object{p}, Z.miss{p}, A(Z.structSize{p}), Znormsqr{p});
        else        
            [fp{p},Gp{p}] = pca_fg(Z.object{p}, A(Z.structSize{p}), Znormsqr{p});        
        end
    end
end

%% Compute overall gradient
G = cell(size(A));
for n = 1:numel(G)
    G{n} = zeros(size(A{n}));
end
for p = 1:P
    for i = 1:length(Z.structSize{p})
        j = Z.structSize{p}(i);
        G{j} = Gp{p}{i};
    end
end
R = length(A{1}(1,:));
C = Z.sameDim;
r = Z.r;
for m = 1: length(C)
    Gtemp = G{C{m}(1)}+G{C{m}(2)};
    G{C{m}(1)}(:,1:r) = Gtemp(:,1:r);
    G{C{m}(1)}(:,r+1:R) =0;
    G{C{m}(2)}(:) =0;
end

%% Compute overall function value
f = sum(cell2mat(fp));

return;
