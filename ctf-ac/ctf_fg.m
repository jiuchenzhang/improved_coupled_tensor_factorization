function [f,G] = ctf_fg(Z,A,Znormsqr)


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

dim = length(Z.csize)/Z.N;

%% compute object f and g
for p = 1:P
    if length(size(Z.object{p}))>=3
        % Tensor
        if isfield(Z,'miss') && ~isempty(Z.miss{p})
            [fp{p},Gp{p}] = tt_cp_wfg(Z.object{p}, Z.miss{p}, A{p}, Znormsqr{p});
        else
            [fp{p},Gp{p}] = tt_cp_fg(Z.object{p},  A{p}, Znormsqr{p});
        end
    elseif length(size(Z.object{p}))==2
        % Matrix
        if isfield(Z,'miss') && ~isempty(Z.miss{p})
            [fp{p},Gp{p}] = pca_wfg(Z.object{p}, Z.miss{p}, A{p}, Znormsqr{p});
        else
            [fp{p},Gp{p}] = pca_fg(Z.object{p}, A{p}, Znormsqr{p});
        end
    end
end

%% Compute overall function value and gradient
%% compute f
f = 0;

% tensor item
for n = 1:Z.N
    weightsOfT(n) = (numel(double(Z.initData{n})))*( 1-Z.m(n) );
    f = f + 2 * fp{n} / weightsOfT(n);
end
    for cf = 1:Z.pCf
        weightsOfF(cf)  = numel(A{1}{cf});
        for n = 1:Z.N-1
             f = f + ( norm(A{n}{cf} - A{n+1}{cf},'fro')^2 ) / weightsOfF(cf);
        end 
        f = f + ( norm(A{1}{cf} - A{Z.N}{cf},'fro')^2 ) /  weightsOfF(cf);
    end



%% compute g
    for c = 1:dim
        if(c <= Z.pCf)
            for n = 1 :Z.N
                if n == 1
                    Gp{n}{c} = 2  * Gp{n}{c} / weightsOfT(n) + 2*(2*A{n}{c} - A{n+1}{c}-A{Z.N}{c}) / weightsOfF(c);
                elseif n == Z.N
                    Gp{n}{c} = 2  * Gp{n}{c} / weightsOfT(n) + 2*(2*A{n}{c} - A{n-1}{c}-A{1}{c})   / weightsOfF(c);
                else
                    Gp{n}{c} = 2  * Gp{n}{c} / weightsOfT(n) + 2*(2*A{n}{c} - A{n-1}{c}-A{n+1}{c}) / weightsOfF(c);
                end
            end
        else
            for n = 1 :Z.N
                Gp{n}{c} = 2  *  Gp{n}{c} / weightsOfT(n);
            end
        end
    end


G = Gp;
