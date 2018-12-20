%% add missing value
function Z = addmv(Z,M,flag_sparse)
X = Z.initData;
P = length(X);
Z.object = cell(P,1);
Z.miss = cell(P,1);
for p=1:P
    W{p}  = tt_create_missing_data_pattern(size(Z.initData{p}), M(p), flag_sparse(p));
    if flag_sparse(p)
        Z.object{p} = W{p}.*sptensor(X{p});
    else
        Z.object{p} = W{p}.*X{p};
    end
    norms(p)    = norm(Z.object{p});
    Z.miss{p}   = W{p};
end