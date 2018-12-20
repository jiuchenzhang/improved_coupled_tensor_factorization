function Z = addMissingVal(Z,M,flag_sparse)
Z.m = M;
sz= Z.size;
modes=  Z.modes;
X = Z.initData;

P = length(X);
for p=1:P
    W{p}  = tt_create_missing_data_pattern(sz(modes{p}), M(p), flag_sparse(p));
    if flag_sparse(p)
        Z.object{p} = W{p}.*sptensor(X{p});
    else
        Z.object{p} = W{p}.*X{p};
    end
    norms(p)    = norm(Z.object{p});
    Z.object{p} = Z.object{p}/norms(p);    
    Z.miss{p}   = W{p};
end