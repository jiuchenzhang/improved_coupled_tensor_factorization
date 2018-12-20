function U = psfCTF_nvecs(Z,n,r)

P = length(Z.object);


A = [];
for p = 1:P
    idx = find(Z.structSize{p} == n);
    for i = idx
        if (isa(Z.object{p},'tensor') && (length(size(Z.object{p}))>=3))
            A = [A, double(tenmat(Z.object{p},i))];
        elseif (isa(Z.object{p}, 'sptensor')  && (length(size(Z.object{p}))>=3))
            A = [A, double(sptenmat(Z.object{p},i))];
        else
            if i == 1
                A = [A double(Z.object{p})];
            elseif i == 2
                A = [A double(Z.object{p})'];
            end
        end
        break
    end
end


Y = A*A';
           
[U,~] = eigs(Y, r, 'LM');
