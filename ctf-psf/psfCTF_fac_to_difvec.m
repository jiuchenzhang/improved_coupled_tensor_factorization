function x = psfCTF_fac_to_difvec(A,Z)
%% Set-up

R = size(A{1},2);

N  = length(Z.object);

for  n = 1:N
    ind = Z.structSize{n};
    sz = Z.eachSize(ind);
    P = sum(sz)*R;
    x{n} = zeros(P,1);
    Mn = length(sz);
    for m = 1:Mn
        idx1 = sum(sz(1:m-1))*R + 1;
        idx2 = sum(sz(1:m))*R;
        x{n}(idx1:idx2) = reshape(A{ind(m)},sz(m)*R,1);
    end
end
