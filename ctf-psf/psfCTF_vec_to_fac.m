function A = psfCTF_vec_to_fac(x,Z)
%% Set-up
P   = length(x);
sz  = Z.eachSize;
N   = length(sz);

R = P / sum(sz);

%% Create A
A = cell(N,1);
for n = 1:N
    idx1 = sum(sz(1:n-1))*R + 1;
    idx2 = sum(sz(1:n))*R;
    A{n} = reshape(x(idx1:idx2),sz(n),R);
end
mz = Z.sameDim;
r = Z.r;
D = length(mz);
for d = 1:D
    A{mz{d}(2)}(:,1:r) = A{mz{d}(1)}(:,1:r); 
end