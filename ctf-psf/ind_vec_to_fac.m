function [x,A] = ind_vec_to_fac(x,xk,Z,index)
%% Set-up
P   = length(x);
sz  = Z.eachSize(Z.structSize{index});
N   = length(sz);

R = P / sum(sz);

%xk  -> x
ind = Z.sameDim{1}(index);
i = find(Z.structSize{index}==ind);
j=  sum(sz(1:i))*R;
x(j-sz(i)*(R-Z.r)+1:j) =xk;

%% Create A
A = cell(N,1);
for n = 1:N
    idx1 = sum(sz(1:n-1))*R + 1;
    idx2 = sum(sz(1:n))*R;
    A{n} = reshape(x(idx1:idx2),sz(n),R);
end