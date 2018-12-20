function [f,g,x] = ind_fun(x, xk, Z, Znormsqr, index)
%% Convert the input vector into a cell array of factor matrices

[x,A]  = ind_vec_to_fac(x,xk,Z,index);

%% Compute the function and gradient values
[f,G] = ind_fg(Z,A,Znormsqr,index);

% Vectorize the cell array of matrices
g0 = psfCTF_fac_to_vec(G);

R = length(x)/sum(Z.eachSize(Z.structSize{index}));
nu =0;
for i = 1:index-1
    nu = nu +sum(Z.eachSize(Z.structSize{i}))*R;
end


ind = Z.sameDim{1}(index);
idx2 = sum(Z.eachSize(1:ind))*R-nu;
idx1 = idx2 -Z.eachSize(ind)*(R-Z.r)+1;
g = g0(idx1:idx2);