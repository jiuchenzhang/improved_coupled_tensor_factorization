function [f,g,xk] = psfCTF_fun(x, Z, Znormsqr)
%% Convert the input vector into a cell array of factor matrices
A  = psfCTF_vec_to_fac(x,Z);

%% Compute the function and gradient values
[f,G] = psfCTF_fg(Z,A,Znormsqr);

%% Vectorize the cell array of matrices
g = psfCTF_fac_to_vec(G);
xk = psfCTF_fac_to_vec(A);