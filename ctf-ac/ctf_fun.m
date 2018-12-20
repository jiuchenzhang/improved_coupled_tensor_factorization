function [f,g] = ctf_fun(x, Z, Znormsqr)
%% Convert the input vector into a cell array of factor matrices
A  = ctf_vec_to_fac(x,Z);

%% Compute the function and gradient values
[f,G] = ctf_fg(Z,A,Znormsqr);

%% Vectorize the cell array of matrices
g = ctf_fac_to_vec(G);