function A = ctf_vec_to_fac(x,Z)
%% Set-up
sz  = Z.csize;
N   = length(sz);
TT = length(Z.object);

R = Z.R;

%% Create A
Ai = cell(N,1);

for n = 1:N
    idx1 = sum(sz(1:n-1))*R + 1;
    idx2 = sum(sz(1:n))*R;
    Ai{n} = reshape(x(idx1:idx2),sz(n),R);
end
A = cell(TT,1);
Ai = reshape(Ai,N/TT,TT)';
for t = 1:TT
    A{t} = Ai(t,:)';
end