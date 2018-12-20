function x =psfCTF_fac_to_vec(A)
%% Set-up
N = length(A);

%% Get sizes
sz = zeros(N,1);
for n = 1:N
    sz(n) = size(A{n},1);
end

R = size(A{1},2);
P = sum(sz)*R;

%% Create x
x = zeros(P,1);
for n = 1:N
    idx1 = sum(sz(1:n-1))*R + 1;
    idx2 = sum(sz(1:n))*R;
    x(idx1:idx2) = reshape(A{n},sz(n)*R,1);
end