function x =ctf_fac_to_vec(A)
%% Set-up
N = length(A{1});
TT = length(A);
%% Get sizes
sz = zeros(N,1);


R = size(A{1}{1},2);

%% Create x
xi = cell(TT,1);
for t = 1: TT
    for n = 1:N
        sz(n) = size(A{t}{n},1);
    end
    P = sum(sz)*R;
    xi{t} = zeros(P,1);
    for n = 1:N
        idx1 = sum(sz(1:n-1))*R + 1;
        idx2 = sum(sz(1:n))*R;
        xi{t}(idx1:idx2) = reshape(A{t}{n},sz(n)*R,1);
    end
end
x = [];
for t = 1:TT
    x =[x;xi{t}];
end
