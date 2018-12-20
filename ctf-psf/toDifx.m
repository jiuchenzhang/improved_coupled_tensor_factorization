function [x, xn] = toDifx(x0, Z, index)
%% get different part x
R = length(x0)/sum(Z.eachSize);
s = 0;
for i = 1: index
    ind = Z.structSize{i};
    sz = Z.eachSize(ind);
    s =s+sum(sz)*R;
end
idx2 = s;

ind = Z.structSize{index};
sz = Z.eachSize(ind);
idx1 =idx2 -sum(sz)*R+1;
x = x0(idx1:idx2);

ind = Z.sameDim{1}(index);
idx2 = sum(Z.eachSize(1:ind))*R;
idx1 = idx2 -Z.eachSize(ind)*(R-Z.r)+1;
xn = x0(idx1:idx2);