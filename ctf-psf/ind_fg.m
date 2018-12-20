function [f,G] = ind_fg(Z,A,Znormsqr,index)

if isa(A,'ktensor')
    A = tocell(A);
end

if ~iscell(A)
    error('A must be a cell array');    
end

 
if length(size(Z.object{index}))>=3
    % Tensor
    if isfield(Z,'miss') && ~isempty(Z.miss{index})            
        [fp,Gp] = tt_cp_wfg(Z.object{index}, Z.miss{index}, A, Znormsqr{index});                
    else            
        [fp,Gp] = tt_cp_fg(Z.object{index}, A, Znormsqr{index});
    end       
elseif length(size(Z.object{index}))==2
    % Matrix
    if isfield(Z,'miss') && ~isempty(Z.miss{index})
        [fp,Gp] = pca_wfg(Z.object{index}, Z.miss{index}, A, Znormsqr{index});
    else        
        [fp,Gp] = pca_fg(Z.object{index}, A, Znormsqr{index});        
    end
end


%% Compute overall gradient

Gi = Gp;
for i = 1:length(Gp)
    Gi{i}(:)= 0;
end

C = Z.sameDim;
r = Z.r;
for m = 1: length(C)
    for t = 1:2
        if ismember(C{m}(t),Z.structSize{index})
            ti = find(Z.structSize{index} == C{m}(t));
            Gi{ti}(:,r+1:end) =Gp{ti}(:,r+1:end);
%             Gp{ti}(:,1:r) =0;
            break;
        end
    end
end

%% Compute overall function value
f = fp;
G = Gi;
return;
