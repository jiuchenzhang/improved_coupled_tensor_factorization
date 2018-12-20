function [out omega]  = conver(x,f,params,out)
%% Check if this is the first call
if ~exist('out','var')
    out.Params = params;
    out.ExitFlag = -1;
    out.F = f;
    out.Iters = 0;
    out.Fun(1) = f;
else
    oldf = out.F;
    if (f <= out.F)
        out.X = x;
        out.F = f;
    end
    relfit = abs((f - oldf)/oldf);
    out.Iters = out.Iters + 1;
end

if out.Iters >= params.Results.MaxIters
    % maximum iterations exceeded
    out.ExitFlag = 1;
    out.ExitDescription = 'Maximum number of iterations exceeded';
elseif (out.Iters > 0) && (relfit <= params.Results.RelFuncTol)
    out.ExitFlag=2;
    out.ExitDescription = 'Relative change in F < RelFuncTol';
end

%% Display iteration information
% First Iteration
if (out.Iters == 0) && ~strcmp(params.Results.Display, 'off')
    fprintf(1,' Iter         F(X)             \n');
    fprintf(1,'------  ----------------\n');
end

% Iteration info
if strcmp(params.Results.Display, 'iter') && (mod(out.Iters,params.Results.DisplayIters)==0 || out.ExitFlag>=0)
    fprintf(1,'%6d %16.8f \n', out.Iters,f);
end

% Final Iteration
if ((out.ExitFlag >= 0) && strcmp(params.Results.Display, 'final'))
    fprintf(1,'%6d %16.8f\n', out.Iters,out.F);
end

out.Fun(out.Iters+1) = f;
