function out = convergence(x,f,g,nfev,params,out)
%% Check if this is the first call
if ~exist('out','var')
    out.Params = params;
    out.ExitFlag = -1;
    out.ExitDescription = '';
    out.X = x;
    out.F = f;
    out.G = g;
    out.FuncEvals = nfev;
    out.Iters = 0;
    if params.Results.TraceX, out.TraceX=[]; end
    if params.Results.TraceFunc, out.TraceFunc=[]; end
    if params.Results.TraceRelFunc, out.TraceRelFunc=[]; end
    if params.Results.TraceGrad, out.TraceGrad=[]; end
    if params.Results.TraceGradNorm, out.TraceGradNorm=[]; end
    if params.Results.TraceFuncEvals, out.TraceFuncEvals=[]; end
else
    oldf = out.F;
    if (f <= out.F)
        out.X = x;
        out.F = f;
        out.G = g;
    end
    if abs(oldf) < eps
        relfit = abs(f - oldf);
    else
        relfit = abs((f - oldf)/oldf);
    end
    out.FuncEvals = out.FuncEvals + nfev;  
    out.Iters = out.Iters + 1;
end

if params.Results.TraceX
    out.TraceX(:,end+1) = x;
end

if params.Results.TraceFunc 
    out.TraceFunc(end+1) = f; 
end

if  params.Results.TraceFunc && params.Results.TraceRelFunc 
    if (out.Iters > 0)
        out.TraceRelFunc(end+1) = relfit;
    end
end

if params.Results.TraceFuncEvals
    out.TraceFuncEvals(end+1) = nfev;
end

if out.Iters >= params.Results.MaxIters
    % maximum iterations exceeded
    out.ExitFlag = 1;
    out.ExitDescription = 'Maximum number of iterations exceeded';
elseif out.FuncEvals >= params.Results.MaxFuncEvals
    % maximum function evaluations exceeded
    out.ExitFlag = 2; 
    out.ExitDescription = 'Maximum number of function evaluations exceeded';
elseif (out.Iters > 0) && (relfit <= params.Results.RelFuncTol)
%     elseif (out.Iters > 0) && (relfit <= params.Results.RelFuncTol)
    % relative fit tolerance is reached
    out.ExitFlag=3;
    out.ExitDescription = 'Relative change in F < RelFuncTol';
end

%% Display iteration information

% First Iteration
if (out.Iters == 0) && ~strcmp(params.Results.Display, 'off')
    fprintf(1,' Iter     FuncEvals          F(X)             \n');
    fprintf(1,'------   ---------     ----------------\n');
end
    
% Iteration info
if strcmp(params.Results.Display, 'iter') && (mod(out.Iters,params.Results.DisplayIters)==0 || out.ExitFlag>=0)
    fprintf(1,'%6d %9d   %16.16f \n', out.Iters, out.FuncEvals, ...
            f);
end

% Final Iteration
if ((out.ExitFlag >= 0) && strcmp(params.Results.Display, 'final'))
%    if strcmp(params.Results.Display, 'iter')
%        fprintf(1,'------ --------- ---------------- ----------------\n');
%    end
    fprintf(1,'%6d %9d    %16.16f\n', out.Iters, out.FuncEvals, ...
            out.F);
end