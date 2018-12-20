function [xk,f,g,a,lsinfo,lsnfev] = poblano_linesearch_joint(FUN,x0,f0,g0,a0,d0,params)


if nargin < 7
    error('POBLANO_LINESEARCH: too few arguments');
end

switch (params.LineSearch_method)
    case 'more-thuente'
        minfun = 'cvsrch_joint';
end

% Check whether user specified an initial step
if params.LineSearch_initialstep > 0
    a0 = params.LineSearch_initialstep;
end

[xk,f,g,a,lsinfo,lsnfev] = feval(minfun,FUN,x0,f0,g0,a0,d0,params);
