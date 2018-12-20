function out = ctf_ncg(FUN,x0,varargin)
%% Parse parameters

% Create parser
params = inputParser;

% Set Poblano parameters
params = poblano_params(params);

% Set parameters for this method
params.addParamValue('RestartIters',20,@(x) x > 0);
params.addParamValue('Update','HS',@(x) ismember(x,{'FR','PR','HS','SD'}));
params.addParamValue('RestartNW',false,@islogical);
params.addParamValue('RestartNWTol',0.1,@(x) x > 0);

% Parse input
params.parse(varargin{:});

%% Check input arguments
if (nargin == 1) && isequal(FUN,'defaults') && (nargout == 1)
    out = params.Results;
    return;
elseif (nargin < 2)
    error('Error: invalid input arguments');
end

%% Initialize
xk = x0;
[fk,gk] = feval(FUN,xk);
out = convergence(xk,fk,gk,1,params);

%% Main loop
while out.ExitFlag == -1
    if out.Iters == 0
        pk = -gk;
        ak = 1.0;
        gkTgk = gk'*gk;
    else    
        % Compute next direction
        if mod(out.Iters,params.Results.RestartIters) == 0
            % restart to prevent stagnation
            bk = 0;
            pk = -gk;
        else
            % direction update
            switch (params.Results.Update)
                case 'FR'
                    % Fletcher-Reeves
                    gkTgk = gk'*gk;
                    if gkTgkold > 0
                        bk = gkTgk/gkTgkold;
                    else
                        fprintf(1,[mfilename,': warning: bk set to 0\n']);
                        bk = 0;
                    end
                case 'PR'
                    % Polak-Ribiere
                    gkTgk = gk'*gk;
                    gkMgkold = gk-gkold;
                    if gkTgkold > 0
                        bk = (gk'*gkMgkold)/gkTgkold;
                    else
                        fprintf(1,[mfilename,': warning: bk set to 0\n']);
                        bk = 0;
                    end
                case 'HS'
                    % Hestenes-Stiefel
                    gkMgkold = gk-gkold;
                    denom = pkold'*gkMgkold;
                    if denom > 0
                        bk = (gk'*gkMgkold)/denom;
                    else
                        fprintf(1,[mfilename,': warning: bk set to 0\n']);
                        bk = 0;
                    end
                case 'SD'
                    % Steepest Descent
                    bk = 0;
                otherwise
                    error('Error: options.Update is not valid. Choices are {FR, PR, HS}');
            end
            % do not allow negative conjugate direction weights
            if bk < 0
                bk = max(0,bk);
            end
             % restart method from Nocedal and Wright
             if params.Results.RestartNW
                 v = params.Results.RestartNWTol;
                 if ((gk'*gkold)/(gkTgkold^2) >= v)
                     bk = 0;
                 end
             end

            % new direction
            pk = -gk + bk*pkold;
        end
    end
    xkold = xk;
    fkold = fk;
    gkold = gk;
    pkold = pk;
    gkTgkold = gkTgk;
    
    % Compute step length
    [xk,fk,gk,ak,lsinfo,nfev] = poblano_linesearch(FUN,xk,fk,gk,ak,pk,params.Results);
    
    out = convergence(xk,fk,gk,nfev,params,out);    
end