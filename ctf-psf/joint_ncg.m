function out = joint_ncg(FUN,FUNI,Todifx,x0,NumOfObject,shared,varargin)
%Joint-NCG   Nonlinear conjugate gradient minimization.

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

for i = 1:NumOfObject
    [xq{i}, xz{i}]= feval(Todifx,xk,i);
    [fks(i),gks{i}] = feval(FUNI,xq{i},xz{i},i);
end
%joint
[fk,gk] = feval(FUN,xk);
gtemp = gk(find(gk~=0));
g = [];
for i = 1:NumOfObject
    g =[g; gks{i}];
end
g = [g;gtemp];
% out = poblano_out(xk,fk,gk,1,params);
out = poblano_out(xk,fk,g,1,params);

%% Main loop
while out.ExitFlag == -1
    if isequal(shared,0)
        %% individual
        for i = 1:NumOfObject
            [xq{i}, xz{i}]= feval(Todifx,xk,i);
            [fks(i),gks{i}] =  feval(FUNI,xq{i},xz{i},i);
            if out.Iters == 0
                pks{i} = -gks{i};
                aks(i) = 1.0;
                gkTgks(i) = gks{i}'*gks{i};
            else    
            % Compute next direction
            if mod(out.Iters,params.Results.RestartIters) == 0
                % restart to prevent stagnation
                bks(i) = 0;
                pks{i} = -gks{i};
            else
                % direction update
                switch (params.Results.Update)
                    case 'FR'
                        % Fletcher-Reeves
                        gkTgks(i) = gks{i}'*gks{i};
                        if gkTgkolds(i) > 0
                            bks(i) = gkTgks(i)/gkTgkolds(i);
                        else
                            fprintf(1,[mfilename,': warning: bk set to 0\n']);
                            bks(i) = 0;
                        end
                    case 'PR'
                        % Polak-Ribiere
                        gkTgks(i) = gks{i}'*gks{i};
                        gkMgkolds{i} = gks{i}-gkolds{i};
                        if gkTgkolds(i) > 0
                            bks(i) = (gks{i}'*gkMgkolds{i})/gkTgkolds(i);
                        else
                            fprintf(1,[mfilename,': warning: bk set to 0\n']);
                            bks(i) = 0;
                        end
                    case 'HS'
                        % Hestenes-Stiefel
                        gkMgkolds{i} = gks{i}-gkolds{i};
                        denoms(i) = pkolds{i}'*gkMgkolds{i};
                        if denoms(i) > 0
                            bks(i) = (gks{i}'*gkMgkolds{i})/denoms(i);
                        else
%                             fprintf(1,[mfilename,': warning: bk set to 0\n']);
                            bks(i) = 0;
                        end
                    case 'SD'
                        % Steepest Descent
                        bks(i) = 0;
                    otherwise
                        error('Error: options.Update is not valid. Choices are {FR, PR, HS}');
                end
                % do not allow negative conjugate direction weights
                if bks(i) < 0
                    bks(i) = max(0,bks(i));
                end

                % restart method from Nocedal and Wright
                if params.Results.RestartNW
                    v = params.Results.RestartNWTol;
                    if ((gks{i}'*gkolds{i})/(gkTgkolds(i)^2) >= v)
                        bks(i) = 0;
                    end
                end

                % new direction
                pks{i} = -gks{i} + bks(i)*pkolds{i};
              end
            end
            xkolds{i} = xz{i};
            fksold(i) = fks(i);
            gkolds{i} = gks{i};
            pkolds{i} = pks{i};
            gkTgkolds(i) = gkTgks(i);

            % Compute step length
            [xq{i},xz{i},fks(i),gks{i},aks(i),lsinfo(i),nfevs(i)] = poblano_linesearch_ind(FUNI,xq{i},xz{i},fks(i),gks{i},aks(i),pks{i},i,params.Results);
            if (lsinfo(i) ~= 1) && strcmp(params.Results.Display, 'iter')
                fprintf(1,[mfilename,': line search warning = %d\n'],lsinfo(i));
            end
        end

        xk = [];
        for j = 1:NumOfObject
            xk = [xk;xq{j}];
        end
        [fk,gk] = feval(FUN,xk);
    end
    %% joint
    if out.Iters == 0
        pk = -gk;
        ak = 1.0;
        gkTgk = gk'*gk;
    else    
        [fk,gk] = feval(FUN,xk);
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
%                                     fprintf(1,[mfilename,': warning: bk set to 0\n']);
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
    [xk,fk,gk,ak,lsinfo,nfev] = poblano_linesearch_joint(FUN,xk,fk,gk,ak,pk,params.Results);
    if (lsinfo ~= 1) && strcmp(params.Results.Display, 'iter')
         fprintf(1,[mfilename,': line search warning = %d\n'],lsinfo);
    end

    gtemp = gk(find(gk~=0));
    g = [];
    for i = 1:NumOfObject
    g =[g; gks{i}];
    end
    g = [g;gtemp];

    %% Update counts, check exit conditions, etc.
    out = poblano_out(xk,fk,g,nfev,params,out);   
end
