
function obj_logL=max_LSTVAR_mleScov(param0,resids,Z0,Omega_length,T,Nvar,theta0,scale_penalty)

% Inputs
% Y = vector dependent variables
% X = vector of lags 
% Z0 = regime shifter for contemporanous response (matrices OMEGA) and for
%      the dynamic response (coefs in the VAR)
% Omega_length = length of the vech(cov matrix of residuals)
% T = sample size
% Nvar = number of variables in the VAR

% convert vector of parameters into approapriate matrices for eatch regime
[Omega0,Omega1,theta]=vec2matScov(param0,Omega_length);


% check if eigenvalues of Omega0 and Omega1 are greater than zero
% if not exit the function
if  theta<0
    if min(eig(Omega0))>0 & min(eig(Omega1))>0

        % log-likelihood
        obj_logL=LSTVAR_mleScov(resids,Z0,Omega0,Omega1,theta,T,Nvar);
        
        % take the negative of the log-likelihood (because the objective
        % function is minimization rather than maximization) and add
        % penalty for deviation of theta (the curvature in the transition
        % probability) from the prior. 
        obj_logL=-obj_logL+(theta-theta0)^2*scale_penalty;
        
    else
        obj_logL=1e12;
    end
else
    obj_logL=1e12;
end

end




