% this function computes the likelihood for the process the VAR with a
% smooth transtion (observed) state
% Transition for state is given by two functions
%   slopes/intercept     = F(Z) = exp(gamma*Z)/(1+exp(gamma*Z))
%   cov matrix of shocks = M(Z) = exp(theta*Z)/(1+exp(theta*Z))

% identifying assuumtpions are that gamma<0, theta<0

% compute VAR coefficients using GLS

% Inputs
% Y = vector dependent variables
% X = vector of lags
% Z0 = regime shifter for contemporanous response (matrices OMEGA)
% Omega0 = cov matrix of residuals in regime 0
% Omega1 = cov matrix of residuals in regime 1
% gamma = transition speed for VAR coefs
% theta = transition speed for coevariance matrix of VAR residuals

% T = sample size
% Nvar = number of variables in the VAR

function logL=LSTVAR_mleScov(resids,Z0,Omega0,Omega1,theta,T,Nvar)


% weights on the "expansion" regime
M_Z_t=1./(1+exp(theta*Z0(:,1)));

logL=0;

for i1=1:T

    % compute time-varying parameters
    Omega_t = Omega0*(1-M_Z_t(i1)) + Omega1*M_Z_t(i1);

    % error term
    e_t=resids(i1,:);

    logL=logL-0.5*( log(det(Omega_t)) + e_t*inv(Omega_t)*e_t');
end

logL=logL; %-T*Nvar*log(2*3.14);
