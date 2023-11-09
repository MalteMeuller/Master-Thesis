% this function computes the likelihood for the process the VAR with a
% smooth transtion state
% Transition for state is given by two functions
%   slopes/intercept     = F(Z) = exp(gamma*Z)/(1+exp(gamma*Z))
%   cov matrix of shocks = M(Z) = exp(theta*Z)/(1+exp(theta*Z))

% identifying assuumtpions are that gamma<0, theta<0

% compute VAR coefficients using GLS

% Inputs
% Y = vector dependent variables
% X = vector of lags
% Z0 = regime shifter for contemporanous response (matrices OMEGA)
% Z1 = regime shifter for dynamic response (matrices PI)
% Omega0 = cov matrix of residuals in regime 0
% Omega1 = cov matrix of residuals in regime 1
% gamma = transition speed for VAR coefs
% theta = transition speed for coevariance matrix of VAR residuals

% T = sample size
% Nvar = number of variables in the VAR
% FX = forecast series
% Nforecs = number of forecasted periods

% trendON: 1 = include a time trend in the VAR; 
%           0 = do not include time trend in the VAR
% 
% trend_regspec:  1 = allow regime specific time trend
%                 0 = do not allow regime specific time trend
% 
% inter_regspec=1: 1 = allow regime specific intercept
%                  0 = do not allow regime specific intercept
                  
                  
function [logL, beta0, beta1]=LSTVAR_mleS_struct_linears(Omega0,s_data,s_desc,s_exo);

Omega1=Omega0;

Y=s_data.Y0;
X=s_data.X;
Z0=s_data.Z0;
Z1=Z0;
T=s_desc.T;
Nvar=s_desc.Nvar;
FX=s_exo.FX;
Nforecs=s_exo.Nforecs;
trendON=s_exo.trendON;
trend_regspec=s_exo.trend_regspec;
inter_regspec=s_exo.inter_regspec;


% weights on the "expansion" regime
M_Z_t=ones(size(Z0(:,1)));
F_Z_t=ones(size(Z0(:,1)));

% construct vector of regressors
% the last two entries capture constant terms and trends
[rX cX]=size(X);

if trendON==0
    if inter_regspec==0
        XM=[X  ones(length(F_Z_t),1) FX];
    else
        XM=[X  (1-F_Z_t) F_Z_t FX];
    end
else
    if inter_regspec==0 & trend_regspec==0
        XM=[X  ones(length(F_Z_t),1) (1:length(F_Z_t))' FX];
    elseif inter_regspec==1 & trend_regspec==0
        XM=[X  (1-F_Z_t) F_Z_t (1:length(F_Z_t))' FX];
    elseif inter_regspec==0 & trend_regspec==1
        XM=[X  ones(length(F_Z_t),1) (1-F_Z_t).*(1:length(F_Z_t))' F_Z_t.*(1:length(F_Z_t))' FX];
    elseif inter_regspec==1 & trend_regspec==1
        XM=[X (1-F_Z_t) F_Z_t (1-F_Z_t).*(1:length(F_Z_t))' F_Z_t.*(1:length(F_Z_t))' FX];        
    end
end


% XM=[X.*repmat(1-F_Z_t,1,cX) X.*repmat(F_Z_t,1,cX) ones(length(F_Z_t),1) FX];
% estimate VAR coefs using GLS
num0=0;
den0=0;

for t=1:T
    wgtM=inv(Omega0);
    den0 = den0 + kron(wgtM,XM(t,:)'* XM(t,:));
    num0 = num0 + vec(XM(t,:)'*Y(t,:)*wgtM);
end
betaK=inv(den0)*num0;
betaK1=unvec(betaK,cX+trendON*(1+trend_regspec)+1+inter_regspec+Nforecs,Nvar);

% VAR coefficients in regime 0
beta0=betaK1(1:cX,:)';
% VAR coefficients in regime 1
beta1=beta0;

% compute residuals
residsM=Y-XM*betaK1;
% compute log likelihood
logL=0;

for i1=1:T

    % compute time-varying parameters
    Omega_t = Omega0;

    % error term
    e_t=residsM(i1,:);

    logL=logL-0.5*( log(det(Omega_t)) + e_t*inv(Omega_t)*e_t');
end

logL=logL; %-T*Nvar*log(2*3.14);
