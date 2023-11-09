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
                  
                  
function [logL, beta0, beta1, beta_const, beta_other]=...
    LSTVAR_mleS_struct(Omega0,Omega1,gamma,theta,s_data,s_desc,s_exo);



Y=s_data.Y0;
X=s_data.X;
Z0=s_data.Z0;
Z1=Z0; %same 
T=s_desc.T;
Nvar=s_desc.Nvar;
FX=s_exo.FX;
Nforecs=s_exo.Nforecs;
trendON=s_exo.trendON;
trend_regspec=s_exo.trend_regspec;
inter_regspec=s_exo.inter_regspec;


% weights on the "expansion" regime
%M_Z_t=exp(theta*Z0(:,1))./(1+exp(theta*Z0(:,1)));
%F_Z_t=exp(gamma*Z1(:,1))./(1+exp(gamma*Z1(:,1)));

%testar logistics istället
M_Z_t=1./(1+exp(theta*Z0(:,1)));
F_Z_t=1./(1+exp(gamma*Z1(:,1)));

% construct vector of regressors
% the last two entries capture constant terms and trends
[rX cX]=size(X);

if trendON==0
    if inter_regspec==0
        if s_desc.VARlags>0
        XM=[X.*repmat(1-F_Z_t,1,cX) X.*repmat(F_Z_t,1,cX) ones(length(F_Z_t),1) FX];
        else
            XM=ones(length(F_Z_t),1);
        end
    else
        XM=[X.*repmat(1-F_Z_t,1,cX) X.*repmat(F_Z_t,1,cX) (1-F_Z_t) F_Z_t FX];
    end
else
    if inter_regspec==0 & trend_regspec==0
        XM=[X.*repmat(1-F_Z_t,1,cX) X.*repmat(F_Z_t,1,cX) ones(length(F_Z_t),1) (1:length(F_Z_t))' FX];
    elseif inter_regspec==1 & trend_regspec==0
        XM=[X.*repmat(1-F_Z_t,1,cX) X.*repmat(F_Z_t,1,cX) (1-F_Z_t) F_Z_t (1:length(F_Z_t))' FX];
    elseif inter_regspec==0 & trend_regspec==1
        XM=[X.*repmat(1-F_Z_t,1,cX) X.*repmat(F_Z_t,1,cX) ones(length(F_Z_t),1) (1-F_Z_t).*(1:length(F_Z_t))' F_Z_t.*(1:length(F_Z_t))' FX];
    elseif inter_regspec==1 & trend_regspec==1
        XM=[X.*repmat(1-F_Z_t,1,cX) X.*repmat(F_Z_t,1,cX) (1-F_Z_t) F_Z_t (1-F_Z_t).*(1:length(F_Z_t))' F_Z_t.*(1:length(F_Z_t))' FX];        
    end
end


% XM=[X.*repmat(1-F_Z_t,1,cX) X.*repmat(F_Z_t,1,cX) ones(length(F_Z_t),1) FX];
% estimate VAR coefs using GLS
num0=0;
den0=0;

for t=1:T
    wgtM=inv(Omega0*(1-M_Z_t(t)) + Omega1*M_Z_t(t));
    den0 = den0 + kron(wgtM,XM(t,:)'* XM(t,:));
    num0 = num0 + vec(XM(t,:)'*Y(t,:)*wgtM);
end
betaK=inv(den0)*num0;
betaK1=unvec(betaK,cX*2+trendON*(1+trend_regspec)+1+inter_regspec+Nforecs,Nvar);

% VAR coefficients in regime 0
beta0=betaK1(1:cX,:)';
% VAR coefficients in regime 1
beta1=betaK1(cX+1:2*cX,:)';

% compute coefs on the constant term and any other control variables
beta_const=betaK1(2*cX+1:2*cX+1,:)';
beta_other=betaK1(2*cX+2:end,:)';  %se här om variablerna är med

% compute residuals
residsM=Y-XM*betaK1;
% compute log likelihood
logL=0;

for i1=1:T

    % compute time-varying parameters, här skapas parametrarna
    Omega_t = Omega0*(1-M_Z_t(i1)) + Omega1*M_Z_t(i1);

    % error term
    e_t=residsM(i1,:);

    logL=logL-0.5*( log(det(Omega_t)) + e_t*inv(Omega_t)*e_t');
end

logL=logL; %-T*Nvar*log(2*3.14);
