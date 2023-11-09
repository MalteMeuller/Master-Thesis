
% convert vector of parameters into approapriate matrices

function [beta0,beta1,Omega0,Omega1,theta,gamma,residsM,beta_const, beta_other,logL]=...
    vec2matSF_struct(param0,s_data,s_desc,s_prior,s_exo,s_MCMC)

% Inputs
% Y = vector dependent variables
% X = vector of lags
% Z0 = regime shifter for contemporanous response (matrices OMEGA)
% Z1 = regime shifter for dynamic response (matrices PI)
% Omega_length = length of the vech(cov matrix of residuals)
% T = sample size
% Nvar = number of variables in the VAR
% trendON = include trend: yes = 1 ; no = 0
% inter_regspec = allow for regime specific intercept: yes = 1; no = 0
% trend_regspec = allow for regime specific trend: yes = 1; no = 0
% FX = list of exogenous variables
% Nforecs= number of exogenous variables

% Outputs:
% beta0 = slopes in VAR in regime 0
% beta1 = slopes in VAR in regime 1
% Omega0 = cov matrix of residuals in VAR in regime 0
% Omega1 = cov matrix of residuals in VAR in regime 1
% theta = curvature in transition probability (for slopes)
% gamma = curvature in transition probability (for cov matrix of resids)
% residsM = residuals
% const0 = constant term in regime 0
% const1 = constant term in regime 1
% logL = log-likelihood

%-------------------------------------------------------------------------

Y=s_data.Y0;
X=s_data.X;
Z0=s_data.Z0;
Z1=Z0;

Omega_length=s_desc.Omega_length_cov;
T=s_desc.T;
Nvar = s_desc.Nvar;

trendON = s_exo.trendON;
inter_regspec = s_exo.inter_regspec;
trend_regspec = s_exo.trend_regspec;
FX=s_exo.FX;
Nforecs=s_exo.Nforecs;

%-------------------------------------------------------------------------

[Omega0,Omega1,theta]=vec2matScov(param0,Omega_length);

gamma=theta;

M_Z_t=1./(1+exp(theta*Z0(:,1)));
F_Z_t=1./(1+exp(gamma*Z1(:,1)));

[rX cX]=size(X);

if trendON==0
    if inter_regspec==0
        if s_desc.VARlags>0
            XM=[X.*repmat(1-F_Z_t,1,cX) X.*repmat(F_Z_t,1,cX) ones(length(F_Z_t),1) FX]; %only using this one
        else
            XM= ones(length(F_Z_t),1);
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

num0=0;
den0=0;

for t=1:T
    wgtM=inv(Omega0*(1-M_Z_t(t)) + Omega1*M_Z_t(t));
    den0 = den0 + kron(wgtM,XM(t,:)'* XM(t,:));
    num0 = num0 + vec(XM(t,:)'*Y(t,:)*wgtM);
end
betaK=inv(den0)*num0;

betaK1=unvec(betaK,cX*2+trendON*(1+trend_regspec)+1+inter_regspec+Nforecs,Nvar);

% compute residuals
residsM=Y-XM*betaK1;

% VAR coefficients in regime 0
beta0=betaK1(1:cX,:)';
% VAR coefficients in regime 1
beta1=betaK1(cX+1:2*cX,:)';

% compute coefs on the constant term and any other control variables
beta_const=betaK1(2*cX+1:2*cX+1,:)';
beta_other=betaK1(2*cX+2:end,:)';

logL=0;

for i1=1:T
    
    % compute time-varying parameters
    Omega_t = Omega0*(1-M_Z_t(i1)) + Omega1*M_Z_t(i1);
    
    % error term
    e_t=residsM(i1,:);
    
    logL=logL-0.5*( log(det(Omega_t)) + e_t*inv(Omega_t)*e_t');
end


end

