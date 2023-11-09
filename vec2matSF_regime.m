
% convert vector of parameters into approapriate matrices

function [beta0,beta1,Omega0,Omega1,theta,gamma,residsM,const0,const1,logL]=...
    vec2matSF_regime(param0,Omega_length,X,Y,Z0,Z1,T,Nvar,FX,Nforecs,...
    trendON,trend_regspec,inter_regspec)

% Inputs
% Y = vector dependent variables
% X = vector of lags
% Z0 = regime shifter for contemporanous response (matrices OMEGA)
% Z1 = regime shifter for dynamic response (matrices PI)
% Omega_length = length of the vech(cov matrix of residuals)
% T = sample size
% Nvar = number of variables in the VAR

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

[Omega0,Omega1,theta]=vec2matScov(param0,Omega_length);

gamma=theta;

M_Z_t=1./(1+exp(theta*Z0(:,1)));
F_Z_t=1./(1+exp(gamma*Z1(:,1)));

[rX cX]=size(X);

if trendON==0
    if inter_regspec==0
        XM=[X.*repmat(1-F_Z_t,1,cX) X.*repmat(F_Z_t,1,cX) ones(length(F_Z_t),1) FX];
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

% constant term in  regime 0
const0=betaK1(2*cX+1,:)';

% constant term in  regime 1
if inter_regspec==1
    const1=betaK1(2*cX+2,:)';
else
    const1=const0;
end

logL=0;

for i1=1:T

    % compute time-varying parameters
    Omega_t = Omega0*(1-M_Z_t(i1)) + Omega1*M_Z_t(i1);

    % error term
    e_t=residsM(i1,:);

    logL=logL-0.5*( log(det(Omega_t)) + e_t*inv(Omega_t)*e_t');
end


end

