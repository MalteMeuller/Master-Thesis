% this code computes the MLE estimate with time-varying volatility and
% slopes/interpect in the VAR. The time variation is given by variable Z
% (which can later be extended to include a vector of variables)
% Assume for now that there is no feedback from X to Z.

% this code compute the slopes in the VAR using GLS isntead of non-linear
% iterations

% this code also allows for itnerecept which can vary across regimes

% this code computes state dependent impulse response
function []=code_fullsample_linear(VARlags000)
% clear all
% close all
warning off all
tic
randn('seed',1234567)	% set the seed for random generator

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           		PRELIMINARY STEPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---------------------------------------------------------------------
%           Sample
%---------------------------------------------------------------------
s_desc.startS  = 1;
s_desc.endS    = 113;
s_desc.sample0 = s_desc.startS:1:s_desc.endS;

%---------------------------------------------------------------------
%            Conventions
%---------------------------------------------------------------------
% Regimes
%  region 0 = recession
%  region 1 = expansion

%---------------------------------------------------------------------
%            Parameters governing the model
%---------------------------------------------------------------------
% VAR parameters
s_desc.VARlags=VARlags000;      % number of lags in the VAR
s_desc.Nvar=4;         % number of variables in the VAR
s_desc.irf_hor=20;     % number of periods in the impulse response
percentile = 60



% Regime specific intercept and trend
s_exo.trendON=0;      % 1 = include a time trend in the VAR;
% 0 = do not include time trend in the VAR

s_exo.trend_regspec=0;  % 1 = allow regime specific time trend
% 0 = do not allow regime specific time trend

s_exo.inter_regspec=0;  % 1 = allow regime specific intercept
% 0 = do not allow regime specific intercept

%---------------------------------------------------------------------
%                   Bayesian VAR parameters
%---------------------------------------------------------------------
s_prior.flagON=0;   % 1 = use a prior in estimation of VAR slopes
% 0 = do not use a prior in estimation of VAR slopes

s_prior.flagOmegaON=0;    % 1 = use a prior in estimation in cov matrix of VAR residuals
% 0 = do not use a prior in cov matrix of VAR
% residuals

s_prior.flagOLS=0;  % 1 = use OLS estimates of the VAR as a prior
% 0 = use Minnesota prior for the VAR (note that MN
% prior does not have a prior for the covariance matrix
% of the residuals)



% penalty for theta which governs the curvature in the transition
% probability function
s_prior.theta0              = -2;    % the slope in the switching equation
s_prior.scale_penalty       = 1e6;
s_prior.scale_penalty_beta  = 0;

% Parameters of the Minnesota prior: random walk for each variable in the VAR
s_prior.lambda_prior=1;   % controls cross-variable responses (smaller values mean stronger prior
s_prior.theta_prior=0.5;  % controls own response of variables

%---------------------------------------------------------------------
%                       MCMC step
%---------------------------------------------------------------------

s_MCMC.runMCMCflag=1; % = 0 do not use MCMC
% = 1 use MCMC

s_MCMC.MCMCn=100000; % number of MCMC draws
s_MCMC.dropfirstT=min(round(s_MCMC.MCMCn*0.1),20000); % burn-in period for MCMC

%---------------------------------------------------------------------
%           			Import data
%---------------------------------------------------------------------


%download variables
data0= readtable("C:\Users\malte\OneDrive\Desktop\thesis\Data\Clean\df_qoqld.csv");

xo = data0(1:113, {'ld_kixcpi', 'ld_kixgdp', 'kixpol_hp'});
xo = table2array(xo);

% switching series: output in data_SV
data_SV = data0(1:113,"inflation");

%---------------------------------------------------------------------
%                    create series for VAR
%---------------------------------------------------------------------

%this is where our data should be created. 
% Core series for the VAR
Y = data0(1:113, {'ld_gdp',  'ld_kpif','polrate_hp', 'ld_kix'});
Y = table2array(Y);
%change here the position of KIX
s_desc.KIX=4 ;
s_desc.KPIF=2;

% Switching variable for the VAR
Z=data_SV; % threshold
Z = table2array(Z);
% adjust the switching variable to have suffiently long recession periods
Th = prctile(Z, percentile), %we do a 70/30 split

Z0=Z-Th;


% vector of exogenous variables
s_exo.FX=xo;
s_exo.FX=s_exo.FX(s_desc.VARlags+1:end,:);

% number of variables in the vector of forecasts
s_exo.Nforecs=3;
% vector of exogenous variables

% NEED: introduce the flag for a trend in the VAR
% create the set of basic regressors
X=[];
for i=1:s_desc.VARlags
    X=[X Y(s_desc.VARlags+1-i:end-i,1:s_desc.Nvar)];
end
[rX cX]=size(X);
Y0=Y(s_desc.VARlags+1:end,1:s_desc.Nvar);

Z0=Z0(s_desc.VARlags+1:end,:);

% this is the first order expasion given the guess about theta (the
% curvature in the transition probability)
Xm=[X X.*kron(Z0,ones(1,cX)) X.*kron(Z0.^2.*sign(Z0),ones(1,cX)) ones(length(Y0),1) s_exo.FX];


s_desc.T=length(Y0);    % sample size

s_data.Y0=Y0;
s_data.X=X;
s_data.Z0=Z0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           	run VARs for initial values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---------------------------------------------------------------------
%            non-linear VAR to generate residuals
%---------------------------------------------------------------------
% estimate parameters
beta=inv(Xm'*Xm)*(Xm'*Y0);
% compute residuals
resids=Y0-Xm*beta;
% compute the covariance matrix of residuals
Omega=cov(resids);
beta_var=kron(Omega,inv(Xm'*Xm));

s_lag.AICnonlin=log(det(Omega))+2/length(resids)*length(beta);
s_lag.BICnonlin=log(det(Omega))+log(length(resids))/length(resids)*length(beta);

%---------------------------------------------------------------------
%            linear VAR to generate residuals
%---------------------------------------------------------------------
Xm1=[X ones(length(Y0),1) s_exo.FX];
% estimate parameters
betaVARLIN=inv(Xm1'*Xm1)*(Xm1'*Y0);
residsVARLIN=Y0-Xm1*betaVARLIN;
% compute the covariance matrix of residuals
OmegaVARLIN=cov(residsVARLIN);
betaVARLIN=betaVARLIN(1:end-1-s_exo.Nforecs,:);

VARLIN.betaVARLIN=betaVARLIN;
VARLIN.OmegaVARLIN=OmegaVARLIN;
VARLIN.residsVARLIN=residsVARLIN;

s_lag.AIClin=log(det(OmegaVARLIN))+2/length(resids)*length(betaVARLIN);
s_lag.BIClin=log(det(OmegaVARLIN))+log(length(resids))/length(resids)*length(betaVARLIN);

disp('Lag length criteria')
disp(['    Lag = ' num2str(s_desc.VARlags)])
disp('     Non-linear model')
disp(['         AIC =' num2str(s_lag.AICnonlin)])
disp(['         BIC =' num2str(s_lag.BICnonlin)])
disp('     Linear model')
disp(['         AIC =' num2str(s_lag.AIClin)])
disp(['         BIC =' num2str(s_lag.BIClin)])


%---------------------------------------------------------------------
%            generate Minnesota prior
%---------------------------------------------------------------------
if s_prior.flagOLS==0
    [s_prior.prior_mean,s_prior.prior_var]=genMNprior(s_desc.Nvar,s_desc.VARlags,VARLIN.OmegaVARLIN,s_prior.lambda_prior,s_prior.theta_prior);
    
    % note that the MN prior does not provide any information for the
    % covariance matrix of residuals. Hence, impose a flat prior "centered"
    % around the OLS estimate
    
    % covariance matrix of residuals in the VAR
    % compute the covariance matrix for the estimate of OmegaNL
    Dn=dupmat(s_desc.Nvar);
    Dplus=inv(Dn'*Dn)*Dn';
    OmegaVARLIN_var=2*Dplus*kron(OmegaVARLIN,OmegaVARLIN)*Dplus'/s_desc.T;
    
    s_prior.OmegaVAR=1e6*ones(size(diag(OmegaVARLIN_var)));
    s_prior.Omega=OmegaVARLIN;
else
end



% store the initial values
Omega0E_1=VARLIN.OmegaVARLIN;
Omega0chol=chol(VARLIN.OmegaVARLIN)';
param01cov=[vech(Omega0chol)];
s_desc.Omega_length_cov = length(vech(Omega0chol));

%---------------------------------------------------------------------
%  Second step: ESTIMATE USING MCMC ROUTINES
%  Run the full linear model
% 
%---------------------------------------------------------------------

[Amat,beta0mat,acceptrate,valJ]=...
    optimizeMCMC_struct_linear(param01cov,s_data,s_desc,s_prior,s_exo,s_MCMC);

%---------------------------------------------------------------------
%           Plot impulse responses with confidence intervals
%---------------------------------------------------------------------

% mean values of parameters in the estimated VAR 
beta0E=unvec(mean(beta0mat(s_MCMC.dropfirstT:end-1,:))',s_desc.Nvar,s_desc.Nvar*s_desc.VARlags);

[Omega0E]=vec2matScov_linear(mean(Amat(s_MCMC.dropfirstT:end-1,:))',s_desc.Omega_length_cov);

Omega_length=length(vech(Omega0E));

Omega0mat=Amat(s_MCMC.dropfirstT:end-1,1:Omega_length);
beta0matF=beta0mat(s_MCMC.dropfirstT:end-1,:);


strA=['code_fullsample_linear.mat'];
save(strA)


