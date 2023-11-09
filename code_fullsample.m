% this code computes the MLE estimate with time-varying volatility and
% slopes/interpect in the VAR.


%Estimation for figure 2 (baseline specification)
%%
function []=code_fullsample(theta000,VARlags000)
warning off all


tic
randn('seed',1234567)	% set the seed for random generator

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           		PRELIMINARY STEPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---------------------------------------------------------------------
%           Sample
%---------------------------------------------------------------------
s_desc.startS  = 1; %start
s_desc.endS    = 113; %end
s_desc.sample0 = s_desc.startS:1:s_desc.endS;

%---------------------------------------------------------------------
%            Conventions
%---------------------------------------------------------------------



%---------------------------------------------------------------------
%            Parameters governing the model
%---------------------------------------------------------------------
% VAR parameters
s_desc.VARlags=VARlags000;      % number of lags in the VAR
s_desc.Nvar=4;                  % number of variables in the VAR
s_desc.irf_hor=20;              % number of periods in the impulse response
percentile = 60;


% Regime specific intercept and trend
s_exo.trendON=0;      % 1 = include a time trend in the VAR;
% 0 = do not include time trend in the VAR

s_exo.trend_regspec=0;  % 1 = allow regime specific time trend
% 0 = do not allow regime specific time trend

s_exo.inter_regspec=0;  % 1 = allow regime specific intercept
% 0 = do not allow regime specific intercept

%
%---------------------------------------------------------------------
%                   Bayesian VAR parameters
%currently turned of
%---------------------------------------------------------------------

s_prior.flagON=0;   % 1 = use a prior in estimation of VAR slopes
% 0 = do not use a prior in estimation of VAR slopes

s_prior.flagOmegaON=0;
% 1 = use a prior in estimation in cov matrix of VAR residuals
% 0 = do not use a prior in cov matrix of VAR residuals

s_prior.flagOLS=0;
% 1 = use OLS estimates of the VAR as a prior
% 0 = use Minnesota prior for the VAR (note that MN
% prior does not have a prior for the covariance matrix
% of the residuals)

% penalty for theta which governs the curvature in the transition
% probability function
s_prior.theta0              = theta000;    % the slope in the switching equation
s_prior.scale_penalty       = 1e6;
s_prior.scale_penalty_beta  = 0;

% Parameters of the Minnesota prior: random walk for each variable in the VAR
s_prior.lambda_prior=1;   % controls cross-variable responses (smaller values mean stronger prior
s_prior.theta_prior=0.5;  % controls own response of variables
%%
%---------------------------------------------------------------------
%                       MCMC step, %ändra här för långsammare MCMC
%---------------------------------------------------------------------

s_MCMC.runMCMCflag=1; % = 1 use MCMC in estimation

s_MCMC.MCMCn=100000;                         % number of MCMC draws, ädndrade här, änra tillbaka
s_MCMC.dropfirstT=min(round(s_MCMC.MCMCn*0.3),20000); % burn-in period for MCMC %ädnrade här, ändra tillbaka

%---------------------------------------------------------------------
%           			Import data
%---------------------------------------------------------------------
%download variables to use
data0= readtable("C:\Users\malte\OneDrive - Handelshögskolan i Stockholm\Thesis\Data\Data_Quarterly\Clean\df_qoqld.csv", 'ReadVariableNames',true);

% switching series: output in data_SV
data_SV = data0(1:s_desc.endS,"ld_kpif_lag2");

%exogenous varaibles
xo = data0(1:s_desc.endS, {'ld_kixcpi', 'ld_kixgdp', 'kixpol_hp'});
xo = table2array(xo);

time=1:113

%---------------------------------------------------------------------
%                    create series for VAR
%---------------------------------------------------------------------

%this is where our data should be created. 
Y = data0(1:s_desc.endS, {'ld_gdp',  'polrate_hp', 'ld_kix', 'ld_kpif' });
Y = table2array(Y);
%change here the position of KIX
s_desc.KIX=3 ;
s_desc.KPIF=4;

%switching variable
Z=data_SV; % growth rate of inf
Z = table2array(Z);

% adjust the switching variable to have suffiently long recession periods
Th = prctile(Z, percentile); %we do a 60/40 split


%moving the threshold
Z0=Z-Th;



% vector of exogenous variables
s_exo.FX=xo;
s_exo.FX=s_exo.FX(s_desc.VARlags+1:end,:);

% number of variables in the vector of forecasts
s_exo.Nforecs=3;

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
%sign takes 1 if sign positive, -1 if sign negative
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

%aik and bic for residuals
s_lag.AICnonlin=log(det(Omega))+2/length(resids)*length(beta);
s_lag.BICnonlin=log(det(Omega))+log(length(resids))/length(resids)*length(beta);

%---------------------------------------------------------------------
%            linear VAR to generate residuals
%---------------------------------------------------------------------
Xm1=[X ones(length(Y0),1) s_exo.FX]; % exogenous here also
% estimate parameters
betaVARLIN=inv(Xm1'*Xm1)*(Xm1'*Y0);
residsVARLIN=Y0-Xm1*betaVARLIN;
% compute the covariance matrix of residuals
OmegaVARLIN=cov(residsVARLIN);
betaVARLIN=betaVARLIN(1:end-1-s_exo.Nforecs,:); %11 parametrar, 4*2 + 2 exo + 1 constant. 

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
    disp('Error: use MN prior!')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           NON-Linear VAR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---------------------------------------------------------------------
%  First step:
%  estimate the process for the distribution of the error term
%  Omega0 in low inflation and Omega1 in high inflation
%---------------------------------------------------------------------

Omega0=Omega;

% take a small departure for one of the regimes to make sure that starting
% values are slightly different
D=randn(length(Omega));
Omega1=Omega+min(eig(Omega))*D*D';

% Apply Cholesky decomposition since we need to look only for a lower
% triangular matrix
Omega0chol=chol(Omega0)';
Omega1chol=chol(Omega1)';
s_desc.Omega_length_cov = length(vech(Omega0chol));

% initial parameter values, to fit the algo
param0cov=[vech(Omega0chol);  vech(Omega1chol);  s_prior.theta0];

itermax=1000;
disp('fitting first step...')
% I need to pass the prior here
%MAXIMUM LIKELYHOOD ESTIMATIONS 
options=optimset('MaxFunEvals',itermax,'MaxIter',itermax,'Display','off','TolFun',1e-6,'TolX',1e-6);
[param01cov,mleLcov,flagquasi]=fminunc(@max_LSTVAR_mleScov,param0cov,options,resids,Z0,s_desc.Omega_length_cov,s_desc.T,s_desc.Nvar,s_prior.theta0,s_prior.scale_penalty);
[param01cov,mleLcov,flagquasi]=fminsearch(@max_LSTVAR_mleScov,param01cov,options,resids,Z0,s_desc.Omega_length_cov,s_desc.T,s_desc.Nvar,s_prior.theta0,s_prior.scale_penalty);
disp(['logL = ' num2str(mleLcov)])

% convert vectorized matrices with parameters into appropriate matrices
[Omega0E,Omega1E,thetaE]=vec2matScov(param01cov,s_desc.Omega_length_cov);

% Omega0E = contemporanoeous covariance in the expansion
% Omega1E = contemporanoeous covariance in the recession

% store the initial values
Omega0E_1=Omega0E;
Omega1E_1=Omega1E;
Omega0chol=chol(Omega0E)';
Omega1chol=chol(Omega1E)';

%this creates the covaraiance matrix from the Cholesky
Omega1E=Omega1chol*Omega1chol';
Omega1E_1=Omega1E;

param01cov=[vech(Omega0chol);  vech(Omega1chol);  thetaE];


%---------------------------------------------------------------------
%  Second step: ESTIMATE USING MCMC ROUTINES
%  Run the full non-linear model
%  - VAR for each regime
%  - covariance matrix of residuals can vary across regimes.
%            non-linear VAR to generate residuals
%---------------------------------------------------------------------

%gives output 250 by some reason
[Amat,beta0mat,beta1mat,acceptrate,valJ]=...
    optimizeMCMC_struct(param01cov,s_data,s_desc,s_prior,s_exo,s_MCMC);

%---------------------------------------------------------------------
%           Plot impulse responses with confidence intervals
%---------------------------------------------------------------------

% mean values of parameters in the estimated VAR with two regimes
beta1E=unvec(mean(beta1mat(s_MCMC.dropfirstT:end-1,:))',s_desc.Nvar,s_desc.Nvar*s_desc.VARlags); % parameters in high inflation
beta0E=unvec(mean(beta0mat(s_MCMC.dropfirstT:end-1,:))',s_desc.Nvar,s_desc.Nvar*s_desc.VARlags); %parameters in low inflation

[Omega0E,Omega1E,thetaE]=vec2matScov(mean(Amat(s_MCMC.dropfirstT:end-1,:))',s_desc.Omega_length_cov);

Omega_length=length(vech(Omega0E));

Omega1mat=Amat(s_MCMC.dropfirstT:end-1,Omega_length+1:2*Omega_length);
Omega0mat=Amat(s_MCMC.dropfirstT:end-1,1:Omega_length);
beta1matF=beta1mat(s_MCMC.dropfirstT:end-1,:);
beta0matF=beta0mat(s_MCMC.dropfirstT:end-1,:);


strA=['code_fullsample_KPIF.mat'];
save(strA)

%high and low regime


toc
end 