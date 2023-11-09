
function obj_logL=max_LSTVAR_mleS_struct(param0,s_data,s_desc,s_prior,s_exo,s_MCMC)

% Y,X,Z0,Z1,Omega_length,T,...
%     Nvar,theta0,scale_penalty,prior_mean,prior_var,FX,Nforecs,...
%     trendON,trend_regspec,inter_regspec)

% Inputs
% Y = vector dependent variables
% X = vector of lags 
% Z0 = regime shifter for contemporanous response (matrices OMEGA) and for
%      dynamic response (matrices PI)
% Omega_length = length of the vech(cov matrix of residuals)
% T = sample size
% Nvar = number of variables in the VAR
% Penalty terms and priors: 
%    theta0 = guess for the curvature in the transition probability
%    scale_penalty = size of the penalty for departures of theta from theta0 
%    prior_mean = prior for slopes in the VAR
%    prior_var  = tightness of the prior for VAR slopes,
%    FX,Nforecs

Y=s_data.Y0;
X=s_data.X;
Z0=s_data.Z0;
Omega_length=s_desc.Omega_length_cov;
T=s_desc.T;
Nvar=s_desc.Nvar;
theta0=s_prior.theta0;
scale_penalty=s_prior.scale_penalty;
prior_mean=s_prior.prior_mean;
prior_var=s_prior.prior_var;
FX=s_exo.FX;
Nforecs=s_exo.Nforecs;

% convert vector of parameters into approapriate matrices for eatch regime
[Omega0,Omega1,theta]=vec2matScov(param0,Omega_length);


% check if eigenvalues of Omega0 and Omega1 are greater than zero
% if not exit the function
if  theta<0
    if min(eig(Omega0))>0 & min(eig(Omega1))>0
        % impose the restriction that gamma = theta (that is the curvature
        % of the transition probability is the same for contemporanous 
        % response (matrices OMEGA) and for dynamic response (matrices PI)
        % and it's governed by the same switching variable
        gamma=theta;
        
        [obj_logL,beta0,beta1,beta_const, beta_other]=LSTVAR_mleS_struct(Omega0,Omega1,gamma,theta,s_data,s_desc,s_exo);


        % restriction of the speed of change in the regime
        penaltry_term0=(theta-theta0)^2*scale_penalty;
        
        % bayesian restriction: VAR slopes
        if s_desc.VARlags>0 & s_prior.flagON==1
        penaltry_term1=0.5*sum(sum((beta1-prior_mean).^2./prior_var))*s_prior.flagON;
        penaltry_term2=0.5*sum(sum((beta0-prior_mean).^2./prior_var))*s_prior.flagON;
        else
            penaltry_term1=0;
            penaltry_term2=0;
        end
        
        % bayesian restriction: VAR covariance matrix of residuals
        penaltry_term3=0.5*sum(sum( vech(Omega1-s_prior.Omega).^2./s_prior.OmegaVAR) )*s_prior.flagON*s_prior.flagOmegaON;
        penaltry_term4=0.5*sum(sum( vech(Omega0-s_prior.Omega).^2./s_prior.OmegaVAR) )*s_prior.flagON*s_prior.flagOmegaON;
        
        obj_logL=-obj_logL + penaltry_term0 + ...
            penaltry_term1 + penaltry_term2 + ...
            penaltry_term3 + penaltry_term4; 
        
    else
        obj_logL=1e12;
    end
else
    obj_logL=1e12;
end

end




