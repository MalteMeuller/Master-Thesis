% this function CIiterscomputes IRFs with a shock and without a shock
% for a non-linear VAR where RHS variables include polynomials of lagged variables
% the shocks is 1% shock to government spending

function [irf_Dimp,irf_CImean,CIup,CIlow,IRFse,...
    cum_irf_Dimp,cum_irf_CImean,cum_CIup,cum_CIlow,cum_IRFse,...
    max_irf_Dimp,max_irf_CImean,max_CIup,max_CIlow,max_IRFse]=...
    irf_VAR_MCMC_struct_alt_dynamic_feedback3(betaNL_0,OmegaNL_0,...
    betaNL_1,OmegaNL_1,...
    Omegamat_0,betamat_0,Omegamat_1,betamat_1,s_desc,s_prior,...
    Y, Z, X,pos000)

% Outputs:
%       irf_imp   = impulse response with a shock
%       irf_CImean = mean value of the IRF across draws
%       CIup      = upper bound on the confidence interval
%       CIlow     = lower bound on the confidence interval
%       IRFse     = standard error for the confidence interval
% Inputs:
%       betaNL      = matrix of estimated VAR coefficients
%       betaNL_var  = variance of estimated VAR coefficients
%       Nvar        = number of variables in the VAR
%       VARlags     = number of lags in the VAR
%       OmegaNL     = covariance matrix of estimated error terms
%       irf_hor     = horizon of the VAR
%       T           = sample size
%       Y_series    = output series
%       Z_series    = swtiching variable based on the growth rate of output
%       pos000      = position for which IRF is computed
% Defintion of regimes
%       0 = dynamics in expansions
%       1 = dynamics in recessions

Nvar=s_desc.Nvar;
VARlags=s_desc.VARlags;
T=s_desc.T;
irf_hor=s_desc.irf_hor;

% position for which to compute the IRF
Z_run0    = Z(pos000-1)-0.8; %manuel input, 

% number of iterations to compute confidence bands and st.error for the
% impulse response function
randn('seed',1234567)
CIiters=1000;

% percentiles for confidence intervals
up_percentile  = 95;
low_percentile = 5;


% companion matrix: regime 0
A_NL_0=betaNL_0;
A_NL_0=[A_NL_0; eye(Nvar*(VARlags-1)) zeros(Nvar*(VARlags-1),length(betaNL_0)-Nvar*(VARlags-1))];

% companion matrix: regime 1
A_NL_1=betaNL_1;
A_NL_1=[A_NL_1; eye(Nvar*(VARlags-1)) zeros(Nvar*(VARlags-1),length(betaNL_0)-Nvar*(VARlags-1))];

%=========================================================================
%                   compute the impulse response
%=========================================================================
f
F_Z_t0=exp(s_prior.theta0*Z_run)./(1+exp(s_prior.theta0*Z_run));
Z_run=Z_run0;
Z_ser=[Z_run];
% impulse response with a shock
% starting value of Z in the impulse response: s_desc.startZir
% compute the covariance matrix of the shock
cov_shock=OmegaNL_0*(1-F_Z_t0)+OmegaNL_1*(F_Z_t0);
cov_shock_chol=chol(cov_shock)';

% compute the size of the G shock = 1%
shockvec_wgt=cov_shock_chol(:,s_desc.Gpos)/cov_shock_chol(s_desc.Gpos,s_desc.Gpos);
shockvec_wgt=shockvec_wgt/100;

% stack together shocks for the VAR
shock_lin=[shockvec_wgt; zeros(Nvar*(VARlags-1),1)];

% compute the lag polynomial which governs the VAR
A_NL_wgt=A_NL_0*(1-F_Z_t0)+A_NL_1*F_Z_t0;

% compute the dynamics at time t
irf_imp=[shock_lin];

% update the state (take into account that dY is serially correlated
% and therefore Z may stay above average for a little longer). 
Z_run=((6/7)+1/7*(0.3))*Z_run+((irf_imp(3,1)-irf_imp(6,1))*100)/7;
Z_ser=[Z_ser; Z_run];
for i=2:irf_hor
    % impulse response with a shock
    % starting value of Z in the impulse response: s_desc.startZirf
    F_Z_t0=exp(s_prior.theta0*Z_run)./(1+exp(s_prior.theta0*Z_run));
    
    % compute the lag polynomial which governs the VAR
    A_NL_wgt=A_NL_0*(1-F_Z_t0)+A_NL_1*F_Z_t0;
    
    % compute the dynamics at time t
    irf_imp=[irf_imp A_NL_wgt*irf_imp(:,i-1)];
    
    % update the state
    Z_run=((6/7)+1/7*(0.3))*Z_run+((irf_imp(3,i)-irf_imp(6,i))*100)/7;
    Z_ser=[Z_ser; Z_run];
end

irf_Dimp=irf_imp;
cum_irf_Dimp=sum(irf_imp(3,:))/sum(irf_imp(1,:));
max_irf_Dimp=max(irf_imp(3,:));

%=========================================================================
%                   compute confidence intervals for the IRF
%                   using MCMC draws
%=========================================================================
% compute confidence intervals for the shock
irf_Dimp0mat=[];

betaNL_vec=vec(betaNL_0);
rand('seed',1234567)

for i1=1:CIiters
    % draw of parameter values in the VAR
    drawposition=max(1,min(round(rand(1,1)*length(betamat_0)),length(betamat_0)));
    
    % draw the covariance matrix of residuals (contemporanous relationship)
    % regime 0
    OmegaNLdraw0=Omegamat_0(drawposition,:)';
    OmegaNLdraw1=unvech(OmegaNLdraw0);
    OmegaNLdraw2=OmegaNLdraw1-triu(OmegaNLdraw1,1);
    OmegaNLdraw3=OmegaNLdraw2*OmegaNLdraw2';
    
    OmegaNL_0=OmegaNLdraw3;
    
    % regime 1
    OmegaNLdraw0=Omegamat_1(drawposition,:)';
    OmegaNLdraw1=unvech(OmegaNLdraw0);
    OmegaNLdraw2=OmegaNLdraw1-triu(OmegaNLdraw1,1);
    OmegaNLdraw3=OmegaNLdraw2*OmegaNLdraw2';
    OmegaNL_1=OmegaNLdraw3;
    
    % draw the dynamic responses
    % regime 0
    paramdraw0=betamat_0(drawposition,:)';
    paramdraw=unvec(paramdraw0,Nvar,Nvar*VARlags);
    A_NL_0=paramdraw;
    A_NL_0=[A_NL_0; eye(Nvar*(VARlags-1)) zeros(Nvar*(VARlags-1),length(betaNL_0)-Nvar*(VARlags-1))];
    
    % regime 1
    paramdraw0=betamat_1(drawposition,:)';
    paramdraw=unvec(paramdraw0,Nvar,Nvar*VARlags);
    A_NL_1=paramdraw;
    A_NL_1=[A_NL_1; eye(Nvar*(VARlags-1)) zeros(Nvar*(VARlags-1),length(betaNL_1)-Nvar*(VARlags-1))];
    
    % start computing the impulse response
    Z_run=Z_run0;
    Z_ser=[Z_run];
    % impulse response with a shock
    % starting value of Z in the impulse response: s_desc.startZirf
    F_Z_t0=exp(s_prior.theta0*Z_run)./(1+exp(s_prior.theta0*Z_run));
    
    % compute the covariance matrix of the shock
    cov_shock=OmegaNL_0*(1-F_Z_t0)+OmegaNL_1*(F_Z_t0);
    cov_shock_chol=chol(cov_shock)';
    
    % compute the size of the G shock = 1%
    shockvec_wgt=cov_shock_chol(:,s_desc.Gpos)/cov_shock_chol(s_desc.Gpos,s_desc.Gpos);
    shockvec_wgt=shockvec_wgt/100;
    
    % stack together shocks for the VAR
    shock_lin=[shockvec_wgt; zeros(Nvar*(VARlags-1),1)];
    
    % compute the lag polynomial which governs the VAR
    A_NL_wgt=A_NL_0*(1-F_Z_t0)+A_NL_1*F_Z_t0;
    
    % compute the dynamics at time t
    irf_imp=[shock_lin];
    
    % update the state
    Z_run=((6/7)+1/7*(0.3))*Z_run+((irf_imp(3,1)-irf_imp(6,1))*100)/7;
    Z_ser=[Z_ser; Z_run];
    for i=2:irf_hor
        % impulse response with a shock
        % starting value of Z in the impulse response: s_desc.startZirf
        F_Z_t0=exp(s_prior.theta0*Z_run)./(1+exp(s_prior.theta0*Z_run));
        
        % compute the lag polynomial which governs the VAR
        A_NL_wgt=A_NL_0*(1-F_Z_t0)+A_NL_1*F_Z_t0;
        
        % compute the dynamics at time t
        irf_imp=[irf_imp A_NL_wgt*irf_imp(:,i-1)];
        
        % update the state
        Z_run=((6/7)+1/7*(0.3))*Z_run+((irf_imp(3,i)-irf_imp(6,i))*100)/7;
        Z_ser=[Z_ser; Z_run];
    end
    cum_irf_imp=sum(irf_imp(3,:))/sum(irf_imp(1,:));

    irf_Dimp0mat(i1,:,:)=irf_imp;
    cum_irf_Dimp0mat(i1)=cum_irf_imp;
    max_irf_Dimp0mat(i1)=max(irf_imp(3,:));
    
end

% compute confidence intervals and standard errors
CIup=[];
CIlow=[];
IRFse=[];
irf_CImean=[];
for i=1:Nvar
    CIup   = [CIup; prctile(squeeze(irf_Dimp0mat(:,i,:)),up_percentile,1)];
    CIlow  = [CIlow; prctile(squeeze(irf_Dimp0mat(:,i,:)),low_percentile,1)];
    IRFse  = [IRFse; std(squeeze(irf_Dimp0mat(:,i,:)),1,1)];
    irf_CImean = [irf_CImean; median(squeeze(irf_Dimp0mat(:,i,:)),1)];
end

cum_CIup=prctile(cum_irf_Dimp0mat,up_percentile,1);
cum_CIlow=prctile(cum_irf_Dimp0mat,low_percentile,1);
cum_IRFse=std(cum_irf_Dimp0mat);
cum_irf_CImean=median(cum_irf_Dimp0mat);

max_CIup=prctile(max_irf_Dimp0mat,up_percentile,1);
max_CIlow=prctile(max_irf_Dimp0mat,low_percentile,1);
max_IRFse=std(max_irf_Dimp0mat);
max_irf_CImean=median(max_irf_Dimp0mat);



end
