% this function computes IRFs with a shock and without a shock
% for a non-linear VAR where RHS variables include polynomials of lagged variables
% the shocks is 1% shock to government spending

% Alternative way to compute standard errors:
% use duplication matrix to compute the variance for the covariance matrix
% of residuals rather than draw the moments from the MCMC chain for the
% covariance matrix of residuals

function [irf_Dimp,irf_CImean,CIup,CIlow,IRFse]=irf_VAR_MCMC_struct_alt(betaNL,OmegaNL,Omegamat,betamat,s_desc)

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

Nvar=s_desc.Nvar;
VARlags=s_desc.VARlags;
T=s_desc.T;
irf_hor=s_desc.irf_hor;

% number of iterations to compute confidence bands and st.error for the
% impulse response function
CIiters=1000;
display ("irf_var_mcmc_struct_alt")

% percentiles for confidence intervals
up_percentile  = 84;
low_percentile = 16;


% companion matrix
A_NL=betaNL;
A_NL=[A_NL; eye(Nvar*(VARlags-1)) zeros(Nvar*(VARlags-1),length(betaNL)-Nvar*(VARlags-1))];

%=========================================================================
%                   Point estimate of the IRF
%=========================================================================
% use cholesky decomposition to compute the shock
CholOmegaNL=chol(OmegaNL)';
shockvec=CholOmegaNL(:,s_desc.KIX)/CholOmegaNL(s_desc.KIX,s_desc.KIX);
% shockvec
% impulse response with a shock
irf_lin=[shockvec; zeros(Nvar*(VARlags-1),1)];
irf_imp=irf_lin;

for i=2:irf_hor
    irf_lin=A_NL*irf_lin;
    irf_imp=[irf_imp irf_lin];
end

% compute the point estimate of the impulse response
irf_Dimp=irf_imp(1:Nvar,:);

%=========================================================================
%                   compute confidence intervals for the IRF
%                   using MCMC draws
%=========================================================================
% compute confidence intervals for the shock
irf_Dimp0mat=[];

betaNL_vec=vec(betaNL);
rand('seed',1234567)

%HERE IS THE DIFFERENCE
OmegaVARLIN=OmegaNL;
Dn=dupmat(s_desc.Nvar);
Dplus=inv(Dn'*Dn)*Dn';
OmegaVARLIN_var=2*Dplus*kron(OmegaNL,OmegaNL)*Dplus'/s_desc.T; %page 22, variance covariance matrix


for i1=1:CIiters
    % draw of parameter values in the VAR
    drawposition=max(1,min(round(rand(1,1)*length(betamat)),length(betamat)));
    % draw the dynamic responses
    paramdraw0=betamat(drawposition,:)';
    paramdraw=unvec(paramdraw0,Nvar,Nvar*VARlags);

    % draw the covariance matrix of residuals (contemporanous relationship)
    OmegaNLdraw2=vech(OmegaVARLIN)+chol(OmegaVARLIN_var)'*randn(size(vech(OmegaVARLIN)));
    OmegaNLdraw3=unvech(OmegaNLdraw2);
    
    CholOmegaNL0=chol(OmegaNLdraw3)';
    shockvec0=CholOmegaNL0(:,s_desc.KIX)/CholOmegaNL0(s_desc.KIX,s_desc.KIX);

    A_NL=paramdraw;
    A_NL=[A_NL; eye(Nvar*(VARlags-1)) zeros(Nvar*(VARlags-1),length(betaNL)-Nvar*(VARlags-1))];


    irf_lin=[shockvec0; zeros(Nvar*(VARlags-1),1)];
    irf_imp=irf_lin;

    for i=2:irf_hor
        irf_lin=A_NL*irf_lin;
        irf_imp=[irf_imp irf_lin];
    end

    irf_Dimp0mat(i1,:,:)=irf_imp;

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


