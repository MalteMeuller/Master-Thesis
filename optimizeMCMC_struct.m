% this function is the MCMC optimizer

function [Amat,beta0mat,beta1mat,acceptrate,valJ,...
    beta_const_mat, beta_other_mat]=...
    optimizeMCMC_struct(param01cov,s_data,s_desc,s_prior,s_exo,s_MCMC)



% Omega_length_cov,T,Nvar,theta0,...
%     scale_penalty,prior_mean,prior_var,FX,Nforecs,VARlags,MCMCn,...
%     trendON,trend_regspec,inter_regspec)

% Inputs:
% param01cov = initial value of parameters
% Y0 = LHS of VAR equations
% X  = RHS of VAR equations
% Z0 = regime switching variable
% Omega_length_cov = length of the
% T  = sample size
% Nvar = number of variables in the VAR
% theta0 = guess about curvature in the regime switching probability
% scale_penalty = scale of penalty for the departure of estimated curvature in the transition probability from the guess
% prior_mean = mean of the Minnesota prior
% prior_var  = tightness of the Minnesota prior
% FX         = vector of forecasts used to augment the VAR
% Nforecs    = number of variables in the vector of forecasts
% VARlags    = # of lags in the VAR
% MCMCn      = # of MCMC draws

% Outputs:
% Amat = matrix of accepted MCMC draws
% beta0mat = Amat converted into slopes of the VAR in regime 0
% beta1mat = Amat converted into slopes of the VAR in regime 1
% acceptrate = rate at which candidate values are accepted in the MCMC
% valJ = value of the objective function for the
% beta_const_mat = matrix of intercepts
% beta_other_mat = matrix of coefs on other controls, SEARCH HERE

% read arguements from imported structures
Omega_length_cov=s_desc.Omega_length_cov;
T=s_desc.T;
Nvar=s_desc.Nvar;

theta0=s_prior.theta0;
scale_penalty=s_prior.scale_penalty;
prior_mean=s_prior.prior_mean;
prior_var=s_prior.prior_var;

FX=s_exo.FX;
Nforecs=s_exo.Nforecs;

VARlags=s_desc.VARlags;
MCMCn=s_MCMC.MCMCn;

Y0=s_data.Y0;
X=s_data.X;
Z0=s_data.Z0;

% initial value and the fit
A=param01cov;
J1=max_LSTVAR_mleS_struct(A,s_data,s_desc,s_prior,s_exo,s_MCMC); %går till ny funktion

% matrices to store draws and the fit
Amat=zeros(MCMCn,length(A));
beta0mat=zeros(MCMCn,Nvar^2*VARlags);
beta1mat=zeros(MCMCn,Nvar^2*VARlags);
acceptrate=zeros(MCMCn,1);
valJ=zeros(MCMCn,1);

% size of the draws for the MCMC procedures, myabe something here
sigmaMH=0.0005*ones(size(param01cov)); %ist för 0.0005
if scale_penalty==0
    sigmaMH(end)=0.05;
else
    sigmaMH(end)=0; %sigmamh last index
end


scaleSIGMA=1;

[beta0E,beta1E,Omega0E,Omega1E,thetaE,gammaE,residsM,beta_const, beta_other,logL]=...
    vec2matSF_struct(A,s_data,s_desc,s_prior,s_exo,s_MCMC);

tic
i=1;
randn('seed',1234567);	% set the seed for random generator
[A_rows, A_cols] = size(A); %ändrade till detta för att få en output i row och columner. %kan det gå fel här?

% Create shocks0 with the same size as A
shocks0 = randn(A_rows, A_cols, MCMCn * 5); %två input istället för en vector input (input must be svalar)
counterB=1;
counterX=0;
while i<MCMCn+1
    
    A1=A+sigmaMH.*shocks0(:,counterB);  % draw new stock to paramter value ###here comes shocks that changes things
    counterB=counterB+1;
    %     A1=A+sigmaMH.*randn(size(A));            % draw new stock to paramter value
    
    [Omega0_candidate,Omega1_candidate,theta_candidate]=vec2matScov(A1,Omega_length_cov);
    
    % check if eigenvalues of Omega0 and Omega1 are greater than zero
    % if not, discard the candidate vector of parameters
    if min(eig(Omega0_candidate))>0 & min(eig(Omega1_candidate))>0
        
        % objective value for a new draw candidate value
        J2=max_LSTVAR_mleS_struct(A1,s_data,s_desc,s_prior,s_exo,s_MCMC);
        
        % stochastically accept/reject candidate values
        v0=exp(J1-J2);                      % compute the exp of ration of loss functions
        v0 = v0 ;
        v1=min(v0,1);                       % determine accept probability,, THIS IS STEP 2
        U=rand(1,1);                        % draw rv from uniform to decide whether accept
        if U<v1                             % accept if the U is less than acceptance rate
            A=A1;
            acceptrate(i,1)=1;
            valJ(i,1)= J2;
            J1=J2;
            [beta0E,beta1E,Omega0E,Omega1E,...
                thetaE,gammaE,residsM,...
                beta_const, beta_other,...
                logL]=...
                vec2matSF_struct(A1,s_data,s_desc,s_prior,s_exo,s_MCMC);
        else
            acceptrate(i,1)=0;
            valJ(i)= J1;
        end
        
        Amat(i,:)= A';
        if s_desc.VARlags>0
            beta0mat(i,:)=vec(beta0E)';
            beta1mat(i,:)=vec(beta1E)';
            
            beta_const_mat(i,:)=vec(beta_const)';
            beta_other_mat(i,:)=vec(beta_other)';
            
        end
        
        i=i+1;
        counterX=0;
        % adjust size of the shocks on the fly.
        if i<400000
            if mod(i,500)==0
                
               
                if mean(acceptrate(max(1,i-500):i))>0.35
                    sigmaMH=1.1*sigmaMH;
                    scaleSIGMA=1.1*scaleSIGMA;
                end
                if mean(acceptrate(max(1,i-500):i))<0.25
                    sigmaMH=0.9*sigmaMH;
                    scaleSIGMA=0.9*scaleSIGMA;
                end
            end
        end
        
        
        if i==250
            disp(i)
        end
        if mod(i,10000)==0
            disp(i)
        end
    else
        counterX=counterX+1;
    end
    
    if counterX>200
        disp('too many draws...  Go back to previous successful draw')
        A=Amat(i-1,:)';
        J1=valJ(i-1);
        counterX=0;
    end

   
end

