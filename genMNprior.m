
% generate Minnesota prior
function [prior_mean,prior_var]=genMNprior(Nvar,VARlags,OmegaVARLIN,lambda_prior,theta_prior)

% Inputs:
% lambda_prior = controls cross-variable responses (smaller values mean stronger prior 
% theta_prior = controls own response of variables
% Nvar = # of variables in the VAR
% VARlags = # of lags in the VAR
% OmegaVARLIN = covariance matrix of residuals to make restrictions scale
% invariant (typically from "first stage" unrestricted VAR)

% Output: matrix of tighness for restrictions

prior_mean=zeros(Nvar,Nvar*VARlags);
for i=1:Nvar
    prior_mean(i,i)=1;
end    

prior_var=zeros(Nvar,Nvar*VARlags);
for L=1:VARlags
    for i=1:Nvar
        for j=1:Nvar
            if i==j
                prior_var(i,j+(L-1)*Nvar)=lambda_prior^2/(L^2);
            else
                prior_var(i,j+(L-1)*Nvar)=(lambda_prior*theta_prior)^2/(L^2)*(OmegaVARLIN(i,i)/OmegaVARLIN(j,j));                
            end
        end
    end
end   



 
