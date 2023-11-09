
% convert vector of parameters into approapriate matrices

function [Omega0,Omega1,theta]=vec2matScov(param0,Omega_length)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                   the curvature in the transition function
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
theta=param0(2*Omega_length+1,1);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                   covariance matrix of the error terms
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% pick vectorized matrices from the vector of parameters
Omega0_vec=param0(1:Omega_length,1);
Omega1_vec=param0(Omega_length+1:2*Omega_length,1);

% unvectorize matrices
Omega0M=unvech(Omega0_vec);
Omega1M=unvech(Omega1_vec);

% since the initially vectorized matrix is lower diagonal remove the
% elements above the diagonal
Omega0M=Omega0M-triu(Omega0M,1);
Omega1M=Omega1M-triu(Omega1M,1);

% ensure the covariance matrices are positive definite (this step is
% effectively using the nature of the Choleski decomposition)
Omega0=Omega0M*Omega0M';
Omega1=Omega1M*Omega1M';

end

