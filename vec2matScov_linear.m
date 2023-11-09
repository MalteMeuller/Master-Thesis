
% convert vector of parameters into approapriate matrices

function [Omega0]=vec2matScov_linear(param0,Omega_length)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                   covariance matrix of the error terms
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% pick vectorized matrices from the vector of parameters
Omega0_vec=param0(1:Omega_length,1);

% unvectorize matrices
Omega0M=unvech(Omega0_vec);

% since the initially vectorized matrix is lower diagonal remove the
% elements above the diagonal
Omega0M=Omega0M-triu(Omega0M,1);

% ensure the covariance matrices are positive definite (this step is
% effectively using the nature of the Choleski decomposition)
Omega0=Omega0M*Omega0M';

end

