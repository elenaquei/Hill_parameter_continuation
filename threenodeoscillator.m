function hillEvaluation = threenodeoscillator(x, lambda, n)
%THREENODEOSCILLATOR - Evaluation function for 3 dimensional toy Hill model for a oscillatory network.  
%
%   THREENODEOSCILLATOR - Evaluates the vector field for the system of ODEs:
%               x1' = -x1 + H-(x3)  corresponds to x3 ---| x1
%               x2' = -x2 + H-(x1)  corresponds to x1 ---| x2
%               x3' = -x3 + H-(x2)  corresponds to x2 ---| x3
%
%   Inputs:
%       x - a 3-dimensional vector of state variables 
%       lambda - a 3-by-3 array of parameters. Each row corresponds to an edge of the network model as follows:
%                lambda(j,:) = [theta_i, l_i, delta_i]
%       n - a 3-dimensional vector of Hill coefficients
%
%   Outputs:
%       f(x, lambda, n) - a 3-dimensional vector
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 09-Jan-2020; Last revision: 09-Jan-2020

f1 = -x(1) + lambda(1,2) + (lambda(1,3).*lambda(1,1).^n(1))./(lambda(1,1).^n(1) + x(3).^n(1));
f2 = -x(2) + lambda(2,2) + (lambda(2,3).*lambda(2,1).^n(2))./(lambda(2,1).^n(2) + x(1).^n(2));
f3 = -x(3) + lambda(3,2) + (lambda(3,3).*lambda(3,1).^n(1))./(lambda(3,1).^n(3) + x(2).^n(3));
hillEvaluation = [f1; f2; f3];
end % end threenodeoscillator

% Revision History:
%{
Author: Elena Queirolo
 - Changed to repressilator; 14-Jan-2020
%}
