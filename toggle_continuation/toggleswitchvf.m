function hillPolynomialEvaluation = toggleswitchvf(x, lambda, n)
%TOGGLESWITCHVF - Evaluation function for 3 dimensional toy Hill model for a oscillatory network.  
%
%   TOGGLESWITCHVF - Evaluates the vector field for the system of ODEs:
%               x1' = -x1 + H-(x2)  corresponds to x2 ---| x1
%               x2' = -x2 + H-(x1)  corresponds to x1 ---| x2
%
%   Inputs:
%       x - a 2-dimensional vector of state variables 
%       lambda - a 2-by-4 array of parameters. Each row corresponds to an edge of the network model as follows:
%                lambda(j,:) = [gamma_j, theta_ij, ell_j, delta_j]
%       n - a 2-dimensional vector of Hill coefficients
%
%   Outputs:
%       f(x, lambda, n) - a 2-dimensional vector
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 09-Jan-2020; Last revision: 09-Jan-2020

% unpack parameters
gamma = lambda(:,1);
theta = lambda(:,2);
ell = lambda(:,3);
delta = lambda(:,4);

% evaluate toggle switch polynomial map
denom1 = theta(2).^n(1) + x(2).^n(1);
denom2 = theta(1).^n(2) + x(1).^n(2);
p1 = (-gamma(1).*x(1) + ell(1)).*denom1 + delta(1).*theta(2).^n(1);
p2 = (-gamma(2).*x(2) + ell(2)).*denom2 + delta(2).*theta(1).^n(2);
hillPolynomialEvaluation = [p1; p2];
end % end TOGGLESWITCHVF

% Revision History:
%{

%}

% % ==================== OLD VERSION  % ====================
% % evaluate toggle switch vector field instead of the Hill polynomials
% f1 = -gamma(1).*x(1) + ell(1) + (delta(1).*theta(2).^n(1))./(theta(2).^n(1) + x(2).^n(1));
% f2 = -gamma(2).*x(2) + ell(2) + (delta(2).*theta(1).^n(2))./(theta(1).^n(2) + x(1).^n(2));
% hillEvaluation = [f1; f2];