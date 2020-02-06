function hillCoefficientGradient = toggleswitchgrad(x, lambda, n)
%TOGGLESWITCHGRAD - Evaluate derivative of equilibrium zero finding map with respect to the Hill coefficient
%
%   TOGGLESWITCHGRAD() - A more detailed description of the function
%
%   Syntax:
%       output = TOGGLESWITCHGRAD(input1, input2)
%       [output1, output2] = TOGGLESWITCHGRAD(input1, input2, input3)
%    
%   Inputs:
%       input1 - Description
%       input2 - Description
%       input3 - Description
%
%   Outputs:
%       output1 - Description
%       output2 - Description
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 15-Jan-2020; Last revision: 15-Jan-2020

% unpack parameters
gamma = lambda(:,1);
theta = lambda(:,2);
ell = lambda(:,3);
delta = lambda(:,4);

if length(n) > 1
    n = n(1);
end

D1n = (-gamma(1).*x(1) + ell(1)).*(log(theta(2)).*theta(2).^n + log(x(2)).*x(2).^n) + delta(1).*log(theta(2)).*theta(2).^n;
D2n = (-gamma(2).*x(2) + ell(2)).*(log(theta(1)).*theta(1).^n + log(x(1)).*x(1).^n) + delta(2).*log(theta(1)).*theta(1).^n;
hillCoefficientGradient = [D1n; D2n];
end % end toggleswitchgrad

% Revision History:
%{

%}
