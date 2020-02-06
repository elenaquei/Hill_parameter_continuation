function dnDf = toggleswitchjacgrad(x, lambda, n)
%TOGGLESWITCHJACGRAD - Evaluate derivative of the Jacobian of the equilibrium zero finding map with respect to the Hill coefficient.
%
%   TOGGLESWITCHJACGRAD() - A more detailed description of the function
%
%   Syntax:
%       output = TOGGLESWITCHJACGRAD(input1, input2)
%       [output1, output2] = TOGGLESWITCHJACGRAD(input1, input2, input3)
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

dn_D1f1 = -gamma(1).*(log(theta(2)).*theta(2).^n(1) + log(x(2)).*x(2).^n(1));
dn_D2f1 = (ell(1) - gamma(1).*x(1)).*x(2).^(n(1)-1).*(1 + n(1)*log(x(2)));
dn_D1f2 = (ell(2) - gamma(2).*x(2)).*x(1).^(n(2)-1).*(1 + n(2)*log(x(1)));
dn_D2f2 = -gamma(2).*(log(theta(1)).*theta(1).^n(2) + log(x(1)).*x(1).^n(2));
dnDf = [dn_D1f1, dn_D2f1; dn_D1f2, dn_D2f2];
end % end toggleswitchjacgrad

% Revision History:
%{

%}
% dn_D2f1 = n(1).*(n(1)-1).*(-gamma(1).*x(1) + ell(1)).*(x(2)^(n(1)-2));
% dn_D1f2 = n(2).*(n(2)-1).*(-gamma(2).*x(2) + ell(2)).*(x(1)^(n(2)-2));