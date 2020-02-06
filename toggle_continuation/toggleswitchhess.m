function hillHessianEvaluation = toggleswitchhess(x, lambda, n)
%TOGGLESWITCHHESS - Evaluate the Hessian of the equilibrium zero finding map.
%
%   TOGGLESWITCHHESS() - A more detailed description of the function
%
%   Syntax:
%       output = TOGGLESWITCHHESS(input1, input2)
%       [output1, output2] = TOGGLESWITCHHESS(input1, input2, input3)
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
ell = lambda(:,3);

H1 = [0, -n(1).*gamma(1).*x(2).^(n(1)-1); -n(1).*gamma(1).*x(2).^(n(1)-1), n(1).*(n(1)-1).*(-gamma(1).*x(1) + ell(1)).*x(2).^(n(1)-2)];
H2 = [0, -n(2).*gamma(2).*x(1).^(n(2)-1); -n(2).*gamma(2).*x(1).^(n(2)-1), n(2).*(n(2)-1).*(-gamma(2).*x(2) + ell(2)).*x(1).^(n(2)-2)];
hillHessianEvaluation = {H1, H2};
end % end toggleswitchhess

% Revision History:
%{

%}
