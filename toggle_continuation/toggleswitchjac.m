function hillDiffEvaluation = toggleswitchjac(x, lambda, n)
%TOGGLESWITCHJAC - Evaluate the Jacobian of the equilibrium zero finding map.
%
%   TOGGLESWITCHJAC() - A more detailed description of the function
%
%   Syntax:
%       output = TOGGLESWITCHJAC(input1, input2)
%       [output1, output2] = TOGGLESWITCHJAC(input1, input2, input3)
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

D1p1 = -gamma(1).*(theta(2).^n(1) + x(2).^n(1));
D2p1 = n(1).*(-gamma(1).*x(1) + ell(1)).*x(2).^(n(1)-1);
D1p2 = n(2).*(-gamma(2).*x(2) + ell(2)).*x(1).^(n(2)-1);
D2p2 = -gamma(2).*(theta(1).^n(2) + x(1).^n(2));
hillDiffEvaluation = [D1p1, D2p1; D1p2, D2p2];
end % end toggleswitchjac

% Revision History:
%{

%}

% % ==================== OLD VERSION  % ====================
% % evaluate toggle switch vector field jacobian instead of the Hill polynomial jacobian
% D1f1 = -gamma(1);
% D2f1 = -(n(1).*delta(1).*x(2).^(n(1)-1).*theta(2).^n(1))./(theta(2).^n(1) + x(2).^n(1)).^2;
% D1f2 = -(n(2).*delta(2).*x(1).^(n(2)-1).*theta(1).^n(2))./(theta(1).^n(2) + x(1).^n(2)).^2;
% D2f2 = -gamma(2);
% hillDiffEvaluation = [D1f1, D2f1; D1f2, D2f2];