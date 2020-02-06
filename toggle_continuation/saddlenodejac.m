function saddleNodeJacobianEvaluation = saddlenodejac(x, lambda, n, v)
%SADDLENODEJAC - Evaluate the jacobian of the saddle node zero finding map
%
%   SADDLENODEJAC() - A more detailed description of the function
%
%   Syntax:
%       output = SADDLENODEJAC(input1, input2)
%       [output1, output2] = SADDLENODEJAC(input1, input2, input3)
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

% % unpack parameters
% gamma = lambda(:,1);
% theta = lambda(:,2);
% ell = lambda(:,3);
% delta = lambda(:,4);

Df = toggleswitchjac(x, lambda, [n,n]);
gradf = toggleswitchgrad(x, lambda, [n,n]);
D2f = toggleswitchhess(x, lambda, [n,n]); 
D2fv = [D2f{1}*v, D2f{2}*v].';
gradDf = toggleswitchjacgrad(x, lambda, [n,n]);
saddleNodeJacobianEvaluation = [Df,gradf,zeros(2); zeros(1,2), 0, [1, 0]; D2fv, gradDf*v, Df]; 
end % end saddlenodejac

% Revision History:
%{

%}
