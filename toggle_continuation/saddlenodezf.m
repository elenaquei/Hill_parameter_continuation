function saddleNodeMapEvaluation = saddlenodezf(x, lambda, n, v)
%SADDLENODEZF - Evaluate the saddle-node zero finding map.
%
%   SADDLENODEZF() - A more detailed description of the function
%
%   Syntax:
%       output = SADDLENODEZF(input1, input2)
%       [output1, output2] = SADDLENODEZF(input1, input2, input3)
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

Df = toggleswitchjac(x, lambda, [n,n]); % evaluate derivative of equilibrium map
g1 = toggleswitchvf(x, lambda, [n,n]); % evaluate equilibrium map
g2 = v(1) - 1; % evaluate eigenvector scaling condition
g3 = Df*v; % evaluate kernel map
saddleNodeMapEvaluation = [g1;g2;g3];
end % end saddlenodezf

% Revision History:
%{

%}
