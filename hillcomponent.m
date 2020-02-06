function [H, dnH, dxH, dxxH, dxnH] = hillcomponent(interactionType, theta, ell, delta, varargin)
%HILLCOMPONENT - Returns a function handle to evaluate a Hill function for positive parameters {theta, ell, delta}
%   HILLCOMPONENT() - A more detailed description of the function
%
%   Syntax:
%       H = HILLCOMPONENT(+, theta, ell, delta) returns the cooperative (monotone increasing) version of the Hill function with specified parameters.
%       H = HILLCOMPONENT(-, theta, ell, delta) returns the destructive (monotone decreasing) version of the Hill function.
%
%   Inputs:
%       interactionType - A string in {'+', '-'} which specifies whether the interaction is cooperation or inhibition.
%       theta - threshold parameter
%       ell - lower bound of dynamic range
%       delta - dynamic range
%
%   Outputs:
%       H - A handle for evaluating a scalar hill function of the form H(x, n) where n is the Hill coefficient
%
%   Subfunctions: none
%   Classes required: none
%   Other m-files required: none
%   MAT-files required: none

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 03-Feb-2020; Last revision: 03-Feb-2020

if strcmp(interactionType, '-')
    H = @(x,n)ell + delta*theta.^n./(theta.^n + x.^n);
    dnH = @(x,n)delta.*theta.^n.*x.^n.*(log(theta) - log(x))./(theta.^n + x.^n).^2; % use difference of logs to avoid complex valued log evaluations due to roundoff error
    dxH = @(x, n)-n.*delta.*theta.^n.*x^(n-1)./(theta.^n + x.^n).^2;
    dxxH = @(x, n)-n.*delta.*theta.^n.*x.^(n-2).*((n-1).*theta.^n - (n+1).*x.^n)./(theta.^n+x.^n).^3;
    dxnH = @(x, n)delta.*theta.^n.*x.^(n-1).*(n.*(theta.^n - x.^n).*(log(theta) - log(x)) - theta.^n - x.^n)./(theta.^n + x.^n).^3;
    
elseif strcmp(interactionType, '+')
    H = @(x,n)ell + delta*x.^n./(theta.^n + x.^n);
    dnH = @(x,n)delta.*theta.^n.*x.^n.*(log(x) - log(theta))./(theta.^n + x.^n).^2; % use difference of logs to avoid complex valued log evaluations due to roundoff error
    dxH = @(x, n)n.*delta.*theta.^n.*x^(n-1)./(theta.^n + x.^n).^2;
    dxxH = @(x, n)n.*delta.*theta.^n.*x.^(n-2).*((n-1).*theta.^n - (n+1).*x.^n)./(theta.^n+x.^n).^3;
    dxnH = @(x, n)delta.*theta.^n.*x.^(n-1).*(-n.*(theta.^n - x.^n).*(log(theta) - log(x)) + theta.^n + x.^n)./(theta.^n + x.^n).^3;
    
else
    error('Interaction type should be + or -')
end
end % end hillcomponent

% Revision History:
%{

%}
