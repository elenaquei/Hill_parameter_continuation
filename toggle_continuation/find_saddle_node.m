%FIND_SADDLE_NODE - One line description of what the script performs (H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%
%   Description:
%       FIND_SADDLE_NODE description
%
%   Output:
%       FIND_SADDLE_NODE output
%
%   Other m-files required: none
%   MAT-files required: none
%
%   See also: OTHER_SCRIPT_NAME,  OTHER_FUNCTION_NAME

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 15-Jan-2020; Last revision: 15-Jan-2020

%% ================================================== SECTION 1 ==================================================
close all
% % PARAMETER SET 1 WHERE THE EQUILIBRIUM SEEMS TO BE HARD TO FIND
% gamma = [1;1];
% ell = [1;1];
% theta = [3;4];
% delta = [5;6];
% n = [3;4];

% % PARAMETER SET 2 WHERE A SADDLE NODE HAPPENS AT SOME CHOICE OF 3 < N < 4
gamma = [1;1];
ell = [1;1];
theta = [3;3];
delta = [5;6];
n = 3.5*[1;1]; 
x0 = [5.7;1.4];


% % PARAMETER SET 3 
% gamma = [1;1.25];
% ell = [1.5;2];
% theta = [2.9;1.7];
% delta = [6.6;8.4];
% n = 10*[1;1]; 
% x0 = [5.6;1.6];



% set parameters
lambda = [gamma, theta, ell, delta];

% nullcline plots
X1Range = [0, 10];
X2Range = [0, 10];
X1 = linspace(X1Range(1), X1Range(2), 200);
X2 = linspace(X2Range(1), X2Range(2), 200);
N1 = ell(1)./gamma(1) + (delta(1).*theta(2).^n(1))./(gamma(1).*(theta(2).^n(1) + X2.^n(1)));   % f1 = 0 nullcline
N2 = ell(2)./gamma(2) + (delta(2).*theta(1).^n(2))./(gamma(2).*(theta(1).^n(2) + X1.^n(2)));   % f1 = 0 nullcline

% find equilibria
f = @(x)toggleswitchvf(x, lambda, n);
Df = @(x)toggleswitchjac(x, lambda, n);
xInit = findroot(f, Df, x0);
xFull = findroot(f, Df, x0, 'FullOrbit', true)
% plot nullclines and equilibrium

figure
hold on
plot(N1, X2, 'b')
plot(X1, N2, 'r')
legend('f1=0','f2=0')
scatter(xInit(1), xInit(2), 'g*')


%% ================================================== FIND SADDLE NODE BIFURCATION ==================================================
v = [1;1];
g = @(u)saddlenodezf(u(1:2), lambda, u(3), u(4:5));
Dg = @(u)saddlenodejac(u(1:2), lambda, u(3), u(4:5));
u0 = [xInit; n(1); v];
u = findroot(g, Dg, u0);
xSol = u(1:2);
nSol = u(3);
vSol = u(4:5);
n = [nSol; nSol];

% find equilibria
f = @(x)toggleswitchvf(x, lambda, n);
Df = @(x)toggleswitchjac(x, lambda, n);


X1Range = [0, 10];
X2Range = [0, 10];
X1 = linspace(X1Range(1), X1Range(2), 200);
X2 = linspace(X2Range(1), X2Range(2), 200);
SNparm = u(3);
N1 = ell(1)./gamma(1) + (delta(1).*theta(2).^SNparm)./(gamma(1).*(theta(2).^SNparm + X2.^SNparm));   % f1 = 0 nullcline
N2 = ell(2)./gamma(2) + (delta(2).*theta(1).^SNparm)./(gamma(2).*(theta(1).^SNparm + X1.^SNparm));   % f1 = 0 nullcline

figure
hold on
plot(N1, X2, 'b')
plot(X1, N2, 'r')
legend('f1=0','f2=0')
scatter(u(1), u(2), 'g*')


dealfig()

uFull = findroot(g, Dg, u0, 'FullOrbit', true)


