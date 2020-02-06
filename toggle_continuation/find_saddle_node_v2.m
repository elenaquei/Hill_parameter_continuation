%FIND_SADDLE_NODE_V2 - One line description of what the script performs (H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%   Optional file header info (to give more details about the function than in the H1 line)
%
%   Description:
%       FIND_SADDLE_NODE_V2 description
%
%   Output:
%       FIND_SADDLE_NODE_V2 output
%
%   Other m-files required: none
%   MAT-files required: none
%
%   See also: OTHER_SCRIPT_NAME,  OTHER_FUNCTION_NAME

%   Author: Shane Kepley
%   email: shane.kepley@rutgers.edu
%   Date: 31-Jan-2020; Last revision: 31-Jan-2020

close all
clear all
clc
addpath('/Users/shane/Hill_parameter_continuation')

%% ================================================== SECTION 1 ==================================================
clc
% % PARAMETER SET 2 WHERE A SADDLE NODE HAPPENS AT SOME CHOICE OF 3 < N < 4
gamma = [1.2;.87];
ell = [1.4;1];
% theta = [3;3];
theta = [2;3.6];

delta = [5;6];
n0 = 4.1;

eq1 = [ell(1); ell(2) + delta(2)]./gamma;
eq2 = theta;
eq3 = [ell(1) + delta(1); ell(2)]./gamma;
x0 = eq2;


% % PARAMETER SET 2 WHERE EVEN THE INITIAL ROOT BECOMES COMPLEX SUDDENLY IN THE 13TH NEWTON STEP
% gamma = [1;1];
% ell = [1;1];
% theta = [3;3];
% delta = [5;6];
% n0 = 3.15;
% x0 = [4.4;3.7];

% set parameters
lambda = [gamma, theta, ell, delta];

% nullcline plots
X1Range = [0, 10];
X2Range = [0, 10];
X1 = linspace(X1Range(1), X1Range(2), 500);
X2 = linspace(X2Range(1), X2Range(2), 500);
N1 = ell(1)./gamma(1) + (delta(1).*theta(2).^n0)./(gamma(1).*(theta(2).^n0 + X2.^n0));   % f1 = 0 nullcline
N2 = ell(2)./gamma(2) + (delta(2).*theta(1).^n0)./(gamma(2).*(theta(1).^n0 + X1.^n0));   % f2 = 0 nullcline

% find equilibria
[H1, gradH1, DH1, D2H1, gradDH1] = hillcomponent('-', theta(2), ell(1), delta(1));
[H2, gradH2, DH2, D2H2, gradDH2] = hillcomponent('-', theta(1), ell(2), delta(2));
f0 = @(x)-lambda(:,1).*x + [H1(x(2), n0); H2(x(1), n0)];
Df0 = @(x)[-gamma(1), DH1(x(2),n0); DH2(x(1),n0), -gamma(2)];

%
xInit = findroot(f0, Df0, x0);
xFull = findroot(f0, Df0, x0, 'FullOrbit', true);


% plot nullclines and equilibrium
figure
hold on
plot(N1, X2, 'b', 'LineWidth', 3)
plot(X1, N2, 'r', 'LineWidth', 3)
legend('f1=0','f2=0')
scatter(xInit(1), xInit(2), 'g*')
title(sprintf('n = %0.4f', n0))
dealfig()

%% ================================================== FIND SADDLE NODE BIFURCATION ==================================================
f = @(x, n)-lambda(:,1).*x + [H1(x(2), n); H2(x(1), n)];
Df = @(x, n)[-gamma(1), DH1(x(2),n); DH2(x(1),n), -gamma(2)];
D2fv = @(x, n, v)[0, D2H1(x(2), n)*v(2); D2H2(x(1), n)*v(1), 0];
gradf = @(x, n)[gradH1(x(2), n); gradH2(x(1), n)];
gradDfv = @(x, n, v)[gradDH1(x(2), n)*v(2); gradDH2(x(1), n)*v(1)];

g = @(u)[f(u(1:2), u(3)); u(4)-1; Df(u(1:2), u(3))*u(4:5)];
Dg = @(u)[Df(u(1:2), u(3)), gradf(u(1:2), u(3)), zeros(2); zeros(1,2), 0, [1,0]; D2fv(u(1:2), u(3), u(4:5)), gradDfv(u(1:2), u(3), u(4:5)), Df(u(1:2), u(3))];
v0= [1;-.7];
% u0 = [x0; n0; v0];
u0 = [xInit; n0; v0];



% ground truth from version 1
% xTT = findroot(@(x)f(x,nTrue), @(x)Df(x, nTrue), xTrue, 'FullOrbit', true);

u = findroot(g, Dg, u0);
xSol = u(1:2);
nSol = u(3);
vSol = u(4:5);
n0 = nSol;


% find equilibria
f0 = @(x)toggleswitchvf(x, lambda, n0);
Df0 = @(x)toggleswitchjac(x, lambda, n0);


X1Range = [0, 10];
X2Range = [0, 10];
X1 = linspace(X1Range(1), X1Range(2), 200);
X2 = linspace(X2Range(1), X2Range(2), 200);
SNparm = nSol;
N1 = ell(1)./gamma(1) + (delta(1).*theta(2).^SNparm)./(gamma(1).*(theta(2).^SNparm + X2.^SNparm));   % f1 = 0 nullcline
N2 = ell(2)./gamma(2) + (delta(2).*theta(1).^SNparm)./(gamma(2).*(theta(1).^SNparm + X1.^SNparm));   % f2 = 0 nullcline

figure
hold on
plot(N1, X2, 'b', 'LineWidth', 3)
plot(X1, N2, 'r', 'LineWidth', 3)
legend('f1=0','f2=0')
scatter(u(1), u(2), 'g*')
title(sprintf('n = %0.4f', n0))

dealfig()
uFull = findroot(g, Dg, u0, 'FullOrbit', true)


