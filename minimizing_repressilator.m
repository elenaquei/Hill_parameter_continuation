% test for minimizing the Hill coefficient of the  Hopf bifurcation in the repressilator system
% 26th Feb

clear all

main_repressilator;
% gamma = [1 1 1];
% theta = [1 1 1]; 
% l = [0 0 0];
% u = [2 2 2];
%
% true_sol(1:3) = x;
% true_sol(4) = n_temp;
% true_sol(5:7) = real(v_comp);
% true_sol(8:10) = imag(v_comp);
% true_sol(11) = imag(lambda);



x = true_sol(1:3);
n = true_sol(4);
v1 = true_sol(5:7);
v2 = true_sol(8:10);
beta = true_sol(11);



parameters_mat = [theta; l; u];
all_parameters = [gamma', parameters_mat(:).'];


vector_data = [n,x.',all_parameters,v1.',v2.',beta];


H = @(n,x,theta,ell,u)ell + (u-ell)*theta.^n./(theta.^n + x.^n);
dnH = @(n,x,theta,ell,u)(u-ell).*theta.^n.*x.^n.*(log(theta) - log(x))./(theta.^n + x.^n).^2; % use difference of logs to avoid complex valued log evaluations due to roundoff error
dxH = @(n,x,theta,ell,u)-n.*(u-ell).*theta.^n.*x^(n-1)./(theta.^n + x.^n).^2;
dxxH = @(n,x,theta,ell,u)-n.*(u-ell).*theta.^n.*x.^(n-2).*((n-1).*theta.^n - (n+1).*x.^n)./(theta.^n+x.^n).^3;
dxnH = @(n,x,theta,ell,u)(u-ell).*theta.^n.*x.^(n-1).*(n.*(theta.^n - x.^n).*(log(theta) - log(x)) - theta.^n - x.^n)./(theta.^n + x.^n).^3;
dthetaH = @(n,x,theta,ell,u)n.*(u-ell).*x.^n.*theta^(n-1)./(theta.^n + x.^n).^2;
dlH = @(n,x,theta,ell,u) 1 - theta.^n./(theta.^n + x.^n);
duH = @(n,x,theta,ell,u) theta.^n./(theta.^n + x.^n);
d_thetalu_H = @(n,x,theta,ell,u) [dthetaH(n,x,theta,ell,u), dlH(n,x,theta,ell,u), duH(n,x,theta,ell,u)];

vector_field = @(n,x,lambda) - diag(lambda(1:3))*x' + [ H(n,x(3), lambda(10),lambda(11),lambda(12))
    H(n,x(2), lambda(7),lambda(8),lambda(9))
    H(n,x(1), lambda(4),lambda(5),lambda(6))];

DxVF = @(n,x,lambda) - diag(lambda(1:3)) + [0 0 dxH(n,x(3), lambda(10),lambda(11),lambda(12))
0 dxH(n,x(2), lambda(7),lambda(8),lambda(9))  0
dxH(n,x(1), lambda(4),lambda(5),lambda(6)) 0 0];

DnVF = @(n,x,lambda)[ dnH(n,x(3), lambda(10),lambda(11),lambda(12))
    dnH(n,x(2), lambda(7),lambda(8),lambda(9))
    dnH(n,x(1), lambda(4),lambda(5),lambda(6))];


DparVFH = @(n,x,lambda) [0 0 0 0 0 0 d_thetalu_H(n,x(3), lambda(10),lambda(11),lambda(12))
    0 0 0 d_thetalu_H(n,x(2), lambda(7),lambda(8),lambda(9)) 0 0 0 
    d_thetalu_H(n,x(1), lambda(4),lambda(5),lambda(6)) 0 0 0 0 0 0];

DlambdaVF = @(n,x,lambda)[ -x(1) 0 0 
    0 -x(2) 0 
    0 0 -x(3)] + DparVFH(n,x,lambda);

eigen_prob = @(n,x,lambda,v1,v2,beta) [ DxVF(n,x,lambda) * v1' + beta * v2' 
    DxVF(n,x,lambda) * v2' - beta * v1'];

% complex amplitude being 1
amplitude = @(v1,v2) [sum(v1.^2+v2.^2)-1
    sum(v1.*v2)];


full_problem = @(n,x,lambda,v1,v2,beta) [vector_field(n,x,lambda)
    eigen_prob(n,x,lambda,v1,v2,beta)
    amplitude(v1,v2)];

full_problem_vector = @(x) full_problem(x(1),x(2:4),x(5:16),x(17:19),x(20:22),x(23));

% construct the function to minimize
minimization_func = @(n,x,lambda,v1,v2,beta,l) n + l * full_problem(n,x,lambda,v1,v2,beta);
minimization_func_vector = @(x) minimization_func(x(1),x(2:4),x(5:16),x(17:19),x(20:22),x(23),x(24:34));
% l needs to be of length 11

e4 = zeros(34,1);
e4(4) = 1;

full_problem_t = @(t,x) full_problem_vector(x);

% ERROR - THIS IS A WRONG USAGE OF THE FUNCTION numjac (needs a square
% derivative)
Ds = @(x) numjac(full_problem_t,0,x,full_problem_t(0,x),10^-12);

gradient_min_func= @(y) [y(24:34)*Ds(y(1:23))
    Zero_finding_problem_vector(y(1:23))]+e4;

starting_vec_for_minimization = [vector_data,rand(1,11)];

% actual minimization w.r.t. lambda and n
% xn = gradient_descent1D(minimization_func_vector,gradient_min_func,starting_vec_for_minimization);
% FMINSEARCH!

starting_vec_for_minimization = [vector_data,rand(1,11)];
xn = fminsearch(minimization_func_vector,starting_vec_for_minimization);

