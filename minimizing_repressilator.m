% minimizing the Hill coefficient of the  Hopf bifurcation in the repressilator system
% 26th Feb


% in the repressilator, to have periodic behavior we need
% % ["{
% % L(Z,X) < THETA(X,Y),
% % THETA(X,Y) < U(Z,X)
% % },{
% % THETA(X,Y)
% % }",
% % "{
% % L(X,Y) < THETA(Y,Z),
% % THETA(Y,Z) < U(X,Y)
% % },{
% % THETA(Y,Z)
% % }",
% % "{
% % L(Y,Z) < THETA(Z,X),
% % THETA(Z,X) < U(Y,Z)
% % },{
% % THETA(Z,X)
% % }"]
% that I think corresponds to l < theta < u
gamma = [1 1 1];
theta = [1 1 1];
l = [0 0 0]+0.5; 
u = [2 2 2];

% find the Hopf position at given values of the parameters
% [x,n,v1,v2,beta] = Hopf_repressilator(gamma, theta, l, u);
% there is probably a bug in there
%
%return

% structuring the problem: defining the Hill function w.r.t all the
% parameters and its derivative w.r.t. anything I could think of

H = @(n,x,theta,ell,u)ell + (u-ell)*theta.^n./(theta.^n + x.^n);
dnH = @(n,x,theta,ell,u)(u-ell).*theta.^n.*x.^n.*(log(theta) - log(x))./(theta.^n + x.^n).^2; % use difference of logs to avoid complex valued log evaluations due to roundoff error
dxH = @(n,x,theta,ell,u)-n.*(u-ell).*theta.^n.*x^(n-1)./(theta.^n + x.^n).^2;
dxxH = @(n,x,theta,ell,u)-n.*(u-ell).*theta.^n.*x.^(n-2).*((n-1).*theta.^n - (n+1).*x.^n)./(theta.^n+x.^n).^3;
dxnH = @(n,x,theta,ell,u)(u-ell).*theta.^n.*x.^(n-1).*(n.*(theta.^n - x.^n).*(log(theta) - log(x)) - theta.^n - x.^n)./(theta.^n + x.^n).^3;
dthetaH = @(n,x,theta,ell,u)n.*(u-ell).*x.^n.*theta^(n-1)./(theta.^n + x.^n).^2;
dlH = @(n,x,theta,ell,u) 1 - theta.^n./(theta.^n + x.^n);
duH = @(n,x,theta,ell,u) theta.^n./(theta.^n + x.^n);
d_thetalu_H = @(n,x,theta,ell,u) [dthetaH(n,x,theta,ell,u), dlH(n,x,theta,ell,u), duH(n,x,theta,ell,u)];

% the repressilator vector field
vector_field = @(n,x,lambda) - diag(lambda(1:3))*x.' + 0*[ H(n,x(3), lambda(10),lambda(11),lambda(12))
    H(n,x(2), lambda(7),lambda(8),lambda(9))
    H(n,x(1), lambda(4),lambda(5),lambda(6))];

% useful derivative
DxVF = @(n,x,lambda) - diag(lambda(1:3)) + ...
    [0 0 dxH(n,x(3), lambda(10),lambda(11),lambda(12))
    0 dxH(n,x(2), lambda(7),lambda(8),lambda(9))  0
    dxH(n,x(1), lambda(4),lambda(5),lambda(6)) 0 0];

% MORE derivatives! %%
DnVF = @(n,x,lambda)[ dnH(n,x(3), lambda(10),lambda(11),lambda(12))
    dnH(n,x(2), lambda(7),lambda(8),lambda(9))
    dnH(n,x(1), lambda(4),lambda(5),lambda(6))];

DparVFH = @(n,x,lambda) [0 0 0 0 0 0 d_thetalu_H(n,x(3),lambda(10),lambda(11),lambda(12))
    0 0 0 d_thetalu_H(n,x(2), lambda(7), lambda(8), lambda(9)) 0 0 0
    d_thetalu_H(n,x(1), lambda(4), lambda(5), lambda(6)) 0 0 0 0 0 0];

DlambdaVF = @(n,x,lambda)cat(2,...
    [-x(1) 0 0
    0 -x(2) 0
    0 0 -x(3)],DparVFH(n,x,lambda));

full_der_VF = @(n,x,lambda)...
    [DnVF(n,x,lambda), DxVF(n,x,lambda), DlambdaVF(n,x,lambda) zeros(3,3*2+1)];

full_der_eigen_prob = @(n,x,lambda,v1,v2,beta)...
    [DxnVF(n,x,lambda)* v1' DxxVF(n,x,lambda,v1)  DlambdaxVF(n,x,lambda, v1) DxVF(n,x,lambda) beta*(1+0*v1')   v2'
    DxnVF(n,x,lambda)* v2' DxxVF(n,x,lambda,v2)  DlambdaxVF(n,x,lambda, v2) -beta*(1+0*v1')  DxVF(n,x,lambda) -v1'];

full_der_amplitude = @(n,x,lambda,v1,v2,beta)[0 0 0 0 0*lambda 2*v1 2*v2 0
    0 0 0 0 0*lambda   v2   v1 0];

e1 = zeros(34,1);
e1(1) = 1;

der_full_problem = @(n,x,lambda,v1,v2,beta) [full_der_vector_field(n,x,lambda)
    full_der_eigen_prob(n,x,lambda,v1,v2,beta)
    full_der_amplitude(v1,v2)];

der_full_problem_vector = @(x) der_full_problem(x(1),x(2:4),x(5:16),x(17:19),x(20:22),x(23));

gradient_min_func = @(y) [y(24:34)*Ds(y(1:23))
    full_problem_vector(y(1:23))]+e1;
% end of unused derivatives

% eigenvalue problem for Hopf 
eigen_prob = @(n,x,lambda,v1,v2,beta) [ DxVF(n,x,lambda) * v1' + beta * v2'
    DxVF(n,x,lambda) * v2' - beta * v1'];

% complex amplitude of the eigenvalue being 1
amplitude = @(v1,v2) [sum(v1.^2+v2.^2)-1
    sum(v1.*v2)];

% putting all the constraints together
full_problem = @(n,x,lambda,v1,v2,beta) ...
    [vector_field(n,x,lambda)
    eigen_prob(n,x,lambda,v1,v2,beta)
    amplitude(v1,v2)];

% turning the input into a vector
full_problem_vector = @(x) full_problem(x(1),x(2:4),x(5:16),x(17:19),x(20:22),x(23));

% lagrange multipliers, with square to ensure positivity
minimization_func = @(n,x,lambda,v1,v2,beta,l) n + l * full_problem(n,x,lambda,v1,v2,beta);
minimization_func_vector = @(x) minimization_func(x(1),x(2:4),x(5:16),x(17:19),x(20:22),x(23),x(24:34)).^2;
% the lagrange multiplier l needs to be of length 11

for i = 1:500
    % set up random initial points 
    gamma = rand(1,3);
    l = rand(1,3);
    theta = l+rand(1,3);
    u = theta+rand(1,3);
    
    parameters_mat = [gamma; theta; l; u];
    all_parameters = parameters_mat(:).';
    
    n = rand*3;
    x = rand(3,1);
    v1 = rand(3,1);
    v2 = rand(3,1);
    beta = rand;
    vector_data = [n,x.',all_parameters,v1.',v2.',beta];
    
    % do the minimization of the langrange function 
    starting_vec_for_minimization = [vector_data,ones(1,11)];
    xn = fminsearch(minimization_func_vector,starting_vec_for_minimization);
    
    % check the data
    n = xn(1);
    x = xn(2:4);
    all_parameters = xn(4+(1:4*3));
    parameters_mat = reshape(all_parameters,4,3);
    gammas = parameters_mat(1,:);
    thetas = parameters_mat(2,:);
    ls = parameters_mat(3,:);
    us = parameters_mat(4,:);
    if any(all_parameters<0)
        warning('Negative parameters at iteration %i', i)
        continue
    elseif any(x<0)
        warning('Negative fixed point at iteration %i', i)
        continue
    elseif any(n<0)
        warning('Negative exponent at iteration %i', i)
        continue
    end
    
    % show off the data
    plot(sum(xn(2:end)),n,'b.','MarkerSize',13)
    hold on;
end
