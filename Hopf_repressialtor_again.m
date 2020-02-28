% another try at finding the Hopf bifurcation
% 26th Feb


% in the repressilator, to have periodic behavior we need
% what I think corresponds to l < theta < u
gamma = [1 1 1];
theta = [1 1 1]-0.2;
l = [0 0 0]+0.5; 
u = [2 2 2];

H = @(n,x,theta,ell,u) ell + (u-ell)*theta.^n./(theta.^n + x.^n);
dxH = @(n,x,theta,ell,u)-n.*(u-ell).*theta.^n.*x^(n-1)./(theta.^n + x.^n).^2;

% the repressilator vector field
vector_field = @(n,x,lambda) -diag(lambda(1:3))*x + ...
    [H(n,x(3), lambda(10),lambda(11),lambda(12))
    H(n,x(2), lambda(7),lambda(8),lambda(9))
    H(n,x(1), lambda(4),lambda(5),lambda(6))];

% useful derivative
DxVF = @(n,x,lambda) - diag(lambda(1:3)) + ...
    [0 0 dxH(n,x(3), lambda(10),lambda(11),lambda(12))
    0 dxH(n,x(2), lambda(7),lambda(8),lambda(9))  0
    dxH(n,x(1), lambda(4),lambda(5),lambda(6)) 0 0];

% eigenvalue problem for Hopf 
eigen_prob = @(n,x,lambda,v1,v2,beta) [ DxVF(n,x,lambda) * v1 + beta * v2
    DxVF(n,x,lambda) * v2 - beta * v1];

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

full_problem_parameters = @(y, lambda) full_problem(y(1),y(2:4),lambda, y(5:7),y(8:10),y(11));

par = [1,1,1,1.1,1.2,1.3,0.5,0.6,0.7,2.1,2.4,2.7];

vector_field_par = @(y) vector_field(y(1),y(2:4),par);
H_par = @(y) H(y(1),y(2), theta(1), l(1), u(1));

full_problem_fixed_par = @(y) full_problem_parameters(y,par);


for i = 2 : 20
approx_sol = [i,[1.1,1.2,1.3],[1,2,3],[4,5,6],3]';

approx_x = Newton_handle(@(x) vector_field(approx_sol(1),x, par),approx_sol(2:4));
approx_sol(2:4) = approx_x;


[v, eigenval] = eigs( DxVF(approx_sol(1),approx_x,par) )
end 

% true_sol = Newton_handle(full_problem_fixed_par,approx_sol);