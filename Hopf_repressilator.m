function [x,n,v1,v2,beta] = Hopf_repressilator(gamma, theta, l, u)
% function [x,n,v1,v2,beta] = Hopf_repressilator(gamma, theta, l, u)
% 
% find Hopf bifurcations in the repressilator system for fixed parameters
% 27th Feb 2020 - Elena Queirolo
%
% INPUT
% gamma, theta, l, u        fixed parameters, each a 3-vector
% OUTPUT
% x                         fixed point, 3-vector
% n                         Hill coefficient
% v1                        real part of the eigenvector associated with the
%                           Hopf bifurcation
% v2                        imag part of the eigenvector associated with the
%                           Hopf bifurcation
% beta                      eigenvalue associated with the Hopf bifurcation

% building the vector field takes quite a while (lines 5-30)
dim = 3;
if length(gamma)~=dim 
    error('Given gamma values do not match the repressilator.')
end
if length(theta)~=dim 
    error('Given theta values do not match the repressilator.')
end
if length(l)~=dim 
    error('Given l values do not match the repressilator.')
end
if length(u)~=dim 
    error('Given u values do not match the repressilator.')
end

Gamma = diag(gamma);

% using the data to construct Hill functions
[H_minus1, ~, dxH_minus1] = ...
    hillcomponent('-', theta(2),l(2), u(2)-l(2));
[H_minus2, ~, dxH_minus2] = ...
    hillcomponent('-', theta(3),l(3), u(3)-l(3));
[H_minus3, ~, dxH_minus3] = ...
    hillcomponent('-', theta(1),l(1), u(1)-l(1));

% assembly the Hill functions
H_minus= @(x,n) [H_minus1(x(2),n);
                H_minus2(x(3),n);
                H_minus3(x(1),n)];
            
dxH_minus = @(x,n)[ 0                   dxH_minus1(x(2),n)  0
                    0                   0                   dxH_minus2(x(3),n)
                    dxH_minus3(x(1),n)  0                   0];

% writing the full vector field
vector_field = @(x,n) -Gamma*x + H_minus(x,n);
Dx_vector_field = @(x,n) -Gamma + dxH_minus(x,n);

eigen_cond_real = @(x,n,w,z,beta) Dx_vector_field(x,n)*w + z*beta;
eigen_cond_imag = @(x,n,w,z,beta) Dx_vector_field(x,n)*z - w*beta;

amplitude_cond = @(w,z) [sum(w.^2+z.^2) - 1;
            2*sum(w.*z)];

Zero_finding_problem = @(x,n,w,z,beta) [vector_field(x,n);
    eigen_cond_real(x,n,w,z,beta);
    eigen_cond_imag(x,n,w,z,beta);
    amplitude_cond(w,z)];

Zero_finding_problem_vector = @(vec) Zero_finding_problem...
    (vec(1:dim),vec(dim+1),vec(dim+1+(1:dim)),vec(2*dim+1+(1:dim)),vec(3*dim+2));

% finding a first approximation for the fixed point with a "low" n
n_temp=4.2;%n_temp=4.8;%n_temp=4.2;%n_temp = 3;
x = Newton_handle(@(x)vector_field(x,n_temp),[1,1,1]',@(x) Dx_vector_field(x,n_temp));

% use the first approximation to find the eigenvalue and eigenvector match
n_temp=4.2;%n_temp=4.8;
[v,eig_loc]= eig(Dx_vector_field(x,n_temp));
all_eigs= diag(eig_loc);
[~,index] = min(abs(real(all_eigs)));
beta = all_eigs(index); 
v_comp = v(:,index);

approx_vec = [x;n_temp;real(v_comp);imag(v_comp);imag(beta)];

% look for a better solution
true_sol = Newton_handle(Zero_finding_problem_vector,approx_vec);

x = true_sol(1:3);
n = true_sol(4);
v1 = true_sol(5:7);
v2 = true_sol(8:10);
beta = true_sol(11);
