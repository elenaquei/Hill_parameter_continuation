% test for finding Hopf bifurcations in the repressilator system
% 30th Jan 2020

% building the vector field takes quite a while (lines 5-30)
dim = 3;
Gamma = diag(ones(dim,1));
theta = [1 1 1]; 
l = [0 0 0];
u = [2 2 2];
% using the data to construct Hill functions
[H_minus1, dnH_minus1, dxH_minus1, dxxH_minus1, dxnH_minus1] = ...
    Hill_minus(theta(2),l(2),u(2));
[H_minus2, dnH_minus2, dxH_minus2, dxxH_minus2, dxnH_minus2] = ...
    Hill_minus(theta(3),l(3),u(3));
[H_minus3, dnH_minus3, dxH_minus3, dxxH_minus3, dxnH_minus3] = ...
    Hill_minus(theta(1),l(1),u(1));

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
Dn_vector_field = @(x,n) dnH_minus(x,n);

eigen_cond_real = @(x,n,w,z,lambda) Dx_vector_field(x,n)*w + z*lambda;
eigen_cond_imag = @(x,n,w,z,lambda) Dx_vector_field(x,n)*z - w*lambda;

amplitude_cond = @(x,n,w,z,lambda) [sum(w.^2+z.^2) - 1;
            sum(w.*z)];

Zero_finding_problem = @(x,n,w,z,lambda) [vector_field(x,n);
    eigen_cond_real(x,n,w,z,lambda);
    eigen_cond_imag(x,n,w,z,lambda);
    amplitude_cond(x,n,w,z,lambda)];

Zero_finding_problem_vector = @(vec) Zero_finding_problem(vec(1:dim),vec(dim+1),vec(dim+1+(1:dim)),vec(2*dim+1+(1:dim)),vec(3*dim+2));

approx_vec = rand(3*dim+2,1);

% finding a first approximation for the fixed point with a "low" n
n_temp = 3;
x = Newton_handle(@(x)vector_field(x,n_temp),[1,1,1]',@(x) Dx_vector_field(x,n_temp));

% use the first approximation to find the eigenvalue and eigenvector match
n_temp=4.2;
[v,l]= eig(Dx_vector_field(x,n_temp));
lambda= l(end);
v_comp = v(:,end);

approx_vec(1:3) = x;
approx_vec(4) = n_temp;
approx_vec(5:7) = real(v_comp);
approx_vec(8:10) = imag(v_comp);
approx_vec(11) = imag(lambda);

% look for a better solution
true_sol = Newton_handle(Zero_finding_problem_vector,approx_vec);

 

