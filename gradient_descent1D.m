function xn = gradient_descent1D(f,gradf,x0)
% function xn = gradient_descent1D(f, gradf,x0)
%
% apply the gradient descent method (or steepest descent) to F with
% initial point x0
%
% INPUT
% f             function from R^n to R
% gradf         gradient of f
% x0            vector in R^n, first approx of minimum
%
% OUTPUT
% xn            local minimum of f

max_alpha = 1;
tol = 10^-4;

Dx0 = -gradf(x0);
alpha0 = line_search(f,x0,gradf(x0),Dx0,max_alpha);
x_j = x0 + alpha0*Dx0;
x_j(1:length(x_j)/2) = max(10^-6,x_j(1:length(x_j)/2)); 

iter = 1000;
old_Dx = Dx0;
old_s = Dx0;

% possibly find better iterate number!

for j = 2:iter
    Dx = - gradf(x_j);
    if any(imag(Dx))
        error('Complex values in gradient')
    end
    betaPR = Dx*(Dx-old_Dx)'/(old_Dx*old_Dx');
    beta = max(0,betaPR);
    s = Dx + beta*old_s;
    alpha = line_search(f,x_j,gradf(x_j),Dx,max_alpha);
    
    if any(imag(alpha*s))
        error('Complex values')
    end
    x_j = x_j + alpha*s;
    
    x_j(1:length(x_j)/2) = max(10^-6,x_j(1:length(x_j)/2)); 
    % reset to zero if negative values
    
    if norm(alpha*s)<tol
        break
    end
end

xn = x_j;