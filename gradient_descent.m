function xn = gradient_descent(F, dF,x0)
% function xn = gradient_descent(F, dF,x0)
%
% apply the gradient descent method (or steepest descent) to F with
% initial point x0
%
% INPUT
% F         function from R^n to R^n, column vector
% dF        d_xF 
% x0        vector in R^n, first approx of minimum
%
% OUTPUT
% xn        local minimum of F

max_alpha = 1;
tol = 10^-6;

% f(x) = norm(F(x))^2
% grad f(x) = F(x)' * dF(x)
f = @(x) sum(F(x).^2);
gradf = @(x) 2*F(x)'*dF(x);

Dx0 = -gradf(x0);
alpha0 = line_search(f,x0,gradf(x0),Dx0,max_alpha);
x_j = x0 + alpha0*Dx0;

iter = 1000;
old_Dx = Dx0;
old_s = Dx0;

% possibly find better iterate number!

for j = 2:iter
    Dx = - gradf(x_j);
    betaPR = Dx*(Dx-old_Dx)'/(old_Dx*old_Dx');
    beta = max(0,betaPR);
    s = Dx + beta*old_s;
    alpha = line_search(f,x_j,gradf(x_j),Dx,max_alpha);
    x_j = x_j + alpha*s;
    if norm(alpha*s)<tol*10^-2
        break
    end
end

xn = x_j;