function alpha_argmin = line_search(f,x, grad_f_x, p, max_alpha)
% function alpha_argmin = line_search(f,x, grad_f_x, p, max_alpha)
% 
% cheap search for 
%     alpha = argmin f(x + alpha p)
%
% INPUT
% f             function to minimize
% x             center of minimization
% grad_f_x      grad f(x)
% p             direction of search
% max_alpha     maximum accepted value of alpha
%
% OUTPUT
% alpha_argmin  approximation of argmin f(x + alpha p)

m = grad_f_x * p';

if numel(m)~=1
    error('m bust be a float')
end

tau = 1/2; % parameters can change
c = 1/2;

t = - c*m;
alpha = max_alpha;

while f(x) - f(x + alpha * p) < alpha*t
    alpha = tau * alpha;
end

alpha_argmin = alpha;