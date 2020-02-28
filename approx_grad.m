function gradf = approx_grad(f,x,dim_x)
% function gradf = approx_grad(f,x,dim_x)
%
% returns 
h = 10^-6;
gradf =zeros(dim_x,1);
for i = 1:dim_x
    e_i = 0*x;
    e_i(i) = 1;
    gradf(i) = ( f(x + h * e_i) - f(x) )/h;
end

if any(isnan(gradf))
    error('Here start the NaNs')
end
if any(imag(gradf))
    error('Here start the complex')
end