function DF = numjac_loc(f,x,h)
% function DF = numjac_loc(f,x,h)
%
% 
if nargin<3 || isempty(h)
    h = 10^-3;
end
dim_x = length(x);
DF =zeros(dim_x,dim_x);
for i = 1:dim_x
    e_i = 0*x;
    e_i(i) = 1;
    DF(:,i) = ( f(x + h * e_i) - f(x) )/h;
end

if any(isnan(DF))
    error('Here start the NaNs')
end
if any(imag(DF))
    error('Here start the complex')
end