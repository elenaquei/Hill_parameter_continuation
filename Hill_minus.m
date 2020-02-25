function [H, dnH, dxH, dxxH, dxnH] = Hill_minus(theta,l,u)
% function [H, dnH, dxH, dxxH, dxnH] = Hill_minus(theta,l,u)
% 
% given l u and theta, returns the Hill function with appropriate parameter
% as a function of theta, x and n, together with plenty of derivatives
%
% INPUT
% theta,l,u           positive floats
% OUTPUT
% H             Hill- function handle, function of x and n
% dnH           a lot of derivatives
% dxH
% dxnH
% dxxH
% dnnH

%[H_plus, dnH_plus, dxH_plus, dxxH_plus, dxnH_plus] = Hill_plus(theta,l,u);
H = @(x,n) l + (u-l)*(theta^n)/(x^n+theta^n);
dnH = @(x,n) (u-l)*theta^n*x^n*(log(theta)-log(x))./(x.^n+theta.^n).^2;
dxnH = @(x, n) (u-l)*(theta^n*x^(n-1)*(n*log(theta/x)*(theta^n-x^n)+...
    (theta^n+x^n)))/(x^n+theta^n)^3;
dxH = @(x, n)(u-l)* n*theta^n*x^(n-1)/(theta^n+x^n)^2;
dxxH = @(x, n)(u-l)* n*theta^n*x^(n-2)*( theta^n*(n-1)-x^n*(n+1))/(theta^n+x^n)^3;
