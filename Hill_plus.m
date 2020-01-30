function [H, dnH, dxH, dxxH, dxnH] = Hill_plus(theta,l,u)
% function [H, dnH, dxH, dxxH, dxnH] = Hill_plus(theta,l,u)
% 
% given l u and theta, returns the Hill function with appropriate parameter
% as a function of theta, x and n, together with plenty of derivatives
%
% INPUT
% theta,l,u           positive floats
% OUTPUT
% H             Hill+ function handle, function of x and n
% dnH           a lot of derivatives
% dxH
% dxnH
% dxxH
% dnnH


H = @(x,n) l + (u-l)*x.^n./(theta.^n+x.^n);
dxH = @(x,n) (u-l)*(n*x.^(n-1)*theta.^n) ./(theta.^n+x.^n).^2;
dnH = @(x,n) (u-l)*(x.^n*ln(x) .*theta.^n - theta.^n*ln(theta))./(theta.^n+x.^n).^2;
dxnH = @(x,n) - (u-l)*x.^(n-1)*theta.^n.*( (theta.^n -x.^n)*n*log(theta./x)-theta.^n-x.^n)./(theta.^n+x.^n).^3;
dxxH = @(x,n) -n*(u-l)*theta.^n*x.^(n-2)*((n-1)*theta.^n-(n+1)*x.^n)./(theta.^n+x.^n).^3;
