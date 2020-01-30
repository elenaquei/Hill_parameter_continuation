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


[H_plus, dnH_plus, dxH_plus, dxxH_plus, dxnH_plus] = Hill_plus(theta,l,u);
H = @(x,n) H_plus(x, n);
dnH = @(x,n) - dnH_plus(x,n);
dxnH = @(x, n) - dxnH_plus(x,n);
dxH = @(x, n) - dxH_plus(x,n);
dxxH = @(x, n) - dxxH_plus(x,n);