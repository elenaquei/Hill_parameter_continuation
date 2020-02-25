function zero_approx = Newton_handle(function_handle, approx, DF_func)

max_iter = 50;
tol = 10^-9;

if nargin<3 || isempty(DF_func)
    function_handle_t = @(t,x) function_handle(x);
     DF_func = @(x) numjac(function_handle_t,0,x,function_handle_t(0,x),tol*10^-12);
end

for i = 1:max_iter
    DF = DF_func(approx);
    F = function_handle(approx);
    rcond_DF = rcond(DF);
    if abs(rcond_DF)>10^12 || abs(rcond_DF)<10^-12 || isnan(rcond_DF) || isinf(rcond_DF)
        error('RCOND not good, %e',rcond_DF)
    end
    step = - DF \ F ;
    fprintf('Norm of stepsize %e at iteration %i\n',norm(step),i)
    if i>1 && norm(step)<tol 
        return
    end
    zero_approx = approx + step;
    %zero_approx(1:3) = max(zero_approx(1:3),0);
    approx = zero_approx;
end
error('Newton did not converge')