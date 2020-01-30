function zero_approx = Newton_handle(function_handle, approx, DF_func)

max_iter = 100;
tol = 10^-7;

if nargin<3 || isempty(DF_func)
    function_handle_t = @(t,x) function_handle(x);
     DF_func = @(x) numjac(function_handle_t,0,x,function_handle_t(0,x),tol*10^-2);
end

for i = 1:max_iter
    DF = DF_func(approx);
    F = function_handle(approx);
    rcond_DF = rcond(DF);
    if abs(rcond_DF)>10^6 || abs(rcond_DF)<10^-6 || isnan(rcond_DF) || isinf(rcond_DF)
        error('RCOND not good')
    end
    step = - DF \ F ;
    if i>1 && norm(step)<tol 
        return
    else
        %disp(norm(step))
    end
    zero_approx = approx + step;
    approx = zero_approx;
end
error('Newton did not converge')