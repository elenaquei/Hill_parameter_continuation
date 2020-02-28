function zero_approx = Newton_handle(function_handle, approx, DF_func)

max_iter = 50;
tol = 10^-9;

if nargin<3 || isempty(DF_func)
    DF_func = @(x) numjac_loc(function_handle,x,10^-10);
end

for i = 1:max_iter
    DF = DF_func(approx);
    F = function_handle(approx);
    rcond_DF = rcond(DF);
    if abs(rcond_DF)>10^12 || abs(rcond_DF)<10^-12 || isnan(rcond_DF) || isinf(rcond_DF)
        error('At iteration %i, RCOND not good, %e ',i,rcond_DF)
    end
    step = - DF \ F ;
    %fprintf('Norm of stepsize %e at iteration %i\n',norm(step),i)
    if i>1 && norm(step)<tol 
        return
    end
    zero_approx = approx + step;
    zero_approx = max(zero_approx,10^-6);
    approx = zero_approx;
end
error('Newton did not converge')