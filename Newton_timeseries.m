function x_new = Newton_timeseries(vector_field, x_in, h)
% function Newton_timeseries(vector_field, x, h)
% 
% look for a periodic orbit near the approximate periodic orbit x provided,
% with the given vector field.
%
% INPUT
% vector_field      handle, \dot x = f(x), takes input in R^n and return in
%                   the same size
% x                 approximate periodic orbit, x\in R^(nxm) x(:,j)
% approximate solution at time hj
% h                 time stepping used for x
% OUTPUT
% x_new             better approximation

% size management
reshape_loc = @(x) reshape(x,size(x_in));
squeeze_function = @(x) x(:);
index_minus1 = @(x) [x(:,end),x(:,1:end-1)];

% Crank Nicolson RK method 
%RK = @(x1,x2) x1 + h/2*(vector_field(x1)+vector_field(x2)); 
RK = @(x1,x2) x1 + h/2*(vector_field(x1)+vector_field(x2)); 

% periodic solution must satisfy this system of equations
zero_finding = @(x) x - RK(index_minus1(x),x);

% required shape
zero_finding_shape = @(x) squeeze_function(zero_finding(reshape_loc(x)));

x_new_reshape = Newton_handle(zero_finding_shape, x_in(:));

x_new = reshape(x_new_reshape,size(x_in));

end