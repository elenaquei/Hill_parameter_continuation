function y_i = flow_time_series(start_point,sink, index, theta)
% function y_i = flow_time_series(start_point,sink, index, theta)
%
% letting the flow starting at start and flowing linearly to end_point
% intersect the hyperplane defined by x(index) = theta
% INPUT
% start         element of R^(MxN)+
% end_point     sink of the linear ODE system 
% index         defining the hyperplane we need to end at
% theta         value of x(index) at the hyperplane
% 
% OUTPUT
% y             points of intersection between the flow started at start and
%               the specified hyperplane


tau = (theta -sink(index) ) /  (-sink(index) +start_point(index));
y_i = repmat(sink,length((0:0.001:-log(tau))),1) + exp(-(0:0.001:-log(tau)))' * (start_point - sink);

end