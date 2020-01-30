function start = backward_flow_projection(y,sink, index, theta)
% function start = backward_flow_projection(y,sink, index, theta)
%
% letting the flow starting at start and flowing linearly to sink
% intersect the hyperplane defined by x(index) = theta
% INPUT
% y             element of R^N+
% sink          sink of the linear ODE system 
% index         defining the hyperplane we need to end at
% theta         value of x(index) at the hyperplane
% 
% OUTPUT
% start         point of intersection between the flow started at start and
%               flown backward to the specified hyperplane
start = 0*y;
for i = 1:size(y,1)
    tau = (theta - y(i,index))/(sink(index) - y(i,index));
    start(i,:) = y(i,:) + tau * (sink - y(i,:));
end