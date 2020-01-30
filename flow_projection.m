function y = flow_projection(start,end_point, index, theta)
% function y = flow_projection(start,end_point, index, theta)
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

y = 0*start;
for i = 1:size(start,1)
    tau = (theta -start(i,index) ) /  (end_point(index) - start(i,index));
    if tau<0
        error('Need negative time to continue')
    end
    y(i,:) = start(i,:) + tau * ( end_point(:)' -start(i,:));
end
end