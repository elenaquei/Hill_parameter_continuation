function [f3,cancelled_indeces] = intersection(f1,f2)
% function intersection(f1,f2)
% 
% given the hypersurface defined by f1 and f2 (a hyperplane in R^3+),
% consider their intersection 
%
% INPUT 
% f1        random points
% f2        2xN trivial hyper-rectangles (along the axes)
% OUTPUT
% f3        selected points of f1 in face f2
% cancelled_indeces     indeces of all the elements outside of f2

size_face = abs(f2(1,:) - f2(2,:))/2;
center_face = (f2(1,:) + f2(2,:))/2;

Center_face = repmat(center_face, size(f1,1),1);
Size_face = repmat(size_face, size(f1,1),1);

cancelled_indeces = find( sum(abs(f1 - Center_face)>Size_face,2),size(f1,1));

f1(cancelled_indeces,:) = [];
f3 = f1;
