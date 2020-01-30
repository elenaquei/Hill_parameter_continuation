function temp_points = cover_face(selected_face,N_rand)
% function temp_points = cover_face(selected_face,N_rand)
%
% given a bunch of random points in the given face
%
% INPUT
% selected_face         2xN matrix with the two opposite edges (product of
%                       rectangles)
% N_rand                number of random points requested (DEFAULT 1000)
% OUTPUT
% temp_points           N_rand x N points on the given face

if nargin == 1
    N_rand = 10^4;
end

if any(any(selected_face ==Inf))
    selected_face = min(selected_face, 10^3);
end

size_face = abs(selected_face(1,:) - selected_face(2,:));
center_face = (selected_face(1,:) + selected_face(2,:))/2;

temp_points = (rand(N_rand, size(selected_face,2))-0.5) .* repmat(size_face,N_rand,1) + center_face;
end