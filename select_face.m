function full_face = select_face(Sections_1and2,Theta)
% function full_face = select_face(Sections_1and2,Theta)
%
% INPUT 
% Sections_1and2    2x1 integer vector indicating two sections in the phase
%                   space
% Theta             thresholds
% OUTPUT
% full_face         2xN matrix defining the face between sections 1 and 2

Section1 = Sections_1and2(:,1);
Section2 = Sections_1and2(:,2);

% compatibility
if sum(abs(Section1-Section2))~=1
    error('The two sectors are not touching')
end

N = length(Theta);
dim = find(Section1-Section2);

lowest_point = zeros(1,N);
highest_point = zeros(1,N);

for i = 1:N
    if i == dim
        continue
    end
    if Section1(i)-1==0
        inf_theta_i = 0;
    else
        inf_theta_i = Theta{i}(Section1(i)-1);
    end
    if Section1(i) > length(Theta{i})
        sup_theta_i = Inf;
    else
        sup_theta_i = Theta{i}(Section1(i));
    end
    lowest_point(i) = inf_theta_i;
    highest_point(i) = sup_theta_i;
end

lowest_point(dim) = Theta{dim}(min(Section1(dim),Section2(dim)));
highest_point(dim) = lowest_point(dim);


full_face = [lowest_point;highest_point];

