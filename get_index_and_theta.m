function [index,theta] = get_index_and_theta(Theta,Sections_1and2)
% function [index,theta] = get_index_and_theta(Theta,Sections_1and2)
% 
% INPUT 
% Theta             cell(N) with thresholds
% Sections_1and2    2xN indication of two sections
% OUTPUT
% index             integer, which dimension is fixed at the common face
% theta             value of threshold at the face

Section1 = Sections_1and2(:,1);
Section2 = Sections_1and2(:,2);

% compatibility
if sum(abs(Section1-Section2))~=1
    error('The two sectors are not touching')
end

index = find(Section1-Section2);

theta = Theta{index}(min(Section1(index),Section2(index)));