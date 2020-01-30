function [start_face] = backward_flow(Sections,Theta,V)
% function [start_face] = backward_flow(Sections,Theta,V)
% 
% INPUT
% Sections      vector with list of visited sections
% Theta         threshold
% V             sinks
% OUTPUT
% start_face    2xN matrix defining the hyper-rectangle with initial points

% define the full initial face
temp_face = select_face(Sections(:,[end,1]),Theta);

for i = 1:size(Sections,2)
    if i>1 && i <size(Sections,2)
        [index_Section(Sections(:,end-i+1),Theta)
        index_Section(Sections(:,end-i+2),Theta)]
        temp_sections = Sections(:,[end-i+2,end-i+1]);
    elseif i ==1
        [index_Section(Sections(:,end),Theta)
        index_Section(Sections(:,1),Theta)]
        temp_sections = Sections(:,[end,1]);
    else
        [index_Section(Sections(:,1),Theta)
        index_Section(Sections(:,2),Theta)]
        temp_sections = Sections(:,[2,1]);
    end
    [index,theta] = get_index_and_theta(Theta,temp_sections);
    sink = V(index_Section(Sections(:,end-i+2),Theta),:); 
    backward_flow_face = backward_flow_projection(temp_face,sink, index, theta);
    temp_face = intersection(backward_flow_face,select_face(temp_sections,Theta));
    if temp_face == 0 
        break
    end
end
start_face = temp_face;