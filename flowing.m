function [start_face,end_face] = flowing(Sections,Theta,V)
% function [start_face,end_face] = flowing(Sections,Theta,V)
% 
% INPUT
% Sections      vector with list of visited sections
% Theta         threshold
% V             sinks
% OUTPUT
% start_face    random point on the initial face that got until the end of
% the periodic orbit
% end_face      forward flow of start_face along the periodic orbit

selected_face = select_face(Sections(:,[end,1]),Theta);
temp_points = cover_face(selected_face);
start_points = temp_points;

for i = 1:size(Sections,2)
    if i<size(Sections,2)
        temp_sections = Sections(:,[i,i+1]);
    else
        temp_sections = Sections(:,[i,1]);
    end
    [index,theta] = get_index_and_theta(Theta,temp_sections);
    sink = V(index_Section(Sections(:,i),Theta),:); 
    forward_flow_face = flow_projection(temp_points,sink, index, theta);
    [temp_points,index_void] = intersection(forward_flow_face,select_face(temp_sections,Theta));
    start_points(index_void,:)=[];
    if size(temp_points,1) <= 0 
        break
    end
end
start_face = start_points;
end_face = temp_points;