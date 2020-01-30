function y = time_series(start_point, Sections,Theta,V)
% function y = time_series(start_point, Sections,Theta,V)
%
% given a starting point, returns a time series of the orbit
%
% INPUT
% start_point   1xN starting point for the orbit (needs to belong to the
%               appropriate face)
% Sections      vector with list of visited sections
% Theta         threshold
% V             sinks
% OUTPUT
% start_face    random point on the initial face that got until the end of
% the periodic orbit
% end_face      forward flow of start_face along the periodic orbit
y = start_point;

for i = 1:size(Sections,2)
    if i<size(Sections,2)
        temp_sections = Sections(:,[i,i+1]);
    else
        temp_sections = Sections(:,[i,1]);
    end
    [index,theta] = get_index_and_theta(Theta,temp_sections);
    sink = V(index_Section(Sections(:,i),Theta),:); 
    y_i = flow_time_series(y(end,:),sink, index, theta);
    y = cat(1,y,y_i);
end