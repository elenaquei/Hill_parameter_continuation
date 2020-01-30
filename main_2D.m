% DATA NEEDED
% play with the 2D with 2,1 thresholds to start with

N = 2;
% thresholds
Theta  = cell(N,1);
Theta{1} = [1,2];
Theta{2} = 1;

Dim_graph = 1;
for i = 1:N
    Dim_graph = Dim_graph*(length(Theta{i})+1);
end

% sinks
V = zeros(Dim_graph,N);
V(1,:) = [2,-10];
V(2,:) = [3,0.5];
V(3,:) = [2.5,5];
V(4,:) = [-1.7,0];
V(5,:) = [0,1];
V(6,:) = [1,2];

Section1 = [1,2,2,1
    1,1,2,2];
Section2 = [1,2,3,3,2,1
            1,1,1,2,2,2];

plot(V(:,1),V(:,2),'*')
hold on

for j= 1:length(Theta{1})
    plot([Theta{1}(j),Theta{1}(j)],[0,3])
end
for j= 1:length(Theta{2})
    plot([0,3],[Theta{2}(j),Theta{2}(j)])
end

[start_face1,end_face1] = flowing(Section2,Theta,V);
%[start_face2,end_face2] = flowing(Section2,Theta,V);

[B,I] = sort(abs(start_face1(:,1) - end_face1(:,1)));
start_point = start_face1(I(1),:);
end_point = end_face1(I(1),:);

y = time_series(start_point, Section2,Theta,V);
% creat the whole length time series 
plot(y(:,1),y(:,2),'*-')


