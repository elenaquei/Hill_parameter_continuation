% main 3D
% after the success of main_2D it is time to go to a real world example,
% and consider the repressilator (oscillator based on one repression)

% parameters
lambda = zeros(3,3);
%       lambda - a 3-by-3 array of parameters. Each row corresponds to an edge of the network model as follows:
%                lambda(j,:) = [theta_i, l_i, delta_i]
% to have oscillations, we need l_i < theta_i+1 < u_i
lambda = [1 0.5 1.5
    1 0.5 1.5
    1 0.5 1.5];

% with lambda and knowing the system, we can recreate the vectors of sink

% we follow the construction of main_2D
N = 3;
% thresholds
Theta  = cell(N,1);
for i = 1:N
    Theta{i} = lambda(i,1);
end

Dim_graph = 1;
for i = 1:N
    Dim_graph = Dim_graph*(length(Theta{i})+1);
end


% sinks
V = zeros(Dim_graph,N);
% construction of the sinks
x = [lambda(3,2)+lambda(3,3), lambda(3,2)];
y = [lambda(2,2)+lambda(2,3), lambda(2,2)];
z = [lambda(1,2)+lambda(1,3), lambda(1,2)];

V(1,:) = [x(1),y(1),z(1)];
V(2,:) = [x(2),y(1),z(1)];
V(3,:) = [x(1),y(2),z(1)];
V(4,:) = [x(2),y(2),z(1)];
V(5,:) = [x(1),y(1),z(2)];
V(6,:) = [x(2),y(1),z(2)];
V(7,:) = [x(1),y(2),z(2)];
V(8,:) = [x(2),y(2),z(2)];

%Section = [1 1 1
%    2 1 1 
%    2 2 1
%    2 2 2
%    1 2 2
%    1 2 1]';
Section = [2 2 1
    2 1 1
    2 1 2
    1 1 2
    1 2 2
    1 2 1]';
    

[start_face1,end_face1] = flowing(Section,Theta,V);
if isempty(start_face1)
    error('no orbit found')
end

[B,I] = sort(sum(abs(start_face1 - end_face1),2));

start_point = start_face1(I(1),:);
end_point = end_face1(I(1),:);

y = time_series(start_point, Section,Theta,V);


