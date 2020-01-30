function index = index_Section(Section_vec,cell_Theta)
% function index = index_Section(Section_vec,cell_Theta)
%
% transfor vector indices into float indices
% INPUT
% Section_vec       vector indicator of a section
% cell_Theta        cell of the dimension of the threshold
% OUTPUT
% index             index of hte Section at hand

if length(Section_vec)~= length(cell_Theta)
    error('Dimensions do not correspond')
end
N = length(Section_vec);
a = cell(N,1);

Dim_graph = zeros(N,1);
for i = 1:N
    Dim_graph(i) = length(cell_Theta{i})+1;
    a{i} = Section_vec(i);
end

[index] = sub2ind(Dim_graph,a{:});
