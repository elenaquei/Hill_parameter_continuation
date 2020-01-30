function bool = is_empty_face(face)
% function bool = is_empty_face(face)
%
% INPUT
% face      2xN vector, indicating the two end points of a N-rectangle
% OUTPUT
% bool      1 if face is empty, 0 otherwise

if any(face(1,:)>face(2,:))
    bool = 1;
else
    bool = 0;
end