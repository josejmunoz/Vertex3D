function [V_new,T_new] = flipnm_chatGPT(V,T,n,m)
% V: Nx3 matrix of vertex coordinates
% T: Mx4 matrix of tetrahedron connectivity
% n,m: indices of the two tetrahedra to flip
% V_new: updated vertex coordinates
% T_new: updated tetrahedron connectivity

% Extract the vertex indices for the two tetrahedra to flip
n_vertices = T(n,:);
m_vertices = T(m,:);

% Find the common face between the two tetrahedra
common_vertices = intersect(n_vertices, m_vertices);

% Find the two additional vertices that form the quadrilateral
additional_vertices = setdiff(union(n_vertices, m_vertices), common_vertices);

% Make sure there are only two additional vertices
assert(length(additional_vertices) == 2, 'Error: Invalid tetrahedron connectivity.');

% Find the opposite vertices of the quadrilateral
opposite_vertices = setdiff(common_vertices, additional_vertices);

% Create the two new tetrahedra
T_new = [n_vertices(opposite_vertices), m_vertices(additional_vertices(1)), m_vertices(additional_vertices(2)), n_vertices(additional_vertices(1))];
T_new = [T_new; n_vertices(opposite_vertices), m_vertices(additional_vertices(1)), m_vertices(additional_vertices(2)), n_vertices(additional_vertices(2))];

% Update the vertex coordinates
V_new = V;
for i = 1:2
    tet_vertices = T_new(i,:);
    V_new(tet_vertices,:) = circumcenter(V(tet_vertices,:));
end

end

function center = circumcenter(vertices)
% Compute the circumcenter of a tetrahedron given its vertex coordinates
d = size(vertices, 2);
if d == 2
    center = mean(vertices);
elseif d == 3
    A = [vertices(2,:) - vertices(1,:); vertices(3,:) - vertices(1,:); vertices(4,:) - vertices(1,:)];
    b = 0.5 * [sum(vertices(2,:).^2 - vertices(1,:).^2); sum(vertices(3,:).^2 - vertices(1,:).^2); sum(vertices(4,:).^2 - vertices(1,:).^2)];
    center = A \ b;
else
    error('Error: Invalid vertex coordinates.');
end
end

