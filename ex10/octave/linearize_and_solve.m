% performs one iteration of the Gauss-Newton algorithm
% each constraint is linearized and added to the Hessian

function dx = linearize_and_solve(g)

nnz = nnz_of_graph(g);
N = length(g.x);

% allocate the sparse H and the vector b
H = spalloc(N, N, nnz);
b = zeros(1,N);

needToAddPrior = true;

% compute the addend term to H and b for each of our constraints
disp('linearize and build system');
for eid = 1:length(g.edges)
  edge = g.edges(eid);

  % pose-pose constraint
  if (strcmp(edge.type, 'P') != 0)
    % edge.fromIdx and edge.toIdx describe the location of
    % the first element of the pose in the state vector
    % You should use also this index when updating the elements
    % of the H matrix and the vector b.
    % edge.measurement is the measurement
    % edge.information is the information matrix
    i = edge.fromIdx;
    j = edge.toIdx;
    xi = g.x(i:i+2);  % the first robot pose
    xj = g.x(j:j+2);  % the second robot pose
    z  = edge.measurement;
    Omega = edge.information;


    % Computing the error and the Jacobians
    % e the error vector
    % A Jacobian wrt x1
    % B Jacobian wrt x2
    [e, A, B] = linearize_pose_pose_constraint(xi, xj, z);

    % compute and add the term to H and b
    % Hii
    H(i:i+2,i:i+2) += A' * Omega * A;
    % Hij
    H(i:i+2,j:j+2) += A' * Omega * B;
    % Hji
    H(j:j+2,i:i+2) += B' * Omega * A;
    % Hjj
    H(j:j+2,j:j+2) += B' * Omega * B;

    % biT
    b(i:i+2) += e' * Omega * A;
    % bjT
    b(j:j+2) += e' * Omega * B;

    if (needToAddPrior)
        % This fixes one node to remain at its current location
        H(1:3,1:3) += eye(3);
        needToAddPrior = false;
    end

  % pose-landmark constraint
    elseif (strcmp(edge.type, 'L') != 0)
        % edge.fromIdx and edge.toIdx describe the location of
        % the first element of the pose and the landmark in the state vector
        % You should use also this index when updating the elements
        % of the H matrix and the vector b.
        % edge.measurement is the measurement
        % edge.information is the information matrix
        x1 = g.x(edge.fromIdx:edge.fromIdx+2);  % the robot pose
        x2 = g.x(edge.toIdx:edge.toIdx+1);      % the landmark
        Omega = edge.information;
        i = edge.fromIdx;
        j = edge.toIdx;

        % Computing the error and the Jacobians
        % e the error vector
        % A Jacobian wrt x1
        % B Jacobian wrt x2
        [e, A, B] = linearize_pose_landmark_constraint(x1, x2, edge.measurement);

        % TODO: compute and add the term to H and b


    end
end

disp('solving system');

% TODO: solve the linear system, whereas the solution should be stored in dx
% Remember to use the backslash operator instead of inverting H

b = b';
dx = H\-b;

end
