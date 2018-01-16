% Computes the total error of the graph
function Fx = compute_global_error(g)

Fx = 0;

% Loop over all edges
for eid = 1:length(g.edges)
    edge = g.edges(eid);

    % pose-pose constraint
    if (strcmp(edge.type, 'P') != 0)

        x1 = g.x(edge.fromIdx:edge.fromIdx+2);  % the first robot pose
        x2 = g.x(edge.toIdx:edge.toIdx+2);      % the second robot pose
        z  = edge.measurement;                  % measurement of the edge

        % calc error between measurement and poses
        [err, _, _] = linearize_pose_pose_constraint(x1, x2, z);

        % add error to global error
        Fx = Fx + err' * edge.information * err;
    % pose-landmark constraint
    elseif (strcmp(edge.type, 'L') != 0)
        x = g.x(edge.fromIdx:edge.fromIdx+2);  % the robot pose
        l = g.x(edge.toIdx:edge.toIdx+1);      % the landmark
        z = edge.measurement;                  % measurement of lm position

        % calculate error
        [err, _, _] = linearize_pose_landmark_constraint(x, l, z);
        % add error to global error
        Fx = Fx + err' * edge.information * err;
    end

end
