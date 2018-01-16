% Computes the total error of the graph
function Fx = compute_global_error(g)

Fx = 0;

% Loop over all edges
for eid = 1:length(g.edges)
  edge = g.edges(eid);

  % pose-pose constraint
  if (strcmp(edge.type, 'P') != 0)

    x1 = v2t(g.x(edge.fromIdx:edge.fromIdx+2));  % the first robot pose
    x2 = v2t(g.x(edge.toIdx:edge.toIdx+2));      % the second robot pose
    z = v2t(edge.measurement);                   % measurement of the edge

    % calc error between measurement and poses
    % same as (x2 - (x1 + z)) in non-homog space
    err = t2v(inv(z) * (inv(x1) * x2));
    % add error to global error
    Fx = Fx + err' * edge.information * err;
    disp('pose-pose')
  % pose-landmark constraint
  elseif (strcmp(edge.type, 'L') != 0)
    x = g.x(edge.fromIdx:edge.fromIdx+2);  % the robot pose
    l = g.x(edge.toIdx:edge.toIdx+1);      % the landmark
    z = edge.measurement;

    % create diff vector between landmark and robot
    d = l - x(1:2);
    % calculate expecte measurement with range and bearing
    zexp = [sqrt(d(1)^2 + d(2)^2);
            normalize_angle(atan2(d(2), d(1)))];
    % calculate error
    err = z - zexp;
    % add error to global error
    Fx = Fx + err' * edge.information * err;
    disp('pose-landmark')
  end

end
