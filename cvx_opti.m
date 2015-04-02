function [ cost ] = cvx_opti( cmd_end_times )

global p_xis;
global p_yis;
global ts;
global dt;
global G;
global z;
% global initCmdTimes;

cmd_end_times = cumsum(footplan(:,1));
pxs = footplan(:, 2);
pys = footplan(:, 3);

% Set up discrete time nodes
t_s = 0;
t_f = 6;
N = 249;

p_xis = interp1(cmd_end_times, pxs, ts, 'next')';
p_yis = interp1(cmd_end_times, pys, ts, 'next')';
  
end

% Set up disciplined convex optimization problem
cvx_begin
    variables x(N) u_x(N) y(N) u_y(N)

    G = 10.0;
    z = 0.5;
    
    vels_x   = (x(2:end) - x(1:end-1)) / dt;
    vels_y   = (y(2:end) - y(1:end-1)) / dt;
    dists_x  = x - p_xis;
    dists_y  = y - p_yis;
    %Change = sum((ts - tInit').^2);
    
    J_x = sum(vels_x.^2) + sum(dists_x.^2) + 30*sum(u_x.^2);
    J_y = sum(vels_y.^2) + sum(dists_y.^2) + 30*sum(u_y.^2);
    minimize( J_x + J_y)
    
    subject to
        for i=3:N
            % Enforce discrete time dynamics as equality constraints
            (x(i) - 2*x(i-1) + x(i-2))/dt^2 == (G/z) * (x(i) - p_xis(i) + u_x(i));
            (y(i) - 2*y(i-1) + y(i-2))/dt^2 == (G/z) * (y(i) - p_yis(i) + u_y(i));
        end
cvx_end

cost = cvx_optval;

end

