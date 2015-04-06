% Objective function

function [ cost, x, y, u_x, u_y, p_xis, p_yis ] = cvx_opti( cmd_end_times )

cmd_end_times

global pxs;
global pys;
global ts;
global dt;
global t_s;
global t_f;
global N;
global G;
global z;

% Make fmincon obey hard constraints
ceq = [cmd_end_times(1) - t_s;
       cmd_end_times(end) - t_f];
c = diff( cmd_end_times);
if(any(c < 0) || any(ceq)) 
    cost = 500;
    return;
end

Fx = griddedInterpolant(cmd_end_times, pxs, 'next');
p_xis = Fx(ts)';
Fy = griddedInterpolant(cmd_end_times, pys, 'next');
p_yis = Fy(ts)';

% Set up disciplined convex optimization problem
cvx_begin quiet
    variables x(N) u_x(N) y(N) u_y(N)

    vels_x   = (x(2:end) - x(1:end-1)) / dt;
    vels_y   = (y(2:end) - y(1:end-1)) / dt;
    dists_x  = x - p_xis;
    dists_y  = y - p_yis;
    %Change = sum((ts - tInit').^2);
    
    J_x = sum(vels_x.^2) + sum(dists_x.^2) + 30*sum(u_x.^2);
    J_y = sum(vels_y.^2) + sum(dists_y.^2) + 30*sum(u_y.^2);
    acc_x = (vels_x(2:end) - vels_x(1:end-1)) / dt;
    J_const_vel_x = sum(acc_x.^2);
    acc_y = (vels_y(2:end) - vels_y(1:end-1)) / dt;
    J_const_vel_y = sum(acc_y.^2);
    minimize( J_x + J_y + J_const_vel_x + J_const_vel_y)
    
    subject to
        for i=3:N
            % Enforce discrete time dynamics as equality constraints
            (x(i) - 2*x(i-1) + x(i-2))/dt^2 == (G/z) * (x(i) - p_xis(i) + u_x(i));
            (y(i) - 2*y(i-1) + y(i-2))/dt^2 == (G/z) * (y(i) - p_yis(i) + u_y(i));
        end
        % Velocity endpoint constraints
        vels_x(1) == 0;
        vels_x(end) == 0;
        vels_y(1) == 0;
        vels_y(end) == 0;
cvx_end

cost = cvx_optval

end

