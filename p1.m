% Part 1

footplan = dlmread('plan001.txt');
footplan = [0, footplan(1,2:end); footplan]; % Dummy start hold command

cmd_end_times = cumsum(footplan(:,1));
pxs = footplan(:, 2);
pys = footplan(:, 3);

% Set up discrete time nodes
t_s = 0;
t_f = 6.0;
N = 240;

ts = linspace(t_s, t_f, N);
dt = ts(2) - ts(1);
p_xis = interp1(cmd_end_times, pxs, ts, 'next')';
p_yis = interp1(cmd_end_times, pys, ts, 'next')';

G = 10.0;
z = 0.5;

% Set up disciplined convex optimization problem
cvx_begin
    variables x(N) u_x(N) y(N) u_y(N)
    
    vels_x   = (x(2:end) - x(1:end-1)) / dt;
    vels_y   = (y(2:end) - y(1:end-1)) / dt;
    dists_x  = x - p_xis;
    dists_y  = y - p_yis;
    
    J_x = sum(vels_x.^2) + sum(dists_x.^2) + 30*sum(u_x.^2);
    J_y = sum(vels_y.^2) + sum(dists_y.^2) + 30*sum(u_y.^2);
    minimize( J_x + J_y )
    
    subject to
        for i=3:N
            % Enforce discrete time dynamics as equality constraints
            (x(i) - 2*x(i-1) + x(i-2))/dt^2 == (G/z) * (x(i) - p_xis(i) + u_x(i));
            (y(i) - 2*y(i-1) + y(i-2))/dt^2 == (G/z) * (y(i) - p_yis(i) + u_y(i));
        end
cvx_end

% Plots
COP_x = p_xis - u_x;
COP_y = p_yis - u_y;

figure;
plot(ts, p_xis, ts, COP_x, ts, x);
title('X Trajectory');
xlabel('Time (s)');
legend('Foot Location', 'COP', 'COM');
legend('Location', 'SouthEast');

figure;
plot(ts, p_yis, ts, COP_y, ts, y);
title('Y Trajectory');
xlabel('Time (s)');
legend('Foot Location', 'COP', 'COM');
legend('Location', 'SouthEast');
