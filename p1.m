footplan = dlmread('plan001.txt');
footplan = [0, footplan(1,2:end); footplan]; % Dummy start hold command

cmd_end_times = cumsum(footplan(:,1));
pxs = footplan(:, 2);
pys = footplan(:, 3);

% Set up discrete time nodes
t_s = 0;
t_f = 6.0;
N = 100;

ts = linspace(t_s, t_f, N);
dt = ts(2) - ts(1);
p_xis = interp1(cmd_end_times, pxs, ts, 'next')';
p_yis = interp1(cmd_end_times, pys, ts, 'next')';

G = 10.0;
z = 5;

cvx_begin
    variables x(N) u_x(N)
    
    vels   = (x(2:end) - x(1:end-1)) / dt;
    dists  = x - p_xis;
    
    minimize( sum(vels.^2) + sum(dists.^2) + 30*sum(u_x.^2))
    
    subject to
        for i=3:N
            (x(i) - 2*x(i-1) + x(i-2))/dt.^2 == (G/z) * (x(i) - p_xis(i) + u_x(i));
        end
cvx_end

COP = u_x + p_xis;

figure;
title('X Trajectory');
xlabel('Time (s)');
plot(ts, p_xis, ts, COP, ts, x);
legend('Foot Location', 'COP', 'COM');
