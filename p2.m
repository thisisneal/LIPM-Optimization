% Part 2 : Optimize ?
% Neal Bhasin & Rick Shanor
%
% Note: Requires CVX installation

footplan = dlmread('plan001.txt');
footplan = [0, footplan(1,2:end); footplan]; % Dummy start hold command

timings_0 = cumsum(footplan(:,1));

% Global constants for use in criterion function
global pxs; 
global pys;
global ts;
global dt;
global t_s;
global t_f;
global N;
global G;
global z;

pxs = footplan(:, 2);
pys = footplan(:, 3);

% Set up discrete time nodes
t_s = 0;
t_f = 6.0;
N = 300;

ts = linspace(t_s, t_f, N);
dt = ts(2) - ts(1);

G = 10.0;
z = 0.75;

algo = 'interior-point';
lb = zeros(size(timings_0,1), 1);
lb(end) = t_f;
ub = t_f * ones(size(timings_0,1), 1);
ub(1) = 0;
options = optimoptions(@fmincon,'Algorithm',algo,'Display','iter', ...
                       'MaxFunEvals', 1000, 'MaxIter', 2000, 'DiffMinChange', dt);
[opt_times,fval,exitflag]=fmincon(@cvx_opti,timings_0,[],[],[],[],lb,ub,@cvx_cstrs,options);

[ ~, x, y, u_x, u_y, p_xis, p_yis ] = cvx_opti( opt_times );

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

figure;
subplot(1,2,1), plot(ts, u_x);
title('X Input');
xlabel('Time (s)');
subplot(1,2,2), plot(ts, u_y);
title('Y Input');
xlabel('Time (s)');

figure;
subplot(1,2,1), plot(ts(2:end), diff(x));
title('X Velocity');
xlabel('Time (s)');
subplot(1,2,2), plot(ts(2:end), diff(y));
title('Y Velocity');
xlabel('Time (s)');
