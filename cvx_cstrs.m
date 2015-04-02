%
function[c, ceq] = cvx_cstrs( cmd_end_times )

global t_s;
global t_f;

ceq = [cmd_end_times(1) - t_s;
       cmd_end_times(end) - t_f];
c = -diff( cmd_end_times);

end