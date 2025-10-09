% constraint_def.m
x_1min=x_min;
x_1max=x_max;
x_2min=x_min;
x_2max=x_max;
x_3min=x_min;
x_3max=x_max;
u_min=u_min;
u_max=u_max;
b_x = [x_1max; x_2max; x_3max; x_1min; x_2min; x_3min];
b_u = [u_max; u_min];