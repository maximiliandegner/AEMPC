% constraint_def.m
x_1min=0.03;
x_1max=0.25;%4;%;
x_2min=0.03;
x_2max=0.25;
x_3min=0.03;
x_3max=0.25;
u_min=0.049;
u_max=0.449;
b_x = [x_1max; x_2max; x_3max; x_1min; x_2min; x_3min];
b_u = [u_max; u_min];