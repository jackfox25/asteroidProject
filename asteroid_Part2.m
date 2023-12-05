asteroid_Part1;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% PART 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
Q = [dt_int^3/3*eye(3)    dt_int^2/2*eye(3);...
     dt_int^2/2*eye(3)    dt_int*eye(3)];

P0 = 10000*eye(6);

R = diag([sigma_u^2 sigma_v^2]);


