mu_A = 4.892e-9;                      % [km3/s2]   (Asteroid gravitational parameter)
T_A = 4.296057 * 3600;                % [s]        (Asteroid rotation period)
w_A = 2*pi/T_A;                       % [rad/s]    (Asteroid rotation rate)
r_S = [1.5e8 0 0]';                   % [km]       (Inertial Sun position w.r.t. asteroid center (1 AU)
SPC = 1e14;                           % [kgkm/s2]  (Solar pressure constant)
rho = 0.4;                            %            (Coefficient of reflectivity)
Aoverm = (1/62)*1e-6;                 % [km2/kg]   (Area-to-mass ratio)
sigma_w = 1e-9;                       % [km/s2]    (Process noise standard deviation)
f = 2.089795918367347e3;              % [pixels]   (Camera focal length)
u0 = 512; v0 = 512;                   % [pixels]   (Camera principal point coordinates)
umin = 0; vmin = 0;                   % [pixels]   (Min. pixel coordinates)
umax = 1024; vmax = 1024;             % [pixels]   (Max. pixel coordinates)
sigma_u = 0.25; sigma_v = 0.25;       % [pixels]   (Measurement standard deviation)
t0 = 0.0;                             % [s]        (Initial epoch)
tf_int = 432000.0;                    % [s]        (Final epoch for integration)
dt_int = 60.0;                        % [s]        (Integration time step)
tf_obs = 259200;                      % [s]        (Final epoch for observation)
dt_obs = 600.0;                       % [s]        (Observation time step)
r0 = [0 -1 0]';                       % [km]       (Nominal initial position)
rdot0 = [0 0 sqrt(mu_A/norm(r0))]';   % [km/s]     (Nominal initial velocity)