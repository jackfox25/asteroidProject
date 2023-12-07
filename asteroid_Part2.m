asteroid_Part1;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% PART 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

NMC = 1; % # of Monte-Carlo simulations


% -- Process noise covariance matrix 
 Q = sigma_w^2*[dt_int^3/3*eye(3)    dt_int^2/2*eye(3);...
                dt_int^2/2*eye(3)    dt_int*eye(3)];
% -- Initial state covariance matrix
 P0 = 0.1^2*eye(6);
% -- Measurement noise covariance matrix (for each landmark)
 R = diag([sigma_u^2 sigma_v^2]);


% Pre-compute Fbar matrices for DT state propagation
Fbar_k = zeros(6,6,length(tspan)-1);
for k = 1:length(tspan)-1    
    Abar_k = linearizedAmat(mu_A, x_nom(k,:)');      % CT
    Fbar_k(:,:,k) = linearizedFmat(dt_int, Abar_k);  % ... to DT, save
end


% ========================= MONTE-CARLO LOOP =========================
for mc=1:NMC

   % Configure tolerances
    tol=1.e-8;
    OPTIONS = odeset('RelTol',3*tol,'AbsTol',tol);

   % ---------------------------------------------------------
   % Generate noisy truth state
   % ---------------------------------------------------------
    x_noisy = zeros(length(tspan),6);
    x_noisy(1,:) = x0';
    tic;
    for t_k = t0+dt_int:dt_int:tf_int
        processNoise = sigma_w^2*randn(3,1);
        [t_noisy_k,x_noisy_k] = ode45(@(t,x) satDynamicModel(t,x,[],processNoise),[t_k-dt_int t_k],x_noisy(t_k/dt_int,:)',OPTIONS); 
        x_noisy(t_k/dt_int+1,:) = x_noisy_k(end,:);
    end
    toc;
   % --------------------------------------------------------- 
   % Generate noisy measurements from noisy truth state
   % ---------------------------------------------------------
    y_noisy_tbl = [];
    for t_k = obsvTimes % For each observation timestep
        
       % Recover s/c position vector, define khatC vector (cam nadir pointing)
        satpos_k_N = x_noisy((tspan==t_k)~=0,1:3)';  % satellite position in inertial frame
        
       % Get rotation matrix R_CtoN and extract unit vectors
        R_CtoN_k = R_CtoN(:,:,(t_k/dt_obs)+1);
        ihatC_k_N = R_CtoN_k(:,1); jhatC_k_N = R_CtoN_k(:,2); khatC_k_N = R_CtoN_k(:,3);
    
       % Get rotation matrix R_AtoN
        theta = w_A*t_k;
        R_AtoN_k = rotZ(theta);
    
        for i = 1:Nlmks % For each landmark
           
           % Get landmark position in intertial frame
            lmkipos_k_A = pos_lmks_A(:,i);
            lmkipos_k_N = R_AtoN_k*lmkipos_k_A;
            disttolmk_k_N = lmkipos_k_N - satpos_k_N;
    
           % Simulate ideal measurement (WITH NOISE)
            u_i = f * (disttolmk_k_N'*ihatC_k_N) / (disttolmk_k_N'*khatC_k_N) + u0 + sigma_u^2*randn(1);
            v_i = f * (disttolmk_k_N'*jhatC_k_N) / (disttolmk_k_N'*khatC_k_N) + v0 + sigma_v^2*randn(1);
    
           % Check whether landmark is in camera FOV
            if u_i >= 0 && u_i <= umax && v_i >= 0 && v_i <= vmax && (disttolmk_k_N'*khatC_k_N) > 0
               % Check whether landmark is facing satellite
                if (lmkipos_k_N'*khatC_k_N) < 0
                   % If so, append [timestamp  lmk_id  u   v] to y_table_ideal
                    y_noisy_tbl = [y_noisy_tbl;...
                                   t_k i u_i v_i]; %#ok<*AGROW>   <--- suppresses size change warning
                end
            end
        end
    end
    


   % ---------------------------------------------------------------------------------------------------------------- 
   % LINEARIZED KALMAN FILTER
   % ----------------------------------------------------------------------------------------------------------------
    
   % Containers 
    x_bar_LKF = zeros(length(tspan)-1,6); 
    P = zeros(6,6,length(tspan)-1);
    
   % Set initial conditions
    x_bar_LKF(1,:) = x_bar(1,:); 
    P(:,:,1) = P0; % Use the provided initial covariance matrix P0

    for k = 1:length(tspan)-1
        t_k = tspan(k);
    
       % Time Update / Prediction Step occurs at every k (dt = 60s)
        x_bar_LKF(k+1,:) = Fbar_k(:,:,k) * x_bar_LKF(k,:)';
        P(:,:,k+1) = Fbar_k(:,:,k) * P(:,:,k) * Fbar_k(:,:,k)' + Q;
        
       % The measurement update occurs at each observation epoch (dt = 600s) 
        if mod(t_k,dt_obs) == 0 && t_k <= y_noisy_tbl(end,1)
           % Get ideal measurements at time k
            y_full_k = y_full((y_full(:,1)==t_k)~=0,:);
    
           % Get actual measurements at time k
            y_noisy_k = y_noisy_tbl((y_noisy_tbl(:,1)==t_k)~=0,:);
    
           % Truncate ideal measurements to only include landmarks which appear in the actual measurements
            y_trunc_k = []; 
            for row = 1:size(y_noisy_k,1)
                lmkid = y_noisy_k(row,2);
                y_trunc_k =[y_trunc_k; y_full_k((y_full_k(:,2) == lmkid),:)];   
            end 
            
           % Calculate and reshape y_bar_k
            y_bar_k = y_noisy_k(:,3:4) - y_trunc_k(:,3:4);  
            y_bar_k = reshape(y_bar_k',[numel(y_bar_k),1]);

           % Get rotation matrix R_CtoN
            R_CtoN_k = R_CtoN(:,:,(t_k/dt_obs)+1);
        
           % Get rotation matrix R_AtoN
            theta = w_A*t_k;
            R_AtoN_k = rotZ(theta);
            pos_lmks_N = R_AtoN_k*pos_lmks_A; 
    
           % Calculate measurement Jacobian and discretize
            Cbar_k = linearizedCmat(f, R_CtoN_k, pos_lmks_N, y_trunc_k, x_nom_k);
            Hbar_k = linearizedHmat(Cbar_k);
    
           % Compute Kalman gain 
            Big_R = diagTile(R,size(Hbar_k,1)/2); 
            Kkp1 = P(:,:,k+1) * Hbar_k' / (Hbar_k * P(:,:,k+1) * Hbar_k' + Big_R);
                    
           % Measurement Update / Correction Step       
            x_bar_LKF(k+1,:) = x_bar_LKF(k+1,:)' + Kkp1 * (y_bar_k - Hbar_k*x_bar_LKF(k+1,:)'); % Update state estimate
            P(:,:,k+1) = (eye(6) - Kkp1 * Hbar_k) * P(:,:,k+1); % Update covariance matrix
            
        end 
    end

end


if PLOTFLAG
    % Position error and +-2sigma (single plot)
    f=figure;
    subplot(3,1,1);
    plot(tspan,x_bar_LKF(:,1),'DisplayName','$\bar{x}$','Color',mlc(1),'LineWidth',2); hold on;
    plot(tspan,2*sqrt(squeeze(P(1,1,:))),'k--','DisplayName','$\pm2\sigma$','LineWidth',1.5);
    plot(tspan,-2*sqrt(squeeze(P(1,1,:))),'k--','LineWidth',1.5,'HandleVisibility','off');
    labels(gca,{'Time [s]','$\bar{x}$ [km]'},'');
    set(gca,'YLim',[-0.002 0.002]); legend('interpreter','latex');
    
    subplot(3,1,2);
    plot(tspan,x_bar_LKF(:,2),'DisplayName','$\bar{y}$','Color',mlc(2),'LineWidth',2); hold on;
    plot(tspan,2*sqrt(squeeze(P(2,2,:))),'k--','DisplayName','$\pm2\sigma$','LineWidth',1.5);
    plot(tspan,-2*sqrt(squeeze(P(2,2,:))),'k--','LineWidth',1.5,'HandleVisibility','off');
    labels(gca,{'Time [s]','$\bar{y}$ [km]'},'');
    set(gca,'YLim',[-0.006 0.006]); legend('interpreter','latex');
    
    subplot(3,1,3);
    plot(tspan,x_bar_LKF(:,3),'DisplayName','$\bar{z}$','Color',mlc(3),'LineWidth',2); hold on;
    plot(tspan,2*sqrt(squeeze(P(3,3,:))),'k--','DisplayName','$\pm2\sigma$','LineWidth',1.5);
    plot(tspan,-2*sqrt(squeeze(P(3,3,:))),'k--','LineWidth',1.5,'HandleVisibility','off');
    labels(gca,{'Time [s]','$\bar{z}$ [km]'},'');
    set(gca,'YLim',[-0.006 0.006]); legend('interpreter','latex');
    fixfig(f,false);
    
    % Velocity error and +-2sigma (single plot)
    f=figure;
    subplot(3,1,1);
    plot(tspan,x_bar_LKF(:,4),'DisplayName','$\bar{\dot{x}}$','Color',mlc(1),'LineWidth',2); hold on;
    plot(tspan,2*sqrt(squeeze(P(4,4,:))),'k--','DisplayName','$\pm2\sigma$','LineWidth',1.5);
    plot(tspan,-2*sqrt(squeeze(P(4,4,:))),'k--','LineWidth',1.5,'HandleVisibility','off');
    labels(gca,{'Time [s]','$\bar{\dot{x}}$ [km/s]'},'');
    set(gca,'YLim',[-0.000001 0.000001]); legend('interpreter','latex');
    
    subplot(3,1,2);
    plot(tspan,x_bar_LKF(:,5),'DisplayName','$\bar{\dot{y}}$','Color',mlc(2),'LineWidth',2); hold on;
    plot(tspan,2*sqrt(squeeze(P(5,5,:))),'k--','DisplayName','$\pm2\sigma$','LineWidth',1.5);
    plot(tspan,-2*sqrt(squeeze(P(5,5,:))),'k--','LineWidth',1.5,'HandleVisibility','off');
    labels(gca,{'Time [s]','$\bar{\dot{y}}$ [km/s]'},'');
    set(gca,'YLim',[-0.000001 0.000001]); legend('interpreter','latex');
    
    subplot(3,1,3);
    plot(tspan,x_bar_LKF(:,6),'DisplayName','$\bar{\dot{z}}$','Color',mlc(3),'LineWidth',2); hold on;
    plot(tspan,2*sqrt(squeeze(P(6,6,:))),'k--','DisplayName','$\pm2\sigma$','LineWidth',1.5);
    plot(tspan,-2*sqrt(squeeze(P(6,6,:))),'k--','LineWidth',1.5,'HandleVisibility','off');
    labels(gca,{'Time [s]','$\bar{\dot{z}}$ [km/s]'},'');
    set(gca,'YLim',[-0.000001 0.000001]); legend('interpreter','latex');
    fixfig(f,false);
end
