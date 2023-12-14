% asteroid_Part1;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% PART 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

    NMC = 1;            % # of Monte-Carlo simulations (MAX 25)
    PLOTFLAG = true;    % decide whether to produce plots or not
    load('x_noisy_MC25.mat')
    load('x_noisy_MC3.mat');
    
   % -- Process noise covariance matrix
     Qfactor = 5e5;
     Q = Qfactor * sigma_w^2*[dt_int^3/3*eye(3)    dt_int^2/2*eye(3);...
                              dt_int^2/2*eye(3)    dt_int*eye(3)];

   % -- Initial state covariance matrix
     P0pos = (0.1)^2*eye(3);
     P0vel = (1.e-3)^2*eye(3);
     P0 = blkdiag(P0pos,P0vel);

   % -- Measurement noise covariance matrix (for each landmark)
     R = 100*diag([sigma_u^2 sigma_v^2]);
    
    
   % Pre-compute Fbar matrices for DT state propagation
    Fbar_k = zeros(6,6,length(tspan));
    for k = 1:length(tspan)    
        Abar_k = linearizedAmat(mu_A, x_nom(k,:)');      % CT
        Fbar_k(:,:,k) = linearizedFmat(dt_int, Abar_k);  % ... to DT, save
    end
    
   % Container for NEES/NIS
    NEES = zeros(length(tspan)-1,NMC);
    NIS = zeros(length(tspan)-1,NMC);

    NEESr1 = zeros(length(tspan)-1,NMC);
    NEESr2 = zeros(length(tspan)-1,NMC);
    NISr1 = zeros(length(tspan)-1,NMC);
    NISr2 = zeros(length(tspan)-1,NMC);
    
    alpha = 0.05;


   % ========================= MONTE-CARLO LOOP =========================
    for mc=1:NMC
    
       % Configure tolerances
        tol=1.e-8;
        OPTIONS = odeset('RelTol',tol,'AbsTol',tol);
    
       % ---------------------------------------------------------
       % Generate noisy truth state
       % ---------------------------------------------------------
        x_noisy = x_noisy_MC(:,:,mc);
        %x_noisy = x_nom;

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
                u_i = f * (disttolmk_k_N'*ihatC_k_N) / (disttolmk_k_N'*khatC_k_N) + u0 + sigma_u*randn(1);
                v_i = f * (disttolmk_k_N'*jhatC_k_N) / (disttolmk_k_N'*khatC_k_N) + v0 + sigma_v*randn(1);
        
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
        x_bar_LKF = zeros(length(tspan),6); 
        P = zeros(6,6,length(tspan));
        
       % Set initial conditions
        x_bar_LKF(1,:) = x_noisy(1,:) - x_nom(1,:); 

        P(:,:,1) = P0; % Use the provided initial covariance matrix P0
    
        for k = -1:length(tspan)-2

            if k >= 0
               % Time Update / Prediction Step occurs at every k (dt = 60s)
                x_bar_LKF((k+1)+1,:) = (Fbar_k(:,:,k+1) * x_bar_LKF(k+1,:)')';
                P(:,:,(k+1)+1) = Fbar_k(:,:,k+1) * P(:,:,k+1) * Fbar_k(:,:,k+1)' + Q;
            end

           % Extract time t(k+1)
            t_kp1 = tspan((k+1)+1);

           % The measurement update occurs at each observation epoch (dt = 600s)
            y_noisy_kp1 = y_noisy_tbl((y_noisy_tbl(:,1)==t_kp1)~=0,:);
            if size(y_noisy_kp1,1) > 0

               % Get nominal state at time k+1
                x_nom_kp1 = x_nom((k+1)+1,:)';

               % Get ideal measurements at time k
                y_full_kp1 = y_full((y_full(:,1)==t_kp1)~=0,:);
        
               % Truncate ideal measurements to only include landmarks which appear in the actual measurements
                y_trunc_kp1 = []; 
                for row = 1:size(y_noisy_kp1,1)
                    lmkid = y_noisy_kp1(row,2);
                    y_trunc_kp1 =[y_trunc_kp1; y_full_kp1((y_full_kp1(:,2) == lmkid),:)];   
                end 
                
               % Calculate and reshape y_bar_k
                y_bar_kp1 = y_noisy_kp1(:,3:4) - y_trunc_kp1(:,3:4);  
                y_bar_kp1 = reshape(y_bar_kp1',[numel(y_bar_kp1),1]);
    
               % Get rotation matrix R_CtoN
                R_CtoN_kp1 = R_CtoN(:,:,(t_kp1/dt_obs)+1);
            
               % Get rotation matrix R_AtoN
                theta = w_A*t_kp1;
                R_AtoN_kp1 = rotZ(theta);
                pos_lmks_N = R_AtoN_kp1*pos_lmks_A; 
        
               % Calculate measurement Jacobian and discretize
                Cbar_kp1 = linearizedCmat(f, R_CtoN_kp1, pos_lmks_N, y_trunc_kp1, x_nom_kp1);
                Hbar_kp1 = linearizedHmat(Cbar_kp1);
        
               % Compute Kalman gain 
                Big_R = diagTile(R,size(Hbar_kp1,1)/2); 
                Kkp1 = P(:,:,(k+1)+1) * Hbar_kp1' / (Hbar_kp1 * P(:,:,(k+1)+1) * Hbar_kp1' + Big_R);
                        
               % Measurement Update / Correction Step       
                x_bar_LKF((k+1)+1,:) = (x_bar_LKF((k+1)+1,:)' + Kkp1 * (y_bar_kp1 - Hbar_kp1*x_bar_LKF((k+1)+1,:)'))'; % Update state estimate
                P(:,:,(k+1)+1) = (eye(6) - Kkp1 * Hbar_kp1) * P(:,:,(k+1)+1); % Update covariance matrix
                
            end 
    
           % Calculate and save NEES/NIS
            S_kp1 = (Hbar_kp1 * P(:,:,(k+1)+1) * Hbar_kp1' + Big_R);
            [NEES_kp1, NIS_kp1] = calculateNEESNIS(x_bar_LKF((k+1)+1,:)', y_bar_kp1, P(:,:,(k+1)+1), S_kp1, alpha, NMC);
            
            NEES((k+1)+1,NMC) = NEES_kp1.NEES; NEESr1((k+1)+1,NMC) = NEES_kp1.r1; NEESr2((k+1)+1,NMC) = NEES_kp1.r2;
            if size(y_noisy_kp1,1) > 0
                NIS((k+1)+1,NMC)  = NIS_kp1.NIS;   NISr1((k+1)+1,NMC)  = NIS_kp1.r1;  NISr2((k+1)+1,NMC)  = NIS_kp1.r2; 
            else
                NIS((k+1)+1,NMC)  = nan;           NISr1((k+1)+1,NMC)  = nan;         NISr2((k+1)+1,NMC)  = nan; 
            end

        end
    
    end
    
    
    figure(53);
    subplot(2,1,1); hold off;
    plot(tspan,NEES','o','DisplayName','NEES'); hold on;
    plot(tspan,NEESr1','r.','DisplayName','$r_{1_{NEES}}$');
    plot(tspan,NEESr2','r.','DisplayName','$r_{2_{NEES}}$');
    legend('interpreter','latex');

    subplot(2,1,2); hold off;
    plot(tspan,NIS','o','DisplayName','NIS'); hold on;
    plot(tspan,NISr1','r.','DisplayName','$r_{1_{NIS}}$');
    plot(tspan,NISr2','r.','DisplayName','$r_{2_{NIS}}$');
    legend('interpreter','latex')


    if PLOTFLAG
        % Position error and +-2sigma (single plot)
        figure(51); 
        % movefig(gcf,'r');
        subplot(3,1,1); hold off;
        plot(tspan,x_bar_LKF(:,1),'DisplayName','$\bar{x}$','Color',mlc(1),'LineWidth',2); hold on;
        plot(tspan,2*sqrt(squeeze(P(1,1,:))),'k--','DisplayName','$\pm2\sigma$','LineWidth',1.5);
        plot(tspan,-2*sqrt(squeeze(P(1,1,:))),'k--','LineWidth',1.5,'HandleVisibility','off');
        labels(gca,{'Time [s]','$\bar{x}$ [km]'},'');
        %set(gca,'YLim',[-0.002 0.002]); legend('interpreter','latex');
        
        subplot(3,1,2); hold off;
        plot(tspan,x_bar_LKF(:,2),'DisplayName','$\bar{y}$','Color',mlc(2),'LineWidth',2); hold on;
        plot(tspan,2*sqrt(squeeze(P(2,2,:))),'k--','DisplayName','$\pm2\sigma$','LineWidth',1.5);
        plot(tspan,-2*sqrt(squeeze(P(2,2,:))),'k--','LineWidth',1.5,'HandleVisibility','off');
        labels(gca,{'Time [s]','$\bar{y}$ [km]'},'');
        %set(gca,'YLim',[-0.006 0.006]); legend('interpreter','latex');
        
        subplot(3,1,3); hold off;
        plot(tspan,x_bar_LKF(:,3),'DisplayName','$\bar{z}$','Color',mlc(3),'LineWidth',2); hold on;
        plot(tspan,2*sqrt(squeeze(P(3,3,:))),'k--','DisplayName','$\pm2\sigma$','LineWidth',1.5);
        plot(tspan,-2*sqrt(squeeze(P(3,3,:))),'k--','LineWidth',1.5,'HandleVisibility','off');
        labels(gca,{'Time [s]','$\bar{z}$ [km]'},'');
        %set(gca,'YLim',[-0.006 0.006]); legend('interpreter','latex');
        fixfig(gcf,false);
        
        % Velocity error and +-2sigma (single plot)
        figure(52); 
        % movefig(gcf,'r');
        subplot(3,1,1); hold off;
        plot(tspan,x_bar_LKF(:,4),'DisplayName','$\bar{\dot{x}}$','Color',mlc(1),'LineWidth',2); hold on;
        plot(tspan,2*sqrt(squeeze(P(4,4,:))),'k--','DisplayName','$\pm2\sigma$','LineWidth',1.5);
        plot(tspan,-2*sqrt(squeeze(P(4,4,:))),'k--','LineWidth',1.5,'HandleVisibility','off');
        labels(gca,{'Time [s]','$\bar{\dot{x}}$ [km/s]'},'');
        %set(gca,'YLim',[-0.000001 0.000001]); legend('interpreter','latex');
        
        subplot(3,1,2); hold off;
        plot(tspan,x_bar_LKF(:,5),'DisplayName','$\bar{\dot{y}}$','Color',mlc(2),'LineWidth',2); hold on;
        plot(tspan,2*sqrt(squeeze(P(5,5,:))),'k--','DisplayName','$\pm2\sigma$','LineWidth',1.5);
        plot(tspan,-2*sqrt(squeeze(P(5,5,:))),'k--','LineWidth',1.5,'HandleVisibility','off');
        labels(gca,{'Time [s]','$\bar{\dot{y}}$ [km/s]'},'');
        %set(gca,'YLim',[-0.000001 0.000001]); legend('interpreter','latex');
        
        subplot(3,1,3); hold off;
        plot(tspan,x_bar_LKF(:,6),'DisplayName','$\bar{\dot{z}}$','Color',mlc(3),'LineWidth',2); hold on;
        plot(tspan,2*sqrt(squeeze(P(6,6,:))),'k--','DisplayName','$\pm2\sigma$','LineWidth',1.5);
        plot(tspan,-2*sqrt(squeeze(P(6,6,:))),'k--','LineWidth',1.5,'HandleVisibility','off');
        labels(gca,{'Time [s]','$\bar{\dot{z}}$ [km/s]'},'');
        %set(gca,'YLim',[-0.000001 0.000001]); legend('interpreter','latex');
        fixfig(gcf,false);
    end


    % ---------------------------------------------------------------------------------------------------------------- 
    % EXTENDED KALMAN FILTER (first draft)
    % ----------------------------------------------------------------------------------------------------------------
    
    % Paste this code into asteroid_Part2
    
    % -- Process noise covariance matrix 
         Qfactor = 5e5;
         Q = Qfactor * sigma_w^2*[dt_int^3/3*eye(3)    dt_int^2/2*eye(3);...
                                  dt_int^2/2*eye(3)    dt_int*eye(3)];
     % -- Measurement noise covariance matrix (for each landmark)
         R = 100*diag([sigma_u^2 sigma_v^2]);
    
    % Initialize
    % mk_ekf = x_bar(1,:)+x_nom(1,:); % Total State 
    mk_ekf = [r0; rdot0]; %initial state
    
    P0_pos = (.05^2)*eye(3);
    P0_vel = (1.e-3^2)*eye(3);
    P0 = blkdiag(P0_pos,P0_vel); 
    P_ekf = P0; %initial covariance 
    % Container
    mu_ekf_hist = zeros(length(tspan),6); 
    % mu_ekf_hist(1,:)= mk_ekf; 
    Pk_ekf_filt_hist = zeros(6,6,length(tspan));
    
    
    % Perform Initial Measurment Update (since we have a measurement @ t=0)
    t_k = 0;
            y_pred_k = y_pred_tbl((y_pred_tbl(:,1)==t_k)~=0,:);
                
           % Get actual measurements at time k
            y_noisy_k = y_noisy_tbl((y_noisy_tbl(:,1)==t_k)~=0,:);
    
           % Truncate predicted measurements to only include landmarks which appear in the actual measurements
            y_trunc_pred_k = []; 
            for row = 1:size(y_noisy_k,1)
                lmkid = y_noisy_k(row,2);
                y_trunc_pred_k =[y_trunc_pred_k; y_pred_k((y_pred_k(:,2) == lmkid),:)];   
            end 

        % Calculate and reshape NL measurement residual/innovation
        innov_kp1 = y_noisy_k(:,3:4) - y_trunc_pred_k(:,3:4);  
        innov_kp1 = reshape(innov_kp1',[numel(innov_kp1),1]);

        % Calculate Jacobian at Predicted State
        Cbar_k = linearizedCmat(f, R_CtoN_k, pos_lmks_N, y_trunc_pred_k,mk_ekf);
        Hbar_k = linearizedHmat(Cbar_k);
    
        % Compute Kalman gain 
        Big_R = diagTile(R,size(Hbar_k,1)/2); 
        Kkp1 = P_ekf * Hbar_k' / (Hbar_k * P_ekf * Hbar_k' + Big_R);
            
       % Measurement Update / Correction Step       
        mk_ekf_p1plus = mk_ekf + Kkp1 * innov_kp1; % Update state estimate

        Pk_ekf_p1plus = (eye(6) - Kkp1 * Hbar_k) * P_ekf; % Update covariance matrix

        %%store results and cycle for next iteration
        mk_ekf = mk_ekf_p1plus'; 
        mu_ekf_hist(1,:) = mk_ekf_p1plus';
        P_ekf = Pk_ekf_p1plus;
        Pk_ekf_filt_hist(:,:,1)=Pk_ekf_p1plus;

%length(tspan)-1
for k = 1:length(tspan)-1
        t_k = tspan(k+1);
        if (k+2) <= length(tspan)
            t_meas = tspan((k+1)+1);
        end
        fprintf("Performing time update")
        fprintf("\n")
        disp(['t_k = ', num2str(t_k)]);
        fprintf("time update for k= %d", k)
        fprintf("\n")
        fprintf("\n")
        disp(['t_meas = ', num2str(t_meas)]);

       % Time Update / Prediction Step occurs at every k (dt = 60s)
        tol=1.e-14;
        OPTIONS = odeset('RelTol',3*tol,'AbsTol',tol);
        % Define the time span for the current integration step
        tspan_current = [tspan(k), tspan(k+1)];
        % Call the nonlinear function using ode45
        [t_NL,x_NL] = ode45(@(t,x) satDynamicModel(t,x,[],processNoise),tspan_current,mk_ekf,OPTIONS);
    
        % Extract the predicted state from the solution
        mk_ekf_p1minus = x_NL(end,:);
        
        % Accordind to the sudo algorithm, we have to linearize about the state estimate
        % so F_EFK would become I + dt_int*A_tilde(state_estimate) 
        F_EKF = eye(6) + dt_int*linearizedAmat(mu_A,mk_ekf); 
        Pk_ekf_p1minus = F_EKF*P_ekf*F_EKF'+Q;



       % The measurement update occurs at each observation epoch (dt = 600s)
       % The measurement update occurs at each observation epoch (dt = 600s)
         y_noisy_kp1 = y_noisy_tbl((y_noisy_tbl(:,1)==t_meas)~=0,:);
        if mod(k,10) == 9 && t_meas <= y_noisy_tbl(end,1) && size(y_noisy_kp1,1) > 0
            t_k=tspan(k+1);
            fprintf("if statement triggered at k=")
            fprintf("\n")
            disp(['t_k = ', num2str(t_k)]);
            fprintf("time update for k= %d", k)
            fprintf("\n")
           % Get ideal measurements at time k
           % t_k = tspan(k+1); 
            y_pred_k = y_pred_tbl((y_pred_tbl(:,1)==t_meas)~=0,:);
            
           % Get actual measurements at time k
            y_noisy_k = y_noisy_tbl((y_noisy_tbl(:,1)==t_meas)~=0,:);
    
           % Truncate predicted measurements to only include landmarks which appear in the actual measurements
            y_trunc_pred_k = []; 
            for row = 1:size(y_noisy_k,1)
                lmkid = y_noisy_k(row,2);
                y_trunc_pred_k =[y_trunc_pred_k; y_pred_k((y_pred_k(:,2) == lmkid),:)];   
            end 

            % Calculate and reshape NL measurement residual/innovation
            innov_kp1 = y_noisy_k(:,3:4) - y_trunc_pred_k(:,3:4);  
            innov_kp1 = reshape(innov_kp1',[numel(innov_kp1),1]);

            % Calculate Jacobian at Predicted State
            Cbar_k = linearizedCmat(f, R_CtoN_k, pos_lmks_N, y_trunc_pred_k,mk_ekf_p1minus');
            Hbar_k = linearizedHmat(Cbar_k);
        
            % Compute Kalman gain 
            Big_R = diagTile(R,size(Hbar_k,1)/2); 
            Kkp1 = Pk_ekf_p1minus * Hbar_k' / (Hbar_k * Pk_ekf_p1minus * Hbar_k' + Big_R);
                
           % Measurement Update / Correction Step       
            mk_ekf_p1plus = mk_ekf_p1minus' + Kkp1 * innov_kp1; % Update state estimate
            Pk_ekf_p1plus = (eye(6) - Kkp1 * Hbar_k) * Pk_ekf_p1minus; % Update covariance matrix

            %%store results and cycle for next iteration
            mk_ekf = mk_ekf_p1plus'; 
            mu_ekf_hist(k+1,:) = mk_ekf_p1plus';
            P_ekf = Pk_ekf_p1plus;
            Pk_ekf_filt_hist(:,:,k+1)=Pk_ekf_p1plus;
            
        else
        % Update state estimate and covariance matrix without measurement update
        mk_ekf = mk_ekf_p1minus'; 
        mu_ekf_hist(k+1, :) = mk_ekf_p1minus'; % Store results for every iteration
        P_ekf = Pk_ekf_p1minus;
        Pk_ekf_filt_hist(:, :, k+1) = Pk_ekf_p1minus;
        end


    end
     

load('x_noisy_MC1.mat')
    % Position error and +-2sigma (single plot)
    figure;
    subplot(3,1,1);
    plot(tspan,mu_ekf_hist(:,1),'DisplayName','$\bar{x}$','Color',mlc(1),'LineWidth',2); hold on;
    plot(tspan,2*sqrt(squeeze(Pk_ekf_filt_hist(1,1,:))),'k--','DisplayName','$\pm2\sigma$','LineWidth',1.5);
    plot(tspan,-2*sqrt(squeeze(Pk_ekf_filt_hist(1,1,:))),'k--','LineWidth',1.5,'HandleVisibility','off');
    plot(tspan,x_noisy_MC(:,1),'DisplayName','$\bar{x}$','Color','r','LineWidth',2);
    labels(gca,{'Time [s]','x [km]'},'');
    legend('total state (ekf)', '+/- 2\sigma','noisy state')
    subtitle('EKF State Position Results (First Code)')
    % set(gca,'YLim',[-0.03 0.002]); 
    
    subplot(3,1,2);
    plot(tspan,mu_ekf_hist(:,2),'DisplayName','$\bar{y}$','Color',mlc(2),'LineWidth',2); hold on;
    plot(tspan,2*sqrt(squeeze(Pk_ekf_filt_hist(2,2,:))),'k--','DisplayName','$\pm2\sigma$','LineWidth',1.5);
    plot(tspan,-2*sqrt(squeeze(Pk_ekf_filt_hist(2,2,:))),'k--','LineWidth',1.5,'HandleVisibility','off');
    plot(tspan,x_noisy_MC(:,2),'DisplayName','$\bar{y}$','Color','r','LineWidth',2);
    labels(gca,{'Time [s]','y [km]'},'');
    % set(gca,'YLim',[-1 1]); legend('interpreter','latex');
    
    subplot(3,1,3);
    plot(tspan,mu_ekf_hist(:,3),'DisplayName','$\bar{z}$','Color',mlc(3),'LineWidth',2); hold on;
    plot(tspan,2*sqrt(squeeze(Pk_ekf_filt_hist(3,3,:))),'k--','DisplayName','$\pm2\sigma$','LineWidth',1.5);
    plot(tspan,-2*sqrt(squeeze(Pk_ekf_filt_hist(3,3,:))),'k--','LineWidth',1.5,'HandleVisibility','off');
    plot(tspan,x_noisy_MC(:,3),'DisplayName','$\bar{z}$','Color','r','LineWidth',2); hold on;
    labels(gca,{'Time [s]','z [km]'},'');
    % set(gca,'YLim',[-1 1]); legend('interpreter','latex');
    % fixfig(gcf,false);


    
    % Velocity error and +-2sigma (single plot)
    figure;
    subplot(3,1,1);
    plot(tspan,mu_ekf_hist(:,4),'DisplayName','\dotx','Color',mlc(1),'LineWidth',2); hold on;
    plot(tspan,2*sqrt(squeeze(Pk_ekf_filt_hist(4,4,:))),'k--','DisplayName','$\pm2\sigma$','LineWidth',1.5);
    plot(tspan,-2*sqrt(squeeze(Pk_ekf_filt_hist(4,4,:))),'k--','LineWidth',1.5,'HandleVisibility','off');
    labels(gca,{'Time [s]','$\bar{\dot{z}}$ [km/s]'},'');
    subtitle('EKF State Velocity Results (First Code)')
    % set(gca,'YLim',[-1e-6 1e-6]); legend('interpreter','latex');
    
    subplot(3,1,2);
    plot(tspan,mu_ekf_hist(:,5),'DisplayName','\doty','Color',mlc(2),'LineWidth',2); hold on;
    plot(tspan,2*sqrt(squeeze(Pk_ekf_filt_hist(5,5,:))),'k--','DisplayName','$\pm2\sigma$','LineWidth',1.5);
    plot(tspan,-2*sqrt(squeeze(Pk_ekf_filt_hist(5,5,:))),'k--','LineWidth',1.5,'HandleVisibility','off');
    labels(gca,{'Time [s]','$\bar{\dot{y}}$ [km/s]'},'');
    % ylim([min(x_EKF(:,5)), max(x_EKF(:,5))]);
    % set(gca,'YLim',[-1e-4 1e-4]); legend('interpreter','latex');
    
    subplot(3,1,3);
    plot(tspan,mu_ekf_hist(:,6),'DisplayName','$\dot{z}$','Color',mlc(3),'LineWidth',2); hold on;
    plot(tspan,2*sqrt(squeeze(Pk_ekf_filt_hist(6,6,:))),'k--','DisplayName','$\pm2\sigma$','LineWidth',1.5);
    plot(tspan,-2*sqrt(squeeze(Pk_ekf_filt_hist(6,6,:))),'k--','LineWidth',1.5,'HandleVisibility','off');
    labels(gca,{'Time [s]','$\bar{\dot{z}}$ [km/s]'},'');
    % set(gca,'YLim',[-1e-4 1e-4]); legend('interpreter','latex');
    % fixfig(gcf,false);



 % ========================= MONTE-CARLO LOOP =========================
 NMC = 1;            % # of Monte-Carlo simulations (MAX 25)

    for mc=1:NMC
    
       % Configure tolerances
        tol=1.e-8;
        OPTIONS = odeset('RelTol',tol,'AbsTol',tol);
        load('x_noisy_MC30.mat')
       % ---------------------------------------------------------
       % Generate noisy truth state
       % ---------------------------------------------------------
        % x_noisy = x_noisy_MC(:,:,mc);
        x_noisy = x_nom;

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
                u_i = f * (disttolmk_k_N'*ihatC_k_N) / (disttolmk_k_N'*khatC_k_N) + u0 + sigma_u*randn(1);
                v_i = f * (disttolmk_k_N'*jhatC_k_N) / (disttolmk_k_N'*khatC_k_N) + v0 + sigma_v*randn(1);
        
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
       % EXTENDED KALMAN FILTER (second draft)
       % ----------------------------------------------------------------------------------------------------------------
        % -- Process noise covariance matrix 
         Qfactor = 5e5;
         Q = Qfactor * sigma_w^2*[dt_int^3/3*eye(3)    dt_int^2/2*eye(3);...
                                  dt_int^2/2*eye(3)    dt_int*eye(3)];
         
        % -- Measurement noise covariance matrix (for each landmark)
        % R = 100*diag([sigma_u^2 sigma_v^2]);
        R = diag([0 0]);
       % Containers 
        x_ekf = zeros(length(tspan),6); 
        P_ekf = zeros(6,6,length(tspan));
        
       % Set initial conditions
        % x_ekf(1,:) = x_noisy(1,:) ; 

        P_ekf(:,:,1) = P0; % Use the provided initial covariance matrix P0
    
        for k = -1:length(tspan)-2

            if k >= 0
                % Time Update / Prediction Step occurs at every k (dt = 60s)
                tol=1.e-14;
                OPTIONS = odeset('RelTol',3*tol,'AbsTol',tol);
                % Define the time span for the current integration step
                tspan_current = [0, 60];
                % Call the nonlinear function using ode45
                [t_NL,x_NL] = ode45(@(t,x) satDynamicModel(t,x,[],[0;0;0]),tspan_current,x_ekf((k+1),:),OPTIONS);
                
                % Set the ode45 solution equal to the predicted state 
                 x_ekf((k+1)+1,:) = x_NL(end,:);
               % Time Update / Prediction Step occurs at every k (dt = 60s)
                F_EKF = eye(6) + dt_int*linearizedAmat(mu_A,x_ekf((k+1),:)); 
                P_ekf(:,:,(k+1)+1) = F_EKF * P_ekf(:,:,k+1) * F_EKF' + Q;
            end

           % Extract time t(k+1)
            t_kp1 = tspan((k+1)+1);

           % The measurement update occurs at each observation epoch (dt = 600s)
            y_noisy_kp1 = y_noisy_tbl((y_noisy_tbl(:,1)==t_kp1)~=0,:);
            if size(y_noisy_kp1,1) > 0


               % Get ideal measurements at time k
                y_full_kp1 = y_full((y_full(:,1)==t_kp1)~=0,:);
        
               % Truncate predicted measurements to only include landmarks which appear in the actual measurements
                y_trunc_pred_k = []; 
                for row = 1:size(y_noisy_k,1)
                    lmkid = y_noisy_k(row,2);
                    y_trunc_pred_k =[y_trunc_pred_k; y_pred_k((y_pred_k(:,2) == lmkid),:)];   
                end  
                
               % Calculate and reshape NL measurement residual/innovation
                innov_kp1 = y_noisy_k(:,3:4) - y_trunc_pred_k(:,3:4);  
                innov_kp1 = reshape(innov_kp1',[numel(innov_kp1),1]);
    
               % Get rotation matrix R_CtoN
                R_CtoN_kp1 = R_CtoN(:,:,(t_kp1/dt_obs)+1);
            
               % Get rotation matrix R_AtoN
                theta = w_A*t_kp1;
                R_AtoN_kp1 = rotZ(theta);
                pos_lmks_N = R_AtoN_kp1*pos_lmks_A; 
        
               % Calculate measurement Jacobian and discretize at the
               % predicted state 
               
                Cbar_kp1 = linearizedCmat(f, R_CtoN_kp1, pos_lmks_N, y_trunc_pred_k, x_ekf((k+1)+1,:)');
                Hbar_kp1 = linearizedHmat(Cbar_kp1);
        
               % Compute Kalman gain 
                Big_R = diagTile(R,size(Hbar_kp1,1)/2); 
                Kkp1 = P_ekf(:,:,(k+1)+1) * Hbar_kp1' / (Hbar_kp1 * P_ekf(:,:,(k+1)+1) * Hbar_kp1' + Big_R);
                        
               % Measurement Update / Correction Step       
                x_ekf((k+1)+1,:) = (x_ekf((k+1)+1,:))' + Kkp1 * innov_kp1; % Update state estimate
                P_ekf(:,:,(k+1)+1) = (eye(6) - Kkp1 * Hbar_kp1) * P_ekf(:,:,(k+1)+1); % Update covariance matrix
                
            end 

                       % Calculate and save NEES/NIS
            S_kp1 = (Hbar_kp1 * P_ekf(:,:,(k+1)+1) * Hbar_kp1' + Big_R);
            [NEES_kp1, NIS_kp1] = calculateNEESNIS(x_ekf((k+1)+1,:)', innov_kp1, P_ekf(:,:,(k+1)+1), S_kp1, alpha, NMC);
            
            NEES((k+1)+1,NMC) = NEES_kp1.NEES; NEESr1((k+1)+1,NMC) = NEES_kp1.r1; NEESr2((k+1)+1,NMC) = NEES_kp1.r2;
            if size(y_noisy_kp1,1) > 0
                NIS((k+1)+1,NMC)  = NIS_kp1.NIS;   NISr1((k+1)+1,NMC)  = NIS_kp1.r1;  NISr2((k+1)+1,NMC)  = NIS_kp1.r2; 
            else
                NIS((k+1)+1,NMC)  = nan;           NISr1((k+1)+1,NMC)  = nan;         NISr2((k+1)+1,NMC)  = nan; 
            end
   

        end
    end 
    figure(35);
    subplot(2,1,1); hold off;
    plot(tspan,NEES','o','DisplayName','NEES'); hold on;
    plot(tspan,NEESr1','r.','DisplayName','$r_{1_{NEES}}$');
    plot(tspan,NEESr2','r.','DisplayName','$r_{2_{NEES}}$');
    legend('interpreter','latex');

    subplot(2,1,2); hold off;
    plot(tspan,NIS','o','DisplayName','NIS'); hold on;
    plot(tspan,NISr1','r.','DisplayName','$r_{1_{NIS}}$');
    plot(tspan,NISr2','r.','DisplayName','$r_{2_{NIS}}$');
    legend('interpreter','latex')


    % Position and +-2sigma (single plot)
    figure;
    subplot(3,1,1);
    plot(tspan,x_ekf(:,1),'DisplayName','$\bar{x}$','Color',mlc(1),'LineWidth',2); hold on;
    plot(tspan,2*sqrt(squeeze(P_ekf(1,1,:))),'k--','DisplayName','$\pm2\sigma$','LineWidth',1.5);
    plot(tspan,-2*sqrt(squeeze(P_ekf(1,1,:))),'k--','LineWidth',1.5,'HandleVisibility','off');
    % plot(tspan,x_noisy_MC(:,1),'DisplayName','$\bar{x}$','Color','r','LineWidth',2);
    labels(gca,{'Time [s]','x [km]'},'');
    legend('total state (ekf)', '+/- 2\sigma','noisy state')
    subtitle('EKF State Position Results (second draft)')
    set(gca,'YLim',[-5 5]); 
    
    subplot(3,1,2);
    plot(tspan,x_ekf(:,2),'DisplayName','$\bar{y}$','Color',mlc(2),'LineWidth',2); hold on;
    plot(tspan,2*sqrt(squeeze(P_ekf(2,2,:))),'k--','DisplayName','$\pm2\sigma$','LineWidth',1.5);
    plot(tspan,-2*sqrt(squeeze(P_ekf(2,2,:))),'k--','LineWidth',1.5,'HandleVisibility','off');
    % plot(tspan,x_noisy_MC(:,2),'DisplayName','$\bar{y}$','Color','r','LineWidth',2);
    labels(gca,{'Time [s]','y [km]'},'');
    set(gca,'YLim',[-5 5]); 
    
    subplot(3,1,3);
    plot(tspan,x_ekf(:,3),'DisplayName','$\bar{z}$','Color',mlc(3),'LineWidth',2); hold on;
    plot(tspan,2*sqrt(squeeze(P_ekf(3,3,:))),'k--','DisplayName','$\pm2\sigma$','LineWidth',1.5);
    plot(tspan,-2*sqrt(squeeze(P_ekf(3,3,:))),'k--','LineWidth',1.5,'HandleVisibility','off');
    % plot(tspan,x_noisy_MC(:,3),'DisplayName','$\bar{z}$','Color','r','LineWidth',2); hold on;
    labels(gca,{'Time [s]','z [km]'},'');
    set(gca,'YLim',[-5 5]); 
    % fixfig(gcf,false);


    
    % Velocity and +-2sigma (single plot)
    figure;
    subplot(3,1,1);
    plot(tspan,x_ekf(:,4),'DisplayName','\dotx','Color',mlc(1),'LineWidth',2); hold on;
    plot(tspan,2*sqrt(squeeze(P_ekf(4,4,:))),'k--','DisplayName','$\pm2\sigma$','LineWidth',1.5);
    plot(tspan,-2*sqrt(squeeze(P_ekf(4,4,:))),'k--','LineWidth',1.5,'HandleVisibility','off');
    labels(gca,{'Time [s]','$\bar{\dot{z}}$ [km/s]'},'');
    subtitle('EKF State Velocity Results (second draft)')

    % set(gca,'YLim',[-1e-6 1e-6]); legend('interpreter','latex');
    
    subplot(3,1,2);
    plot(tspan,x_ekf(:,5),'DisplayName','\doty','Color',mlc(2),'LineWidth',2); hold on;
    plot(tspan,2*sqrt(squeeze(P_ekf(5,5,:))),'k--','DisplayName','$\pm2\sigma$','LineWidth',1.5);
    plot(tspan,-2*sqrt(squeeze(P_ekf(5,5,:))),'k--','LineWidth',1.5,'HandleVisibility','off');
    labels(gca,{'Time [s]','$\bar{\dot{y}}$ [km/s]'},'');
    % ylim([min(x_EKF(:,5)), max(x_EKF(:,5))]);
    % set(gca,'YLim',[-1e-4 1e-4]); legend('interpreter','latex');
    
    subplot(3,1,3);
    plot(tspan,x_ekf(:,6),'DisplayName','$\dot{z}$','Color',mlc(3),'LineWidth',2); hold on;
    plot(tspan,2*sqrt(squeeze(P_ekf(6,6,:))),'k--','DisplayName','$\pm2\sigma$','LineWidth',1.5);
    plot(tspan,-2*sqrt(squeeze(P_ekf(6,6,:))),'k--','LineWidth',1.5,'HandleVisibility','off');
    labels(gca,{'Time [s]','$\bar{\dot{z}}$ [km/s]'},'');
    % set(gca,'YLim',[-1e-4 1e-4]); legend('interpreter','latex');
    % fixfig(gcf,false);

     % Position Errors (single plot)
    figure;
    subplot(3,1,1);
    plot(tspan,x_ekf(:,1)-x_noisy_MC(:,1),'DisplayName','$\bar{x}$','Color',mlc(1),'LineWidth',2); hold on;
    % plot(tspan,2*sqrt(squeeze(P_ekf(1,1,:))),'k--','DisplayName','$\pm2\sigma$','LineWidth',1.5);
    % plot(tspan,-2*sqrt(squeeze(P_ekf(1,1,:))),'k--','LineWidth',1.5,'HandleVisibility','off');
    % plot(tspan,x_noisy_MC(:,1),'DisplayName','$\bar{x}$','Color','r','LineWidth',2);
    labels(gca,{'Time [s]','x [km]'},'');
    legend('total state (ekf)', '+/- 2\sigma','noisy state')
    subtitle('Position Errors (EKF) (second code)')
    % set(gca,'YLim',[-0.03 0.002]); 
    
    subplot(3,1,2);
    plot(tspan,x_ekf(:,2)-x_noisy_MC(:,2),'DisplayName','$\bar{y}$','Color',mlc(2),'LineWidth',2); hold on;
    % plot(tspan,2*sqrt(squeeze(P_ekf(2,2,:))),'k--','DisplayName','$\pm2\sigma$','LineWidth',1.5);
    % plot(tspan,-2*sqrt(squeeze(P_ekf(2,2,:))),'k--','LineWidth',1.5,'HandleVisibility','off');
    % plot(tspan,x_noisy_MC(:,2),'DisplayName','$\bar{y}$','Color','r','LineWidth',2);
    labels(gca,{'Time [s]','y [km]'},'');
    % set(gca,'YLim',[-1 1]); legend('interpreter','latex');
    
    subplot(3,1,3);
    plot(tspan,x_ekf(:,3)-x_noisy_MC(:,3),'DisplayName','$\bar{z}$','Color',mlc(3),'LineWidth',2); hold on;
    % plot(tspan,2*sqrt(squeeze(P_ekf(3,3,:))),'k--','DisplayName','$\pm2\sigma$','LineWidth',1.5);
    % plot(tspan,-2*sqrt(squeeze(P_ekf(3,3,:))),'k--','LineWidth',1.5,'HandleVisibility','off');
    % plot(tspan,x_noisy_MC(:,3),'DisplayName','$\bar{z}$','Color','r','LineWidth',2); hold on;
    labels(gca,{'Time [s]','z [km]'},'');
    % set(gca,'YLim',[-1 1]); legend('interpreter','latex');
    % fixfig(gcf,false);


    
    % Velocity Errors (single plot)
    figure;
    subplot(3,1,1);
    plot(tspan,x_ekf(:,4)-x_noisy_MC(:,4),'DisplayName','\dotx','Color',mlc(1),'LineWidth',2); hold on;
    % plot(tspan,2*sqrt(squeeze(P_ekf(4,4,:))),'k--','DisplayName','$\pm2\sigma$','LineWidth',1.5);
    % plot(tspan,-2*sqrt(squeeze(P_ekf(4,4,:))),'k--','LineWidth',1.5,'HandleVisibility','off');
    labels(gca,{'Time [s]','$\bar{\dot{z}}$ [km/s]'},'');
    subtitle('Velocity Errors (EKF) (second code)')

    % set(gca,'YLim',[-1e-6 1e-6]); legend('interpreter','latex');
    
    subplot(3,1,2);
    plot(tspan,x_ekf(:,5)-x_noisy_MC(:,5),'DisplayName','\doty','Color',mlc(2),'LineWidth',2); hold on;
    % plot(tspan,2*sqrt(squeeze(P_ekf(5,5,:))),'k--','DisplayName','$\pm2\sigma$','LineWidth',1.5);
    % plot(tspan,-2*sqrt(squeeze(P_ekf(5,5,:))),'k--','LineWidth',1.5,'HandleVisibility','off');
    labels(gca,{'Time [s]','$\bar{\dot{y}}$ [km/s]'},'');
    % ylim([min(x_EKF(:,5)), max(x_EKF(:,5))]);
    % set(gca,'YLim',[-1e-4 1e-4]); legend('interpreter','latex');
    
    subplot(3,1,3);
    plot(tspan,x_ekf(:,6)-x_noisy_MC(:,15),'DisplayName','$\dot{z}$','Color',mlc(3),'LineWidth',2); hold on;
    % plot(tspan,2*sqrt(squeeze(P_ekf(6,6,:))),'k--','DisplayName','$\pm2\sigma$','LineWidth',1.5);
    % plot(tspan,-2*sqrt(squeeze(P_ekf(6,6,:))),'k--','LineWidth',1.5,'HandleVisibility','off');
    labels(gca,{'Time [s]','$\bar{\dot{z}}$ [km/s]'},'');
    % set(gca,'YLim',[-1e-4 1e-4]); legend('interpreter','latex');
    % fixfig(gcf,false);


% ---------------------------------------------------------------------------------------------------------------- 
% UNSCENTED KALMAN FILTER
% ----------------------------------------------------------------------------------------------------------------
% AQ1: Unscented Kalman Filter (UKF) Implementation

function ukf_main()
    % Loading observation data log
    load('orbitdetermination-finalproj_data_2023_11_14.mat'); 
    
    % Configure integration time and tolerances (adding here for
    % convenience)
    tspan = t0:dt_int:tf_int;
    tol=1.e-14;
    OPTIONS = odeset('RelTol',3*tol,'AbsTol',tol);
    
    % Parameters
    num_states = 6;  % Number of states in the system
    num_measurements = 3;  % Number of measurements in the system
    num_steps = length(tspan);


    % Initial state estimate and covariance (tuning needed)
    x_hat = [r0; rdot0];  % Initial state estimate
    P_hat = eye(num_states);  % Initial state covariance

    % UKF Parameters (tuning needed- working before but inaccurate)
    alpha = 0.1;  % UKF scaling parameter
    beta = 2;  % UKF scaling parameter
    kappa = 0;  % UKF scaling parameter

    % Time parameters (remember- adjustment might be needed)
    dt = 1;  % Time step

    % UKF Initialization
    ukf_params = struct('alpha', alpha, 'beta', beta, 'kappa', kappa);
    ukf = init_ukf(x_hat, P_hat, Q, R, ukf_params);

    % Containers 
    x_ekf = zeros(length(tspan),6); 
    P_ekf = zeros(6,6,length(tspan));

    % Containers- block placement check
    % x_ukf = zeros(length(orbitdetermination_finalproj_data_2023_11_14), num_states);
    % P_ukf = zeros(num_states, num_states, length(orbitdetermination_finalproj_data_2023_11_14));

    % Set initial conditions- block placement check
    x_ukf(1,:) = x_hat';
    P_ukf(:,:,1) = P_hat;

    % Initializing arrays to store filter estimates
    Nlmks = size(pos_lmks_A,2);
    x_hat_history = zeros(num_states, num_steps);
    P_hat_history = zeros(num_states, num_states, num_steps);

    %%%% UKF Main Loop %%%%
    for k = 1:num_steps
        % Prediction Step
        ukf = ukf_predict(ukf, dt);

        % Update Step
        measurement = orbitdetermination_finalproj_data_2023_11_14(:, k);
        ukf = ukf_update(ukf, measurement);

        % Storing results
        x_hat_history(:, k) = ukf.x_hat;
        P_hat_history(:, :, k) = ukf.P_hat;
    end

    % Plot Results
    plot_results(x_hat_history, P_hat_history, orbitdetermination_finalproj_data_2023_11_14);
end

%%%% UKF Initialisation
function ukf = init_ukf(x_hat, P_hat, Q, R, params)
    ukf.x_hat = x_hat;
    ukf.P_hat = P_hat;
    ukf.Q = Q;
    ukf.R = R;
    ukf.params = params;
end


%%%% UKF Prediction
function ukf = ukf_predict(ukf, dt)
    % UKF Prediction Step starts here- - SPF Sigma Point Filtering starts**
    % Implementing the prediction step using the unscented transformation
    
    % Dynamics Model
    ukf.x_hat = satDynamicModel(dt, ukf.x_hat, zeros(size(ukf.x_hat)), ukf.Q);
    Abar_k = linearizedAmat(mu_A, ukf.x_hat);
    Fbar_k = linearizedFmat(dt, Abar_k);
    
    % Prediction
    ukf.x_hat = Fbar_k * ukf.x_hat;
    ukf.P_hat = Fbar_k * ukf.P_hat * Fbar_k' + ukf.Q;
end


%%%% UKF Update
function ukf = ukf_update(ukf, measurement)
    % UKF Update Step
    % Implementing the update step using the unscented transformation
    % Might need to update ukf.x_hat and ukf.P_hat but not sure

    % Linearization and Measurement Model
    Cbar_k = linearizedCmat(f, R_CtoN_k, pos_lmks_N, measurement, ukf.x_hat);
    Hbar_k = linearizedHmat(Cbar_k);
    
    % Kalman Gain
    innovation = measurement - Hbar_k * ukf.x_hat;
    innovation_cov = Hbar_k * ukf.P_hat * Hbar_k' + ukf.R;
    kalman_gain = ukf.P_hat * Hbar_k' / innovation_cov;

    % Update
    ukf.x_hat = ukf.x_hat + kalman_gain * innovation;
    ukf.P_hat = ukf.P_hat - kalman_gain * innovation_cov * kalman_gain';
end

%%%% UKF Results
function plot_results(x_hat_history, P_hat_history, orbitdetermination_finalproj_data_2023_11_14)
    % Plotting Results
    time = 1:length(orbitdetermination_finalproj_data_2023_11_14);
    figure;

    subplot(2, 1, 1);
    plot(time, x_hat_history(1:3, :), 'b-', 'LineWidth', 2);
    hold on;
    plot(time, orbitdetermination_finalproj_data_2023_11_14(1:3, :), 'r.', 'MarkerSize', 10);
    xlabel('Time');
    ylabel('Position (km)');
    legend('Filter Estimate', 'Ground Truth');

    subplot(2, 1, 2);
    plot(time, x_hat_history(4:6, :), 'b-', 'LineWidth', 2);
    hold on;
    plot(time, orbitdetermination_finalproj_data_2023_11_14(4:6, :), 'r.', 'MarkerSize', 10);
    xlabel('Time');
    ylabel('Velocity (km/s)');
    legend('Filter Estimate', 'Ground Truth');
end