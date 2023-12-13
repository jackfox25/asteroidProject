asteroid_Part1;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% PART 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

    NMC = 1;            % # of Monte-Carlo simulations (MAX 25)
    PLOTFLAG = true;    % decide whether to produce plots or not
    
    load('x_noisy_MC40.mat');
    
   % -- Process noise covariance matrix
     Qfactor = 1;
     Q = Qfactor * sigma_w^2*[dt_int^3/3*eye(3)    dt_int^2/2*eye(3);...
                              dt_int^2/2*eye(3)    dt_int*eye(3)];

   % -- Initial state covariance matrix
     P0pos = (0.01)^2*eye(3);
     P0vel = (1.e-6)^2*eye(3);
     P0 = blkdiag(P0pos,P0vel);

   % -- Measurement noise covariance matrix (for each landmark)
     R = 1*diag([sigma_u^2 sigma_v^2]);
    
    
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
            if size(y_noisy_kp1,1) > 0 && true

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
            
            NEES((k+1)+1,mc) = NEES_kp1.NEES; NEESr1((k+1)+1,mc) = NEES_kp1.r1; NEESr2((k+1)+1,mc) = NEES_kp1.r2;
            if size(y_noisy_kp1,1) > 0
                NIS((k+1)+1,mc)  = NIS_kp1.NIS;   NISr1((k+1)+1,mc)  = NIS_kp1.r1;  NISr2((k+1)+1,mc)  = NIS_kp1.r2; 
            else
                NIS((k+1)+1,mc)  = nan;           NISr1((k+1)+1,mc)  = nan;         NISr2((k+1)+1,mc)  = nan; 
            end

        end
    
    end
    
    
    figure(53); movefig(gcf,'r');
    subplot(2,1,1); hold off;
    plot(tspan,mean(NEES,2)','o','DisplayName','NEES'); hold on;
    plot(tspan,mean(NEESr1,2)','r.','DisplayName','$r_{1_{NEES}}$');
    plot(tspan,mean(NEESr2,2)','r.','DisplayName','$r_{2_{NEES}}$');
    legend('interpreter','latex');

    subplot(2,1,2); hold off;
    plot(tspan,mean(NIS,2)','o','DisplayName','NIS'); hold on;
    plot(tspan,mean(NISr1,2)','r.','DisplayName','$r_{1_{NIS}}$');
    plot(tspan,mean(NISr2,2)','r.','DisplayName','$r_{2_{NIS}}$');
    legend('interpreter','latex')


    if PLOTFLAG
        % Position error and +-2sigma (single plot)
        figure(51); movefig(gcf,'r');
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
        figure(52); movefig(gcf,'r');
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
