% =================================================================================================
% Function: linearizedCmat
% Calculates the CT measurement Jacobian at time k.
% =================================================================================================
% Inputs:
%            f   camera focal length (px)
%     R_CtoN_k   (3 x 3) rotation matrix from camera frame to inertial frame, at time k
%   pos_lmks_N   (3 x Nlmks_k) matrix of landmark position vectors in inertial frame [km]
%      y_nom_k   (Nlmks_k x 4) nominal output at time k, each row consisting of [time lmkid u v]
%      x_nom_k   (6 x 1) nominal state at time k
% Outputs:
%       Cbar_k   (2*Nlmks_k x 6) linearized C matrix, corresponding to ybar_k = Cbar_k*xbar_k where
%                                ybar_k = y_k - y_nom_k , xbar_k = x_k - x_nom_k
% =================================================================================================
function Cbar_k = linearizedCmat(f, R_CtoN_k, pos_lmks_N, y_nom_k, x_nom_k)

    % Parse camera unit vectors, and components of khat
    ihat = R_CtoN_k(:,1); jhat = R_CtoN_k(:,2); khat = R_CtoN_k(:,3); 
    k1 = khat(1); k2 = khat(2); k3 = khat(3);

    % Extract nominal position and components from state vector
    xpos_k_nom = x_nom_k(1:3);
    x1 = xpos_k_nom(1); x2 = xpos_k_nom(2); x3 = xpos_k_nom(3);

    % Define number of nominally visible landmarks, instantiate Cbar
    Nlmks_k = size(y_nom_k,1);
    Cbar_k = zeros(2*Nlmks_k,6); 

    for j=1:Nlmks_k % For each landmark that SHOULD be visible at a given time:
        
        % Get landmark position vector in inertial frame, extract components
        lmkj_id = y_nom_k(j,2);
        lmkj_pos_k_N = pos_lmks_N(:,lmkj_id);
        l1 = lmkj_pos_k_N(1); l2 = lmkj_pos_k_N(2); l3 = lmkj_pos_k_N(3);
        
        % Define common factor for each component Cbar_k_j of Cbar_k
        cf = f / ( (lmkj_pos_k_N - xpos_k_nom)'*khat )^2;

        % Create each component of Cbar_k_j
        Cbar_k_j_11 = cf*[-k2*(l2-x2)-k3*(l3-x3); k1*(l2-x2); k1*(l3-x3)]*ihat;
        Cbar_k_j_21 = cf*[-k2*(l2-x2)-k3*(l3-x3); k1*(l2-x2); k1*(l3-x3)]*jhat;
        Cbar_k_j_12 = cf*[k2*(l1-x1); -k1*(l1-x1)-k3*(l3-x3); k2*(l3-x3)]*ihat;
        Cbar_k_j_22 = cf*[k2*(l1-x1); -k1*(l1-x1)-k3*(l3-x3); k2*(l3-x3)]*jhat;
        Cbar_k_j_13 = cf*[k3*(l1-x1); k3*(l2-x2); -k1*(l1-x1)-k2*(l1-x2)]*ihat;
        Cbar_k_j_23 = cf*[k3*(l1-x1); k3*(l2-x2); -k1*(l1-x1)-k2*(l2-x2)]*jhat;

        % Assemble
        Cbar_k_j = [Cbar_k_j_11 Cbar_k_j_12 Cbar_k_j_13 0 0 0;...
                    Cbar_k_j_21 Cbar_k_j_22 Cbar_k_j_23 0 0 0];

        % Insert into Cbar_k
        Cbar_k(2*j-1:2*j,:) = Cbar_k_j;

    end

end