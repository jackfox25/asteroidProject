function xdot = satDynamicModel(t,x,u,processNoise)

    simParameters;
    r = x(1:3); rmag = norm(r);
    
    a2b = -mu_A/rmag^3 * r;
    aSRP = -SPC/(norm(r_S))^2 * (1 + 4*rho/9) * Aoverm * r_S/norm(r_S);
    
    xdot = [x(4);...
            x(5);...
            x(6);...
            a2b + aSRP + processNoise];

end

