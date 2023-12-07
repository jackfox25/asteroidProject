asteroid_Part1;

%%
Nsets = 25;

x_noisy_MC25 = zeros(length(tspan),6,Nsets);

% Configure tolerances
tol=1.e-8;
OPTIONS = odeset('RelTol',3*tol,'AbsTol',tol);

for i=1:Nsets
tic;

fprintf('Generating set  %d\n',i);

    x_noisy = zeros(length(tspan),6);
    x_noisy(1,:) = x0';
    tic;
    for t_k = t0+dt_int:dt_int:tf_int
        processNoise = sigma_w^2*randn(3,1);
        [t_noisy_k,x_noisy_k] = ode45(@(t,x) satDynamicModel(t,x,[],processNoise),[t_k-dt_int t_k],x_noisy(t_k/dt_int,:)',OPTIONS); 
        x_noisy(t_k/dt_int+1,:) = x_noisy_k(end,:);
    end

    x_noisy_MC25(:,:,i) = x_noisy;

toc;

end

