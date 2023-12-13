% =================================================================================================
% Script: generateNoisyTruth
% Pre-generates realizations of noisy state date for asteroid project.
% =================================================================================================
% Parameters:
%        Nsets   number of sets to generate
%     saveFile   .mat filename to save out datasets
% =================================================================================================

asteroid_Part1;

Nsets = 40;
saveFile = "x_noisy_MC40"; % do not put .mat here, gets added later

x_noisy_MC = zeros(length(tspan),6,Nsets);

% Configure tolerances
tol=1.e-8;
OPTIONS = odeset('RelTol',3*tol,'AbsTol',tol);


P0pos = (0.01)^2*eye(3);
P0vel = (1.e-6)^2*eye(3);
P0 = blkdiag(P0pos,P0vel);

for i=1:Nsets
tic;

fprintf('Generating set  %d\n',i);

    x_noisy = zeros(length(tspan),6);
    x_noisy(1,:) = [r0;rdot0]' + mvnrnd(zeros(6,1),P0,1);
    tic;
    for t_k = t0+dt_int:dt_int:tf_int
        processNoise = sigma_w*randn(3,1);
        [t_noisy_k,x_noisy_k] = ode45(@(t,x) satDynamicModel(t,x,[],processNoise),[t_k-dt_int t_k],x_noisy(t_k/dt_int,:)',OPTIONS); 
        x_noisy(t_k/dt_int+1,:) = x_noisy_k(end,:);
    end

    x_noisy_MC(:,:,i) = x_noisy;

toc;

end

% make sure existing datasets are not overwritten
counter = 1;
saveFileOrig = saveFile;
while isfile(saveFile)
    saveFile = saveFileOrig + sprintf("_%d",counter) + ".mat";
    counter = counter + 1;
end

save(saveFile,'x_noisy_MC');

