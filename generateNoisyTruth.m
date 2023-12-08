% =================================================================================================
% Script: generateNoisyTruth
% Pre-generates realizations of noisy state date for asteroid project.
% =================================================================================================
% Parameters:
%        Nsets   number of sets to generate
%     saveFile   .mat filename to save out datasets
% =================================================================================================

asteroid_Part1;

Nsets = 25;
saveFile = "x_noisy_MC25"; % do not put .mat here, gets added later

x_noisy_MC = zeros(length(tspan),6,Nsets);

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