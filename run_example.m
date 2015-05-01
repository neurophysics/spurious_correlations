%% EXAMPLE: HOW TO COMPUTE SIGNIFICANCE WITH SURROGATE DATA PERMUTATION TEST
% N Schaworonkow, DAJ Blythe, J Kegeles, G Curio, VV Nikulin: 
% Power-law dynamics in neuronal and behavioral data introduce spurious 
% correlations. Human Brain Mapping. 2015.
% http://doi.org/10.1002/hbm.22816

%% example data with known alpha
nr_samples = 1000;

alpha1 = 1;
alpha2 = 0.9;

x1 = simulate_powerlaw(nr_samples,alpha1);
x2 = simulate_powerlaw(nr_samples,alpha2);

%% compute significance
nr_repetitions = 1e4; 
[p_val_new, p_val_org] = get_significance(x1, x2, nr_repetitions)

%% plot result
figure; hold on; 

subplot(2,1,1)
plot(x1); 
title(['original corr: ' num2str(p_val_org), ...
         ', corrected: ' num2str(p_val_new)])
subplot(2,1,2)
plot(x2,'r')