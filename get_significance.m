%% COMPUTE SIGNIFICANCE WITH SURROGATE DATA PERMUTATION TEST
% N Schaworonkow, DAJ Blythe, J Kegeles, G Curio, VV Nikulin: 
% Power-law dynamics in neuronal and behavioral data introduce spurious 
% correlations. Human Brain Mapping. 2015.
% http://doi.org/10.1002/hbm.22816

% Surrogates via AAFT: Adjusted Amplitude Fourier Transform
% James Theiler, Stephen Eubank, Andre; Longtin, Bryan Galdrikian, 
% and J. Doyne Farmer. 1992. Testing for nonlinearity in time series: 
% the method of surrogate data. Phys. D 58, 1-4 (September 1992), 77-94. 
% DOI=10.1016/0167-2789(92)90102-S 
% http://dx.doi.org/10.1016/0167-2789(92)90102-S

% IN:
%     TS      [nr_samples x 1] : first time series
%     TS2     [nr_samples x 1] : second time series
%     nr_repetitions [integer] : permutation test iterations
% OUT:
%              p_val_new  [float] : corrected p_value
%              p_val_org  [float] : p_value obtained by standard testing

function[p_val_new, p_val_org] = get_significance(TS1, TS2, nr_repetitions)

TS1 = TS1(:); TS2 = TS2(:);
nr_samples = numel(TS1);

% original correlation
[r_org, p_val_org] = corr(TS1, TS2, 'type', 'Spearman');

% AAFT
r_surrogate = zeros(nr_repetitions,1);
for i = 1:nr_repetitions
    
    if mod(i,1000)==0
        display(['run iteration: ' num2str(i) '/' num2str(nr_repetitions)])
    end
    % create white noise vector with n entries
    white_noise = sort(randn(nr_samples,1));
    % Sort z and extract the ranks
    [sorted_z, ranks] = sort(TS1);
    [~, idx_ranks] = sort(ranks);
    
    % random phase surrogate on white noise
    y = fft(white_noise(idx_ranks)');
    Z_amps = abs(y);
    rand_phases = rand(1,floor(nr_samples/2))*2*pi;
    start = length(rand_phases);
    
    if mod(nr_samples,2) == 0
        start = start-1;
    end
    
    rand_phases = [0, rand_phases, -rand_phases(start:-1:1)];
    % put amps and phases together for complex Fourier spectrum
    white_noise_rand_phase = Z_amps .* exp(1i*rand_phases);
    % project the complex spectrum back to the time domain
    white_noise_rand_phase = real(ifft(white_noise_rand_phase));
    
    % extract the ranks of the phase randomized white noise
    [~, ranks] = sort(white_noise_rand_phase);
    [~, idx_ranks] = sort(ranks);
    
    % assign ranks of phase randomized normal deviates to sorted_z,
    % obtain AAFT surrogates
    z_surrogate = sorted_z(idx_ranks);
    z_surrogate = z_surrogate(:);
    
    r_surrogate(i) = corr(z_surrogate, TS2, 'type', 'Spearman');

end

p_val_new = sum(abs(r_surrogate)>repmat(abs(r_org),nr_repetitions,1)) ...
                /nr_repetitions;

end