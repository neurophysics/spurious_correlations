%% PLOT SIMULATION FIGURES

% plot figures 2a, 2b, 2c from:
% N Schaworonkow, DAJ Blythe, J Kegeles, G Curio, VV Nikulin: 
% Power-law dynamics in neuronal and behavioral data introduce spurious 
% correlations. Human Brain Mapping. 2015.
% http://doi.org/10.1002/hbm.22816

addpath('helper')

nr_repetitions = 100;
nr_samples = 800;
p_alpha = 0.05;     % set significance level

file_name = [ 'sim_iter_' num2str(nr_repetitions) ...
                 '_length_' num2str(nr_samples) ];

nr_timeseries = 251;
alpha = linspace(0.5, 1.5, nr_timeseries);

timeseries = zeros(nr_timeseries,nr_samples);
correlations = zeros(nr_repetitions,nr_timeseries,nr_timeseries);

for j = 1:nr_repetitions
    display(['simulating iteration ' num2str(j) ' ...'])
    for i=1:nr_timeseries
       timeseries(i,:) = simulate_powerlaw(nr_samples,alpha(i));
    end

    correlations(j,:,:) = 1-squareform(pdist(timeseries,'spearman'));
  
end

% calculate significant correlation boundary value with t-approximation
t = tinv(1-p_alpha/2,nr_samples-2);
significant_threshold = sqrt(1/((nr_samples-2)/t^2 +1));

% count number of correlations > boundary value for each alpha-combo
percentages = squeeze(sum(abs(correlations)>significant_threshold,1))...
                /nr_repetitions;

save(file_name, 'correlations', 'percentages', 'alpha')

%% PLOT: COLORMAP
map = div_colormap;
nancolor = [0.285, 0.732, 0.401];
map2 = map(floor(end/2)+1:end,:);

%% FIG: 2A
set(0,'defaulttextinterpreter','latex')

correlation_1try = squeeze(correlations(1,:,:));

figure; 
imagesc(alpha,alpha,correlation_1try)

xlabel('$\alpha$-exponent (time series 1)', 'FontSize', 18); 
ylabel('$\alpha$-exponent (time series 2)', 'FontSize', 18);
axis square; axis xy
cm = colormap(map);

set(gca, ...
            'FontSize',                               14 , ...
            'FontName',                       'CMU Serif', ...
            'YTick',                          [0.5 1 1.5], ...
            'XTick',                          [0.5 1 1.5]);
    
% colorbar
cBar = colorbar;
caxis([-0.61 0.61])
labels = -0.6:0.2:0.6;
set(cBar,'FontName','CMU Serif','YTickLabel',labels, 'FontSize', 14);
ylabel(cBar,'correlation')

%% FIG: 2C
percentages(percentages < 0.05) = NaN;

figure;
imagesc(alpha,alpha,percentages);

xlabel('$\alpha$-exponent (time series 1)', 'FontSize', 18); 
ylabel('$\alpha$-exponent (time series 2)', 'FontSize', 18) 
axis square; axis xy

set(gca, ...
            'FontSize',                               14 , ...
            'FontName',                       'CMU Serif', ...
            'YTick',                          [0.5 1 1.5], ...
            'XTick',                          [0.5 1 1.5]);

% colorbar
amin=0.05; amax=1;
dmap=(amax-amin)/size(map2,1);
colormap([nancolor; map2]);
caxis([amin-dmap amax])

cBar = colorbar;
ylim(cBar,[amin+0.005 amax])
set(cBar,'YTick',[0.05 0.2:0.2:1])
set(cBar,'YTickLabel',[0.05 0.2:0.2:1], 'FontSize', 14);
ylabel(cBar,'fraction of sign. correlation')



%% FIG: 2B
set(0,'defaulttextinterpreter','none')

%example values for alpha-combinations
alpha1 = [.76 .9  1.06 1.2];  
alpha2 = [.8  .96 1.1  1.26];

counter = 1;
h = figure; 

for x = 1:4
    for y = 1:4     
        hhh = subplot(4,4,counter); 
        i = find(alpha==alpha1(x)); 
        j = find(alpha==alpha2(y));
        
        hist(correlations(:,i,j),50)
        
        axis([-0.5 0.5 0 nr_repetitions/10])
        title(['$\alpha_1=' num2str(alpha(i)) ...
               ', \alpha_2=$' num2str(alpha(j))], ...
               'Interpreter', 'latex', 'FontSize', 14)
        
        if counter ~= 1
           set(hhh,'XTick', [])
           set(hhh,'YTick', [])
        end
         
        counter = counter+1;
    end
end

set(h, 'Position', [0 0 750 650])