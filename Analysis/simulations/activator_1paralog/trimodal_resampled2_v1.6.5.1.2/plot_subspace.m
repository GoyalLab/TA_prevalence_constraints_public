% plot 1.6.5.1, 1.6.5.1.2, 1.6.5.1.3, 1.6.5.1.4, and 1.6.5.2 paramsets
% note that the constraints are on 
% basal_nitc_on_ratio
% A1_Aprime1_addon_ratio
% A1_Aprime_prodon_ratio

psc_sampdir51 = '~/Documents/grn_nitc_data/simulations/v1.6.5.1/';
psc_sampdir52 = '~/Documents/grn_nitc_data/simulations/v1.6.5.2/';
psc_sampdir512 = '~/Documents/grn_nitc_data/simulations/v1.6.5.1.2/samples/';
psc_sampdir513 = '~/Documents/grn_nitc_data/simulations/v1.6.5.1.3/samples/';
psc_sampdir514 = '~/Documents/grn_nitc_data/simulations/v1.6.5.1.4/samples/';

psc_outdir512 = '~/Documents/grn_nitc_data/simulations/v1.6.5.1.2/exploratory_analysis/';

% psc_sampdir52 = '/Volumes/IAMYG1/grn_nitc_data/v1.6.5.2/samples/';

% get parameter hypercube ranges
range_table = readtable('~/Documents/grn_nitc_data/simulations/v1.6.5.1.2/more_larger_resamp_bounds.csv');
range_table.custom_radius = [0.1;0.01;0.1;0.1;1;0.1;0.1;0.1];
range_table.min_val_custom_log10 = log10(range_table.center - range_table.custom_radius);
range_table.custom_range = log10(range_table.center + range_table.custom_radius) - range_table.min_val_custom_log10;

%% load lhs

lhs52 = readtable([psc_sampdir52, 'latinhyp_sampledSets_log10.csv']);
lhs51 = readtable([psc_sampdir51, 'latinhyp_sampledSets.csv']);
lhs512 = readtable([psc_sampdir512, 'latinhyp_sampledSets_log10.csv']);
lhs513 = readtable([psc_sampdir513, 'latinhyp_sampledSets_log10.csv']);
lhs514 = readtable([psc_sampdir514, 'latinhyp_sampledSets_log10.csv']);

%% downsample for legibility (500 from 1.6.5.2 and 100 from 1.6.5.3)

rng(253);
lhs52_sub = lhs52(randsample(size(lhs52,1),500),:);
lhs51_sub = lhs51(randsample(size(lhs51,1),100),:);
lhs512_sub = lhs512(randsample(size(lhs512,1),100),:);
lhs513_sub = lhs513(randsample(size(lhs513,1),100),:);
lhs514_sub = lhs514(randsample(size(lhs514,1),100),:);

colors = [repmat([0,0,0], 500, 1);repmat([0,0,1], 100, 1)];
lhs_subs = [lhs52_sub; lhs512_sub];

%% plot
fig1 = figure;
scatter3(lhs52_sub.basal_nitc_on_ratio, lhs52_sub.A1_Aprime1_addon_ratio, lhs52_sub.A1_Aprime_prodon_ratio, 'filled', 'black')
hold on
scatter3(log10(lhs51_sub.basal_nitc_on_ratio), log10(lhs51_sub.A1_Aprime1_addon_ratio), log10(lhs51_sub.A1_Aprime_prodon_ratio), 'filled', 'blue')
plotcube([range_table.custom_range(1), range_table.custom_range(3), range_table.custom_range(4)], [range_table.min_val_custom_log10(1), range_table.min_val_custom_log10(3), range_table.min_val_custom_log10(4)], 0.1, [0,0,1]);
scatter3(lhs512_sub.basal_nitc_on_ratio, lhs512_sub.A1_Aprime1_addon_ratio, lhs512_sub.A1_Aprime_prodon_ratio, 'filled', 'blue')
plotcube([range_table.min_range(1), range_table.min_range(3), range_table.min_range(4)], [range_table.min_val_min_log10(1), range_table.min_val_min_log10(3), range_table.min_val_min_log10(4)], 0.1, [0,0,1]);
scatter3(lhs513_sub.basal_nitc_on_ratio, lhs513_sub.A1_Aprime1_addon_ratio, lhs513_sub.A1_Aprime_prodon_ratio, 'filled', 'blue')
plotcube([range_table.perc10_range(1), range_table.perc10_range(3), range_table.perc10_range(4)], [range_table.min_val_10_log10(1), range_table.min_val_10_log10(3), range_table.min_val_10_log10(4)], 0.1, [0,0,1]);
scatter3(lhs514_sub.basal_nitc_on_ratio, lhs514_sub.A1_Aprime1_addon_ratio, lhs514_sub.A1_Aprime_prodon_ratio, 'filled', 'blue')
plotcube([range_table.perc50_range(1), range_table.perc50_range(3), range_table.perc50_range(4)], [range_table.min_val_50_log10(1), range_table.min_val_50_log10(3), range_table.min_val_50_log10(4)], 0.1, [0,0,1]);
title({'Parameter spaces for trimodal resampling analysis','Downsampled for legibility','Full space (black), Nested decision tree subspaces (blue)'})
xlabel('basal nitc on ratio')
ylabel('A1 Aprime1 addon ratio')
zlabel('A1 Aprime prodon ratio')
hold off
%%
exportgraphics(fig1, [psc_outdir512, 'subspace_parametersets.pdf'])
exportgraphics(fig1, [psc_outdir512, 'subspace_parametersets.eps'])
%%
saveas(fig1, [psc_outdir512, 'subspace_parametersets.svg'],'svg')
saveas(fig1, [psc_outdir512, 'subspace_parametersets.pdf'],'pdf')

%%
saveas(fig1, [psc_outdir512, 'subspace_parametersets_painters.svg'],'svg')