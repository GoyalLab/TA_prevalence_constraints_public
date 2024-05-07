% plot 1.6.5.1, 1.6.5.1.2, 1.6.5.1.3, 1.6.5.1.4, and 1.6.5.2 paramsets
% note that the constraints are on 
% basal_nitc_on_ratio
% A1_Aprime1_addon_ratio
% A1_Aprime_prodon_ratio

psc_sampdir521 = '~/Documents/grn_nitc_data/simulations/v1.6.5.2.1/samples/';
psc_tracedir521 = '~/Documents/grn_nitc_data/simulations/v1.6.5.2.1/fullTraces/';
psc_sampdir52 = '~/Documents/grn_nitc_data/simulations/v1.6.5.2/';

psc_outdir521 = '~/Documents/grn_nitc_data/simulations/v1.6.5.2.1/exploratory_analysis/';

% psc_sampdir52 = '/Volumes/IAMYG1/grn_nitc_data/v1.6.5.2/samples/';

% get parameter hypercube ranges
% range_table = readtable('~/Documents/grn_nitc_data/simulations/v1.6.5.1.2/more_larger_resamp_bounds.csv');
% range_table.custom_radius = [0.1;0.01;0.1;0.1;1;0.1;0.1;0.1];
% range_table.min_val_custom_log10 = log10(range_table.center - range_table.custom_radius);
% range_table.custom_range = log10(range_table.center + range_table.custom_radius) - range_table.min_val_custom_log10;

%% load lhs

lhs52 = readtable([psc_sampdir52, 'latinhyp_sampledSets_log10.csv']);
lhs521 = readtable([psc_sampdir521, 'latinhyp_sampledSets_log10.csv']);
ranges_t = readtable([psc_tracedir521, 'latinhyp_sampledSets_log10_ranges.csv']);
%% downsample for legibility (500 from 1.6.5.2 and 100 from 1.6.5.3)

rng(253);
lhs52_sub = lhs52(randsample(size(lhs52,1),500),:);
lhs521_sub = lhs521(randsample(size(lhs521,1),100),:);


colors = [repmat([0,0,0], 500, 1);repmat([0,0,1], 100, 1)];
lhs_subs = [lhs52_sub; lhs521_sub];

%% plot
fig1 = figure;
scatter3(lhs52_sub.basal_nitc_on_ratio, lhs52_sub.A1_Aprime1_addon_ratio, lhs52_sub.A1_Aprime_prodon_ratio, 'filled', 'black')
hold on
scatter3(lhs521_sub.basal_nitc_on_ratio, lhs521_sub.A1_Aprime1_addon_ratio, lhs521_sub.A1_Aprime_prodon_ratio, 'filled', 'blue')
plotcube([ranges_t.basal_nitc_on_ratio(1), ranges_t.A1_Aprime1_addon_ratio(1), ranges_t.A1_Aprime_prodon_ratio(1)], [0.324, -1, -0.468], 0.1, [0,0,1]);
title({'Parameter spaces for unimodal symmetric resampling analysis','Downsampled for legibility','Full space (black), Decision tree subspace (blue)'})
xlabel('basal nitc on ratio')
ylabel('A1 Aprime1 addon ratio')
zlabel('A1 Aprime prodon ratio')
hold off
%%
exportgraphics(fig1, [psc_outdir521, 'subspace_parametersets.pdf'])
exportgraphics(fig1, [psc_outdir521, 'subspace_parametersets.eps'])
saveas(fig1, [psc_outdir521, 'subspace_parametersets.svg'],'svg')
saveas(fig1, [psc_outdir521, 'subspace_parametersets_sa.pdf'],'pdf')
%%
saveas(fig1, [psc_outdir521, 'subspace_parametersets_painters.svg'])
