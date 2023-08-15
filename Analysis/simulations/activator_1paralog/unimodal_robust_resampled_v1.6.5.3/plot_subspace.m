% plot 1.6.5.3 and 1.6.5.2 paramsets
% note that the constraints are on 
% basal_nitc_on_ratio
% A1_Aprime1_addon_ratio
% A1_Aprime_prodon_ratio

psc_sampdir52 = '/Volumes/IAMYG1/grn_nitc_data/v1.6.5.2/samples/';
psc_sampdir53 = '/Volumes/IAMYG2/grn_nitc_data/v1.6.5.3/samples/';
psc_outdir53 = '/Volumes/IAMYG2/grn_nitc_data/v1.6.5.3/exploratory_analysis/';

%% load lhs

lhs52 = readtable([psc_sampdir52, 'latinhyp_sampledSets_log10.csv']);
lhs53 = readtable([psc_sampdir53, 'latinhyp_sampledSets_log10.csv']);

%% downsample for legibility (500 from 1.6.5.2 and 100 from 1.6.5.3)

rng(253);
lhs52_sub = lhs52(randsample(size(lhs52,1),500),:);
lhs53_sub = lhs53(randsample(size(lhs53,1),100),:);

colors = [repmat([0,0,0], 500, 1);repmat([0,0,1], 100, 1)];
lhs_subs = [lhs52_sub; lhs53_sub];

%% plot
fig1 = figure;
scatter3(lhs52_sub.basal_nitc_on_ratio, lhs52_sub.A1_Aprime1_addon_ratio, lhs52_sub.A1_Aprime_prodon_ratio, 'filled', 'black')
hold on
scatter3(lhs53_sub.basal_nitc_on_ratio, lhs53_sub.A1_Aprime1_addon_ratio, lhs53_sub.A1_Aprime_prodon_ratio, 'filled', 'blue')
plotcube([1.1192, (1-0.3166), 1.4310], [-1, 0.3166, -0.4310], 0.1, [0,0,1]);
title({'Parameter spaces for unimodal robustness analysis','Downsampled for legibility','Full space (black), Decision tree subspace (blue)'})
xlabel('basal nitc on ratio')
ylabel('A1 Aprime1 addon ratio')
zlabel('A1 Aprime prodon ratio')
hold off
saveas(fig1, [psc_outdir53, 'subspace_parametersets.svg'],'svg')
