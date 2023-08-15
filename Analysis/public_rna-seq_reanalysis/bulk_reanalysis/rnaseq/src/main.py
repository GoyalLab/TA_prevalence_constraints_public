import random
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import sem

random.seed(10)
BOOTSTRAP_REPS = 10000

def prepare_for_plot(df, filter_params, pos_only=False, sig_only=False):
    # add in support for p-value coloring
    p_max = filter_params['p_max']
    df['significant'] = df['padj'].apply(lambda x: 1 if x<=p_max else 0)

    # add in support for percent upreg
    fc2_min = filter_params['log2fc_min']
    df['upreg'] = df[['significant', 'FC2']].apply(lambda x: 1 if (x['FC2'] >= fc2_min and x['significant'] == 1) else 0, axis=1)

    if pos_only:
        # filter out any average negative expression
        summarized = df.groupby('KO_Gene').mean()
        relevant_samples = summarized[summarized['FC2'] > 0.1]
        df = df[df['KO_Gene'].isin(list(relevant_samples.index))]

    if sig_only:
        # filter out any average non-significant expression
        summarized = df.groupby('KO_Gene').sum(numeric_only=True)
        relevant_samples = summarized[summarized['significant'] > 0]
        df = df[df['KO_Gene'].isin(list(relevant_samples.index))]

    return df


def filter_data(df, filter_params):
    # filter by base mean
    min_basemean = filter_params['baseMean_min']
    df = df[df['baseMean'] > min_basemean]

    return df


def gen_pvals(df, boot_dists):
    summarized = df.groupby('KO_Gene').mean(numeric_only=True).reset_index()
    data = []
    for _, row in summarized.iterrows():
        dist = boot_dists[row['KO_Gene']]
        perc_upreg_obs = row['upreg']
        perc_greater_than_obs = len(
            [i for i in dist if i >= perc_upreg_obs]) / len(dist)
        data.append([row['KO_Gene'], perc_greater_than_obs])

    df = pd.DataFrame(data, columns=['KO_Gene', 'p_value'])
        
    return df


def plot_data(df, boot_df, boot_dists, filter_params, geo_id, output_dir):
    if len(df) == 0:
        print('* Warning: '+geo_id+' empty dataframe pre-filtering')
        return
    fig_size = (13,8.5) 
    if 'Base_Mean' in df: # add in significance coloring
        hue = 'significant'
        palette = 'gist_gray_r'
        color = None
    else:
        color= None #'black'
        hue = 'Perc_ID'
        palette = None
    
    df = prepare_for_plot(df, filter_params, pos_only=False, sig_only=True)
    if len(df) == 0:
        print('* Warning: '+geo_id+' empty dataframe post-filtering')
        return
    df = df.sort_values(by="KO_Gene")
    if len(set(df['KO_Gene'])) > 17:
        rotation = 45
        bottom = True
    else:
        rotation = 0
        bottom = False

    plt.figure(figsize=fig_size)
    sns.set_theme(context='notebook', style='whitegrid', palette='deep',
                      font='sans-serif', font_scale=1, color_codes=True, rc=None)
    ax1 = sns.stripplot(x='KO_Gene', y='FC2', data=df, color=color, palette=palette, hue=hue)
    sns.barplot(x="KO_Gene", y="FC2", errorbar=None,color='gray', data=df, ax=ax1, alpha=.3)
    plt.xticks(rotation=rotation)
    ax1.tick_params(bottom=bottom)
    plt.savefig(output_dir + 'FC2/' + geo_id + '.png', dpi=200)
    plt.clf()
    df.to_csv(output_dir + 'FC2/tabular_format/'+geo_id+'.csv')

    # setup df for percent upregulated
    pval_df = gen_pvals(df, boot_dists)
    boot_df = boot_df.groupby(['KO_Gene']).mean(numeric_only=True).reset_index()
    combined_df = boot_df.merge(df.groupby(['KO_Gene'])['upreg'].mean(numeric_only=True), left_on='KO_Gene', right_on='KO_Gene')
    combined_df = combined_df.rename({'upreg': 'Paralogs', 'boot_upreg': 'Bootstrapped Genes'}, axis=1)
    upreg_data = pd.melt(combined_df, id_vars=['KO_Gene'], value_vars=['Paralogs', 'Bootstrapped Genes',])
    upreg_data = upreg_data.rename({'value':'Percent Upregulated'}, axis=1)
    upreg_data = upreg_data.merge(boot_df, how='left', left_on=['KO_Gene', 'Percent Upregulated'], right_on=['KO_Gene', 'boot_upreg'])
    upreg_data = upreg_data.merge(pval_df, how='left', left_on='KO_Gene', right_on='KO_Gene')

    # plot percent upregulated
    plt.figure(figsize=fig_size)
    sns.set_theme(context='notebook', style='whitegrid', palette='deep',
                      font='sans-serif', font_scale=1, color_codes=True, rc=None)
    ax2 = sns.barplot(x="KO_Gene", y="Percent Upregulated", palette='Paired',
                      hue='variable', data=upreg_data)
    x_coords = [p.get_x() + 0.5*p.get_width() for p in ax2.patches]
    y_coords = [p.get_height() for p in ax2.patches]
    plt.errorbar(x=x_coords, y=y_coords, yerr=upreg_data["SE"], fmt="none", c="k")
    plt.ylim(0, 1)
    plt.xticks(rotation=rotation)
    ax2.tick_params(bottom=bottom)
    p_val_df = upreg_data.groupby(['KO_Gene']).max().reset_index()
    for _, row in p_val_df.iterrows():
        text = '-'+str(round(row.p_value, 2))+'-' if row.p_value <=0.05 else ''
        ax2.text(row.name, row['Percent Upregulated']+0.01, text,
                color='black', ha='center')
    plt.savefig(output_dir + 'percent_upregulated/' + geo_id + '.png', dpi=200)
    plt.clf()
    upreg_data.to_csv(
        output_dir + 'percent_upregulated/tabular_format/'+geo_id+'.csv')

    plt.close('all')

    return


def process_deseq_data(deseq_file, paralog_data, ko_genes, filter_params):
    fc_data = []
    deseq_file['baseMean'] = pd.to_numeric(deseq_file['baseMean'], errors='coerce')
    deseq_file.dropna(subset=['baseMean'], inplace=True)
    deseq_file = filter_data(deseq_file, filter_params)
    bootstrap_data = []
    bootstrap_dists = {}
    for gene in ko_genes:
        paralogs = list(paralog_data[paralog_data['Gene'] == gene]['Paralog'])
        perc_ids = list(paralog_data[paralog_data['Gene'] == gene]['Perc_ID'])

        gene_data = deseq_file[deseq_file['sampleKO']==gene]
        gene_data = gene_data.set_index('gene_name')
        gene_data = gene_data.sort_values('baseMean')
        gene_data['rank'] = gene_data['baseMean'].rank()

        # confirm KO
        ko_data = gene_data.filter(regex='(^|[_])'+gene+'([_]|$)', axis=0)
        
        if len(ko_data) > 1:
            ko_data = ko_data.replace("NA", np.nan)
            ko_data = ko_data.dropna()
        if len(ko_data) != 1:
            ko_fc = np.nan
            print('Warning: ' + gene + ' KO gene not found')
        else:
            ko_fc = float(ko_data['log2FoldChange'][0])

        bootstrap_fc2_data = [[] for _ in range(len(paralogs))]
        bootstrap_padj_data = [[] for _ in range(len(paralogs))]
        for x, paralog in enumerate(paralogs):
            try:
                par_data = gene_data.filter(regex='(^|[_])'+paralog+'([_]|$)', axis=0)
                if len(par_data) > 1:
                    par_data = par_data.replace("NA", np.nan)
                    par_data = par_data.dropna()
                assert len(par_data) == 1
            except:
                print('Warning: ' + paralog + ' paralog not found')
                continue

            par_fc2 = float(par_data['log2FoldChange'][0])
            par_basemean = float(par_data['baseMean'][0])
            par_padj = float(par_data['padj'][0])

            rank = int(par_data['rank'][0])
            for _ in range(0,BOOTSTRAP_REPS):
                bootstrap_idx = random.randint(*random.choice([(rank-50, rank-1), (rank+1, rank+50)]))
                bootstrap_idx = 0 if (bootstrap_idx < 0) else bootstrap_idx
                bootstrap_idx = -1 if (bootstrap_idx >= len(gene_data)) else bootstrap_idx
                bootstrap_fc2_data[x].append(float(gene_data.iloc[bootstrap_idx]['log2FoldChange']))
                bootstrap_padj_data[x].append(float(gene_data.iloc[bootstrap_idx]['padj']))

            fc_data.append([gene, paralog, perc_ids[x], par_fc2, par_basemean, par_padj, ko_fc])
        
        if len(bootstrap_fc2_data) > 0:
            bootstrap_stats, bootstrap_dist = process_bootstrap_data(
                gene, gene, bootstrap_fc2_data, bootstrap_padj_data, filter_params)
            bootstrap_data.append(bootstrap_stats)
            bootstrap_dists[gene] = bootstrap_dist

    df = pd.DataFrame(fc_data, columns=['KO_Gene', 'Paralog', 'Perc_ID', 'FC2', 'Base_Mean',
                      'padj', 'KO_FC2'])
    boot_df = pd.DataFrame(bootstrap_data, columns=['KO_Gene', 'Sample', 'boot_upreg', 'SE','StdDev'])

    return df, boot_df, bootstrap_dists


def process_bootstrap_data(gene, sample, fc2_data, padj_data, filter_params):
    # analyze the boostrapping samples
    # first remove paralogs that weren't found
    fc2_data = [l for l in fc2_data if len(l) != 0]
    padj_data = [l for l in padj_data if len(l) != 0]
    if len(fc2_data) == 0:
        return [gene, sample, np.nan, np.nan], []

    # start calculating averages
    perc_upregs = []
    for i in range(0, BOOTSTRAP_REPS):
        upreg_sum = 0
        for j, paralog_samples in enumerate(fc2_data):
            if (paralog_samples[i] > filter_params['log2fc_min']) and (padj_data[j][i] <= filter_params['p_max']):
                upreg_sum += 1

        perc_upreg = upreg_sum / len(fc2_data)
        perc_upregs.append(perc_upreg)

    avg_perc_upreg = sum(perc_upregs) / len(perc_upregs)
    se_perc_upreg = sem(perc_upregs)
    std_perc_upreg = np.std(perc_upregs)
    
    return [gene, sample, avg_perc_upreg, se_perc_upreg, std_perc_upreg], perc_upregs


def run_dataset(row, filters, paralog_directory, norm_dataset_directory, output_directory):
    paralog_data = pd.read_csv(paralog_directory+row['GEO_ID']+'-paralogs.csv')
    ko_genes = row['ko_genes'].split(',')
    geo = row['GEO_ID']

    print('starting ' + geo)
    
    seq_directory = norm_dataset_directory + geo + '/' + 'differentialExpression_DESeq_allTargets.csv'
    try:
        seq_file = pd.read_csv(seq_directory)
    except FileNotFoundError:
        print('Warning: ' + row['GEO_ID'] + ' DESeq not found')
        return
    if len(seq_file.columns) != 8:
        print('Warning: ' + row['GEO_ID'] + ' DESeq file empty')
        print(len(seq_file.columns))
        return
    
    for filter_name, filter_params in filters.items():
        df, boot_df, boot_dists = process_deseq_data(seq_file, paralog_data, ko_genes, filter_params)
        plot_data(df, boot_df, boot_dists, filter_params, row['GEO_ID'], output_directory+filter_name+'/')
    
    return


def main():
    dataset_metadata = pd.read_csv('./annotations/dataset_metadata.csv')
    paralog_directory = './paralog_data/'
    norm_dataset_directory = './deseq_files/'
    output_directory = './de_analysis/'

    filters = {'all_data':      {'p_max': 1,    'log2fc_min': 0,   'baseMean_min': 0},
               'p_only':        {'p_max': 0.05, 'log2fc_min': 0,   'baseMean_min': 0},
               'p_fc':          {'p_max': 0.05, 'log2fc_min': 0.5, 'baseMean_min': 0},
               'p_fc_basemean': {'p_max': 0.05, 'log2fc_min': 0.5, 'baseMean_min': 10}}

    for _, row in dataset_metadata.iterrows():
        run_dataset(row, filters, paralog_directory, norm_dataset_directory, output_directory)
    
    return



main()
