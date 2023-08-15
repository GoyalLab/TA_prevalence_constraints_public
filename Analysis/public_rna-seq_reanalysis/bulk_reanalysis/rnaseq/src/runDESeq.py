import itertools
import os
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


def deseq(counts, targets, controls, conditions, dir):

    # process counts df
    counts = counts.groupby(['gene_name'], as_index=False).agg('sum')
    counts = counts.set_index('gene_name')

    # find control cols
    control_samples = [list(counts.filter(regex='^'+x+'$', axis=1)) for x in controls]
    control_samples = list(itertools.chain.from_iterable(control_samples))

    # filter out samples not matching condition
    samples_to_keep = sorted(list(set([x for x in list(counts) if (all(map(lambda y: y in x, conditions)) or conditions[0] == 'na')]+list(control_samples))))
    counts = counts[samples_to_keep]

    running_df = pd.DataFrame()
    # run deseq for each target
    for ko_gene in targets:
        # filter to correct samples
        ko_samples = counts.filter(regex='(^|[_])'+ko_gene+'([_]|$)', axis=1)
        if ko_samples.shape[1] == 0:
            print(ko_gene + " KO sample not found")
            continue
        ko_counts = counts[list(ko_samples)+list(control_samples)]

        # create clinical df
        clin_df = pd.DataFrame({'sample':list(ko_counts)})
        clin_df['condition'] = clin_df['sample'].apply(lambda x: 'control' if x in control_samples else 'KO')
        clin_df = clin_df.set_index('sample')

        # filter and format counts df
        ko_counts = ko_counts.transpose()
        genes_to_keep = ko_counts.columns[ko_counts.sum(axis=0) >= 10]
        ko_counts = ko_counts[genes_to_keep]

        # run deseq
        print('running '+ ko_gene)
        dds = DeseqDataSet(
            counts=ko_counts,
            clinical=clin_df,
            design_factors="condition",
            refit_cooks=True,
            n_cpus=None,
        )

        dds.deseq2()

        stat_res = DeseqStats(dds, contrast=["condition", "KO", "control"], n_cpus=None)
        stat_res.summary()
        df = stat_res.results_df
        df['sampleKO'] = [ko_gene]*len(df.index)
        
        if len(running_df) == 0:
            running_df = df
        else:
            running_df = pd.concat([running_df, df])
    
    if len(running_df) != 0:
        running_df.to_csv(dir+'differentialExpression_DESeq_allTargets.csv')



def main():
    dataset_metadata = pd.read_csv('./annotations/dataset_metadata.csv')
    dataset_directory = './raw_datasets/'
    deseq_directory = './deseq_files/'
    
    for x, row in dataset_metadata.iterrows():
        geo = row['GEO_ID']
        geo_unique = geo.split('-')[0]
        expression_data = pd.read_csv(dataset_directory+geo_unique+'.csv')

        targets = row['ko_genes'].split(',')
        controls = row['control_samples'].split(',')
        condition = row['condition'].split(',')
        prenorm = row['pre_normalized']
        if prenorm:
            continue
        
        path = deseq_directory+geo+'/'
        if not os.path.exists(path):
            os.mkdir(path)
        
        if os.path.exists(path+'differentialExpression_DESeq_allTargets.csv'):
            print(geo+ ': already run')
            continue
        print("starting " + geo)
        deseq(expression_data, targets, controls, condition, path)
        try:
            print('')
        except:
            print("issue with dataset: " + geo)
            continue

        print(str(x+1)+' of ' + str(len(dataset_metadata)) + ' datasets complete\n')


main()
