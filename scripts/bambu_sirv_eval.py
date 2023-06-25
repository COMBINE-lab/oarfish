# coding: utf-8
import pandas as pd
import numpy as np

def mae(y_true, predictions):
    y_true, predictions = np.array(y_true), np.array(predictions)
    return np.mean(np.abs(y_true - predictions))

def main():
    e = pd.read_excel('../../data/bambu_paper/ground_truth/nanopore_sample_sheet.xlsx')
    e = e.loc[pd.notna(e['Spike in concentration']), :]
    t = pd.read_csv('../../data/bambu_paper/ground_truth/spike_in.csv')
    t['conc_norm'] = (t.conc / t.conc.sum()) 

    xd = {}
    prefix_d = '../../data/bambu_paper/quants/'
    for sn in e.Sample:
        m = 'oarfish'
        n = '/'.join([prefix_d, m, sn])
        k = (sn, m)
        xd[k] = pd.read_csv(f'{n}.tsv', sep='\t')
        xd[k]['num_reads_norm'] = (xd[k].num_reads / xd[k].num_reads.sum()) 
        m = 'nanocount'
        n = '/'.join([prefix_d, m, sn])
        k = (sn, m)
        xd[k] = pd.read_csv(f'{n}.tsv', sep='\t')
        xd[k]['num_reads_norm'] = (xd[k].est_count / xd[k].est_count.sum()) 
 
    res = []
    for k,v in xd.items():
        if k[1] == 'oarfish':
            mt = pd.merge(v, t, left_on='tname', right_on='tx_name', how='inner')
        else:
            mt = pd.merge(v, t, left_on='transcript_name', right_on='tx_name', how='inner')
        c = mt.loc[(mt.tx_name.str.startswith('R') & (mt.spike_in == 'sequinMixAV1')) ,:].corr(method='spearman').loc[ 'conc_norm', 'num_reads_norm' ]
        subset = mt.loc[(mt.tx_name.str.startswith('R') & (mt.spike_in == 'sequinMixAV1')), :].fillna(0.0)
        maev = mae(subset['conc_norm'].values, subset['num_reads_norm'].values)
        res.append((k, c, maev))
    res = pd.DataFrame.from_records(res, columns=['Sample', 'correlation', 'mae'])
    print(res)

if __name__ == "__main__":
    main()
