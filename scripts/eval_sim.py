#!/usr/bin/env python3

import pandas as pd
import argparse

def eval_sim(args):

    x = pd.read_csv(args.true_abundances, names=['txp_name', 'num_reads'], sep='\t')

    y1 = pd.read_csv(args.no_coverage, sep='\t')
    y1['txp_name'] = y1.tname.str.split('.').str.get(0)

    y2 = pd.read_csv(args.uniform_abundances, sep='\t')
    y2['txp_name'] = y2.tname.str.split('.').str.get(0)

    y3 = pd.read_csv(args.multinomial_abundances, sep='\t')
    y3['txp_name'] = y3.tname.str.split('.').str.get(0)

    y4 = pd.read_csv(args.Poisson_dis_abundances, sep='\t')
    y4['txp_name'] = y4.tname.str.split('.').str.get(0)

    y5 = pd.read_csv(args.Poisson_Con_abundances, sep='\t')
    y5['txp_name'] = y5.tname.str.split('.').str.get(0)

    m = pd.merge(x, y1[['txp_name', 'num_reads']], on='txp_name', how='outer', suffixes=('_true', '_NoCoverage'))
    m = pd.merge(m, y2[['txp_name', 'num_reads']], on='txp_name', how='outer', suffixes=('', '_Uniform'))
    m = pd.merge(m, y3[['txp_name', 'num_reads']], on='txp_name', how='outer', suffixes=('', '_multinomial'))
    m = pd.merge(m, y4[['txp_name', 'num_reads']], on='txp_name', how='outer', suffixes=('', '_dis_Poisson'))
    m = pd.merge(m, y5[['txp_name', 'num_reads']], on='txp_name', how='outer', suffixes=('', '_con_Poisson'))
    m = m.fillna(0.0)
    m.columns = ['txp_name', 'Ground_Truth', 'No_Coverage', 'Uniform', 'Multinomial', 'Dis_Poisson', 'Con_Poisson']

    numeric_cols = m.select_dtypes(include='number').columns
    correlation = m[numeric_cols].corr(method='spearman')
    print(correlation)
    #print(m.corr(method='spearman'))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            prog='eval-sim',
            description='evaluate accuracy of predicted abundances')
    

    parser.add_argument('true_abundances', type=argparse.FileType('r'), help='true abundances')
    parser.add_argument('no_coverage', type=argparse.FileType('r'), help='No Coverage')
    parser.add_argument('uniform_abundances', type=argparse.FileType('r'), help='uniform abundances')
    parser.add_argument('multinomial_abundances', type=argparse.FileType('r'), help='multinomial abundances')
    parser.add_argument('Poisson_dis_abundances', type=argparse.FileType('r'), help='discrete Poisson abundances')
    parser.add_argument('Poisson_Con_abundances', type=argparse.FileType('r'), help='continuous Poisson abundances')

    args = parser.parse_args()
    eval_sim(args)

