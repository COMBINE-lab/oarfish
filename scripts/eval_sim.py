#!/usr/bin/env python3

import pandas as pd
import argparse

def eval_sim(args):
    x = pd.read_csv(args.true_abundances, names=['txp_name', 'num_reads'], sep='\t')
    y = pd.read_csv(args.predicted_abundances, sep='\t')
    y['txp_name'] = y.tname.str.split('.').str.get(0)

    m = pd.merge(x, y[['txp_name', 'num_reads']], on='txp_name', how='outer')
    m = m.fillna(0.0)

    numeric_cols = m.select_dtypes(include='number').columns
    correlation = m[numeric_cols].corr(method='spearman')
    print(correlation)
    #print(m.corr(method='spearman'))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            prog='eval-sim',
            description='evaluate accuracy of predicted abundances')
    

    parser.add_argument('true_abundances', type=argparse.FileType('r'), help='true abundances')
    parser.add_argument('predicted_abundances', type=argparse.FileType('r'), help='predicted abundances')

    args = parser.parse_args()
    eval_sim(args)

