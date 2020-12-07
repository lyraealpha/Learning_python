#!/usr/bin/env python3.8
# -*- coding: utf-8 -*-

import os
import sys
import math
import time
from scipy import stats
import argparse

def timer(defined_function):  # time decorator
    def wrapper(*args, **kwargs):
        time_start = time.time()
        wrapper_out = defined_function(*args, **kwargs)
        time_end = time.time()
        print('Done !\nRuning script {} finished: {}s elapsed.'.format(sys.argv[0], round((time_end - time_start), 4)))
        return wrapper_out

    return wrapper


class CorrelationClculator:

    def __init__(self, lncRNA_exp, mRNA_exp, method, outpath, key, cor, pvalue):
        self.mRNA_exp, self.lncRNA_exp, self.outpath =[os.path.abspath(x) for x in [mRNA_exp, lncRNA_exp, outpath]]
        self.key = key
        self.cor = cor
        self.pvalue = pvalue
        self.cal_cor = 0
        if method == 'pearson':
            self.cal_cor = stats.pearsonr
        elif method == 'spearman':
            self.cal_cor = stats.spearmanr
        elif method == 'kendall':
            self.cal_cor = stats.kendalltau
        else:
            raise Exception('method must be pearson or spearman or kendall')

    @timer
    def cal_cor_pval(self):
        ofile = f'{self.outpath}/co_expression.{self.key}.xls'
        with open(self.mRNA_exp, 'r') as mRNA,\
            open(self.lncRNA_exp, 'r') as lnc,\
            open(ofile, 'w') as out:
            out.write(f'RNA1\tRNA2\tcoefficient\tpvalue\n')
            mRNA_lines = [x for x in mRNA.readlines() if not x.startswith('#')]
            lnc_lines = [x for x in lnc.readlines() if not x.startswith('#')]
            for lnc_line in lnc_lines:
                lnc_line = lnc_line.strip().split()
                lnc_id, lnc_exp = lnc_line[0], [float(x) for x in lnc_line[1:]]
                for mRNA_line in mRNA_lines:
                    mRNA_line = mRNA_line.strip().split()
                    mRNA_id, mRNA_exp =mRNA_line[0], [float(x) for x in mRNA_line[1:]]
                    if mRNA_id == lnc_id:
                        continue
                    correlation, pval = self.cal_cor(lnc_exp, mRNA_exp)
                    if any([math.isnan(correlation), math.isnan(pval)]):
                        continue
                    if abs(correlation) >= self.cor and pval <= self.pvalue:
                        out.write(f'{lnc_id}\t{mRNA_id}\t{correlation}\t{pval}\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description='Script')
    parser.add_argument('--rna1',
                        type=str,
                        dest='lncRNA_exp',
                        nargs='?',
                        help='lncRNA_exp file')
    parser.add_argument('--rna2',
                        type=str,
                        dest='mRNA_exp',
                        nargs='?',
                        help='mRNA_exp file')
    parser.add_argument('--outpath',
                        type=str,
                        dest='outpath',
                        nargs='?',
                        help='result dir')
    parser.add_argument('--key',
                        type=str,
                        dest='key',
                        nargs='?',
                        default='relation',
                        help='output keyword,default:relation')
    parser.add_argument('--method',
                        type=str,
                        dest='method',
                        nargs='?',
                        default='pearson',
                        help='cor-relation method(pearson, kendall, spearman),default: pearson')
    parser.add_argument('--cor',
                        type=float,
                        dest='cor',
                        nargs='?',
                        default=0.0,
                        help='cor-relation default: 0.0')
    parser.add_argument('--pvalue',
                        type=float,
                        dest='pvalue',
                        nargs='?',
                        default=1.0,
                        help='pvalue cutoff , default: 1.0')
    args = parser.parse_args()
    if not all([args.lncRNA_exp, args.mRNA_exp, args.outpath]):
        parser.print_help()
        sys.exit(1)
    CorrelationClculator(args.lncRNA_exp, args.mRNA_exp, args.method,
                         args.outpath, args.key, args.cor, args.pvalue).cal_cor_pval()
