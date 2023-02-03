#!/usr/bin/env python 

import time
import resource
import os, sys, pandas, numpy, pathlib, getopt
import CytoSig

from scipy import stats
from statsmodels.stats.multitest import multipletests

err_tol = 1e-8
verbose_flag = False

fpath = pathlib.Path(__file__).parent.absolute()

pivots = ['TGFB1', 'TRAIL', 'PGE2']

def profile_geneset_signature(expression):
    signature = []
    fin = open(os.path.join(fpath, 'Tres.kegg'))
    for l in fin:
        fields = l.strip().split('\t')
        
        s = fields[2:]
        signature.append(pandas.Series(numpy.ones(len(s)), index=s, name=fields[0]))
    fin.close()
    
    signature = pandas.concat(signature, axis=1, join='outer', sort=False)
    signature.fillna(0, inplace=True)
    
    common = expression.index.intersection(signature.index)
    signature, expression = signature.loc[common], expression.loc[common]
    
    background = signature.mean(axis=1)
    background.name = 'study bias'
    
    X = signature.loc[:, ['KEGG_CELL_CYCLE', 'KEGG_DNA_REPLICATION']].mean(axis=1)
    X.name = 'Proliferation'
    
    X = pandas.concat([background, X], axis=1, join='inner')
    
    result = CytoSig.ridge_significance_test(X, expression, alpha=0, alternative="two-sided", nrand=0, verbose=verbose_flag)
    
    return result[2].loc['Proliferation']




def interaction_test(expression, X, y):
    signal = X.loc[:, 'pivot']
    
    failed = []
    merge = []
    
    for gid, arr in expression.iterrows():
        #if gid not in ['FIBP', 'TMEM222', 'IL7R']: continue
        
        X.loc[:,'partner'] = arr
        X.loc[:,'interaction'] = arr * signal
    
        # other covariates have no sufficient variation
        if arr.std() < err_tol or X.loc[:,'interaction'].std() < err_tol: continue
        
        try:
            y = pandas.DataFrame(y)
            result = CytoSig.ridge_significance_test(X, y, alpha=0, alternative="two-sided", nrand=0, flag_normalize=False, verbose=verbose_flag)
        
        except ArithmeticError:
            failed.append(gid)
            continue
        
        tvalue = result[2].loc['interaction'].iloc[0]
        pvalue = result[3].loc['interaction'].iloc[0]
        
        merge.append(pandas.Series([tvalue, pvalue], index=['t', 'p'], name=gid))    
    
    result = pandas.concat(merge, axis=1, join='inner').transpose()
    result['q'] = multipletests(result['p'], method='fdr_bh')[1]
    
    return result


def prediction(expression, output_file, count_thres):
    signature = pandas.read_csv(os.path.join(fpath, 'Tres.prediction_signature.gz'), sep='\t', index_col=0)['Tres']
    
    common = expression.index.intersection(signature.index)
    if common.shape[0] < count_thres:
        sys.stderr.write('Low overlap between input data genes and Tres signature genes.\n')
        sys.exit(1)
    
    print(common.shape[0], 'overlapping genes')
    
    expression = expression.loc[common]
    signature = signature.loc[common]
    
    merge = []
    
    try:
        for title, arr in expression.items():
            r, p = stats.pearsonr(arr, signature)
            merge.append(pandas.Series([r, p], index=['r', 'p'], name=title))
    except:
        sys.stderr.write('Pearson correlation failed.\n')
        sys.exit(1)
    
    result = pandas.concat(merge, axis=1, join='outer').transpose().sort_values('r')
    result.to_csv(output_file, sep='\t', index_label='ID')




def run_Tres(expression, output_file, count_thres):
    """
    run Tres model to get the signature
    """
    
    ###############################################################
    # read signature
    signature = os.path.join(sys.prefix, 'bin', 'signature.centroid.expand')
    
    if not os.path.exists(signature):
        sys.stderr.write('Cannot find signature file %s. Please make sure your CytoSig installation is successful.\n' % signature)
        sys.exit(1)
    else:
        signature = pandas.read_csv(signature, sep='\t', index_col=0)
    
    ###############################################################
    # compute signaling activity
    try:
        result_signaling = CytoSig.ridge_significance_test(signature, expression, alpha=1E4, verbose=verbose_flag)
    
    except ArithmeticError:
        sys.stderr.write('CytoSig regression failed.\n')
        sys.exit(1)
    
    # get the z-scores
    result_signaling = result_signaling[2]
    
    ###############################################################
    # compute proliferation signature
    try:
        result_prolifertion = profile_geneset_signature(expression)
    
    except ArithmeticError:
        sys.stderr.write('KEGG proliferation computation failed.\n')
        sys.exit(1)
    
    ###############################################################
    # Interaction test
    flag_group = ['.'.join(v.split('.')[:2]) for v in expression.columns]
    
    expression_group = expression.groupby(flag_group, axis=1)
    
    merge = []
    
    for title, expression_sub in expression_group:
        N = expression_sub.shape[1]
        if N < count_thres: continue
        
        print('process', title, N)
        
        # remove rows all zeros
        flag_nonzero = (expression_sub == 0).mean(axis=1) < 1
        if sum(~flag_nonzero) > 0: expression_sub = expression_sub.loc[flag_nonzero]
        
        y = result_prolifertion.loc[expression_sub.columns]
        
        # regression scaffold
        X = pandas.DataFrame(numpy.zeros((N,4)),
            columns = ['const', 'pivot', 'partner', 'interaction'],
            index = expression_sub.columns)
        
        X.loc[:, 'const'] = 1
        
        for pivot in pivots:
            X.loc[:,'pivot'] = result_signaling.loc[pivot, expression_sub.columns]
        
            result = interaction_test(expression_sub, X, y)
            result.columns += ('.%s.%s' % (title, pivot))
            merge.append(result)
    
    result = pandas.concat(merge, axis=1, join='outer')
    result.to_csv(output_file, sep='\t', index_label='ID')



def main():
    run_mode = 'run'
    input_file = output_file = None
    count_thres = 100
    normalization = False
        
    prompt_msg = 'Usage:\nTres.py -m <running mode (run|predict), default: %s> -i <input single-cell data> -o <output prefix> -c <count threshold, default: %d> -n <mean-normalization, default: %d>\n' % (run_mode, count_thres, normalization)
    
    try:
        opts, _ = getopt.getopt(sys.argv[1:], "hm:i:o:c:n:", [])
    
    except getopt.GetoptError:
        sys.stderr.write('Error input\n' + prompt_msg)
        sys.exit(2)
    
    if len(opts) == 0:
        sys.stderr.write('Please input some parameters or try: Tres.py -h\n')
        sys.exit(2)
    
    for opt, arg in opts:
        if opt == '-h':
            print(prompt_msg)
            sys.exit()
        
        elif opt in ("-m"):
            run_mode = arg
            
            if run_mode not in ['run', 'predict']:
                sys.stderr.write('Cannot recognize run mode %s. Please only use run or predict.\n' % run_mode)
                sys.exit(1)
        
        elif opt in ("-i"):
            input_file = arg
        
        elif opt in ("-o"):
            output_file = arg
        
        elif opt in ("-n"):
            normalization = (int(arg) != 0)
        
        elif opt in ("-c"):
            try:
                count_thres = int(arg)
            except:
                sys.stderr.write('count threshold %s is not a valid integer number.\n' % arg)
                sys.exit(1)
            
            if count_thres < 10:
                sys.stderr.write('count threshold %d < 10. Please input a number >= 10\n' % count_thres)
                sys.exit(1)
    
    if input_file is None:
        sys.stderr.write('Please provide a input file\n')
        sys.exit(1)
    
    elif not os.path.exists(input_file):
        sys.stderr.write('Cannot find input file %s\n' % input)
        sys.exit(1)
    
    if output_file is None:
        output_file = input_file + '.Tres_output'
        sys.stderr.write('No output file input. Automatically generate one as %s\n' % output_file)
    
    ###############################################################
    # read input
    try:
        f = os.path.basename(input_file)
        
        if f.find('.pickle') >= 0:
            expression = pandas.read_pickle(input_file)
        else:
            expression = pandas.read_csv(input_file, sep='\t', index_col=0)
    except:
        sys.stderr.write('Fail to open input file %s\n' % input_file)
        sys.exit(1)
    
    # gene and sample names must be unique
    assert expression.index.value_counts().max() == 1
    assert expression.columns.value_counts().max() == 1
    
    print('input matrix dimension', expression.shape)
    
    if normalization:
        background_expression = expression.mean(axis=1)
        expression = expression.subtract(background_expression, axis=0)
    
    if run_mode == 'run':
        run_Tres(expression, output_file, count_thres)    
    
    elif run_mode == 'predict':
        prediction(expression, output_file, count_thres)
    
    else: # impossible to get here
        assert False
    
    return 0



if __name__ == '__main__':
    time_start = time.perf_counter()
    main()
    time_elapsed = (time.perf_counter() - time_start)
    memMb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024.0/1024.0
    
    print ("%5.1f secs %5.1f MByte" % (time_elapsed,memMb))
