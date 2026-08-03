import sys, csv, numpy as np
sys.path.insert(0,'.')
from scipy.stats import rankdata
from metrics import read_truth, strip_key
base='/scratch1/rob/long-read-ecosystem/oarfish-evaluation-data/sim-panel'
def load(sample, arm):
    """quant rows in file order: (key, len, reads)"""
    ks,ln,rd=[],[],[]
    with open('panelB/work/%s/%s.quant'%(arm,sample)) as fh:
        r=csv.DictReader(fh,delimiter='\t')
        for row in r:
            ks.append(strip_key(row['tname'])); ln.append(float(row['len'])); rd.append(float(row['num_reads']))
    return ks,np.array(ln),np.array(rd)
def ambig(sample):
    u,a=[],[]
    with open('panelB/work/logistic/%s.ambig_info.tsv'%sample) as fh:
        r=csv.DictReader(fh,delimiter='\t')
        for row in r: u.append(float(row['unique_reads'])); a.append(float(row['ambig_reads']))
    return np.array(u),np.array(a)
def sp(a,b): return np.corrcoef(rankdata(a),rankdata(b))[0,1]
for sample,tp,lbl in [('nanosim-NA12878-cdna','nanosim/ground_truth/cdna_ground_truth.csv','ONT cDNA (logistic helps)'),
                      ('tksm-SQ2','tksm/ground_truth/SQ2_ground_truth.csv','PacBio HiFi (logistic hurts)')]:
    truth=read_truth(base+'/'+tp,'counts')
    ks,ln,_=load(sample,'logistic')
    x=np.array([truth.get(k,0.0) for k in ks])
    uq,am=ambig(sample)
    print('='*78); print('%s  -- %s'%(sample,lbl))
    print('  reads: unique=%.0f  ambiguous=%.0f  (%.1f%% of assigned mass is ambiguous)'
          %(uq.sum(),am.sum(),100*am.sum()/(uq.sum()+am.sum())))
    print('  %-34s spearman'%'estimator')
    print('  %-34s %.4f'%('unique reads only (no EM)',sp(x,uq)))
    print('  %-34s %.4f'%('unique+ambig, naive (no EM)',sp(x,uq+am)))
    for arm in ['none','logistic','adaptive-bare']:
        _,_,y=load(sample,arm); print('  %-34s %.4f'%('EM: '+arm,sp(x,y)))
    # residual structure vs length and ambiguity, logistic vs none
    _,_,y_log=load(sample,'logistic'); _,_,y_non=load(sample,'none')
    expr=x>0
    def resid(y): return np.log2((y[expr]+1)/(x[expr]+1))
    rl,rn=resid(y_log),resid(y_non)
    L=ln[expr]; A=(am/(uq+am+1e-9))[expr]
    print('  residual log2(est+1 / truth+1), mean by transcript-length quintile:')
    qs=np.percentile(L,[20,40,60,80])
    print('    %-10s %-9s %-9s %s'%('len bin','none','logistic','n'))
    for i in range(5):
        lo=-np.inf if i==0 else qs[i-1]; hi=np.inf if i==4 else qs[i]
        m=(L>lo)&(L<=hi)
        print('    %-10s %-+9.4f %-+9.4f %d'%('Q%d'%(i+1),rn[m].mean(),rl[m].mean(),m.sum()))
    print('  residual by ambiguous-fraction quintile:')
    aq=np.percentile(A,[20,40,60,80])
    print('    %-10s %-9s %-9s %s'%('ambig bin','none','logistic','n'))
    for i in range(5):
        lo=-np.inf if i==0 else aq[i-1]; hi=np.inf if i==4 else aq[i]
        m=(A>lo)&(A<=hi)
        if m.sum(): print('    %-10s %-+9.4f %-+9.4f %d'%('Q%d'%(i+1),rn[m].mean(),rl[m].mean(),m.sum()))
