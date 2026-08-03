import sys, numpy as np
sys.path.insert(0,'.')
from scipy.stats import rankdata
from metrics import read_quant, read_truth
base='/scratch1/rob/long-read-ecosystem/oarfish-evaluation-data/sim-panel'
SAMP=[('nanosim-NA12878-cdna','nanosim/ground_truth/cdna_ground_truth.csv'),
      ('nanosim-NA12878-drna','nanosim/ground_truth/drna_ground_truth.csv'),
      ('nanosim-H9-cdna','nanosim/ground_truth/H9_1DcDNA_ground_truth.csv'),
      ('nanosim-H9-drna','nanosim/ground_truth/H9_directRNA_ground_truth.csv'),
      ('tksm-RSII','tksm/ground_truth/RSII_ground_truth.csv'),
      ('tksm-SQ2','tksm/ground_truth/SQ2_ground_truth.csv')]
ARMS=['none','logistic','adaptive-bare','auto-new','+calib','+censor','+prune','+rank','auto-old']
def sp(a,b): return np.corrcoef(rankdata(a),rankdata(b))[0,1]
agg={a:{'low':[],'mid':[],'high':[],'expr':[],'fp':[],'glob':[]} for a in ARMS}
for s,tp in SAMP:
    truth=read_truth(base+'/'+tp,'counts')
    uni=sorted(read_quant('panelB/work/logistic/%s.quant'%s))
    x=np.array([truth.get(k,0.0) for k in uni]); expr=x>0; zero=~expr
    xs=x[expr]; q1,q2=np.percentile(xs,[33,66]); lab=np.where(xs<=q1,0,np.where(xs<=q2,1,2))
    for a in ARMS:
        q=read_quant('panelB/work/%s/%s.quant'%(a,s))
        y=np.array([q.get(k,0.0) for k in uni]); ye=y[expr]
        agg[a]['glob'].append(sp(x,y)); agg[a]['expr'].append(sp(xs,ye))
        for i,g in enumerate(['low','mid','high']):
            m=lab==i; agg[a][g].append(sp(xs[m],ye[m]))
        agg[a]['fp'].append((y[zero]>0).sum())
print('Mean over all 6 Panel B samples (Spearman unless noted)')
print('%-14s %-9s %-9s %-9s %-9s %-9s %s'%('arm','global','expressed','low','mid','high','FP zeros'))
for a in ARMS:
    d=agg[a]
    print('%-14s %-9.4f %-9.4f %-9.4f %-9.4f %-9.4f %.0f'%(a,np.mean(d['glob']),np.mean(d['expr']),
        np.mean(d['low']),np.mean(d['mid']),np.mean(d['high']),np.mean(d['fp'])))
print()
print('Delta vs logistic (expressed-only), per sample — is the kernel ever ahead?')
print('%-14s %s'%('arm',' '.join('%-9s'%s.split("-",1)[1][:8] for s,_ in SAMP)))
base_e={a:agg[a]['expr'] for a in ARMS}
for a in ARMS:
    if a=='logistic': continue
    print('%-14s %s'%(a,' '.join('%-+9.4f'%(base_e[a][i]-base_e['logistic'][i]) for i in range(len(SAMP)))))
