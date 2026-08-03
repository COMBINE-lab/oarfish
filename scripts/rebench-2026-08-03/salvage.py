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
ARMS=['logistic','adaptive-bare','auto-new']
THR=[0.0,0.25,0.5,1.0,2.0,5.0]
def sp(a,b): return np.corrcoef(rankdata(a),rankdata(b))[0,1]
res={a:{t:[] for t in THR} for a in ARMS}
for s,tp in SAMP:
    truth=read_truth(base+'/'+tp,'counts')
    uni=sorted(read_quant('panelB/work/logistic/%s.quant'%s))
    x=np.array([truth.get(k,0.0) for k in uni])
    for a in ARMS:
        q=read_quant('panelB/work/%s/%s.quant'%(a,s))
        y0=np.array([q.get(k,0.0) for k in uni])
        for t in THR:
            y=np.where(y0<t,0.0,y0)
            res[a][t].append(sp(x,y))
print('Global Spearman after zeroing estimates below a threshold')
print('(same post-hoc filter applied to every arm; mean over 6 Panel B samples)')
print('%-14s %s'%('arm',' '.join('%-9s'%('t=%g'%t) for t in THR)))
for a in ARMS:
    print('%-14s %s'%(a,' '.join('%-9.4f'%np.mean(res[a][t]) for t in THR)))
print()
print('delta adaptive-bare - logistic:')
print('%-14s %s'%('',' '.join('%-+9.4f'%(np.mean(res['adaptive-bare'][t])-np.mean(res['logistic'][t])) for t in THR)))
print('delta auto-new     - logistic:')
print('%-14s %s'%('',' '.join('%-+9.4f'%(np.mean(res['auto-new'][t])-np.mean(res['logistic'][t])) for t in THR)))
