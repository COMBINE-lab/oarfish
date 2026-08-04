import sys, csv, numpy as np
sys.path.insert(0,'.')
from scipy.stats import rankdata
from metrics import read_truth, strip_key
base='/scratch1/rob/long-read-ecosystem/oarfish-evaluation-data/sim-panel'
SAMP=[('nanosim-NA12878-cdna','nanosim/ground_truth/cdna_ground_truth.csv'),
      ('nanosim-NA12878-drna','nanosim/ground_truth/drna_ground_truth.csv'),
      ('nanosim-H9-cdna','nanosim/ground_truth/H9_1DcDNA_ground_truth.csv'),
      ('nanosim-H9-drna','nanosim/ground_truth/H9_directRNA_ground_truth.csv'),
      ('tksm-RSII','tksm/ground_truth/RSII_ground_truth.csv'),
      ('tksm-SQ2','tksm/ground_truth/SQ2_ground_truth.csv')]
def q(s,a):
    ks,rd=[],[]
    for row in csv.DictReader(open('panelB/work/%s/%s.quant'%(a,s)),delimiter='\t'):
        ks.append(strip_key(row['tname'])); rd.append(float(row['num_reads']))
    return ks,np.array(rd)
def amb(s):
    u,a=[],[]
    for row in csv.DictReader(open('panelB/work/logistic/%s.ambig_info.tsv'%s),delimiter='\t'):
        u.append(float(row['unique_reads'])); a.append(float(row['ambig_reads']))
    return np.array(u),np.array(a)
def sp(a,b): return np.corrcoef(rankdata(a),rankdata(b))[0,1]
THR=[0.0,0.2,0.4,0.6,0.8,1.01]
print('Probe: use `none` where the transcript is well-determined (ambig frac <= t),')
print('`logistic` where it is ambiguous. t=0 is pure logistic; t=1.01 is pure none.')
print('%-22s %s'%('sample',' '.join('%-8s'%('t=%g'%t) for t in THR)))
agg={t:[] for t in THR}
for s,tp in SAMP:
    truth=read_truth(base+'/'+tp,'counts')
    ks,yl=q(s,'logistic'); _,yn=q(s,'none')
    u,a=amb(s); frac=a/(u+a+1e-9)
    x=np.array([truth.get(k,0.0) for k in ks])
    row=[]
    for t in THR:
        y=np.where(frac<=t,yn,yl); v=sp(x,y); row.append(v); agg[t].append(v)
    print('%-22s %s'%(s,' '.join('%-8.4f'%v for v in row)))
print('%-22s %s'%('MEAN',' '.join('%-8.4f'%np.mean(agg[t]) for t in THR)))
print()
print('pure logistic mean = %.4f ; best gated mean = %.4f (t=%g)'%(
    np.mean(agg[0.0]), max(np.mean(agg[t]) for t in THR),
    max(THR,key=lambda t:np.mean(agg[t]))))
