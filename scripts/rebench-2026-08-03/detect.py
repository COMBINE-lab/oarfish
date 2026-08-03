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
    ks,ln,rd=[],[],[]
    for row in csv.DictReader(open('panelB/work/%s/%s.quant'%(a,s)),delimiter='\t'):
        ks.append(strip_key(row['tname'])); ln.append(float(row['len'])); rd.append(float(row['num_reads']))
    return ks,np.array(ln),np.array(rd)
def amb(s):
    u,a=[],[]
    for row in csv.DictReader(open('panelB/work/logistic/%s.ambig_info.tsv'%s),delimiter='\t'):
        u.append(float(row['unique_reads'])); a.append(float(row['ambig_reads']))
    return np.array(u),np.array(a)
def sp(a,b): return np.corrcoef(rankdata(a),rankdata(b))[0,1]
RULES=['as shipped','global t=0.5','global t=1.0','no unique read & <1','no unique read & <5',
       'no unique read & <20','reads/kb < 0.5','oracle: perfect detection']
res={r:[] for r in RULES}
for s,tp in SAMP:
    truth=read_truth(base+'/'+tp,'counts')
    ks,ln,y=q(s,'logistic'); u,_=amb(s)
    x=np.array([truth.get(k,0.0) for k in ks]); z=x==0
    def ev(yy): return sp(x,yy)
    res['as shipped'].append(ev(y))
    res['global t=0.5'].append(ev(np.where(y<0.5,0,y)))
    res['global t=1.0'].append(ev(np.where(y<1.0,0,y)))
    for t,lbl in [(1,'no unique read & <1'),(5,'no unique read & <5'),(20,'no unique read & <20')]:
        res[lbl].append(ev(np.where((u==0)&(y<t),0,y)))
    res['reads/kb < 0.5'].append(ev(np.where(y/(ln/1000.0)<0.5,0,y)))
    yz=y.copy(); yz[z]=0; res['oracle: perfect detection'].append(ev(yz))
print('Detection rules applied post-EM to the logistic estimate (mean over 6 Panel B samples)')
print('%-28s %-9s %s'%('rule','spearman','gain'))
b=np.mean(res['as shipped'])
for r in RULES:
    print('%-28s %-9.4f %+.4f'%(r,np.mean(res[r]),np.mean(res[r])-b))
