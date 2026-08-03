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
def sp(a,b): return np.corrcoef(rankdata(a),rankdata(b))[0,1]
def mard(x,y):
    d=np.abs(x)+np.abs(y); return float(np.where(d>0,np.abs(x-y)/d,0.0).mean())
NB=20
print('Oracle length-bias correction, EXPRESSED transcripts only (17%% of universe).')
print('Upper bound on what a perfect length-aware term could buy for abundance accuracy.')
print('%-22s %-9s %-9s %-9s %-9s'%('sample','spearman','->oracle','mard','->oracle'))
A=[]
for s,tp in SAMP:
    truth=read_truth(base+'/'+tp,'counts')
    ks,ln,y=q(s,'logistic')
    x=np.array([truth.get(k,0.0) for k in ks]); e=x>0
    xe,ye,le=x[e],y[e],ln[e]
    edges=np.percentile(le,np.linspace(0,100,NB+1)); edges[0]-=1; edges[-1]+=1
    b=np.clip(np.digitize(le,edges)-1,0,NB-1)
    yc=ye.copy()
    for i in range(NB):
        m=b==i
        if m.sum()<20: continue
        r=np.median(np.log2((ye[m]+1)/(xe[m]+1)))
        yc[m]=np.maximum(0.0,(ye[m]+1)/(2.0**r)-1.0)
    row=(sp(xe,ye),sp(xe,yc),mard(xe,ye),mard(xe,yc)); A.append(row)
    print('%-22s %-9.4f %-9.4f %-9.4f %-9.4f'%(s,*row))
A=np.array(A)
print()
print('MEAN  spearman %.4f -> %.4f (%+.4f)   mard %.4f -> %.4f (%+.4f)'%(
    A[:,0].mean(),A[:,1].mean(),A[:,1].mean()-A[:,0].mean(),
    A[:,2].mean(),A[:,3].mean(),A[:,3].mean()-A[:,2].mean()))
