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
def sp(a,b): return np.corrcoef(rankdata(a),rankdata(b))[0,1]
R={k:[] for k in ['logistic','oracle_zero','oracle_abund','frac_mass_on_zeros','n_fp']}
for s,tp in SAMP:
    truth=read_truth(base+'/'+tp,'counts')
    ks,y=q(s,'logistic')
    x=np.array([truth.get(k,0.0) for k in ks]); z=x==0
    R['logistic'].append(sp(x,y))
    yz=y.copy(); yz[z]=0.0                      # perfect detection, estimated abundance
    R['oracle_zero'].append(sp(x,yz))
    ya=y.copy(); ya[~z]=x[~z]                   # perfect abundance, estimated zeros
    R['oracle_abund'].append(sp(x,ya))
    R['frac_mass_on_zeros'].append(y[z].sum()/y.sum())
    R['n_fp'].append((y[z]>0).sum())
print('Decomposing logistic\'s remaining error on Panel B (mean over 6 samples)')
print('  %-46s %.4f'%('logistic, as shipped',np.mean(R['logistic'])))
print('  %-46s %.4f'%('+ perfect DETECTION (zeros forced to 0)',np.mean(R['oracle_zero'])))
print('  %-46s %.4f'%('+ perfect ABUNDANCE (expressed set to truth)',np.mean(R['oracle_abund'])))
print()
print('  mass misassigned to truth-zero transcripts: %.3f%% of total'%(100*np.mean(R['frac_mass_on_zeros'])))
print('  false-positive transcripts (est>0, truth=0): %.0f mean'%np.mean(R['n_fp']))
print()
print('  headroom from fixing detection only : %+.4f'%(np.mean(R['oracle_zero'])-np.mean(R['logistic'])))
print('  headroom from fixing abundance only : %+.4f'%(np.mean(R['oracle_abund'])-np.mean(R['logistic'])))
