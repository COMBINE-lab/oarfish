import gzip, struct, pickle, sys, bisect
from collections import Counter
BAM=sys.argv[1]; CAP=int(sys.argv[2]) if len(sys.argv)>2 else 8_000_000
junc=pickle.load(open('junc/junctions.pkl','rb'))
def strip(n):
    i=n.rfind('.')
    return n[:i] if i>0 and n[i+1:].isdigit() else n
# unversioned -> (len, junctions)
J={strip(k):v for k,v in junc.items()}
f=gzip.open(BAM,'rb'); assert f.read(4)==b'BAM\x01'
lt,=struct.unpack('<i',f.read(4)); f.read(lt)
nref,=struct.unpack('<i',f.read(4))
refs=[]
for _ in range(nref):
    ln,=struct.unpack('<i',f.read(4)); nm=f.read(ln)[:-1].decode(); f.read(4)
    refs.append(strip(nm))
CONSUME={0,2,3,7,8}
W=[0,2,5,10,20]
hit={'true':Counter(),'decoy':Counter()}
tot={'true':0,'decoy':0}
hit5={'true':Counter(),'decoy':Counter()}
n=0; skipped=0
while n<CAP:
    b=f.read(4)
    if len(b)<4: break
    bs,=struct.unpack('<i',b); rec=f.read(bs)
    if len(rec)<bs: break
    refid,pos=struct.unpack('<ii',rec[0:8])
    lrn=rec[8]; ncig=struct.unpack('<H',rec[12:14])[0]
    flag=struct.unpack('<H',rec[14:16])[0]
    n+=1
    if refid<0 or flag&0x4: continue
    name=rec[32:32+lrn-1].decode('ascii','replace')
    t=refs[refid]
    if t not in J: skipped+=1; continue
    L,js=J[t]
    if not js: continue
    o=32+lrn
    span=0
    for i in range(ncig):
        v,=struct.unpack('<I',rec[o+4*i:o+4*i+4])
        if (v&0xf) in CONSUME: span+=v>>4
    a=pos; e=pos+span
    # true source from NanoSim read name: <acc>_<pos>_aligned_...
    for sep in ('_aligned_','_unaligned_','_perfect_'):
        if sep in name:
            src=name.split(sep)[0].rsplit('_',1)[0]; break
    else:
        continue
    k='true' if src==t else 'decoy'
    tot[k]+=1
    i=bisect.bisect_left(js,e); best=1<<30
    for c in (i-1,i):
        if 0<=c<len(js): best=min(best,abs(js[c]-e))
    i2=bisect.bisect_left(js,a); best5=1<<30
    for c in (i2-1,i2):
        if 0<=c<len(js): best5=min(best5,abs(js[c]-a))
    for w in W:
        if best<=w: hit[k][w]+=1
        if best5<=w: hit5[k][w]+=1
print('records scanned: %d   alignments used: true=%d decoy=%d   (skipped no-struct: %d)'%(n,tot['true'],tot['decoy'],skipped))
print()
print("3' alignment END within w of an INTERNAL junction of t")
print('%-6s %-12s %-12s %-10s'%('w(nt)','P(hit|true)','P(hit|decoy)','LR decoy/true'))
for w in W:
    pt=hit['true'][w]/max(tot['true'],1); pd=hit['decoy'][w]/max(tot['decoy'],1)
    print('%-6d %-12.4f %-12.4f %-10s'%(w,pt,pd,'%.2fx'%(pd/pt) if pt>0 else 'inf'))
print()
print("5' alignment START within w of an INTERNAL junction of t")
print('%-6s %-12s %-12s %-10s'%('w(nt)','P(hit|true)','P(hit|decoy)','LR decoy/true'))
for w in W:
    pt=hit5['true'][w]/max(tot['true'],1); pd=hit5['decoy'][w]/max(tot['decoy'],1)
    print('%-6d %-12.4f %-12.4f %-10s'%(w,pt,pd,'%.2fx'%(pd/pt) if pt>0 else 'inf'))
