# Canonical CRISPRme nuclease/PAM set (reconstructed). Format per file:
#   3' PAM: <guideN><PAMmotif> <pam_len>        (positive = PAM at 3' end)
#   5' PAM: <PAMmotif><guideN> -<pam_len>       (negative = PAM at 5' start)
# filename: <guidelen>bp-<motif>-<enzyme>.txt   (enzyme parsed as split('-')[2:])
PAMS = [
    # (guidelen, motif, enzyme, pam_len_signed)   3' unless pam_len<0
    (20, "NGG",    "SpCas9",       3),
    (20, "NG",     "SpCas9NG",     2),    # SpCas9-NG / xCas9
    (20, "NGA",    "SpCas9VQR",    3),
    (20, "NGCG",   "SpCas9VRER",   4),
    (20, "NAG",    "SpCas9",       3),
    (21, "NNGRRT", "SaCas9",       6),
    (21, "NNNRRT", "SaCas9KKH",    6),
    (23, "TTTV",   "Cas12a",      -4),    # 5' PAM (Cpf1)
    (23, "TTTN",   "Cas12a",      -4),
]
IUPAC=set("ACGTRYSWKMBDHVN")
import os,sys
outdir=sys.argv[1]; os.makedirs(outdir,exist_ok=True)
def build(gl,motif,pl):
    Ns="N"*gl
    if pl>0: return Ns+motif, pl          # 3': guide then PAM
    return motif+Ns, pl                    # 5': PAM then guide
made=[]
for gl,motif,enz,pl in PAMS:
    seq,plen=build(gl,motif,pl)
    assert set(motif)<=IUPAC, motif
    # sanity: total length == guide + |pam|, motif recoverable
    n=abs(plen)
    rec = seq[-n:] if plen>0 else seq[:n]
    assert rec==motif, (rec,motif)
    fn=f"{gl}bp-{motif}-{enz}.txt"
    open(os.path.join(outdir,fn),"w").write(f"{seq} {plen}\n")
    made.append((fn, f"{seq} {plen}"))
for fn,c in made: print(f"{fn:32} {c}")
