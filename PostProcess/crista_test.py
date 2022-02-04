import pandas as pd
from CRISTA_score import CRISTA_predict

crista_1450 = pd.read_csv("./BCL11A_DE_OTs_1450_CRISTA_scores.csv")
sgrna_20nt = list(crista_1450['sgRNA_20nt'])[:20]
dna_seq_29nt = list(crista_1450['dna_seq_29nt'])[:20]
print(CRISTA_predict(sgrna_20nt, dna_seq_29nt))