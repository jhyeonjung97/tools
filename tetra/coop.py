from aloha.cobi_analysis import *
import sys

cobi=Coopout('.')
cobi.pcobi() 
# cobi.pcobi(label=1,lm_orbital={'Fe':('d'),'O':('p')})
# cobi.pcobi(label=1,lm_orbital={'Fe':('d')})
#label=1, lm_orbital={'Ir':('dxy','dz2'),'O':('p')}), summed_spin_channels=False
#d orbital (dxy, dyz, dz2, dxz, dx2) & p orbital (px, py, pz)

# p1=coop.d["1"]["lm_orbital"][1]['pcoop']
#p2=a.d["2"]["lm_orbital"][1]['pcoop']
#p3=a.d["3"]["lm_orbital"][1]['pcoop']

# sys.stdout = open('pcoop1.txt', 'w')
# print(p1)
#sys.stdout = open('pcoop2.txt', 'w')
#print(p2)
#sys.stdout = open('pcoop3.txt', 'w')
#print(p3)