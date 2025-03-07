from aloha.cohp_analysis import *
import sys

cohp=Cohpout('.')
cohp.pcohp() 
# cohp.pcohp(label=1,lm_orbital={'Fe':('d'),'O':('p')})
# cohp.pcohp(label=1,lm_orbital={'Fe':('d')})
#label=1, lm_orbital={'Ir':('dxy','dz2'),'O':('p')}), summed_spin_channels=False
#d orbital (dxy, dyz, dz2, dxz, dx2) & p orbital (px, py, pz)

# p1=cohp.d["1"]["lm_orbital"][1]['pcohp']
#p2=a.d["2"]["lm_orbital"][1]['pcohp']
#p3=a.d["3"]["lm_orbital"][1]['pcohp']

# sys.stdout = open('pcohp1.txt', 'w')
# print(p1)
#sys.stdout = open('pcohp2.txt', 'w')
#print(p2)
#sys.stdout = open('pcohp3.txt', 'w')
#print(p3)