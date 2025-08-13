#!/bin/env python3

data = [
    # [catalyst, e_, e_oh, e_o, e_ooh, e_o2]
    # ['Mn-N-C', -281.65, -292.61, -287.44, -296.92, np.nan],
    ['Fe-N-C', -280.18, -290.95, -285.94, -295.19, -290.62],
    ['Co-N-C', -279.55, -289.62, -284.79, -294.19, -289.90],
    # ['Ni-N-C', -277.21, -287.28, -281.31, -291.72, np.nan],
]

# gas
h2 = -6.77149190
h2o = -14.23091949

zpeh2o = 0.558
cvh2o = 0.103
tsh2o = 0.675

zpeh2 = 0.268
cvh2 = 0.0905
tsh2 = 0.408

gh2o = h2o + zpeh2o - tsh2o + cvh2o
gh2 = h2 + zpeh2 - tsh2 + cvh2
go2 = 2*gh2o - 2*gh2 + 4.92

# ads
zpeoh = 0.376
cvoh = 0.042
tsoh = 0.066

zpeo = 0.064
cvo = 0.034
tso = 0.060

zpeooh = 0.471
cvooh = 0.077
tsooh = 0.134

dgo = zpeo + cvo - tso
dgoh = zpeoh + cvoh - tsoh
dgooh = zpeooh + cvooh - tsooh
dgh = dgoh - dgo
dgo2 = dgooh - dgh

steps = ['O2->OOH*', 'OOH*->O*', 'O*->OH*', 'OH*->H2O']

print('catal\tgoh\tgo\tgooh\tdg0\tdg1\tdg2\tdg3\tdg4\top14\top04')
for row in data:
    catalyst, e_, e_oh, e_o, e_ooh, e_o2 = row
    g_ = e_
    g_oh = e_oh + dgoh
    g_o = e_o + dgo
    g_ooh = e_ooh + dgooh 
    g_o2 = e_o2 + dgo2

    g1 = g_oh + gh2/2 - g_ - gh2o
    g2 = g_o + gh2 - g_ - gh2o
    g3 = g_ooh + gh2*3/2 - g_ - gh2o*2

    dg0 = g_o2 - g_ - go2
    dg1 = g_ooh - g_ - go2 - gh2/2 + 1.23
    dg2 = g_o + gh2o - g_ooh - gh2/2 + 1.23
    dg3 = g_oh - g_o - gh2/2 + 1.23
    dg4 = g_ + gh2o - g_oh - gh2/2 + 1.23

    dg14 = [dg1, dg2, dg3, dg4]
    # dg04 = [dg0, dg1-dg0, dg2, dg3, dg4]
    # dg04 = [dg0, dg1-dg0, dg2, dg3, dg4]
    overpotential14 = max(dg14)
    overpotential04 = max(dg04)
    
    print(f'{catalyst}\t{round(g1, 2)}\t{round(g2, 2)}\t{round(g3, 2)}\t{round(dg0, 2)}\t{round(dg1, 2)}\t{round(dg2, 2)}\t{round(dg3, 2)}\t{round(dg4, 2)}\t{round(overpotential14, 2)}\t{round(overpotential04, 2)}')