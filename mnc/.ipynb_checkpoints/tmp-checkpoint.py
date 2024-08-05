for metal in [2_Ti  3_V  4_Cr  5_Mn  6_Fe  7_Co  8_Ni  9_Cu]:
    for spin in [1_LS 2_IS 2_HS 3_HS]:
        for dz in [1_ 2_ 3_ 4_ 5_ 6_]:
            path = os.path(f'/scratch/x2755a09/3_MNC/3d/{metal}/{spin}/{dz}')
            if path.exist():
                cd path
                qsub submit.sh