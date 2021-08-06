#####################################################################

def make_ham_single_imp_anderson_realspace( NL, NR, Vg, U, t, Vbias, tleads=1.0, Full=False ):
    #subroutine to generate one and two electron integrals for the symmetric single impurity anderson impurity model in real-space
    #ie have hopping terms between the lead states, as is done in RT-DMRG (see Heidrich-Meisner PRB 2009)
    #Returns an (NL+NR+1)x(NL+NR+1) array of the 1 e- terms 
    #Returns if Full=True an (NL+NR+1)^4 array of the 2e- terms or if Full=False simply returns U
    #The sites are ordered as the dot (impurity) site being first, followed in ascending order of the distance from the dot (impurity) with the left lead states coming before the right lead states
    #NL - number of sites in the left lead
    #NR - number of sites in the right lead
    #Vg - gate voltage ie the energy of the dot (impurity) site
    #U  - Hubbard repulsion on the dot (impurity)
    #t - hopping term between the dot (impurity) and the first site on the right/left lead
    #Vbias - constant energy term applied asymetrically to the leads to mimic a bias, Vbias/2 added to the right lead and subtracted from the left lead
    #tleads  - hopping term between lead sites, taken as the energy scale
    #Full - logical stating whether to print out the full 2e- integrals


    #Input error check
    if( np.absolute(NL-NR) != 1 and np.absolute(NL-NR) != 0 ):
        print('ERROR: Difference between number of sites in the left and right leads must be 0 or 1')
        exit()

    #Scale all energy values by tleads
    Vg = Vg * tleads
    U  = U * tleads
    t  = t * tleads
    Vbias = Vbias * tleads


    #Initialize
    N = NL+NR+1
    hmat = np.zeros([N,N])


    #Coupling part of the 1e- terms
    #dot-leads
    hmat[0,1] = -t
    hmat[0,2] = -t
    #left lead
    for lead in range(NL-1):
        indx1 = 2*lead+1
        indx2 = 2*(lead+1)+1
        hmat[indx1,indx2] = -tleads
    #right lead
    for lead in range(NR-1):
        if( NR > NL and lead == NR-2 ):
            indx1 = N-2 
            indx2 = N-1 
        else:
            indx1 = 2*(lead+1)
            indx2 = 2*(lead+2)
        hmat[indx1,indx2] = -tleads
    #complex conjugate
    hmat = hmat + hmat.conjugate().transpose()

    #Diagonal part of the 1e- terms
    #dot (impurity)
    hmat[0,0] = Vg
    #left lead
    for lead in range(NL):
        indx = 2*lead+1
        hmat[indx,indx] = -Vbias/2.0
    #right lead
    for lead in range(NR):
        if( NR > NL and lead == NR-1 ):
            indx = N-1
        else:
            indx = 2*(lead+1)
        hmat[indx,indx] = Vbias/2.0


    #Form the trivial two electron terms
    if( Full ):
        Vmat = np.zeros( [ N, N, N, N ] )
        Vmat[0,0,0,0] = U
    else:
        Vmat = U

    return hmat, Vmat

#####################################################################
