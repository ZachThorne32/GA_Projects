#####################################################################
import numpy as np
def make_ham_multi_imp_anderson_realspace( Nimp, NL, NR, Vg, U, timp, timplead, Vbias, tleads=1.0, Full=False ):
    #subroutine to generate one and two electron integrals for the symmetric multi-impurity anderson impurity model in real-space
    #ie have hopping terms between the lead states, as is done in RT-DMRG (see Heidrich-Meisner PRB 2009)
    #Returns an (NL+NR+Nimp)x(NL+NR+Nimp) array of the 1 e- terms 
    #Returns if Full=True an (NL+NR+Nimp)^4 array of the 2e- terms or if Full=False simply returns U
    #The sites are ordered from left to right, with impurities being in the center of the two leads
    #If the number of impurities equals number of sites you get back Hubbard model with timp=timplead=tleads
    #Nimp - number of impurity sites
    #NL - number of sites in the left lead
    #NR - number of sites in the right lead
    #Vg - gate voltage ie the energy of the impurity sites
    #U  - Hubbard repulsion on the impurities
    #timplead - hopping term between the furthest left and furtherst right impurity and the first site on the left/right lead
    #timp - hopping term between the impurities
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
    timp  = timp * tleads
    timplead  = timplead * tleads
    Vbias = Vbias * tleads

    #Initialize
    N = NL+NR+Nimp
    hmat = np.zeros([N,N])

    #Lists for indices of left-lead, impurities, and right-lead
    left_indx = np.arange(NL)
    imp_indx  = np.arange(NL,NL+Nimp)
    right_indx = np.arange(NL+Nimp,N)

    #Coupling part of the 1e- terms
    #left lead
    for lead in left_indx[:-1]:
        hmat[ lead, lead+1 ] = tleads
        hmat[ lead+1, lead ] = tleads
    #left lead - impurity
    hmat[ left_indx[-1], imp_indx[0] ] = timplead
    hmat[ imp_indx[0], left_indx[-1] ] = timplead      
    #impurities
    for imp in imp_indx[:-1]:
        hmat[ imp, imp+1 ] = timp
        hmat[ imp+1, imp ] = timp
    #impurity - right lead
    hmat[ imp_indx[-1], right_indx[0] ] = timplead
    hmat[ right_indx[0], imp_indx[-1] ] = timplead
    #right lead
    for lead in right_indx[:-1]:
        hmat[ lead, lead+1 ] = tleads
        hmat[ lead+1, lead ] = tleads

    #Diagonal part of the 1e- terms
    #impurities
    for imp in imp_indx:
        hmat[imp,imp] = Vg
    #left lead
    for lead in left_indx:
        hmat[lead,lead] = -Vbias/2.0
    #right lead
    for lead in right_indx:
        hmat[lead,lead] = Vbias/2.0

    #Form the trivial two electron terms
    if( Full ):
        Vmat = np.zeros( [ N, N, N, N ] )
        for imp in imp_indx:
            Vmat[imp,imp,imp,imp] = U
    else:
        Vmat = U

    return hmat,Vmat

#####################################################################
