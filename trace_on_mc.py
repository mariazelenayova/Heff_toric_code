#storing gs and plattice after each measurment
# Claudios imp + hole

import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib import colors as c
from scipy import linalg as la
import math
import time

def sites_diamond():
    sites=np.ones(((Lrow+1),(Lcol+1)))
    N=0
    for r in range(0,int((Lrow-1)/2)):
        N+=1+2*r
        for c in range(0,int(Lcol/2)-r):
            sites[r][c]=0
            sites[Lrow-r][c]=0
            sites[r][Lcol-c]=0
            sites[Lrow-r][Lcol-c]=0
    N=int(2*(N+Lcol+1))
    return sites,N
def defects():
    Ndef=0
    sites[8][1]=0
    sites[9][1]=0
    
    sites[6][7]=0
    sites[6][8]=0
    sites[6][9]=0
    sites[7][6]=0
    sites[7][7]=0
    sites[7][8]=0
    sites[7][9]=0
    sites[8][5]=0
    sites[8][6]=0
    sites[8][8]=0
    sites[9][4]=0
    sites[9][5]=0

    sites[12][5]=0
    sites[13][5]=0

    sites[14][7]=0
    sites[15][7]=0
    Ndef=18
    return sites, Ndef
def create_dictionaries(some_lattice):
    # some_lattice = information about existing sites
    # ijmat holds infromation about position of kth existing site
    # kmat holds information about numebr of spin on the required position
    k=0
    ijmat={}
    kmat={}
    for i in range(0,Lrow+1):
        for j in range(0,Lcol+1):
            if some_lattice[i][j]==1: #if site exist
                ijmat[k]=(i,j) #call it "k" and remember it's position
                kmat[i,j]=k
                k+=1
    return ijmat, kmat

def plaquette_start(sites): 
    lattice = np.zeros((Lrow,Lcol))
    if sites[half][0]==0 or sites[half+1][0]==0:
        lattice[half][0]=0
    else:
        lattice[half][0]=1 if np.random.rand() < 0.5 else -1
    
    if sites[half][Lcol]==0 or sites[half+1][Lcol]==0:
        lattice[half][Lcol-1]=0
    else:
        lattice[half][Lcol-1]=1 if np.random.rand() < 0.5 else -1
    for r in range (1,Lrow):
        for c in range(1,Lcol):
            if sites[r][c]==0:
                lattice[r][c]=0
                lattice[r][c-1]=0
                lattice[r-1][c]=0
                lattice[r-1][c-1]=0
            else:
                lattice[r][c]=1 if np.random.rand() < 0.5 else -1
    for r in range(0,int((Lrow-1)/2)):
        for c in range(0,int(Lcol/2)-r):
            lattice[r][c]=0
            lattice[Lrow-1-r][c]=0
            lattice[r][Lcol-1-c]=0
            lattice[Lrow-1-r][Lcol-1-c]=0
    return lattice
    
def pi_chose(lattice):
    #pick the plaquette, assign pi
    lattice[6][7]=-1
    return lattice 
    
def plaquette_compact(lattice):
    #pmat will contain only existing plaquettes
    pmat={}
    pkmat={}
    k=0
    for i in range(0,Lrow):
        for j in range(0,Lcol):
            if lattice[i][j]==1 or lattice[i][j]==-1:
                pmat[k]=(1,i,j)
                pkmat[i,j]=(k)
                k+=1
    return pmat,k

def random_start(pmat):
    plattice=np.zeros((Lrow,Lcol))
    #will change only existing plaquettes
    for k in range (0,len(pmat)):
        plattice[pmat[k][1]][pmat[k][2]]=1 if np.random.rand() < 0.5 else -1
    return plattice

def no_vison(pmat):
    plattice=np.zeros((Lrow,Lcol))
    #will change only existing plaquettes
    for k in range (0,len(pmat)):
        plattice[pmat[k][1]][pmat[k][2]]=-1
    return plattice    

def bond_sign(lattice):
    #by bond you mean vertical one
    bond=np.ones((Lrow,Lcol+1))
    for i in range (0,Lrow):
        for j in range(0,Lcol):
            if lattice[i][j]==1:
                for c in range (1,Lcol+1-j):
                    bond[i][j+c]*=-1
    return bond

def spins_compact(kmat, plaq_lattice):
    #bond mat contains infabout spins creating the bond
    # spinmat contains inf about plaqu to be flipped
    spinmat={}
    bondmat={}
    l=0
    for i in range(1,Lrow-1):
        for j in range(1,Lcol-1):
            if plaq_lattice[i][j]==1 or plaq_lattice[i][j]==-1:
                spinmat[l]=(i,j,i,j-1) #left
                spinmat[l+1]=(i,j,i-1,j) #down
                bondmat[l]=(kmat[i,j],kmat[(i+1),j]) #left bond
                bondmat[l+1]=(kmat[i,j],kmat[i,(j+1)]) #down bond
                l+=2
                if plaq_lattice[i][j+1]==0: #right bound
                    spinmat[l]=(i,j,i,j+1)
                    bondmat[l]=(kmat[i,(j+1)],kmat[(i+1),(j+1)])
                    l+=1
                if plaq_lattice[i+1][j]==0:#upper bound
                    spinmat[l]=(i,j,i+1,j)
                    bondmat[l]=(kmat[(i+1),j],kmat[(i+1),(j+1)])
                    l+=1
    spinmat[l]=(half,Lcol-1,half,Lcol-2) #0,0 so no plaquet would be flipped, cause plaquete[0,0] is zero
    spinmat[l+1]=(half,Lcol-1,half-1,Lcol-1)
    spinmat[l+2]=(half,Lcol-1,0,0)
    spinmat[l+3]=(half,Lcol-1,half+1,Lcol-1)

    bondmat[l]=(kmat[half,Lcol-1],kmat[half+1,Lcol-1])
    bondmat[l+1]=(kmat[half,Lcol-1],kmat[half,Lcol])
    bondmat[l+2]=(kmat[half,Lcol],kmat[half+1,Lcol])
    bondmat[l+3]=(kmat[half+1,Lcol-1],kmat[half+1,Lcol])
    return spinmat,bondmat, l+4 #(l+4) stores number of spins active - which when flipped change the vison

def MC_step(plaq_lattice, hmat, E1,T):
    for step in range(0, numplaq): #this is number_plaq*microstep
        r1 =np.random.randint(0,numspin)
        hmat[bondmat[r1][0]][bondmat[r1][1]]*=-1
        hmat[bondmat[r1][1]][bondmat[r1][0]]*=-1
        eigval=la.eigh(hmat,eigvals_only=True,eigvals=(0,0))
        E2=eigval
        r2=np.random.rand()
        if (E1-E2)>=0:
            new_energy=E2 #accept and don't forget to change plaquettes
            plaq_lattice[spinmat[r1][0]][spinmat[r1][1]]*=-1
            plaq_lattice[spinmat[r1][2]][spinmat[r1][3]]*=-1
            if r1 in hole_bonds:
                plaq_lattice[6][7]*=-1
        elif math.exp(-(E2-E1)/T)>r2:
            new_energy=E2 #accept and don't forget to change plaquettes
            plaq_lattice[spinmat[r1][0]][spinmat[r1][1]]*=-1
            plaq_lattice[spinmat[r1][2]][spinmat[r1][3]]*=-1
            if r1 in hole_bonds:
                plaq_lattice[6][7]*=-1
        else:
            new_energy=E1 #don't accept
            hmat[bondmat[r1][0]][bondmat[r1][1]]*=-1
            hmat[bondmat[r1][1]][bondmat[r1][0]]*=-1 
            #flip back cause we don't accept
        E1=new_energy
    return plaq_lattice, new_energy, hmat

def thermalization(plaq_lattice, hmat, E1, T):
    for loop in range(0,1000):
        plaq_lattice,E1,hmat=MC_step(plaq_lattice,hmat, E1,T)
    return plaq_lattice, E1, hmat

def hamiltonian(N,ijmat,bond):
    hmat=np.zeros((N-Ndef,N-Ndef))
    for k in range(0,N-Ndef):
        for kprime in range (0,N-Ndef):            
            if k!=kprime:
                if ((ijmat[k][0]-ijmat[kprime][0])*(ijmat[k][0]-ijmat[kprime][0])+\
                    (ijmat[k][1]-ijmat[kprime][1])*(ijmat[k][1]-ijmat[kprime][1]))>1:
                    hmat[k][kprime]=0
                elif ijmat[k][0]==ijmat[kprime][0]: #y-pos same-->differ in x direction byt at max one
                    hmat[k][kprime]=1
                else:
                    if ijmat[k][0]>ijmat[kprime][0]:
                        hmat[k][kprime]=bond[ijmat[kprime][0]][ijmat[k][1]]
                    else:
                        hmat[k][kprime]=bond[ijmat[k][0]][ijmat[k][1]]
    return hmat

#INITIALIZING GLOBAL VARIABLES
T_final=0.0001
print(T_final)
T_initial=0.1
print(T_initial/T_final)
b=(T_initial/T_final)**(1/60)
print(b)
Lrow=15
Lcol=14
N=0
Ndef=0 
half=int((Lrow-1)/2)
numplaq=0
numspin=0
hole_bonds=[35,38,41,44,58,74,87,97,96,115,113,111,110,131,128,125,108,95,86,85,73,72,57,56]

sites,N=sites_diamond() #exsisting sites=1, or 0 if they don't exist.
sites,Ndef=defects() #set missing sites to zero
ijmat,kmat=create_dictionaries(sites) #position of sites which exists

pstart=plaquette_start(sites) #inicializing plaquettes
pmat,numplaq=plaquette_compact(pstart) #pmat keeps track of existing plaqs as ijmat/kmat of sites.

plattice=no_vison(pmat)
spinmat,bondmat,numspin =spins_compact(kmat, plattice) #spinmat keeps track of spins and corresponding plaquettes which they flip

#-----STARTING WITH NO VISONS
bond=bond_sign(plattice)

hmat=hamiltonian(N,ijmat,bond)
eigval,eigvec=la.eigh(hmat,eigvals=(0,0))
T=T_initial
E1=eigval#find ground-state
print(math.exp(-1/(k*T)*E1)

# plattice=pi_chose(plattice)
# #-----------------------

# plattice,E1,hmat=thermalization(plattice,hmat, E1,T)
# idx=0
# temp=np.zeros(61)
    # T=T_initial
    # while T>T_final:
    # print(T)
    # start_time = time.time()
    # energy=np.zeros(50)
    # numvis=np.zeros(50)
    # eigvector = np.array([])
    # vison=np.array([])
    # plattice,E1,hmat=thermalization(plattice,hmat, E1,T)

    # for measurment in range(0,50):
    # for loop in range(0,10):
    # plattice,E1,hmat=MC_step(plattice,hmat, E1, T)
    # eigval,eigvec=la.eigh(hmat,eigvals=(0,0))
    # #eigh gives real things
    # #eigvals=(0,1) gives teh lowest eigvalue ==min(eigval)
    # E1=eigval #=eigval[0]
    # gsvector=eigvec
    # eigvector=np.append(eigvector,gsvector)
    # vison=np.append(vison,plattice)
    # energy[measurment]=E1
    # numvis[measurment]=np.count_nonzero(plattice == 1)
    # temp[idx]=T
    # idx += 1
    # T=T/b
    # print("temperature", idx,"%s seconds ---" % (time.time() - start_time))
    # np.savetxt('temperature.txt',temp)