#implementation of defects
#Lrow=15
#Lcol=16

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors as c

def hot_start(): 
#==random start
    lattice = np.random.random_integers(0,1,(Lrow,Lcol)) #choose from 0 or 1
    lattice[lattice==0] =- 1 #where lattice equal zero, set it to -1
    return lattice

def bond_sign(lattice):
    bond=np.ones((Lrow,Lcol+1))
    for i in range (0,Lrow):
        for j in range(0,Lcol):
            if lattice[i][j]==-1:
                for c in range (1,Lcol+1-j):
                    bond[i][j+c]*=-1
    return bond

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

# def defects():
#     Ndef=0
#     sites[1][2]=0
#     sites[2][2]=0
  
#     Ndef=2
#     return sites, Ndef

def create_dictionaries(some_lattice):
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

def hamiltonian(N,ijmat,bond):
    hmat=np.zeros((N-Ndef,N-Ndef))
    for k in range(0,N-Ndef):
        print("kacko",ijmat[k][0],ijmat[k][1])
        for kprime in range (0,N-Ndef):
            print("ka prime",ijmat[kprime][0],ijmat[kprime][1])
            if k!=kprime:
                if ((ijmat[k][0]-ijmat[kprime][0])*(ijmat[k][0]-ijmat[kprime][0])+\
                    (ijmat[k][1]-ijmat[kprime][1])*(ijmat[k][1]-ijmat[kprime][1]))>1:
                    hmat[k][kprime]=0
                elif ijmat[k][0]==ijmat[kprime][0]: #y-pos same-->differ in x direction byt at max one
                    hmat[k][kprime]=1
                    print("tento poklada za horiz")
                else:
                    if ijmat[k][0]>ijmat[kprime][0]:
                        hmat[k][kprime]=bond[ijmat[kprime][0]][ijmat[k][1]]
                    else:
                        hmat[k][kprime]=bond[ijmat[k][0]][ijmat[k][1]]
    return hmat


Lrow=3
Lcol=2
sites,N=sites_diamond() #exsisting sites=1, or 0 if they don't exist.
# sites,Ndef=defects()
Ndef=0

ijmat,kmat=create_dictionaries(sites) #position of sites which exists



plattice=hot_start()
bond=bond_sign(plattice)
hmat=hamiltonian(N,ijmat,bond)
print(sites)
print(plattice)
print(bond)
print(hmat)


