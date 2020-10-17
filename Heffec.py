import numpy as np

#L- lenght of the side
#Lcol=even
#Lrow=odd

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
                for c in range (1,Lcol-j):
                    bond[i][j+c]*=-1
    return bond

def sites_diamond():
    sites=np.ones(((Lrow+1),(Lcol+1)))
    N=0
    print(sites)
    for r in range(0,int((Lrow-1)/2)):
        N+=1+2*r
        for c in range(0,int(Lcol/2)-r):
            sites[r][c]=0
            sites[Lrow-r][c]=0
            sites[r][Lcol-c]=0
            sites[Lrow-r][Lcol-c]=0
    print(sites)
    N=int(2*(N+Lcol+1))
    print(N)
    return sites,N

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


Lrow=15
Lcol=14

plattice=hot_start()
bond=bond_sign(plattice)
sites,N=sites_diamond()
ijmat,kmat=create_dictionaries(sites)
print(sites)
print(N)

hmat=np.zeros((N,N))

for k in range(0,N):
    for kprime in range (0,N):
        if k!=kprime:
            if ((ijmat[k][0]-ijmat[kprime][0])*(ijmat[k][0]-ijmat[kprime][0])+\
                (ijmat[k][1]-ijmat[kprime][1])*(ijmat[k][1]-ijmat[kprime][1]))>1:
                hmat[k][kprime]=0
            elif ijmat[k][1]==ijmat[kprime][1]: #y-pos same-->differ in x direction byt at max one
                hmat[k][kprime]=1
            else:
                hmat[k][kprime]=bond[ijmat[k][0]][min(ijmat[k][1],ijmat[kprime][1])]

print(hmat)

#compare min vs if...