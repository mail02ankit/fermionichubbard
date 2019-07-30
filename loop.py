import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.interpolate import griddata
import matplotlib.ticker as ticker

#Shell environment
import os
import sys
import subprocess       # to send python variables to shell.

# Permutations
import itertools


## Functions
#### Number of combinations
def choose(n, k):
#     https://stackoverflow.com/questions/3025162/statistics-combinations-in-python/3027128
    if 0 <= k <= n:
        ntok = 1
        ktok = 1
        for t in range(1, min(k, n - k) + 1):
            ntok *= n
            ktok *= t
            n -= 1
        return ntok // ktok
    else:
        return 0

#### Get the sign for the movement of a single fermionic operator  Uk...UjUi...U1U0|0> .
def getFermionSign(_iOpetrator, _iOfBasisOnlyOnes):
    _list_below = np.where( np.array(_iOfBasisOnlyOnes)> _iOpetrator )[0]
    if len(_list_below) >0:
        _num_loop = _list_below[-1] + 1 # This condition will pick the element where we need to move the operator at.
    else:
        _num_loop = 0
    _fermionSign = (-1)**_num_loop
#     print("n", _iTemp, _num_loop, _fermionSign)
    return _fermionSign,  np.insert(_iOfBasisOnlyOnes,_num_loop,_iOpetrator)


#### Basis vectors generation and indexing
def getBasis( _ns, _neU, _neD):
    '''
    * for NS sites, [0,NS-1] exponents are possible for base 2.
    * Out of those exponents, we pick NEU/D exponents to create a basis.
    * Binary representation will give us a readable basis.
    '''
    _site_exponents = np.array( [_ins for _ins in range(_ns)])
    _site_fillingU = list(itertools.combinations(_site_exponents,_neU))
    _site_fillingD = list(itertools.combinations(_site_exponents,_neD))
    
    _basisH=[]
    for _ifillU in range( len(_site_fillingU) ):
        for _ifillD in range( len(_site_fillingD) ):
            # Spin up eletron configuration decimal value.
            _id_basisU = 0
            for _ineU in range(_neU):
                _id_basisU += 2**_site_fillingU[_ifillU][_ineU]
            # Spin down eletron configuration decimal value. 
            _id_basisD = 0
            for _ineD in range(_neD):
                _id_basisD += 2**_site_fillingD[_ifillD][_ineD]    
            
            # Fill the basis
            _basisH.append([_id_basisU, _id_basisD])
    _basisH = np.array( _basisH )
    return _basisH

#### Number conversion
# Get binary
def getBinary(_n):
    return [int(x) for x in bin(_n)[2:]]

def getDecimal(_n): 
    _out = 0
    for _in, _value in enumerate(_n):
        _out += _value*2**(len(_n)-1-_in)
    return _out

#########

######################## Single one

from scipy.sparse import identity
from scipy.sparse.linalg import eigs
from scipy.sparse import lil_matrix

#### Calculate filling
def calc_Hmu(_iBasis, _basis):
#     print(_basis)
    _sumTemp = 0.0
    for _iS in range(NS):
        if _basis[0][_iS] == 1:
            _sumTemp += mu
        if _basis[1][_iS] == 1 :
            _sumTemp += mu
    if np.fabs(_sumTemp) >0.0:
        if np.fabs( H[_iBasis,_iBasis]) < 0.0000000001:
            H[_iBasis,_iBasis] = -_sumTemp
        else:
            H[_iBasis,_iBasis] += -_sumTemp


#### H_intraction act on a basis
def calc_HU(_iBasis, _basis):
#     print(_basis)
    _sumTemp = 0.0
    for _iS in range(NS):
        if _basis[0][_iS] == 1 and _basis[1][_iS] ==1 :
            _sumTemp += u
    if _sumTemp >0.0:
        if np.fabs( H[_iBasis,_iBasis]) < 0.0000000001:
            H[_iBasis,_iBasis] = _sumTemp
        else:
            H[_iBasis,_iBasis] += _sumTemp
#             print(_iS)
#     print("Done with the interaction part")

def calc_HK(_iBasis, _basis):
    _posOnlyOnesU = []
    _posOnlyOnesD = []
    ### Find the operator sequence that create the input basis. 
    for _iS in range(NS): 
        _iSS = NS-1 -_iS # Decrement
        if _basis[0][_iSS] == 1:
            _posOnlyOnesU.insert(0,_iS)
        if _basis[1][_iSS] == 1:
            _posOnlyOnesD.insert(0,_iS)
    _posOnlyOnes = [_posOnlyOnesU, _posOnlyOnesD]
      
    #### Act HK and find the fermionic sign  
    
    ### Forward shift: c*(i+1)c(i)
    _mapped_basis = np.copy(_basis)
    for _iS in range(NS): 
        _iSS = NS-1 -_iS # Decrement
        ## Non-zero kinetic term index
        for _iSpin in range(2): # Loop over spin up and down basis part
            if _basis[_iSpin][_iSS] == 1 and _basis[_iSpin][ (_iSS-1)%NS ] ==0 : 
                # Find the mapped basis index
                _mapped_basis[_iSpin][(_iSS-1)%NS] = 1
                _mapped_basis[_iSpin][_iSS] = 0
                _target_basis = [getDecimal(_mapped_basis[0]), getDecimal(_mapped_basis[1])]
                # Search the index of the target basis
                _iTargetUpSpin = np.where( _target_basis[0] == basisH[:,0])
                _iTarget_basis = _iTargetUpSpin[0][ np.where( basisH[_iTargetUpSpin][:,1] == _target_basis[1] )[0] ] 
                # Get the sign
                _fermionSign1, _newSequence1 = getFermionSign(_iS, _posOnlyOnes[_iSpin] )
                _fermionSign2, _newSequence2 = getFermionSign((_iS+1)%NS, _newSequence1 )
                _fermionSign = _fermionSign1*_fermionSign2
                
                H[_iBasis,_iTarget_basis] = - _fermionSign * t
                _mapped_basis = np.copy(_basis)
         
        ### Backward shift: c*(i)c(i+1)
        for _iSpin in range(2): # Loop over spin up and down basis part
            if _basis[_iSpin][_iSS] == 0 and _basis[_iSpin][ (_iSS-1)%NS ] ==1 : 
                # Find the mapped basis index
                _mapped_basis[_iSpin][(_iSS-1)%NS] = 0
                _mapped_basis[_iSpin][_iSS] = 1
                _target_basis = [getDecimal(_mapped_basis[0]), getDecimal(_mapped_basis[1])]
                # Search the index of the target basis
                _iTargetUpSpin = np.where( _target_basis[0] == basisH[:,0])
                _iTarget_basis = _iTargetUpSpin[0][ np.where( basisH[_iTargetUpSpin][:,1] == _target_basis[1] )[0] ]
                # Get the sign
                _fermionSign1, _newSequence1 = getFermionSign((_iS+1)%NS, _posOnlyOnes[_iSpin] )
                _fermionSign2, _newSequence2 = getFermionSign(_iS, _newSequence1 )
                _fermionSign = _fermionSign1*_fermionSign2
                
                H[_iBasis,_iTarget_basis] = - _fermionSign * t
                _mapped_basis = np.copy(_basis)

# Loop
######################## Single one

#### Parameters
NS = 6  # Number of sites
NSpinU = 3 # Number of spin up
NSpinD = 3 # Number of spin down
u = 0.0
t = 0.25
mu = 0.0
# mu = u/2
####

#### H definition
dimH = choose(NS,NSpinU) * choose(NS,NSpinD)
H = lil_matrix((dimH, dimH))
### convert to csc
# print("Size in bytes: H_lil:",sys.getsizeof(H), "csr:", sys.getsizeof(H.tocsr()), " Full:", sys.getsizeof(H.toarray()) )
# H = H.tocsr()  # I do not think we need to convert to csr to get eigen values.
print("[NS, NU, ND, dim]", [NS, NSpinU, NSpinD, dimH])
# print("[u, t, mu]", [u, t, mu])
####

### Create basis
basisH = getBasis(NS,NSpinU, NSpinD)
####

gs = []
# for t in [ it*0.1 for it in range(20)]:
for u in [ iu*0.5 for iu in range(30)]:
#     print("\n [NS, NU, ND, dim]", [NS, NSpinU, NSpinD, dimH])
    print("[u, t, mu]", [u, t, mu])

    #### Fill H elements
    for iBasis in range(len(basisH)):
        ### Extarct spin resolved basis and pad zeros at the begining.
        basisU = getBinary(basisH[iBasis][0])
        basisU = np.pad(basisU, ( NS - len(basisU),0), 'constant', constant_values=(0, 0))

        basisD = getBinary(basisH[iBasis][1])
        basisD = np.pad(basisD, ( NS - len(basisD),0), 'constant', constant_values=(0, 0))
        ### Combine
        basisComb = np.array( [basisU, basisD])
        ### Calculate H
        calc_HU( iBasis, basisComb )
        calc_HK( iBasis, basisComb )
        calc_Hmu(iBasis, basisComb )
    # print(H)
#     print(H.toarray())

    #### Eigen values using sparse eigs routine
    eig_val, eig_vec = eigs(H, 2, which='SR')
#     print("Eig values Sparse:", np.sort(real(eig_val)))
#     eig_val = eigvals(H.toarray())
#     print("Eig values Full:", np.sort(eig_val))
    
    gs.append( [u, t , mu, np.min(np.real( eig_val))])

gs = np.array(gs)
np.savetxt( ''.join( ('./out.dat') ), gs  ,delimiter = " ", fmt="%0.12f")
print("\n", gs)

# Plots
### Plot
fig, ax = plt.subplots(1,1, figsize=(3,3), dpi=300)
marksize = 3
fontl = 8

# GS energy
ax.plot(gs[:,0],gs[:,3], '-o', markerfacecolor='r', markersize = marksize, label='P(x)' )

# ax.set_ylim(0,10)
ax.set_xlabel('U', fontsize = 1.2*fontl)
ax.set_ylabel(r'$\epsilon_0$', fontsize = 1.2*fontl)
# ax.legend( loc='upper right', frameon = 'True')
ax.tick_params(axis='both', which='major', labelsize=1.0*fontl) 
#ax.tick_params(axis='y', left='off', top='off', right='off', bottom='off', labelleft='off', labeltop='off', labelright='off', labelbottom='off', labelsize=0.4*fontl) 

fig.savefig('test.eps', bbox_inches='tight')
fig.savefig('test.pdf', bbox_inches='tight')
fig.savefig('test.jpg', bbox_inches='tight')

