# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 11:32:20 2018

@author: Shreyansh Agrawal
"""

import numpy as np
from numpy import linalg as ln

class WaveFunction:
    
    def __init__(self, n, xi, yi, zi, l):
        
        """
        initializing the object that would hold the informations
        of one particular wave function
        for a particle in a three dimension box problem.
        The solution for for box with equal edges.
        The potentials beyond the boundaries are assumed to be infinity
        """
        
        # h is the unit interval for discretizing axes
        self.h = l/(n-1)
        # n is the number of points in discretized solution
        self.n = n
        # the starting point of the box
        self.xi = xi
        self.yi = yi
        self.zi = zi
        # ending point of the box
        self.xf = xi + l
        self.yf = yi + l
        self.zf = zi + l
        # A is the Hamiltonian matrix
        
        # initializing A by Kinetic Energy function
        L = np.zeros((n,n),dtype=float)
        np.fill_diagonal(L,1)
        for i in range(n-1):
            L[i][i+1] = -1
        L = np.add(L,L.T)
        I = np.identity(n)
        A = np.kron(I,np.kron(I,L))+np.kron(I,np.kron(L,I))+np.kron(L,np.kron(I,I))
        
        self.A = A
        
        # discretizing the axes
        self.xx = np.linspace(xi, self.xf, self.n)
        self.yy = np.linspace(yi, self.yf, self.n)
        self.zz = np.linspace(zi, self.zf, self.n)
    
        
    def addPotential(self, x1, x2, y1, y2, z1, z2, f):
        
        """
        This function will add a potential energy matrix
        to the pre existing Hamiltonian matrix.
        The potential energy will be calculated from x1,y1,z1 to x2,y2,z2
        and the function giving the potential is f.
        f should be a function of 'x', 'y' and 'z' only.
        """
        
        # calculating the vaue of function f
        # on the pre-existing discretized axes
        xx = self.xx
        yy = self.yy
        zz = self.zz
        
        # calculating the value of the function
        # and adding it to the potential energy matrix
        for z in zz[zz>=z1]:
            if z > z2:
                break
            # k is the index corresponding to values of z
            k = np.where(zz==z)[0][0]*self.n*self.n
            for y in yy[yy>=y1]:
                if y > y2:
                    break
                # i is the index corresponding to values of y
                i = np.where(yy==y)[0][0]*self.n
                for x in xx[xx>=x1]:
                    if x > x2:
                        break
                    
                    # j is the index corresponding to values of x
                    j = np.where(xx==x)[0][0]
                    try:
                        self.A[i+j+k][i+j+k] = self.A[i+j+k][i+j+k]+eval(f)*self.h*self.h
                    except:
                        self.A[i+j+k][i+j+k] = 10**6

        
    
    def evaluate(self):
        
        """
        This evaluates the eigenvalues and eigenvectors of the
        hamiltonian matrix.
        'evals' stores the eigen values which also give the energies of
        different states.
        'evects' stores the eigen vectors. Each column corresponds
        to the eigenvalue in 'evals'. Thus, each column is the value of
        the wave function on the discretized x,y,z axis corresponding
        to one energy state.
        """
        self.evals, self.evects = ln.eigh(self.A)
    
    
    def getFunction(self, state):
        
        # returns the value of energy of state 'state'
        # and the wave function corresponding to energy state 'state'
        return self.evals[state], self.evect[:state]