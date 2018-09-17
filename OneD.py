# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 15:13:49 2018

@author: Shreyansh Agrawal
"""

import numpy as np
from numpy import linalg as ln
from matplotlib import pyplot as plt

class WaveFunction:
    
    def __init__(self, n, xi, xf):
        
        """
        initializing the object that would hold the informations
        of one particular wave function
        for a particle in a one dimension box problem
        the potentials beyond the boundaries are assumed to be infinity
        """
        
        # n is the number of points in discretized solution
        self.n = n
        # xi is the starting point of the box
        self.xi = xi
        # xf ending point of the box
        self.xf = xf 
        # A is the Hamiltonian matrix
        self.A = np.zeros((n,n)) 
        
        # initializing A by Kinetic Energy function
        A = self.A
        A = np.zeros((n,n),dtype=float)
        np.fill_diagonal(A,1)
        for i in range(n-1):
            A[i][i+1] = -1
        A = np.add(A,A.T)
        
        # discretizing the x axis
        self.xx = np.linspace(xi, xf, n)
        # calculating the unit interval 'h'
        self.h = (xf-xi)/(n-1)
        
        
    def addPotential(self, x1, x2, f):
        
        """
        This function will add a potential energy matrix
        to the pre existing Hamiltonian matrix.
        The potential energy will be calculated from x1 to x2
        and the function giving the potential is f.
        f should be a function of 'x' only.
        """
        
        # V is the potential energy matrix that will be added
        V = np.zeros((self.n,self.n))
        # calculating the vaue of function f
        # on the pre-existing discretized x axis
        xx = self.xx
        
        # calculating the value of the function from x1 to x2
        # and adding it to the potential energy matrix
        for x in xx[xx>=x1]:
            if x > x2:
                break
            # j is the index corresponding to values of x
            j = np.where(xx==x)[0][0]
            V[j][j] = eval(f)*self.h*self.h
        
        # adding the potential energy matrix
        self.A = self.A + V
        
    
    def evaluate(self):
        
        """
        This evaluates the eigenvalues and eigenvectors of the
        hamiltonian matrix.
        'evals' stores the eigen values which also give the energies of
        different states.
        'evects' stores the eigen vectors. Each column corresponds
        to the eigenvalue in 'evals'. Thus, each column is the value of
        the wave function on the discretized x axis corresponding
        to one energy state.
        """
        self.evals, self.evects = ln.eigh(self.A)
    
    
    def getFunction(self, state):
        
        # returns the value of energy of state 'state'
        # and the wave function corresponding to energy state 'state'
        return self.evals[state], self.evect[:state]
    
    
    def drawFunction(self, state):
        
        # graphing the wave function
        e, y = self.getFunction(state)
        
        #inserting the end where wave function= 0
        x = np.insert(self.xx,[0],self.xx[0]-self.h)
        x = np.append(x, self.xx[-1]+self.h)
        y = np.insert(y,[0],0)
        y = np.append(y,0)
        plt.plot(x,y)
        plt.show()
    