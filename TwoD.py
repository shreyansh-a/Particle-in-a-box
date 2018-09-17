# -*- coding: utf-8 -*-
"""
Created on Sat Jun  9 15:26:33 2018

@author: Shreyansh Agrawal
"""

import numpy as np
from numpy import linalg as ln
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm


class WaveFunction:
    
    def __init__(self, h, xi, xf, yi, yf):
        
        """
        initializing the object that would hold the informations
        of one particular wave function
        for a particle in a two dimension box problem.
        The potentials beyond the boundaries are assumed to be infinity
        """
        
        # h is the unit interval for discretizing axes
        self.h = h
        # n is the number of points in discretized solution
        self.nx = (int)(((xf-xi)/h)+1)
        self.ny = (int)(((yf-yi)/h)+1)
        # xi and yi is the starting point of the box
        self.xi = xi
        self.yi = yi
        # xf ending point of the box
        self.xf = xf
        self.xf = xf
        # A is the Hamiltonian matrix
        
        # initializing A by Kinetic Energy function
        Ax = np.zeros((self.nx,self.nx))
        np.fill_diagonal(Ax,1)
        for i in range(self.nx-1):
            Ax[i][i+1] = -1
        Ax = np.add(Ax,Ax.T)
        
        Bx = np.identity(self.nx, float)
        
        Ay = np.zeros((self.ny,self.ny))
        np.fill_diagonal(Ay,1)
        for i in range(self.ny-1):
            Ay[i][i+1] = -1
        Ay = np.add(Ay,Ay.T)
        
        By = np.identity(self.ny, float)
    
        C = np.kron(Ay,Bx) + np.kron(By,Ax)
        
        self.A = C
        
        # discretizing the axes
        self.xx = np.linspace(xi, xf, self.nx)
        self.yy = np.linspace(yi, yf, self.ny)
    
        
    def addPotential(self, x1, x2, y1, y2, f):
        
        """
        This function will add a potential energy matrix
        to the pre existing Hamiltonian matrix.
        The potential energy will be calculated from x1,y1 to x2,y2
        and the function giving the potential is f.
        f should be a function of 'x' and 'y' only.
        """
        
        # V is the potential energy matrix that will be added
        V = np.zeros((self.nx*self.ny,self.nx*self.ny))
        # calculating the vaue of function f
        # on the pre-existing discretized x and y axis
        xx = self.xx
        yy = self.yy
        
        # calculating the value of the function from x1,y1 to x2,y2
        # and adding it to the potential energy matrix
        for y in yy[yy>=y1]:
            if y > y2:
                break
            # i is the index corresponding to values of y
            i = np.where(yy==y)[0][0] * self.nx
            for x in xx[xx>=x1]:
                if x > x2:
                    break
                
                # j is the index corresponding to values of x
                j = np.where(xx==x)[0][0]
                try:
                    V[i+j][i+j] = eval(f)*self.h*self.h
                except ZeroDivisionError:
                    V[i+j][i+j] = 10**6
        
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
        the wave function on the discretized x,y axis corresponding
        to one energy state.
        """
        self.evals, self.evects = ln.eigh(self.A)
    
    
    def getFunction(self, state):
        
        # returns the value of energy of state 'state'
        # and the wave function corresponding to energy state 'state'
        return self.evals[state], self.evect[:state]
    
    
    def drawFunction(self, state):
        
        # graphing the wave function
        e, z = self.getFunction(state)
        z = z.reshape((self.nx,self.ny))
        
        # inserting the end points where wave function= 0
        x = np.insert(self.xx,[0],self.xx[0]-self.h)
        x = np.append(x, self.xx[-1]+self.h)
        y = np.insert(self.yy,[0],self.yy[0]-self.h)
        y = np.append(y, self.yy[-1]+self.h)
        XX, YY = np.meshgrid(x,y)
        
        # padding the value of wave function at the end points
        u = np.zeros(self.ny)
        z = np.vstack((u,np.vstack((z,u))))
        u = np.zeros((self.nx+2,1))
        z = np.hstack((u,np.hstack((z,u))))
        
        fig = plt.figure(2)
        axi = fig.add_subplot(111,projection='3d')
        axi.plot_surface(XX,YY,z,cmap=cm.coolwarm)
        plt.show()