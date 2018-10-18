#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Oct 2018

Author: Stefan Tionanda
"""

from scipy import zeros_like 
from scipy.integrate import odeint 
import matplotlib.pyplot as plt 
import numpy as np

def deriv(y,t,k1,nx):
    y = np.insert(y, 0, 0)
    doty = zeros_like(y)
    
    for a in range(1, nx+1): 
        x1 = 0
        x2 = 0

        for i in range(1, a + 1): 
            x1 += y[i]*y[a-i] #first term of the equation
        
        for j in range(1, nx-a+1):           
            x2 += y[a]*y[j]
    
        doty[a] = k1*x1 - k2*x2
    
    return doty[1:] 

# =============================================================================
# Parameters and Calculations
# =============================================================================
k1 = 3.3*(10**-4) # Rate of first term of equation, which is kappa
#k1 = 2.8*(10**-5)
k2 = k1 # Rate of the second term of equation, which is kappa
nx = 20 # Nmax parameter which is the assumed maximum polymer achievable
frames = 300 # Determines the faster frame rate in animation
stoptime = 600 # Stoppage time in seconds
time = np.linspace(0,stoptime,frames) # Time function
yinit = np.array([100] +[0] * nx) # Initial array as 100% polymer is in monomer

y = odeint(deriv,yinit,time, args=(k1,nx))
#minlengthreq = # Minimum threshold length that we require
#result = # Results
# =============================================================================
# Figure function of polymer over time
# =============================================================================
plt.figure()
for i in range(0,nx):
    labl = i + 1
    plt.plot(time,y[:,i], label = labl)#,time,y[:,1],time,y[:,2],time,y[:,18],time,y[:,19]) 
#plt.plot(t,y1,'r-',linewidth=2,label='k=0.1')
#plt.plot(t,y2,'b--',linewidth=2,label='k=0.2')
#plt.plot(t,y3,'g:',linewidth=2,label='k=0.5')

#plt.legend(loc=1)
#plt.xlim(0,1300)
plt.xlabel('Time (seconds)')
plt.ylabel('Percent')
#plt.title('Npu DnaE(N) + Ssp DnaE(C)')
plt.title('Ssp DnaE(C)')
plt.show()

# =============================================================================
# Animation function of ratio in histograms in new window
# =============================================================================

def animate(i):
    yplot = y[i,:]
    N = len(yplot)
    xplot = range(1,N+1)
    plt.xlabel('Polymer Species')
    plt.ylabel('Percent')
#    plt.title('Npu DnaE(N) + Ssp DnaE(C)')
    plt.title('Ssp DnaE')
    return plt.bar(xplot, y[i,:], color="blue")

import matplotlib.animation as animation
fig = plt.figure()
import sys


if __name__ == '__main__':
    # FuncAnimation will call the 'update' function for each frame; here
    # animating over 10 frames, with an interval of 200ms between frames.
    anim = animation.FuncAnimation(fig, animate, interval=70, blit=True)
    if len(sys.argv) > 1 and sys.argv[1] == 'save':
        anim.save('2.gif', dpi=80, writer='stefan')
    else:
        # plt.show() will just loop the animation forever.
        plt.show()