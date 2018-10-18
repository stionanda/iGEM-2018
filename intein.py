import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
#import matplotlib.animation as animation

def model(X,t,a,b,c,d):
    return [(a*x-b*X[0]), (c*X[0]-d*X[1])]


# all time variable is in seconds
a = 0.08 # mRNA formation (related to transcription rate)
b = 0.00231+ 0.45/60 # mRNA degradation (related to half-life and dilution factor)
c = 0.084 # protein translation (related to fluorescence)
d = 0.45/60 # protein degradation (related to dilution factor)

y0 = [0,0]
x = 100 # copy number
time = 2000
t = np.linspace(0,time,time)

m1 = odeint(model,y0,t,args=(a,b,c,d))
mrna = m1[:,0]
protein = m1[:,1]

plt.plot(t, mrna, "-", label="mRNA")
plt.plot(t, protein, "-", label="Protein")

#fig = plt.figure()
#plt.xlim(0,100)
#plt.ylim(0,14000)
plt.xlabel('time (minutes)')
plt.ylabel('molecules')

#def animate(i):
#    if i == 1:
#        plt.plot(m1[0:i, 0], 'b-', label='mRNA')
#        plt.plot(m1[0:i, 1], 'r-', label='protein')
#        plt.legend(loc=2)
#    else:
#        plt.plot(m1[0:i, 0], 'b-')
#        plt.plot(m1[0:i, 1], 'r-')
#
#ani = animation.FuncAnimation(fig, animate, interval=1)
