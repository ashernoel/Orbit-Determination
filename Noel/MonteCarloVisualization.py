#making sexy histograms


import numpy as np
import matplotlib.pyplot as plt

def gaussian(x, a, mean, std, amp):

    return amp/(np.sqrt(np.pi * 2) * std)*np.exp(-(x-mean)**2/(2*std**2))


a = np.load("MCaList.npy")

fig,ax = plt.subplots(1,1)
ax.hist(a, alpha = 1, bins = 100, histtype='step')

'''
xrange = np.linspace(1.4,1.55,50) # 
ax.plot(xrange, gaussian(xrange, None, np.median(a), np.std(a), 180), color="black")
'''


plt.title("Sample Distribution of Semimajor Axis", fontname="Times New Roman", fontsize=14)
plt.ylabel("Count", fontname="Times New Roman", fontsize=14)
plt.xlabel("Semimajor Axis (a) [AU]", fontname="Times New Roman", fontsize=14)

plt.legend()
plt.show()
