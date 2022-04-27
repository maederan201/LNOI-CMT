import numpy as np
import matplotlib.pyplot as plt

result_array = np.load('result_array.npz')['a']

fig, ax = plt.subplots()
n = 1

for j in range(n):
    ax.plot(result_array[:,0],result_array[:,j+1],marker='.')
ax.set_xlabel(r'$\Delta x$')
ax.set_ylabel('Effective Index')
fig.show()