import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from collections import namedtuple


n_groups = 21

means_men = (8.5, 4.25, 6., 3.25, 4.75,1.25,4.75,1.75,4.,2.,2.25,1.,1.5,2.,1.5,1.5,1.5,0.75,0.,0.75,4.75)
std_men = (5.259911279, 1.707825128, 2.449489743, 2.217355783, 1.5,0.5,2.217355783,0.957427108,0.816496581,0.816496581,1.258305739,1.154700538,1.290994449,1.154700538,0.577350269,1.290994449,0.577350269, 0.5,0.,0.957427108,0.5)



fig, ax = plt.subplots()

index = np.arange(n_groups)
bar_width = 0.35

opacity = 0.4
error_config = {'ecolor': '0.3'}

rects1 = ax.bar(index, means_men, bar_width,
                alpha=opacity, color='b',
                yerr=std_men, error_kw=error_config,)

ax.set_xlabel('Polymer lengths')
ax.set_ylabel('# of polymers')
ax.set_title('Distribution of polymer lengths, runs=4')
ax.set_xticks(index + bar_width / 2)
ax.set_xticklabels(('1', '2', '3', '4', '5', '6','7','8','9','10','11','12','13','14','15', '16','17','18','19','20','>20'))

fig.tight_layout()
plt.show()