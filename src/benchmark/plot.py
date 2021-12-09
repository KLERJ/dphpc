#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib.ticker

########### GRAPH PROPERTIES ############

TITLE = 'Runtime deriche. LARGE dataset. Apple M1 Pro, avg. of 10 runs'
Y_LABEL = 'Runtime [seconds]'
X_LABEL = '# cores'
FILENAME = 'test'
AXIS_LABEL_COLOR = (.4, .4, .4)
BACKGROUND_COLOR = (.88, .88, .88)
HORIZONTAL_LINES_COLOR = (1, 1, 1)

################# DATA ##################

procs = ['1', '2', '4', '8']
reference_impl = [0.1356376 for i in range(len(procs))]
baseline = [0.21723179999999997, 0.09163249999999999, 0.0488108, 0.033523599999999994]
rdma = [0.2195272, 0.0859053, 0.044615699999999994, 0.0324671]
data = [reference_impl, baseline, rdma]

#########################################

procs = [int(x) for x in procs]
line_styles = ['bo-', 'ro-', 'go-', 'yo-', 'mo-']

fig, ax = plt.subplots()
for i in range(len(data)):
    plt.plot(procs, data[i], line_styles[i])

# line labels
plt.text(x=5, y=0.14, s='ref impl', color='blue')
plt.text(x=5, y=0.05, s='MPI baseline', color='red')
plt.text(x=3, y=0.03, s='MPI RDMA', color='green')

# axes
# plt.ylim(0.025, 0.25)
# plt.xlim(0, 9)
# plt.yticks(np.arange(0, 2.76, .25))

# title, axes labels
plt.title(TITLE, y=1.07, loc='left', x=-0.1)
plt.ylabel(' ' + Y_LABEL, rotation=0, horizontalalignment='left', y=1.02, color=AXIS_LABEL_COLOR)
plt.xlabel(X_LABEL, color=AXIS_LABEL_COLOR)

# remove y axis ticks
plt.tick_params(axis='y', which='both', left=False, right=False)

# horizontal lines, background
plt.grid(axis='y', color=HORIZONTAL_LINES_COLOR)
ax.set_facecolor(BACKGROUND_COLOR)

# hide box, add digit separator
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.get_xaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))

# save
# plt.savefig(FILENAME + '.svg')
plt.savefig(FILENAME + '.png', dpi=300)