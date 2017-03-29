
import run_analysis as ra
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import solver as slv
from imp import reload
#reload(slv)
#reload(ra)
xfmt = ticker.ScalarFormatter(useMathText=True, useOffset=False)
xfmt.set_powerlimits([-3, 3])

plt.style.reload_library()
plt.style.use('seaborn-whitegrid')
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['savefig.facecolor'] = 'white'
start_t = 0
t = ra.t
t1 = ra.t1
y_nocat = 1000 * ra.y2
y_nocat_chan = slv.tochannel(y_nocat)
y_withcat = 1000 * ra.y
y_withcat_chan = slv.tochannel(y_withcat)

# WARNING WE ARE CONVERTING TO nM

tot_withcat = slv.tototal(y_withcat)
tot_nocat = slv.tototal(y_nocat)

ax = []
t_lim_low = -1
t_lim_hi = 51
plt.figure(1)
plt.close()
plt.figure(1)

###############################################
# E with and without catalysis
ax.append(plt.subplot(2, 2, 3))
plt.plot(t[start_t:], y_nocat_chan[start_t:, 0, 1], label='nocat E')
plt.plot(t[start_t:], y_withcat_chan[start_t:, 0, 1], label='withcat E',
         color='#E41B17')
plt.yticks([0, 100, 200])
plt.xlabel('time (s)')
plt.ylabel('[E] (nM)')
plt.ylim([-10, 230])
plt.xlim([t_lim_low,t_lim_hi])
# plt.tick_params(axis='both', length=4, direction='in')

#############################################
# total E, to show overall effect of catalysis
ax.append(plt.subplot(2, 2, 1))
plt.plot(t[start_t:],
         y_nocat_chan[start_t:, 0, 1]+y_nocat_chan[start_t:, 2, 1],
         label='nocat totE')
plt.plot(t[start_t:],
         y_withcat_chan[start_t:, 0, 1]+y_withcat_chan[start_t:, 2, 1],
         label='withcat totE',
         color='#E41B17')
plt.yticks([200, 210, 220])
plt.ylim([199, 221])
plt.ylabel('[E] + [ES] (nM)')
plt.xlabel('time (s)')
plt.xlim([t_lim_low, t_lim_hi])

#################################################
# ES in central channel, to show the interplay with E
ax.append(plt.subplot(2, 2, 4))
plt.plot(t[start_t:], y_nocat_chan[start_t:, 2, 1], label='No Catalysis')
plt.plot(t[start_t:], y_withcat_chan[start_t:, 2, 1],
         label='With Catalysis',
         color='#E41B17')

plt.xlabel('time (s)')
plt.ylabel('[ES] (nM)')
plt.yticks([0, 100, 200])
plt.xlim([t_lim_low, t_lim_hi])
plt.ylim([-10, 220])
plt.legend(loc='lower center', bbox_to_anchor=(-.3, -.6), ncol=2)

#######################################################
# total ES+P as a proxy for number of forward binding events
ax.append(plt.subplot(2, 2, 2))
plt.plot(t, tot_withcat[:, 3] + tot_withcat[:, 2],
         label='With Catalysis',
         color='#E41B17')
plt.plot(t, tot_nocat[:, 3] + tot_nocat[:, 2], label='No Catalysis')

# plt.yscale('log')
plt.minorticks_off()
plt.xlim([t_lim_low,t_lim_hi])

plt.xlabel('time (s)')
plt.ylabel('[ES] + [P] (nM)')
plt.ylim([-.1e7, 1.8e7])
#plt.yticks([0, 1e9, 2e9, 3e9])
# plt.yticks([1e4, 1e6, 1e8, 1e10])
plt.tight_layout(pad=3)
plt.subplots_adjust(bottom=.2)  # or whatever

plt.show()

plt.figure(2)
plt.close()
fig = plt.figure(2, figsize=(2.3, 1.8))
Ntzoom = t1.size
plt.style.use('seaborn-ticks')
t_multiplier = 1000
ax.append(fig.add_subplot(111))
plt.plot(t1*t_multiplier,
         tot_withcat[:Ntzoom, 3] + tot_withcat[:Ntzoom, 2],
         label='With Catalysis',
         color='#E41B17')
plt.plot(t1*t_multiplier, tot_nocat[:Ntzoom, 3] + tot_nocat[:Ntzoom, 2],
         label='No Catalysis')
plt.yticks([0, 100])
plt.xticks([0, 0.05,0.1])
plt.xlabel('time (ms)')
for axes in ax:
    axes.yaxis.set_major_formatter(xfmt)
    for axis in ['top', 'bottom', 'left', 'right']:
        axes.spines[axis].set_color('black')

ax[4].set_axis_bgcolor('white')
ax[4].tick_params(direction = 'in')
plt.tight_layout()
plt.show()

plt.savefig('../img/figS5bis.png', dpi = 300, frameon  = True,
            transparent = True)
plt.figure(1)

plt.savefig('../img/figS5.png', dpi=300)
