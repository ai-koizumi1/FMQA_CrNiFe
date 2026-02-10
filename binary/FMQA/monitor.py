from statistics import mean, median,variance,stdev
from matplotlib.ticker import MultipleLocator
import numpy as np
import csv
import sys
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

N_ini=5
ys = np.asarray(np.loadtxt('y_tmp.csv',skiprows=0, delimiter=',') )
initial_data = np.asarray(np.loadtxt('../initial_data/initial_TOTEN.csv',skiprows=0, delimiter=',') )
sampling_method_ = np.asarray(np.loadtxt('samplingmethod_tmp.csv',skiprows=0, delimiter=',',dtype='int') )

sampling_method=np.array([[k,0] for k in range(len(ys)-N_ini)])
for s,flag in sampling_method_:
    sampling_method[s][1]=flag

step_random=[i for i in range(len(initial_data))]
step_FMQA=[]
step=[i for i in range(len(ys))]
random_data=list(ys[:N_ini])
FMQA_data=[]
for s,flag in sampling_method:
    if flag!=0:
        random_data.append(ys[int(s+N_ini)])
        step_random.append(int(s+N_ini))
    else:
        FMQA_data.append(ys[int(s+N_ini)])
        step_FMQA.append(s+N_ini)

step_random=np.array(step_random)
step=np.array(step)
current_max=[]
current_min=[]
for i in range(len(ys)):
    current_max.append(max(ys[:i+1]))
    current_min.append(min(ys[:i+1]))

plt.rcParams["font.size"] = 14
fig = plt.figure(figsize=(10,5))
gs = gridspec.GridSpec(1, 4, width_ratios=[3,3,1,1])
ax1 = fig.add_subplot(gs[0, 0:2])

ax1.scatter(step_FMQA,FMQA_data,s=100,c='orange',label='FMQA')
ax1.scatter(step_random,random_data,s=200,marker='*',c='blue',edgecolors='k',label='Random sampling')
ax1.plot(step,current_min,label='current min',c='k',lw=2,linestyle='dashed')
ax1.plot(step,current_max,label='current max',c='k',lw=2,linestyle=':')
ax1.axvspan(-1, N_ini-1, color='blue', alpha=0.2,label='Initial data')
ax1.set_xlim([-1,len(step)+1])
ax1.set_xlabel('Cycle')
ax1.set_ylabel('Total energy [eV]')
ax1.set_xticks(step)
ax1.xaxis.set_major_locator(MultipleLocator(10))
ax2 = fig.add_subplot(gs[0, 2])

ax2.hist(FMQA_data,bins=30,range=[min(ys),max(ys)],alpha=0.8,color='orange',orientation='horizontal',label='FMQA')
ax2.hist(random_data,bins=30,range=[min(ys),max(ys)],alpha=0.8,color='blue',orientation='horizontal',label='Random')
    
ax2.set_xlabel('Frequency')
ax2.set_yticklabels(labels=[])

ax1.legend(
    loc='lower center',
    bbox_to_anchor=(0.5, 1.02),
    ncol=2,
    fontsize=10
)

plt.tight_layout()
plt.savefig('history.png')
plt.show()

