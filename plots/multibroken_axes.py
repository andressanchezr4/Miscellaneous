# -*- coding: utf-8 -*-
"""
Created on 2022

@author: andres.sanchez
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

santos_order = ["7TM1","Nuclear receptor","VGIC"]

xx = pd.DataFrame({"ComCount": [6, 8, 10, 12, 14, 16, 18, 20, 22],
                   "Family": santos_order * 3, 
                  'Group': ["HMDB"]*3  + ["Both HMDB & DrugBank"]*3 + ["DrugBank"]*3})

xx = pd.read_csv('./cpperfampred.csv', sep = ';')
xx['ComCount'] = np.where(pd.isna(xx['ComCount']), 0, xx['ComCount'])

f, (ax1, ax2, ax3) = plt.subplots(ncols=1, nrows=3,
                             sharex=True)
sz = 20
plt.rcParams["figure.figsize"] = (30, 12)
#plt.xticks(fontsize= sz*1.2)

ax3.spines['top'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.spines['bottom'].set_visible(False)
ax1.spines['bottom'].set_visible(False)

ax1.tick_params(axis='both', which='major', labelsize=sz*1.1)
ax2.tick_params(axis='both', which='major', labelsize=sz*1.1)
ax3.tick_params(axis='both', which='major', labelsize=sz*1.1)

# we want the "Test" to appear on the x axis as individual parameters
# "Latency in ms" should be what is shown on the y axis as a value
# hue should be the "Experiment Setup"
# this will result three ticks on the x axis with X1...X3 and each with three bars for T1...T3
# (you could turn this around if you need to, depending on what kind of data you want to show)
ax1 = sns.barplot(x='Family', y='ComCount',
                  hue='Group', data=xx, ax=ax1)

# we basically do the same thing again for the second plot
ax2 = sns.barplot(x='Family', y='ComCount',
                  hue='Group', data=xx, ax=ax2)

ax3 = sns.barplot(x='Family', y='ComCount',
                  hue='Group', data=xx, ax=ax3)

ax1.set_ylim(22000, 42000)
ax2.set_ylim(550, 900)
ax3.set_ylim(0, 300)

ax1.get_xaxis().set_visible(False)
ax2.get_xaxis().set_visible(False)

ax1.set_ylabel("")
ax2.set_ylabel("")
ax3.set_ylabel("")
ax3.set_xlabel("Target Family", size = sz*1.6)
# then, set a new label on the plot (basically just a piece of text) and move it to where it makes sense (requires trial and error)
f.text(0.1, 0.5, '# Compounds', va='center', rotation='vertical', fontsize = sz*1.6)

# by default, seaborn also gives each subplot its own legend, which makes no sense at all
# soe remove both default legends first
ax1.get_legend().remove()
ax2.get_legend().remove()
ax3.get_legend().remove()
# then create a new legend and put it to the side of the figure (also requires trial and error)
ax2.legend(loc=(0.61, 1.3), fontsize = 30, title="")

# let's put some ticks on the top of the upper part and bottom of the lower part for style
ax1.xaxis.tick_top()
ax2.xaxis.tick_bottom()

# finally, adjust everything a bit to make it prettier (this just moves everything, best to try and iterate)
f.subplots_adjust(left=0.15, right=0.85, bottom=0.15, top=0.85)

d = .01  # how big to make the diagonal lines in axes coordinates
# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
ax1.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)
ax2.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
ax2.plot((1 - d, 1 + d), (-d, +d), **kwargs)   # bottom-right diagonal

kwargs.update(transform=ax3.transAxes)  # switch to the bottom axes
ax3.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax3.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

plt.xticks(rotation=90)
ax1.set_title('Compounds by Predicted Target Families', fontsize = sz*2.3,  y=1.05)
# display new plot again

# values display
plt.margins(y = 0.15)
for p in ax3.patches:
    height = p.get_height()
    if height < 350:
        delta = 20
        ax3.text(x = p.get_x()+(p.get_width()/2), y = height+delta, s = "{:.0f}".format(height),  
            ha = "center", rotation = 90, size = 20.0) 
    if (height > 550) & (height < 900):
        delta = 20
        ax2.text(x = p.get_x()+(p.get_width()/2), y = height+delta, s = "{:.0f}".format(height),  
            ha = "center", rotation = 90, size = 20.0)  
    if height > 20000:
        delta = 3000
        ax1.text(x = p.get_x()+(p.get_width()/2), y = height+delta, s = "{:.0f}".format(height),  
            ha = "center", rotation = 90, size = 20.0) 

# https://stackoverflow.com/questions/50452455/plt-show-does-nothing-when-used-for-the-second-time
from IPython.display import display
display(f)  # Shows plot again

plt.show()

