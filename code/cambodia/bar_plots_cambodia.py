# -*- coding: utf-8 -*-
"""
Created on Sun Nov 14 03:43:51 2021

@author: Rowan
"""

import json
import numpy as np
from scipy.integrate import solve_ivp
import scipy.stats as st
import matplotlib
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator
import matplotlib.gridspec as gridspec
from datetime import date
import scipy.optimize
import csv
import xlrd  # for reading excel
from malaria_utils import *
import pandas as pd
import itertools

def barplots(df, place_figs, filename):
    matplotlib.style.use('grayscale')
    matplotlib.rcParams['axes.prop_cycle'] = matplotlib.cycler('color', ['#1F497D', '#558ED5', '#8EB4E3', '#C6D9F1', '#EEEEEE', '#FFFFFF', '#CCCCCC']) * matplotlib.cycler('hatch', ['', '/'])
    matplotlib.rcParams['hatch.linewidth'] = 5.
    matplotlib.rcParams['hatch.color'] = '#FFFFFF'
    # matplotlib.rcParams['grid.color'] = '#DDDDDD'
    # matplotlib.rcParams['lines.color'] = '#DDDDDD'
    # matplotlib.rcParams['axes.edgecolor'] = '#DDDDDD'
    # matplotlib.rcParams['axes.labelcolor'] = '#FFFFFF'
    matplotlib.rcParams['patch.edgecolor'] = '#FFFFFF'


    arr = df.values
    inds = [2,3,4,5]
    for rind in range(len(arr)):
        arr[rind] = list(arr[rind][:2])+[arr[rind][i] - arr[rind][i-1] for i in inds]
    df.loc[:,:] = arr

    
    df = df.groupby(["Province", "Scenario"], sort=False)
    df = df.sum()
        
    clusters = df.index.levels[0]
    inter_graph = 0
    maxi = np.max(np.sum(df, axis=1))
    total_width = len(df)+inter_graph*(len(clusters)-1)
    
    fig = plt.figure(figsize=(total_width,10))
    gridspec.GridSpec(1, total_width)
    axes=[]
    
    ax_position = 0
    for cluster in clusters:
        subset = df.loc[cluster]
        ax = subset.plot(kind="bar", stacked=True, width=0.8, ax=plt.subplot2grid((1,total_width), (0,ax_position), colspan=len(subset.index)))
        axes.append(ax)
        ax.set_title(cluster)
        ax.set_xlabel("")
        ax.set_ylim(0,maxi)
        # ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        ax_position += len(subset.index)+inter_graph
        
        #top section of each
        ax.patches[-1].set_hatch('/')
        ax.patches[-2].set_hatch('/')
        ax.patches[-3].set_hatch('/')
    
    for i in range(1,len(clusters)):
        axes[i].set_yticklabels("")
        axes[i-1].legend().set_visible(False)
    axes[0].set_ylabel("Proportion of outcomes from 300 model iterations")
    
    # fig.suptitle('Timing of elimination of indigenous transmission of Plasmodium vivax malaria by radical cure scenario\nHigh baseline incidence\n', fontsize="x-large")
    handles, labels = axes[-1].get_legend_handles_labels()
    legend = axes[-1].legend(handles=reversed(handles), labels=reversed(labels), loc='center right',
                             fontsize=16, bbox_to_anchor=(2.7, 0.5), framealpha=1).get_frame()
    legend.set_linewidth(1)
    legend.set_edgecolor("black")
    
    
    
    pl.savefig(place_figs + filename, bbox_inches='tight', transparent=True)
    # pl.close('all')

def load_and_plot_bars(place_figs):
    filename = place_figs + 'elimination_comparison.xlsx'
    data = pd.read_excel(filename)
    df = pd.DataFrame(data)
    df.replace('Mondul_Kiri', 'Mondulkiri', inplace=True)
    df.replace('Kampong_Chhnang', 'Kampong Chhnang', inplace=True)
    df.replace('no_prim', 'Status quo', inplace=True)
    df.replace('best prim m15+', 'Best case PQ males 15+', inplace=True)
    df.replace('perfect prim all', 'Perfect radical cure', inplace=True)
    
    df = df.set_index('Province')
    
    df.pop('Baseline incidence')
    df['No elimination by 2040'] = 1.0
    
    df = df.loc[df['Parameter'] == 'indig_cases']
    df.pop('Parameter')
    
    hv = df.loc[df['Baseline scen'] == 'high']
    hv.pop('Baseline scen')
    lv = df.loc[df['Baseline scen'] == 'low']
    lv.pop('Baseline scen')
    
    
    
    barplots(hv, place_figs, 'high_baseline')
    barplots(lv, place_figs, 'low_baseline')

if __name__ == '__main__':
    place_figs = 'D:\\GitHub\\vivax-primaquine-Cambodia\\code\\cambodia\\generated_figures_300_runs_seed_89\\'
    load_and_plot_bars(place_figs)
    
    