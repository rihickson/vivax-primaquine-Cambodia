# -*- coding: utf-8 -*-
"""
Created on Mon Nov 15 01:53:04 2021

@author: Rowan
"""
# -*- coding: utf-8 -*-
"""

@author: RIH
"""

# Written to aggregate the results of the runs in `cambodia.py` to feed into the decision support around surveillance by DJP

import json
import numpy as np
from scipy.integrate import solve_ivp
import scipy.stats as st
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator
import matplotlib.gridspec as gridspec
from datetime import date
import scipy.optimize
import csv
import xlrd  # for reading excel
from malaria_utils import *
import pandas as pd
from bar_plots_cambodia import *


def run_baseline_processing(flag_make_plots=True, annual_aggregation=True, res_rename=''):
    flag_make_plots = True  # if true, make new plots
    flag_elimination_comparison = True # if true, generate output spreadsheet for how frequently elimination occurs in sampled ensemble
    # annual_aggregaton = True # if true aggregate values annually to provide smoother plots of just annual case totals

    prov_list = ['Pursat', 'Mondul_Kiri', 'Kampong_Chhnang', 'Battambang', 'Takeo', 'Pailin'] #['Battambang'] #

    base_incidence = ['high', 'low']
    base_color = ['r', 'b']

    path_pre = f'./Results{res_rename}/'
    path_post = '/Raw_results/'
    object_path_post = '/Objects/'
    date = '20211012' #'20200427'
    user = 'rmh'
    
    scenarios = {'Status quo': {
        'filename': f'export_raw_Calibrated_{date}.xlsx',
        'best_filename': f'Scenario result_Calibrated_{date}',
        'object_filename': f'sampled_results_Calibrated_{date}.obj',
        'scen_name': 'no_prim',
        'scen_line_col': 'r',
        },
        # 'Best case PQ males 15+':{
        # 'filename': f'export_raw_scenario_0.9_0.9_0.88_males_{date}.xlsx',
        # 'object_filename': f'sampled_results_scenario_0.9_0.9_0.88_males_{date}.obj',
        # 'scen_name': 'best prim m15+',
        # 'scen_line_col': 'b',            
        # },
        # 'perfect PQ males 15+':{
        # 'filename': f'export_raw_scenario_1.0_1.0_1.0_males_{date}.xlsx',
        # 'object_filename': f'sampled_results_scenario_1.0_1.0_1.0_males_{date}.obj',
        # 'scen_name': 'perfect prim m15+',
        # 'scen_line_col': 'g',       
            
        # },
        # 'Perfect radical cure':{
        # 'filename': f'export_raw_scenario_1.0_1.0_1.0_all_{date}.xlsx',
        # 'object_filename': f'sampled_results_scenario_1.0_1.0_1.0_all_{date}.obj',
        # 'scen_name': 'perfect prim all',
        # 'scen_line_col': 'y',   
        # },
        }
    # filenames = [scen['filename'] for scen in scenarios.values()]
    object_filenames = [scen['object_filename'] for scen in scenarios.values()]
    best_filenames   = [scen['best_filename']   for scen in scenarios.values()]
    
    place_figs = f'./generated_figures{res_rename}/' #where to save results
    # if the directory we specify to put the files doesn't exist, create it
    import os
    if not os.path.isdir(place_figs):
        os.mkdir(place_figs)

    pop_names = ['M 15+', 'Gen']
    pop_plots = {pop: [pop] for pop in pop_names}
    # pop_plots['Total'] = pop_names
    # pop_plots = {'Total': ['M 15+', 'Gen']}
    
    scen_names = [scen['scen_name'] for scen in scenarios.values()]
    scen_plot_names = list(scenarios.keys())
    scen_line_col = [scen['scen_line_col'] for scen in scenarios.values()]
    time_step = 5.0/365.  # 5 day timestep is the default, update if you changed this
    
    par_plots = {#'indigenous_cases': {'indig_cases': ('New indigenous cases', '-.')},
                 'active_cases':     {'h_cases': ('Incident active cases', '-')},
                 #'diagnosed_cases':  {'h_inci_diag': ('Diagnosed cases', '-')}
                # 'cases': {'indig_cases': ('New indigenous cases', '-'), 'h_cases': ('Incident active cases', '-')},
                 } #CRITICAL these must be SUMMABLE by population.
    pars = list(set([par for parlist in par_plots.values() for par in parlist])) #unique parameters in any of the desired plots

    # import the populations
    population_data = []
    pops = dict()
    for prov in prov_list:
        pops[prov] = {}

    with open(r'../../csvs/cambodia_population_estimates.csv') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        header = next(reader)
        end_year = int(header[-1])
        for row in reader:
            population_data.append(row)
            if row[0] in prov_list:
                pops[row[0]][row[1]] = row[2:]

    """Setup"""
    if flag_elimination_comparison:
        print ('Generating a comparison sheet of how often elimination targets may be met for each sampled scenario')
        elim_target_years = [2025, 2030, 2035, 2040]
        
        elim_detailed_runs = None #store when elimination happens in each run of the model
        num_runs = None
        summary = [['Province', 'Baseline scen', 'Baseline incidence', 'Scenario', 'Parameter'] + [f'Elimination by {year}' for year in elim_target_years]] #store when elimination occurs.



    """Loop through each set of results, only load each result once"""
    
    tex_to_write = ''   
    # base tex needed for the figure
    start_fig = r'\begin{figure}[h!]' + '\n' + r'\centering' + '\n'                     
    
    for prov in prov_list:  # for each province
        
        #only need to get data once - raw data is the same as input in both low and high, only calibration is different
        proj_folder = path_pre + prov + '_' + base_incidence[0] + '_' + date + '_' + user + '//'
        filename = f'cambodia_{prov}_{base_incidence[0]}_{date}.prj'
        P = sc.loadobj(folder=proj_folder, filename=filename)
        # return P
    
        tex_to_write += start_fig
        prov_base_data = dict()
        prov_best_data = dict()
        for baseline in base_incidence:  # for each baseline incidence calibration
            prov_base_data[baseline] = dict()
            prov_best_data[baseline] = dict()
        
            for scen_ind, scen in enumerate(scen_names):  # for each primaquine scenario (no, males, full)
                summary_str = f'{prov}_{baseline}_{scen_plot_names[scen_ind]}'
                
                print (f'Reading... {summary_str}')
                prov_base_data[baseline][scen] = dict()
                prov_best_data[baseline][scen] = dict()
                # read ensemble data
                sr = sc.loadobj(path_pre + prov + '_' + baseline + '_' + date + '_' + user + object_path_post + object_filenames[scen_ind])
                num_runs = len(sr)
                run_years = sr[0][0].t
                plot_years = range(int(min(run_years)), int(max(run_years)), 1) if annual_aggregation else run_years

                #read the best deterministic data
                srb = sc.loadobj(path_pre + prov + '_' + baseline + '_' + date + '_' + user + object_path_post + best_filenames[scen_ind])

                if flag_make_plots:
                    
                    for par in pars: #parameters_of_interest is a list of dictionaries
                        prov_base_data[baseline][scen][par] = dict()
                        prov_best_data[baseline][scen][par] = dict()
                        for pop in pop_names:
                            par_pop_array = np.array([sri[0].get_variable(par, pop)[0].vals for sri in sr])
                            if annual_aggregation:
                                time_steps_per_year = 1/time_step
                                par_pop_array = np.array([[np.average(pa[int(yr*time_steps_per_year):int((yr+1)*time_steps_per_year)], axis=0) for yr in range(len(plot_years))] for pa in par_pop_array])
                            prov_base_data[baseline][scen][par][pop] = par_pop_array
                            
                            par_pop_array_best = srb.get_variable(par, pop)[0].vals
                            if annual_aggregation:
                                time_steps_per_year = 1/time_step
                                par_pop_array_best = np.array([np.average(par_pop_array_best[int(yr*time_steps_per_year):int((yr+1)*time_steps_per_year)], axis=0) for yr in range(len(plot_years))])
                            prov_best_data[baseline][scen][par][pop] = par_pop_array_best
            
        """Now do things with prov_base_data (for a given province and baseline calibration)"""
        if flag_make_plots:
                            
            for par_plot_label in par_plots.keys():
                par_plot = par_plots[par_plot_label]
                parameters_of_interest = list(par_plot.keys())
                par_labels     = {par: par_plot[par][0] for par in parameters_of_interest}
                par_line_style = {par: par_plot[par][1] for par in parameters_of_interest}
                                
                def end_fig(prov):
                    end_fig_str = r'\caption{\csentence{PQ intervention for \pv~in ' + prov + r'.} Vertical dashed line is the current elimination target for \pv.}\label{fig:pq_' + prov + r'}' + '\n' \
                r'\end{figure}' + '\n' + '\n'
                    return end_fig_str
        
                def mid_fig(pop, scenario, figfile):
                    sub_cap_str = r'\subcaptionbox{' + pop + ', ' + scenario + r'.\label{' + pop + '_' + scenario + r'}}{\includegraphics[width=.45\linewidth]{' + figfile + r'}} ' + '\n'
                    return sub_cap_str
        
                for pop_plot_ind, pop_key in enumerate(pop_plots):  # skip the first row since it was column names
                    print (f'Generating figure for {prov}, {baseline}, {pop_key}')
                    pl.figure()
                    
                    for base_ind, baseline in enumerate(base_incidence):
                        for scen_ind, scen in enumerate(scen_names):  # for each primaquine scenario (no, males, full)
                            start_year = 2012.
                            try:
                                start_ind = list(plot_years).index(start_year)
                            except:
                                if start_year < min(plot_years):
                                    start_year = min(plot_years)
                                    start_ind = 0
                                else:
                                    raise Exception(f'Try a different year - {start_year} not found either in or before {list(plot_years)}')
                        
                            for par in parameters_of_interest:
                                style = base_color[base_ind] + par_line_style[par]
                                label = par_labels[par] + ' ' + scen_plot_names[scen_ind] if len(parameters_of_interest)>1 else scen_plot_names[scen_ind]
                                for res_ind in range(len(prov_base_data[baseline][scen][par][pop_plots[pop_key][0]])):
                                    data = np.sum(np.array([prov_base_data[baseline][scen][par][pop_label][res_ind] for pop_label in pop_plots[pop_key]]), axis=0)
                                    pl.plot(plot_years[start_ind:], data[start_ind:], style, label=None, alpha=0.01)
                                
                                data = np.sum(np.array([prov_best_data[baseline][scen][par][pop_label] for pop_label in pop_plots[pop_key]]), axis=0)
                                pl.plot(plot_years[start_ind:], data[start_ind:], style, label=f'Calibration to {baseline} incidence', alpha=1)
                                
                    #Only plot the data once at the end       
                    valid_inds = np.where(np.array(P.parsets[0].pars[par].ts[0].t) > start_year)
                    case_data = np.sum(np.array([P.parsets[0].pars[par].ts[pop_label].vals for pop_label in pop_plots[pop_key]]), axis=0)
                    case_data = np.array(case_data)[valid_inds]
                    time = np.array(P.parsets[0].pars[par].ts[0].t)[valid_inds]
                    pl.plot(time, case_data, 'ko', label='CNM estimated vivax malaria cases')

                    # add labels and legend
                    pl.xlabel('Year')
                    ylabel = 'Number of cases' if len(parameters_of_interest)>1 else par_labels[parameters_of_interest[0]]
                    pl.ylabel(ylabel)
                    pl.legend(loc='best')

                    # add elimination target vertical dashed line
                    axes = pl.gca()
                    ymax = min(max(axes.get_ylim()), 4*max(case_data))
                    # pl.plot([2026, 2026], [0, ymax], ':k', label='Elimination target')
                    pl.ylim((0, ymax))

                    # save figures to file, and close so not too many open in python
                    filename = 'calibration_' + prov + '_' + par_plot_label + '_' + pop_key + '.pdf'
                    pl.savefig(place_figs + filename)
                    pl.close('all')

                    # add info to string to write to tex file later
                    tex_to_write += mid_fig(pop=pop_key, scenario=baseline, figfile=filename)

        # add end of figure text to string to write to tex file later
        tex_to_write += end_fig(prov=prov)
                    
    if flag_make_plots:
        # directories etc
        tex_filename = place_figs + 'pv_pq_figs.tex'
        # write the 4-part figures into tex
        print ('Generating tex output for figures')
        f = open(tex_filename, 'w')
        f.write(tex_to_write)
        f.close()


if __name__ == '__main__':
    P = run_baseline_processing(res_rename='_300_runs_seed_89')


