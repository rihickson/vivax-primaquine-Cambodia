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


def run_post_processing(flag_make_plots=True, flag_elimination_comparison=True, annual_aggregation=True):
    # flag_make_plots = True  # if true, make new plots
    # flag_elimination_comparison = True # if true, generate output spreadsheet for how frequently elimination occurs in sampled ensemble
    # annual_aggregaton = True # if true aggregate values annually to provide smoother plots of just annual case totals

    prov_list = ['Pursat', 'Mondul_Kiri', 'Kampong_Chhnang', 'Battambang', 'Takeo', 'Pailin'] #['Battambang'] #

    base_incidence = ['high', 'low']

    path_pre = './Results/'
    path_post = '/Raw_results/'
    object_path_post = '/Objects/'
    date = '20211012' #'20200427'
    user = 'rmh'
    
    scenarios = {'Status quo': {
        'filename': f'export_raw_Calibrated_{date}.xlsx',
        'object_filename': f'sampled_results_Calibrated_{date}.obj',
        'scen_name': 'no_prim',
        'scen_line_col': 'r',
        },
        'Best case PQ males 15+':{
        'filename': f'export_raw_scenario_0.9_0.9_0.88_males_{date}.xlsx',
        'object_filename': f'sampled_results_scenario_0.9_0.9_0.88_males_{date}.obj',
        'scen_name': 'best prim m15+',
        'scen_line_col': 'b',            
        },
        # 'perfect PQ males 15+':{
        # 'filename': f'export_raw_scenario_1.0_1.0_1.0_males_{date}.xlsx',
        # 'object_filename': f'sampled_results_scenario_1.0_1.0_1.0_males_{date}.obj',
        # 'scen_name': 'perfect prim m15+',
        # 'scen_line_col': 'g',       
            
        # },
        'Perfect radical cure':{
        'filename': f'export_raw_scenario_1.0_1.0_1.0_all_{date}.xlsx',
        'object_filename': f'sampled_results_scenario_1.0_1.0_1.0_all_{date}.obj',
        'scen_name': 'perfect prim all',
        'scen_line_col': 'y',   
        },
        }
    # filenames = [scen['filename'] for scen in scenarios.values()]
    object_filenames = [scen['object_filename'] for scen in scenarios.values()]
    
    elimination_pars = ['indig_cases', 'h_cases', 'h_inci_diag'] #parameters that could be considered stages of elimination to output detailed analysis for

    place_figs = './generated_figures/' #where to save results
    # if the directory we specify to put the files doesn't exist, create it
    import os
    if not os.path.isdir(place_figs):
        os.mkdir(place_figs)

    pop_names = ['M 15+', 'Gen']
    pop_plots = {pop: [pop] for pop in pop_names}
    pop_plots['Total'] = pop_names
    
    scen_names = [scen['scen_name'] for scen in scenarios.values()]
    scen_plot_names = list(scenarios.keys())
    scen_line_col = [scen['scen_line_col'] for scen in scenarios.values()]
    time_step = 5.0/365.  # 5 day timestep is the default, update if you changed this
    
    par_plots = {#'indigenous_cases': {'indig_cases': ('New indigenous cases', '-')},
                 #'active_cases':     {'h_cases': ('Incident active cases', '-')},
                 'diagnosed_cases':  {'h_inci_diag': ('Diagnosed cases', '-')}} #CRITICAL these must be SUMMABLE by population.

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
    
        tex_to_write += start_fig
        for baseline in base_incidence:  # for each baseline incidence calibration
            prov_base_data = dict()
        
            for scen_ind, scen in enumerate(scen_names):  # for each primaquine scenario (no, males, full)
                summary_str = f'{prov}_{baseline}_{scen_plot_names[scen_ind]}'
                
                print (f'Reading... {summary_str}')
                prov_base_data[scen] = dict()
                # read ensemble data
                sr = sc.loadobj(path_pre + prov + '_' + baseline + '_' + date + '_' + user + object_path_post + object_filenames[scen_ind])
                num_runs = len(sr)
                run_years = sr[0][0].t
                plot_years = range(int(min(run_years)), int(max(run_years))+1, 1) if annual_aggregation else run_years
                
                """Save elimination data to spreadsheets"""
                if flag_elimination_comparison:
                    #get a baseline incidence in 2019 for comparison purposes
                    base_case_year = 2019.
                    base_cases = 0.
                    
                    if elim_detailed_runs is None: #first run, set up how many cases there are in the top row
                        elim_detailed_runs = [['Province', 'Baseline scen', 'Baseline incidence', 'Scenario', 'Parameter'] + list(range(num_runs))]
                    
                    for elim_par in elimination_pars:
                        elim_by = []
                        
                        for run, result in enumerate(sr):
                            tot_cases = result[0].get_variable(elim_par)[0].vals + result[0].get_variable(elim_par)[1].vals
                            elimination_starts = None
                            elimination_confirmed = np.inf
                            for ind, t in enumerate(result[0].t):
                                if t >= base_case_year and t - 1. < base_case_year:
                                    base_cases += tot_cases[ind]
                                if tot_cases[ind] == 0.:
                                    if elimination_starts is None:
                                        elimination_starts = t
                                    elif t - 3. == elimination_starts: #3 full years with no malaria diagnoses
                                        # print (f'elmination in run {run} at {t}!')
                                        elimination_confirmed = t
                                        break
                                else:
                                    # if elimination_starts and t - 3. > elimination_starts:
                                    #     print (f'elimination in run {run} lost at {t} :(')
                                    elimination_starts = None
                            elim_by.append(elimination_confirmed)
                        
                        base_inci = base_cases * 1./time_step
                        
                        summary.append([prov, baseline, base_inci, scen, elim_par] + [float(len(np.where(np.array(elim_by)<target_year+1)[0]))/num_runs for target_year in elim_target_years])
                        elim_detailed_runs.append([prov, baseline, base_inci, scen, elim_par] + [eb if np.isfinite(eb) else "" for eb in elim_by])
                
                if flag_make_plots:
                    pars = list(set([par for parlist in par_plots.values() for par in parlist])) #unique parameters in any of the desired plots
                    for par in pars: #parameters_of_interest is a list of dictionaries
                        prov_base_data[scen][par] = dict()
                        for pop in pop_names:
                            par_pop_array = np.array([sri[0].get_variable(par, pop)[0].vals for sri in sr])
                            if annual_aggregation:
                                time_steps_per_year = 1/time_step
                                par_pop_array = np.array([[np.average(pa[int(yr*time_steps_per_year):int((yr+1)*time_steps_per_year)], axis=0) for yr in range(len(plot_years))] for pa in par_pop_array])
                            prov_base_data[scen][par][pop] = np.median(par_pop_array, axis=0)
            
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
    
                        for scen_ind, scen in enumerate(scen_names):  # for each primaquine scenario (no, males, full)
                            start_year = 2015.
                            start_ind = list(plot_years).index(start_year)
                        
                            for par in parameters_of_interest:
                                style = scen_line_col[scen_ind] + par_line_style[par]
                                label = par_labels[par] + ' ' + scen_plot_names[scen_ind]
                                data = np.sum(np.array([prov_base_data[scen][par][pop_label] for pop_label in pop_plots[pop_key]]), axis=0)
                                pl.plot(plot_years[start_ind:], data[start_ind:], style, label=label)
    
                        # add labels and legend
                        pl.xlabel('Year')
                        pl.ylabel('Number of cases')
                        pl.legend(loc='best')
    
                        # add elimination target vertical dashed line
                        axes = pl.gca()
                        ymax = max(axes.get_ylim())
                        pl.plot([2026, 2026], [0, ymax], ':k', label='Elimination target')
                        pl.ylim((0, ymax))
    
                        # save figures to file, and close so not too many open in python
                        filename = prov + '_' + baseline + '_' + par_plot_label + '_' + pop_key + '.png'
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

    if flag_elimination_comparison:
        sc.savespreadsheet(filename='elimination_comparison.xlsx', data = [summary, elim_detailed_runs], folder = place_figs, sheetnames = ['summary', 'details'])
                 
        if flag_make_plots:
            load_and_plot_bars(place_figs)
                   
if __name__ == '__main__':
    run_post_processing(flag_all_prevalences=False, flag_resort_data=True, flag_make_plots=True)


