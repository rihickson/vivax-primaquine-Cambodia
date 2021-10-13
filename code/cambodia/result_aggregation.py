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
from datetime import date
import scipy.optimize
import csv
import xlrd  # for reading excel
from malaria_utils import *


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
    filenames = [f'export_raw_Calibrated_{date}.xlsx', f'export_raw_scenario_1.0_1.0_1.0_males_{date}.xlsx', f'export_raw_scenario_1.0_1.0_1.0_all_{date}.xlsx']
    object_filenames = [f'sampled_results_Calibrated_{date}.obj',  f'sampled_results_scenario_1.0_1.0_1.0_males_{date}.obj', f'sampled_results_scenario_1.0_1.0_1.0_all_{date}.obj']

    place_figs = './generated_figures/' #where to save results
    # if the directory we specify to put the files doesn't exist, create it
    import os
    if not os.path.isdir(place_figs):
        os.mkdir(place_figs)

    pop_names = ['M 15+', 'Gen']
    pop_plots = {pop: [pop] for pop in pop_names}
    pop_plots['Total'] = pop_names
    
    scen_names = ['no_prim', 'prim m15+', 'prim_all']
    scen_plot_names = ['no PQ', 'PQ males 15+', 'PQ all']
    scen_line_col = ['r', 'b', 'g']
    time_step = 5.0/365.  # 5 day timestep is the default, update if you changed this
    
    parameters_of_interest = {'all_par': 'Latent prevalence', 'h_cases': 'Incident cases'} #CRITICAL these must be SUMMABLE by population.
    par_line_style = ['--', '-']

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
        summary = [['Scenario'] + [f'Elimination by {year}' for year in elim_target_years]] #store when elimination occurs.



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
                elim_by = []
                print (f'Reading... {summary_str}')
                prov_base_data[scen] = dict()
                # read ensemble data
                sr = sc.loadobj(path_pre + prov + '_' + baseline + '_' + date + '_' + user + object_path_post + object_filenames[scen_ind])
                num_runs = len(sr)
                run_years = sr[0][0].t
                plot_years = range(int(min(run_years)), int(max(run_years))+1, 1) if annual_aggregation else run_years

                """Save elimination data to spreadsheets"""
                if flag_elimination_comparison:
                    if elim_detailed_runs is None: #first run, set up how many cases there are in the top row
                        
                        elim_detailed_runs = [['Scenario \ Run'] + list(range(num_runs))]
                    for run, result in enumerate(sr):
                        tot_diag = result[0].get_variable('h_inci_diag')[0].vals + result[0].get_variable('h_inci_diag')[1].vals
                        elimination_starts = None
                        elimination_confirmed = np.inf
                        for ind, t in enumerate(result[0].get_variable('h_inci_diag')[0].t):
                            if tot_diag[ind] == 0.:
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
                    
                    summary.append([summary_str] + [float(len(np.where(np.array(elim_by)<target_year+1)[0]))/num_runs for target_year in elim_target_years])
                    elim_detailed_runs.append([summary_str] + [eb if np.isfinite(eb) else "" for eb in elim_by])

                for par in parameters_of_interest.keys():
                    prov_base_data[scen][par] = dict()
                    for pop in pop_names:
                        par_pop_array = np.array([sri[0].get_variable(par, pop)[0].vals for sri in sr])
                        if annual_aggregation:
                            time_steps_per_year = 1/time_step
                            par_pop_array = np.array([[np.average(pa[int(yr*time_steps_per_year):int((yr+1)*time_steps_per_year)], axis=0) for yr in range(len(plot_years))] for pa in par_pop_array])
                        prov_base_data[scen][par][pop] = np.median(par_pop_array, axis=0)
            
            """Now do things with prov_base_data (for a given province and baseline calibration)"""
            if flag_make_plots:
                       
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
                    
                        for p_ind, par in enumerate(parameters_of_interest.keys()):
                            style = scen_line_col[scen_ind] + par_line_style[p_ind]
                            label = parameters_of_interest[par] + ' ' + scen_plot_names[scen_ind]
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
                    filename = prov + '_' + baseline + '_' + pop_key + '.png'
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
                                    
if __name__ == '__main__':
    run_post_processing(flag_all_prevalences=False, flag_resort_data=True, flag_make_plots=True)


