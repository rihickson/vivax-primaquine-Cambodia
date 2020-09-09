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


def run_post_processing(flag_all_prevalences=True, flag_resort_data=True, flag_make_plots=True):

    # flag_all_prevalences = True  # if true, writes out all data for all pops. Otherwise, writes out dataframe as per DJP's request
    # flag_resort_data = False  # if true, make spreadsheets of reordered data
    # flag_make_plots = True  # if true, make new plots

    prov_list = ['Pursat', 'Mondul_Kiri', 'Kampong_Chhnang', 'Battambang', 'Takeo', 'Pailin']

    base_incidence = ['high', 'low']

    path_pre = './Results/'
    path_post = '/Raw_results/'
    date = '20200427'
    user = 'guest'
    filenames = ['export_raw_scenario_0.0_0.0_0.0_' + date + '.xlsx', 'export_raw_scenario_1.0_1.0_1.0_' + date + '.xlsx']

    pop_names = ['M15+', 'Gen']
    scen_names = ['no_prim', 'prim']
    scen_plot_names = ['no PQ', 'PQ']
    time_step = 365.0/5.0  # 5 day timestep is the default, update if you changed this

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

    if flag_all_prevalences:
        file_to_open = 'aggregated_results_prevalences.csv'
        title_line = ['Scenario']
        [title_line.append(str(x)) for x in range(2011, end_year+1)]
        # title_line = ['Scenario', '2011', '2012', '2013', '2014', '2015', '2016', '2017', '2018', '2019', '2020', '2021', '2022', '2023', '2024', '2025']
    else:
        file_to_open = 'aggregated_results.csv'
        title_line = ['Province', 'Baseline incidence', 'Intervention scenario', 'Year', 'Prevalence in M15+']

    if flag_resort_data:
        with open(file_to_open, mode='w') as csvfile:
            csv_writer = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

            if not flag_all_prevalences:  # only need to write once for this format
                csv_writer.writerow(title_line)

            for prov in prov_list:  # for each province

                if flag_all_prevalences:  # write province and repeat title line for this format
                    csv_writer.writerow([prov])
                    csv_writer.writerow(title_line)

                for baseline in base_incidence:  # for each baseline incidence calibration

                    for scen_ind in range(2):  # for each primaquine scenario (no or full)

                        # read data from relevant sheet
                        reader = xlrd.open_workbook(path_pre + prov + '_' + baseline + '_' + date + '_' + user + path_post + filenames[scen_ind])
                        worksheet = reader.sheet_by_index(0)

                        # pull out column names
                        times = []
                        for col in range(worksheet.ncols):
                            times.append(worksheet.cell_value(0, col))

                        # need to pull out `M 15 +: h_cases` = Row 48; `Gen: h_cases` = Row 251
                        # values from col 4
                        rows = [47, 250]  # -1 because python indexes from 0
                        col_data_start = 4
                        data = []

                        for pop_ind in range(2):  # skip the first row since it was column names

                            elm = []
                            elm_check = []
                            # get column indices for start and end of each year for amalgamation
                            year_vals = np.unique([int(x) for x in times[col_data_start:]])
                            start_year_indices = []
                            end_year_indices = []
                            for current_year in year_vals:
                                all_matches = [times.index(i) for i in times if str(current_year) in str(i)]
                                start_year_indices.append(all_matches[0])
                                end_year_indices.append(all_matches[-1])

                            for current_year_ind in range(len(start_year_indices)-1):
                                col_year_start = start_year_indices[current_year_ind]
                                col_year_end = end_year_indices[current_year_ind]
                                year_sum = 0
                                for col in range(col_year_start, col_year_end + 1):  # +1 because open bracket end of range
                                    year_sum += (worksheet.cell_value(rows[pop_ind], col) / time_step)

                                elm.append(year_sum / int(pops[prov][pop_names[pop_ind]][current_year_ind]))  #NB: just use `year_sum` here if want raw numbers of cases
                                if not flag_all_prevalences:
                                    # write desired years to file
                                    if pop_ind==0 and year_vals[current_year_ind] == 2020:
                                        csv_writer.writerow([prov, baseline, scen_names[scen_ind], '2020', elm[-1]])
                                    if pop_ind==0 and year_vals[current_year_ind] == 2025:
                                        csv_writer.writerow([prov, baseline, scen_names[scen_ind], '2025', elm[-1]])
                            data.append(elm)

                        if flag_all_prevalences:
                             # write to file
                            current_scenario = baseline + '_' + scen_names[scen_ind]
                            for scen_ind in range(2):  # for each population type, add the scenario description then write to file
                                data[scen_ind].insert(0, pop_names[scen_ind] + '_' + current_scenario)
                                csv_writer.writerow(data[scen_ind])


                print(prov + ' done')

    # construct plots for paper from raw data files
    if flag_make_plots:
        # # tracking things to plot -- if wanted to store
        # active_cases = dict()
        # latent_cases = dict()
        # need to pull out `M 15 +: h_cases` = Row 48; `Gen: h_cases` = Row 251
        # values from col 4
        active_malaria_rows = [35, 238]  # all_inf; M15+, Gen
        latent_malaria_rows = [5, 208]  # hL; M15+, Gen,           all these are -1 because python indexes from 0
        col_data_start = 4
        line_col = ['b', 'r']

        # directories etc
        place_figs = './generated_figures/'
        tex_filename = place_figs + 'pv_pq_figs.tex'

        # if the directory we specify to put the files doesn't exist, create it
        import os
        if not os.path.isdir(place_figs):
            os.mkdir(place_figs)

        # base tex needed for the figure
        start_fig = r'\begin{figure}[h!]' + '\n' \
        r'\centering' + '\n'

        def end_fig(prov):
            end_fig_str = r'\caption{\csentence{PQ intervention for \pv~in ' + prov + r'.} Vertical dashed line is the current elimination target for \pv.}\label{fig:pq_' + prov + r'}' + '\n' \
        r'\end{figure}' + '\n' + '\n'
            return end_fig_str

        def mid_fig(pop, scenario, figfile):
            sub_cap_str = r'\subcaptionbox{' + pop + ', ' + scenario + r'.\label{' + pop + '_' + scenario + r'}}{\includegraphics[width=.45\linewidth]{' + figfile + r'}} ' + '\n'
            return sub_cap_str

        tex_to_write = ''

        for prov in prov_list:  # for each province

            # # adding prov to things to plot dict
            # active_cases[prov] = dict()
            # latent_cases[prov] = dict()

            # one 4 part figure per prov, start tex
            tex_to_write += start_fig

            for baseline in base_incidence:  # for each baseline incidence calibration
                # # adding prov to things to plot dict
                # active_cases[prov][baseline] = dict()
                # latent_cases[prov][baseline] = dict()

                for pop_ind in range(2):  # skip the first row since it was column names
                    # # adding prov to things to plot dict
                    # active_cases[prov][baseline][pop_names[pop_ind]] = dict()
                    # latent_cases[prov][baseline][pop_names[pop_ind]] = dict()

                    pl.figure()

                    for scen_ind in range(2):  # for each primaquine scenario (no or full)


                        # read data from relevant sheet
                        reader = xlrd.open_workbook(path_pre + prov + '_' + baseline + '_' + date + '_' + user + path_post + filenames[scen_ind])
                        worksheet = reader.sheet_by_index(0)
                        # pull out column names
                        times = []
                        for col in range(worksheet.ncols):
                            times.append(worksheet.cell_value(0, col))


                        all_inf = []
                        latent = []
                        # get column indices for start and end of each year for amalgamation
                        year_vals = list(np.unique([int(x) for x in times[col_data_start:]]))
                        start_year_indices = []
                        end_year_indices = []
                        for current_year in year_vals:
                            all_matches = [times.index(i) for i in times if str(current_year) in str(i)]
                            start_year_indices.append(all_matches[0])
                            end_year_indices.append(all_matches[-1])

                        for current_year_ind in range(len(start_year_indices)):
                            col_year_start = start_year_indices[current_year_ind]
                            col_year_end = end_year_indices[current_year_ind]
                            year_sum_all_inf = 0  # this is `all_ac, hE, haP, haPs, haPt`
                            year_sum_lt = 0
                            for col in range(col_year_start, col_year_end + 1):  # +1 because open bracket end of range
                                year_sum_all_inf += worksheet.cell_value(active_malaria_rows[pop_ind], col) / time_step
                                year_sum_lt += (worksheet.cell_value(latent_malaria_rows[pop_ind], col)) / time_step

                            all_inf.append(year_sum_all_inf)
                            latent.append(year_sum_lt)
                            # latent.append(year_sum_lt * int(pops[prov][pop_names[pop_ind]][current_year_ind]))  # if using `lt_prev`


                        # # adding prov to things to plot dict
                        # active_cases[prov][baseline][scen_names[scen_ind]] = all_inf
                        # latent_cases[prov][baseline][scen_names[scen_ind]] = all_lt

                        # make plot
                        # from 2020
                        year_ind = year_vals.index(2015)
                        pl.plot(year_vals[year_ind:-1], all_inf[year_ind:-1], line_col[scen_ind] + '-', label='Active Pv, ' + scen_plot_names[scen_ind])
                        pl.plot(year_vals[year_ind:-1], latent[year_ind:-1], line_col[scen_ind] + '--', label='Latent Pv, ' + scen_plot_names[scen_ind])

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
                    filename = prov + '_' + baseline + '_' + pop_names[pop_ind] + '.png'
                    pl.savefig(place_figs + filename)
                    pl.close('all')

                    # add info to string to write to tex file later
                    tex_to_write += mid_fig(pop=pop_names[pop_ind], scenario=baseline, figfile=filename)

            # add end of figure text to string to write to tex file later
            tex_to_write += end_fig(prov=prov)

        # write the 4-part figures into tex
        f = open(tex_filename, 'w')
        f.write(tex_to_write)
        f.close()

if __name__ == '__main__':
    run_post_processing(flag_all_prevalences=False, flag_resort_data=True, flag_make_plots=True)


