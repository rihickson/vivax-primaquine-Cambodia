# -*- coding: utf-8 -*-
"""
Created on Mon Aug 20 11:37:03 20118

@author: RMH
"""

from malaria_utils import *
import gc
from scens_cambodia import *

def run_me(book_key):

    start_year = 2011
    end_year = 2041.0
    uncertainty_flag = True  # whether or not to run with uncertainty

    primdict = {'coverage': [0.9, 1.0], 'sensitivity_G6PDd': [0.9, 1.0],'efficacy': [0.88, 1.0]}  # extreme range to determine beyond best case possibility of elimination by 2025

    # date = '20200427'  # Use this line if you want to continue working on results from a previous day, otherwise it will save results to a new folder!
    date = '20211012'
    utils_malaria.date = date

    country = 'cambodia'
    book_str = '' if book_key == '' else book_key + '_'
    currency = 'USD'
    project_name = '%s_%s%s' % (country, book_str, date)

    unfunded_progs = []  # programs that should be removed from the displayed budget and optimization entirely (not part of malaria budget)
    coverage_progs = ['Treat', 'Test']  # programs that have fixed coverage instead of varying cost

    # # # specify how funding might be adjusted over time
    time_adjustments = {
        'LLIN_m': {'start_time': 2011., 'end_time': end_year, 'cover_time': 0.5, 'timespan': 3.0, 'adjust_type': 'flat'}}

    project_folder = './Project/' #get_apps_folder(country=country)

    # %%OPTIONS Adjust these parameters to change the quality of the output
    todo = [  'plots',
        'raw_results'
    ]
    save_figs = True  # Probably leave as True, but False if you want to just see them in the console instead without saving
    load_objects = True  # re-use previously generated optimization results and uncertainty sampling if they exist (otherwise rerun)
    test_run = False  # run optimizations and sampling for a bare minimum number of times to test that everything is working
    key_year = 2017  # for e.g. calibration cascades
    plot_language = 'English'

    confidence_interval = 'quartile'  # options: quartile for 75% range, ci for 95% range, samples for beautiful plots with every line individually that may not have appropriate legends
    plot_against_data_only = False  # if true, only output all the line plots where data exists (useful/faster for calibration!)
    np.random.seed(72)  # set a random seed for reproducible results
    run_parallel = False  # True may or may not work

    # get paths for model framework and data book, expected to be in a `<country>/Projects/` directory structure
    framework_path = get_paths(project_folder, 'framework', version='latest', verbose=True)
    databook_path = get_paths(project_folder, ['Databook', book_key], version='latest', verbose=True)
    if not framework_path or not databook_path:
        raise Exception('Could not find a valid framework or databook, exiting.')

    if test_run:
        plot_quality = 'preview'  # preview = low dpi, include legends, final = high dpi, separate legends
        num_samples = 100
    else:
        plot_quality = 'final'  # preview = low dpi, include legends, final = high dpi, separate legends
        num_samples = 1000  # number of samples for the uncertainty plots - 1000 may be necessary for good results

    results_folder = project_folder + '..%sResults%s%s%s' % (sep, sep, book_str + user_version(),
                                                             sep)  # e.g. project files are in COUNTRY\Project\, results saved in COUNTRY\Results\date\
    objects_folder = results_folder + 'Objects' + sep
    raw_results_folder = results_folder + 'Raw_results' + sep

    start_logging(folder=results_folder)  # create a log file of all the console output

    parset_name = 'default'
    progset_name = ''  # 'default'

    """
    Instructions below are each of the different scenarios that will be run and output
    Each instruction will be run separately with uncertainty and plots will be saved to a different folder in combination with the plot_years
    This will take a lot of disk space, prune unused results!
    The example 'calibrated' instruction is a default run with no programs used
    The example 'Primaquine' instruction is explores the effect of Primaquine being provided to males 15+ from October 2020
    """
    # %%SCENARIOS
    # each scenario has a function defined in opt_country.py
    scen_fns = sc.odict()
    scen_fns['Calibrated'] = scenfn_calibration
    scen_fns['Primaquine males 15+'] = scenfn_primaquine_males
    scen_fns['Primaquine all']       = scenfn_primaquine_all


    plot_sets = {
        # 'Calibration': {'resnames': ['Calibrated'], 'plot_years': [start_year, end_year], 'plot_uncertainty': uncertainty_flag,
        #                 'plots': ['keypars']},
        'Primaquine males 15+': {'resnames': ['Calibrated', 'Primaquine males 15+'], 'plot_years': [start_year, end_year], 'plot_uncertainty': uncertainty_flag,
                       'plots': ['keypars', 'cascade_advanced']},
        'Primaquine all': {'resnames': ['Calibrated', 'Primaquine all'], 'plot_years': [start_year, end_year], 'plot_uncertainty': uncertainty_flag,
                       'plots': ['keypars', 'cascade_advanced']}
    }

# %%RUN EVERYTHING
    results = []

    # %%LOAD FRAMEWORK AND DATABOOK
    F = at.ProjectFramework(framework_path)
    P = at.Project(name=project_name, framework=F, sim_dt=5. / 365., sim_start=2011., sim_end=end_year, stochastic = True, do_run=False)
    P.load_databook(databook_path=databook_path, make_default_parset=True, do_run=False)
    
    
    # %%CREATE SCENARIOS (inc calibration)
    def_prog_args = {'parset': parset_name, 'progset': progset_name,
                     'progset_instructions': at.ProgramInstructions(start_year=2019)}
    def_prog_result = try_loading(function=P.run_sim, fnargs=def_prog_args,
                                  obj_filename='Baseline_with_progs_%s' % (date),
                                  obj_folder=objects_folder, run_if_unfound=True, load_objects=load_objects)

    for scen_name in scen_fns.keys():
        if np.array([scen_name in ps['resnames'] for ps in plot_sets.values()]).any():
            inc_fn = scen_fns[scen_name]
            scen = inc_fn(P=P, scen_name=scen_name, result=def_prog_result, progset_name=progset_name,
                          parset_name=parset_name, time_adjustments=time_adjustments, coverage_progs=coverage_progs, coverage=primdict['coverage'], sensitivity_G6PDd=primdict['sensitivity_G6PDd'], efficacy=primdict['efficacy'])
            if isinstance(scen, at.Scenario):
                P.scens[scen_name] = scen
            elif isinstance(scen,
                            dict):  # if the function returned a dict of scenarios instead of a single scenario, each with their own name
                for detailed_scen_name in scen.keys():
                    print('Adding %s' % detailed_scen_name)
                    P.scens[detailed_scen_name] = scen[detailed_scen_name]
                # also replace the name in any result sets with all the names of the sub-results
                for plot_set_name in plot_sets.keys():
                    rns = plot_sets[plot_set_name]['resnames']
                    if scen_name in rns:
                        ind = rns.index(scen_name)
                        plot_sets[plot_set_name]['resnames'] = rns[:ind] + list(scen.keys()) + rns[ind + 1:]
            else:
                raise Exception('Cannot convert %s to one or more scenarios.' % (type(scen)))

    # %%RUN ALL SCENARIOS (inc optimizations)
    for scen_name in P.scens.keys():
        scen = P.scens[scen_name]
        print('Loading or running %s' % (scen_name))
        run_if_unfound = np.array([scen_name in ps['resnames'] for ps in plot_sets.values()]).any()
        print('Run if unfound: %s' % (run_if_unfound))

        if run_if_unfound:
            parset = scen.get_parset(parset=P.parsets[parset_name], project=P)
            progset = None if isinstance(scen, at.ParameterScenario) else scen.get_progset(
                progset=P.progsets[progset_name], project=P)
            instructions = None if isinstance(scen, at.ParameterScenario) else scen.get_instructions(
                progset=P.progsets[progset_name], project=P)
            run_args = {'parset': parset, 'progset': progset, 'progset_instructions': instructions}
            result = try_loading(function=P.run_sim, fnargs=run_args,
                                 obj_filename='Scenario result_%s_%s' % (scen_name, date),
                                 obj_folder=objects_folder, run_if_unfound=True, load_objects=load_objects)
            result.name = scen_name
            results += [result]

    # %%PLOT RESULTS
    if 'plots' in todo or 'comparisons' in todo:
        for plot_set in plot_sets.keys():
            gc.collect()
            plot_years = plot_sets[plot_set]['plot_years']
            plots = plot_sets[plot_set]['plots']
            uncertainty = plot_sets[plot_set]['plot_uncertainty']
            local_results = [results[[rn.name for rn in results].index(res)] for res in plot_sets[plot_set]['resnames']]

            if 'plots' in todo:
                local_samples = []
                for res_name in [res.name for res in local_results]:
                    if uncertainty and ('plots' in todo or 'raw_results' in todo):
                        print('Loading or running sampled results with uncertainty for %s for use in %s' % (
                        res_name, plot_set))
                        scen = P.scens[res_name]
                        parset = scen.get_parset(parset=P.parsets[parset_name], project=P)
                        progset = None if isinstance(scen, at.ParameterScenario) else scen.get_progset(
                            progset=P.progsets[progset_name], project=P)
                        instructions = None if isinstance(scen, at.ParameterScenario) else scen.get_instructions(
                            progset=P.progsets[progset_name], project=P)
                        run_args = {'parset': parset, 'progset': progset, 'progset_instructions': instructions,
                                    'n_samples': num_samples, 'parallel': run_parallel}
                        sampled_result = try_loading(function=P.run_sampled_sims, fnargs=run_args,
                                                     obj_filename='sampled_results_%s_%s.obj' % (res_name, date),
                                                     obj_folder=objects_folder, run_if_unfound=True,
                                                     load_objects=load_objects)
                    else:
                        sampled_result = None  # don't load potentially large sampled object files unless it's necessary!

                    local_samples += [sampled_result]

                for plot in plots:
                    type_str = 'uncertainty' if uncertainty else 'best'
                    year_str = '%s-%s' % (int(P.settings.sim_start) if plot_years[0] is None else plot_years[0],
                                          int(P.settings.sim_end) if plot_years[1] is None else plot_years[1])
                    plots_folder = results_folder + str.join('_', [plot_set] + [plot] + [year_str] + [type_str]) + sep
                    plot_pars = {'P': P, 'results': local_results, 'results_folder': plots_folder,
                                 'save_figs': save_figs,
                                 'plot_years': plot_years, 'plot_quality': plot_quality,
                                 'confidence_interval': confidence_interval,
                                 'plot_against_data_only': plot_against_data_only, 'sampled_results': local_samples,
                                 'key_year': key_year, 'currency': currency, 'unfunded_progs': unfunded_progs,
                                 'plot_language': plot_language, 'progset_name': progset_name,
                                 'time_adjustments': time_adjustments, }
                    run_plots([plot], plot_pars)

                print('Finished output for  %s, clearing sampled_results from memory.' % (plot_set))
                del local_samples
                gc.collect()

            if 'comparisons' in todo:
                at.export_results(results=local_results,
                                  filename=results_folder + '%s_export_comparison_%s.xlsx' % (plot_set, date))
                export_quick_comparison(results=local_results, variables=comparison_variables,
                                        year_dict=comparison_years,
                                        filename='%s_key_change_comparison_%s.xlsx' % (plot_set, date),
                                        folder=results_folder)


    # %%SAVE USEFUL INFORMATION ABOUT THE RESULTS
    run_info = 'These results generated using %s\n' % __file__.split(sep)[-1]
    run_info += 'Framework: %s\n' % framework_path.split(sep)[-1]
    run_info += 'Databook: %s\n' % databook_path.split(sep)[-1]
    save_run_info(run_info=run_info, folder=results_folder)
    if 'raw_results' in todo:
        export_raw_results(results, raw_results_folder)

    try:
        P.save(filename='%s.prj' % project_name, folder=results_folder)
    except:
        print('WARNING: Could not save project file.')


if __name__ == '__main__':
    book_key_all = ['Pailin_low']
    np.random.seed(10)

    global book_key

    for book_key in book_key_all:
        run_me(book_key)
