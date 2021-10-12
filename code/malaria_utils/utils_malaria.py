# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 12:13:25 2019

@author: rowan.martin-hughes
"""

import atomica as at
import sciris as sc
import numpy as np

from os import sep, path

from datetime import datetime
time = datetime.now().strftime("%H:%M:%S")
date = datetime.now().strftime("%Y%m%d")

def _get_gdrive_folder(country=None):
    import socket
    user = socket.gethostname()
    
    if country is None:
        country_str = ''
    else:
        country_str = '%s%sProject%s'%(country, sep, sep)

    apps_folder = '../'
    user = 'rmh'
    user_initials = 'rmh'
        
    return path.join(path.abspath(apps_folder),''), user_initials  # ensure trailing separator

gdrive, user_initials = _get_gdrive_folder()
def user_version(): #defined as a function as 'date' might be changed to match previous results
    return '%s_%s' % (date, user_initials)

root = path.abspath(path.join(at.parent_dir(),'..'))+path.sep # repository root dir

def get_malaria_analyses_version():
    return at.fast_gitinfo(__file__)

def try_loading(function, fnargs, obj_filename, obj_folder, load_objects=True, run_if_unfound=True):
    """helper function to try loading a pre-run object from a folder, or recreate it if not available"""
    try:
        assert load_objects
        obj = sc.loadobj(obj_filename, folder=obj_folder)
    except:
        if run_if_unfound:
            obj = function(**fnargs)
            sc.saveobj(obj=obj, filename=obj_filename, folder=obj_folder)
        else:
            return None
    return obj



def get_apps_folder(country=None):
 
    if country is None:
        country_str = ''
    else:
        country_str = country #'%s%sProject%s'%(country, sep, sep)
        
    apps_folder = gdrive + country_str
    
    return apps_folder

def get_paths(folder, inclusions=[], extensions=['.xlsx'], exclusions=['~'],
             version='latest', verbose=True, case_sensitive=False, folder_depth=0):
    #find the highest alphabetical databook, framework, and progbook in the project_folder
    #if multiple exist, would pick the one with the highest date
    #ignore anything with a '~' as these are temporary files
    from glob import glob
    contains = sc.promotetolist(inclusions)+sc.promotetolist(extensions)
    excludes = sc.promotetolist(exclusions)
    extra_str = '*%s'%(sep) * folder_depth + '*'
    files = glob(folder+extra_str)
    if case_sensitive:
        valid_files = [file for file in files if np.array([con.lower() in file.lower() for con in contains]).all()]
        valid_files = [file for file in valid_files if np.array([not exc.lower() in file.lower() for exc in excludes]).all()]
    else:
        valid_files = [file for file in files if np.array([con in file for con in contains]).all()]
        valid_files = [file for file in valid_files if np.array([not exc in file for exc in excludes]).all()]
        
    if len(valid_files)==0:
        if verbose: print('Could not find a file containing all of %s and none of %s in folder %s'%(contains, excludes, folder))
        return None
    elif version=='latest':
        valid_files = valid_files[-1]
    
    if verbose: print('Found \'%s\' file(s): %s'%(inclusions, valid_files))
    return valid_files
    
    
def save_run_info(run_info: str = '', folder=None, filename=None):
    if filename is None:
        filename = 'Run info_%s.txt'%(date)
    if folder is None:
        folder = '' #save it into the script folder
    
    filepath = folder+filename
    
    run_info+= 'Date %s, %s\n'%(date, time)
    run_info+= 'Atomica version %s\n'%at.__version__
    run_info+= 'Atomica git details %s\n'%(at.__gitinfo__)
    run_info+= 'Sciris version %s, %s\n'%(sc.__version__, sc.__versiondate__)
    run_info+= 'Malaria analyses git details %s\n'%(get_malaria_analyses_version())
    
    try:
        file = open(filepath, 'w+')
        file.write(run_info)
        file.close()
        return True
    except:
         print('WARNING: Could not save results to %s, file is probably already open.'%(filepath))
         return False
    
def export_raw_results(results, results_folder):
    if isinstance(results, at.Project):
        results = results.results.values()
    
    import os
    os.makedirs(results_folder, exist_ok=True)
    for result in results:
        results_path = results_folder+'export_raw_%s_%s.xlsx'%(result.name, date)
        try:
            result.export_raw(filename=results_path)
        except: 'WARNING: could not save raw results %s.'%(results_path)

def getpops(P, pop='hum'):
    return [key for key, details in P.data.pops.items() if details['type']==pop]

def allequal(x):
    '''return true if all elements of a list/array are equal'''
    return len(set(x)) <=1

def sigfigs(x):
    if type(x)==str: return x
    return sc.sigfig(X=x, sigfigs=2, SI=False, sep=True, keepints=False)

def start_logging(filename=None, folder=None):
    """save an output file of console output"""
    import logging
    import os
    
    if filename is None: filename='script_logging.log'
    
    if not folder is None:
        try:
            os.makedirs(folder, exist_ok=True)
        except:
            print('ERROR: Could not create a results folder at %s'%(folder))
            return False
        
    logger = logging.getLogger()
    h = logging.FileHandler(folder+filename,mode='w')
    logger.addHandler(h)
    return logger





def _filter_pops_by_output(result, output) -> list:
    """
    Helper function for plotting quantities
    Copied from results.py

    With population types, a given output/output aggregation may only be defined
    in a subset of populations. To deal with this when plotting Result objects,
    it's necessary to work out which population the requested output aggregation can be
    plotted in. This function takes in an output definition and returns a list of populations
    matching this.

    :param output: An output aggregation string e.g. 'alive' or ':ddis' or {['lt_inf','lteu']} (supported by PlotData/get_variable)
    :return: A list of population code names

    """

    if sc.isstring(output):
        vars = result.get_variable(output)
    elif isinstance(output, list):
        vars = result.get_variable(output[0])
    elif isinstance(output, dict):
        v = list(output.values())[0]
        if isinstance(v, list):
            vars = result.get_variable(v[0])
        elif sc.isstring(v):
            # It could be a function aggregation or it could be a single one
            _, deps = parse_function(v)
            vars = result.get_variable(deps[0])
    else:
        raise Exception('Could not determine population type')
    return [x.pop.name for x in vars]

def get_data(results,  output, tdict):
    """
    Convert an output to a DataFrame for a group of results
    Copied from results.py and adapted

    This function takes in a list of results, and an output specification recognised by :class:`PlotData`.
    It extracts the outputs from all results and stores them in a 3-level MultiIndexed dataframe, which is
    returned. The index levels are the name of the output, the name of the results, and the populations.

    In addition, this function attempts to aggregate the outputs, if the units of the outputs matches
    known units. If the units lead to anver obvious use of summation or weighted averating, it will be used.
    Otherwise, the output will contain NaNs for the population-aggregated results, which will appear as empty
    cells in the Excel spreadsheet so the user is able to fill them in themselves.

    :param results: List of Results
    :param output_name: The name to use for the output quantity
    :param output: An output specification/aggregation supported by :class:`PlotData`
    :param tdict: Outputs will be interpolated onto the times in this dictionary of lists of tvecs (typically would be annual)
    :return: a PlotData

    """
    output_name = output
    
    pops = _filter_pops_by_output(results[0], output)
    pop_labels = {x: y for x, y in zip(results[0].pop_names, results[0].pop_labels) if x in pops}
    data = dict()
    
    popdata = at.PlotData(results, pops=pops, outputs=output)
    if popdata.series[0].units in {at.FrameworkSettings.QUANTITY_TYPE_NUMBER, results[0].model.pops[0].comps[0].units}:
        action_fn = sum
        pop_aggregation = 'sum'
    elif popdata.series[0].units in {at.FrameworkSettings.QUANTITY_TYPE_FRACTION,
                                         at.FrameworkSettings.QUANTITY_TYPE_PROPORTION,
                                         at.FrameworkSettings.QUANTITY_TYPE_PROBABILITY} or 'prev' in output:
        action_fn = np.mean
        pop_aggregation = 'weighted'
    else:
        print ('Not clear which units to use over time for %s, sum.'%(output_name))
        action_fn = sum

    for tname, tvals in tdict.items():
        popdata = at.PlotData(results, pops=pops, outputs=output)
        assert len(popdata.outputs) == 1, 'Framework plot specification should evaluate to exactly one output series - there were %d' % (len(popdata.outputs))
        popdata.interpolate(tvals)
        
        for result in popdata.results:
            for pop_name in popdata.pops:
                time_vals = popdata[result, pop_name, popdata.outputs[0]].vals
                
                data[(output_name, popdata.results[result], pop_labels[pop_name], tname)] = action_fn(time_vals)
                
                
        agg_popdata = at.PlotData(results, outputs=output, pops={'total': pops}, pop_aggregation=pop_aggregation)
        agg_popdata.interpolate(tvals)
        for result in agg_popdata.results:
            data[(output_name, agg_popdata.results[result], 'Total (sum)', tname)] = action_fn(agg_popdata[result, agg_popdata.pops[0], agg_popdata.outputs[0]].vals)
    
#        df = pd.DataFrame(data, index=tvals)
#        df = df.T
#        df.index = df.index.set_names(['output', 'result', 'pop'])  # Set the index names correctly so they can be reordered easily
    return data

    
def export_quick_comparison(results, variables, pops = None, res_names = None, year_dict={'2016': 2016, '2016 - 2030': list(range(2016, 2031)), '2030': 2030},
                    filename='parameter comparison.xlsx', folder=None):
    data = [['Parameter (display)', 'Parameter (model)', 'Self-comparison', 'Same-year', 'Population', 'Baseline result', 'Baseline year(s)', 'Baseline value',
             'Outcome result', 'Outcome year(s)', 'Outcome value', 'Percentage difference', 'Absolute difference']]
    all_resnames = [res.name for res in results]
    if res_names is None:
        res_names = all_resnames
    if not pops is None:
        comp_pops = pops #same for all
        
    
    for var in variables:
        var_name = results[0].model.framework.get_variable(var)[0]['display name'] #assume results all have the same framework...

        for brn in res_names:
            for byn, bys in year_dict.items():
                bys = sc.promotetolist(bys)
                for orn in res_names:
                    for oyn, oys in year_dict.items():
                        oys = sc.promotetolist(oys)
                        #it only makes sense to compare changes between the same result and different years, or the same year and different results.
                        if bys == oys or (brn == orn and len(oys) == len(bys)): 
                            df = get_data(results=results, output=var, tdict = year_dict)
#                            print (df.keys())
                            if pops is None:
                                comp_pops = list(set([x[2] for x in df.keys()])) #assume pops are the same through all results.......
                            for pop in comp_pops:
                                self_comp = 'Y' if brn == orn else ''
                                same_year = 'Y' if bys == oys else ''
                                base_val = df[var, brn, pop, byn]
                                out_val  = df[var, orn, pop, oyn]
                                
                                percent_diff = (out_val - base_val)/max(base_val, 1e-15)
                                abs_diff     = (out_val - base_val)
        
                                data.append([var_name, var, self_comp, same_year, pop, brn, byn, base_val,
                                             orn, oyn, out_val, percent_diff, abs_diff])
    sc.savespreadsheet(filename=filename, data=data, folder=folder)
