# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 12:47:51 2021

@author: Rowan
"""

from malaria_utils import *
    
date = '20211012'
utils_malaria.date = date

    
project_folder = './Project/' #get_apps_folder(country=country)
res_rename='_300_runs_seed_89'
results_folder = f'./generated_figures{res_rename}/' #where to save results

year_range = range(2011, 2020) #for output

def run_me(book_key):
    country = 'cambodia'
    book_str = '' if book_key == '' else book_key + '_'
    currency = 'USD'
    project_name = '%s_%s%s' % (country, book_str, date)
    

    framework_path = get_paths(project_folder, 'framework', version='latest', verbose=True)
    databook_path = get_paths(project_folder, ['Databook', book_key], version='latest', verbose=True)
    
    F = at.ProjectFramework(framework_path)
    start_year = 2011
    end_year = 2041.0
    P = at.Project(name=project_name, framework=F, sim_dt=5. / 365., sim_start=2011., sim_end=end_year, stochastic = False, do_run=False)
    P.load_databook(databook_path=databook_path, make_default_parset=True, do_run=False)
    
    """
    mls_rate: \item the rate of developing malaria-like symptoms for each person in a given year; 
    p_rdt: \item the daily rate of testing for people with non-severe malaria-like symptoms such as fever; 
    p_rdt_s: \item the daily rate of testing for people with severe malaria-like symptoms; 
    E_act: \item the average duration of the latent period (i.e. until hypnozoite reactivation); 
    inc_clear: \item the proportion incompletely clearing hypnozoites after naturally recovering; and 
    aP_prop: \item the proportion of new malaria cases that are asymptomatic   
    sP_prop: proportion severe
    avg_dur_imm: duration of asymptomatic
    uPdeath_prop: death if uncomplicated malaria left untreated
    sPdeath_prop: death if severe malaria left untreated
    m_pop_variation: seasonality of malaria biting
    """
    
    
    par_keys = ['mls_rate', 'p_rdt', 'p_rdt_s', 'E_act', 'inc_clear', 'aP_prop', 'sP_prop', 'h_bite_rate', 'avg_dur_imm', 'm_pop_variation']
    pops = ['M 15+', 'Gen', 'A. Funestus']
    
    local_data = []
    
    assert len(P.parsets.keys())==1, 'Too many parsets!'
    for par in par_keys:
        for pop in pops:
            if pop in P.parsets[0].pars[par].ts.keys():
                if P.parsets[0].pars[par].ts[pop].has_time_data:
                    vals = list(P.parsets[0].pars[par].ts[pop].interpolate(year_range))
                    best = np.mean(vals)
                    low = min(vals)
                    high = max(vals)
                else:
                    vals = ['' for _ in year_range]
                    best = P.parsets[0].pars[par].ts[pop].assumption
                    low = ''
                    high = ''
                local_data.append([book_key, par, pop, best, low, high] + vals)
    
    return local_data
    


if __name__ == '__main__':
    book_key_all = ['Pursat_low', 'Pursat_high', 'Mondul_Kiri_low', 'Mondul_Kiri_high', 'Kampong_Chhnang_low', 'Kampong_Chhnang_high', 'Battambang_low', 'Battambang_high', 'Takeo_low', 'Takeo_high', 'Pailin_low', 'Pailin_high']

    global book_key

    data = [['Province', 'Parameter', 'Population', 'Best (mean)', 'Min', 'Max'] + list(year_range)]
    
    for book_key in book_key_all:
        data += run_me(book_key)
    
    import sciris as sc
    sc.savespreadsheet(filename='calibration_parameters.xlsx', folder=results_folder, data=data)