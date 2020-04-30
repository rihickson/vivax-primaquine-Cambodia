# -*- coding: utf-8 -*-
"""
Created on Mon Aug 20 11:37:03 2018

@author: User
"""

import atomica as at
import sciris as sc
import numpy as np
from os import sep
import matplotlib.pyplot as plt
import atomica.plotting as aplt

from malaria_utils import getpops, get_optimized_results, allequal

from matplotlib import ticker

defaultblue = (55/256.,126/256.,185/256.,1.)

def addconfidence(P, figs, labels_to_use=[], handlers_to_use=None, colors_to_use=None, pops=None, datapars=None,
                  plot_against_data_only=False, maincolor=None, lightcolor=None):
    '''
    Add confidence intervals and/or data points (possibly with error bars) to a figure after checking years etc are appropriate
    param: labels should be a list of the form ['Data line label', 'Label for confidence interval', 'Label for data estimate']
    param: figs is assumed to be either one figure with a legend, or two figures where the first is a plot and the second is a legend.
            - adds aconfidence intervals to figs[0], replaces the legend for figs[-1]
    params: labels_to_use, handlers_to_use, colors_to_use are for setting up the legend (text label, style of marker in the legend, color of the line in the legend)
    param: pops should be either a single population, or a set of three where they represent best, low, high
    param: datapars should be either a single parameter (either because there is only a 'best' value to use for a single pop, or because best, low, high are specified by different populations within the same parameter),
            or a list of three parameters, used for the best, low, and high values respectively
    '''
    legendsettings = {'loc': 'center left', 'bbox_to_anchor': (1.05, 0.5), 'ncol': 1, 'framealpha':0, 'borderpad':0, }
    
    if maincolor is None:
        maincolor = defaultblue
    if lightcolor is None:
        lightcolor = (maincolor[0],maincolor[1],maincolor[2],0.5) 
#    transpcolors = [(col[0],col[1],col[2],0.5) for col in colors] #50% transparency
    
    
    fig = figs[0]
    if handlers_to_use is None:
        handlers_to_use = ['patch' if 'uncertainty' in ltu else 'line' for ltu in labels_to_use]
    if colors_to_use is None:
        colors_to_use = [lightcolor if 'uncertainty' in ltu else maincolor for ltu in labels_to_use]
    
    keep_plot = False if plot_against_data_only else True
    
#    try:
    for a in [1]:
        if not datapars is None:
            if pops is None:
                pops = ['Total (best)', 'Total (low)', 'Total (high)']
            elif type(pops)==type('string'):
                pops = [pops, pops, pops] #best, low high use the same population, different parameters
            assert len(pops)==3 #must be best, low, high pops
            
            
            if type(datapars)==type('string'):
                if pops[0]==pops[1] and pops[0]==pops[2]:
                    datapars = [datapars, None, None] #there is only one population specified, and only one datapar specified: just use the 'best' value for data and don't plot low/high
                else:
                    datapars = [datapars, datapars, datapars] #getting best low and high from the same parameter, using different 'populations'
            assert len(datapars)==3 #must be best, low, high data pars
            
            datapars = [datapar if datapar in P.data.tdve.keys() else None for datapar in datapars] #remove any datapars if data doesn't actually exist!            
            pops     = [pop if (not datapars[i] is None) and pop in P.data.tdve[datapars[i]].ts.keys() else None for i, pop in enumerate(pops)]
            
            y_limit = fig.axes[0].get_ylim()[1] #current max based on data, may need to extend y_limit upper bound based on data input
            
            bestp, lowp, highp = pops
            bestd, lowd, highd = datapars
            ibest = None if bestd is None else P.data.tdve[bestd]
            years = []
            lows = []
            highs = []
            #for each year that is in both low and high estimates of incidence
            if not (lowd is None or highd is None or lowp is None or highp is None):    
                ilow  = P.data.tdve[lowd]
                ihigh = P.data.tdve[highd]
                
                for yl, year in enumerate(ilow.ts[lowp].t):
                    if year in ihigh.ts[highp].t:
                        yh = ihigh.ts[highp].t.index(year)
                        years.append(year)
                        lows.append(ilow.ts[lowp].vals[yl])
                        highs.append(ihigh.ts[highp].vals[yh])
                if len(years)>0:
                    if len(years)==1:
                        if len(ibest.ts[bestp].t)==0: #unless there's a confidence interval but no best estimate
                            colors_to_use += [lightcolor]
                            handlers_to_use += ['circle']
                            labels_to_use += ['Data (uncertainty range)']
                        #skip label as it should be clear from best value
                        fig.axes[0].errorbar(years, [0], yerr = [[-lows[0]], [highs[0]]], fmt='-')
                    else:
                        colors_to_use += [lightcolor]
                        handlers_to_use += ['patch']
                        labels_to_use+= ['Data (uncertainty range)']
                        fig.axes[0].fill_between(years, lows, highs, facecolor=lightcolor, interpolate=True, alpha = 0.6)
                    y_limit = max(y_limit, max(highs)*1.05)

            #plot the best data points - do this after the range so it appears on top!
            if (not ibest is None) and (not bestp is None) and (len(ibest.ts[bestp].t)>0):
                y_limit = max(y_limit, max(ibest.ts[bestp].vals)*1.05)
                sigma = sc.dcp(P.data.tdve[bestd].ts[bestp].sigma)
#                col = maincolor
                if len(years)==0 and not sigma is None: #uncertainty was entered and there aren't any explicit lows/highs
                    fig.axes[0].errorbar(ibest.ts[bestp].t,ibest.ts[bestp].vals, yerr=sigma, fmt='o', elinewidth=3, ecolor=maincolor)
                else:
                    fig.axes[0].scatter(ibest.ts[bestp].t,ibest.ts[bestp].vals, marker='o', s=40, linewidths=3,edgecolor='face',facecolor=lightcolor)
                
                colors_to_use+=[maincolor] #[(col[0][0], col[0][1],col[0][2],0.5)]
                handlers_to_use+= ['circle']
                labels_to_use+= ['Data (best estimate)']
            
            fig.axes[0].set_ylim(top=y_limit)
#    except Exception as e:
#        print('WARNING: tried to add uncertainty for figure with labels %s, pops %s, datapars %s but failed, skipping uncertainty!\n%s'%(labels_to_use, pops, datapars, e))

#    fig.axes[0].legend(labels_to_use, **legendsettings)
    entries = sc.odict(zip(labels_to_use, colors_to_use))
    figs[-1] = at.plot_legend(entries=entries, plot_type=handlers_to_use, fig=figs[-1], legendsettings=legendsettings)

    if np.array(['Data' in ltu for ltu in labels_to_use]).any():
        keep_plot = True
    if not keep_plot:
        figs = None #get rid of the figure TODO implementation could be cleaner!

    return figs

def pair_plot(P, results, results_folder=None, save_figs=False, plot_years=[2000, 2030], sampled_results = None, **kwargs):
    plot_outputs = ['h_bite_rate', 'h_nets_coverage', 'foi']
    plot_pops = getpops(P, 'hum')
    
    colors = sc.gridcolors(n=len(results))
    
#    legends = ['Model (%s)'%result.name for result in results]
#    
#    
#    fig = None
    figs = None
    
    for sr, sampled_result in enumerate(sampled_results):
        if not sampled_result is None:  #include uncertainty as sampled
            baseline = sc.dcp(results[sr])
            mapping_function = lambda x: at.PlotData(x,outputs=plot_outputs, pops= plot_pops, pop_aggregation='weighted')
            ensemble = at.Ensemble(mapping_function=mapping_function)

#            legends += ['Model (uncertainty range %s)'%baseline.name]  
            baseline.name = sampled_result[0][0].name # baseline name needs to match sampled results name
            ensemble.update(sampled_result)
            for sm in ensemble.samples:
                sm.set_colors(colors[sr], overwrite=True)
            
            ensemble.set_baseline(baseline)
            ensemble.baseline.set_colors(colors[sr], overwrite=True)
            ensemble.name = 'PAIR_TEST_%s'%(sampled_result[0][0].name)
            
            figs = ensemble.pairplot(year=plot_years[1], pops=plot_pops, outputs=plot_outputs)
            
            fnames = ['%s_%s'%(ensemble.name, pop) for pop in plot_pops]
            if save_figs: at.save_figs(figs, path=results_folder, prefix='', fnames=fnames)
    

def standard_plot(P, results, res_pars=None, data_pars=None, pop_aggregation='sum', res_pops=None, data_pops=None, title=None,
                  ylabel=None, results_folder=None, save_figs=False, plot_years=[2000, 2030], sampled_results = None, 
                  t_bins=None, plot_against_data_only=False, confidence_interval = 'ci', **kwargs):
    """
    Standard TB plot covering most situations and including both model and input uncertainty
    'results': a list of results all of which to plot against each other
    'res_pars': parameter name to plot from results (or a list to aggregate),
    'data_pars': parameter name used for data (or a list if using uncertainty in entry in [best, low, high] form)
    'pop_aggregation': for populations weighted (e.g. prevalence) or sum (number of cases),
    'res_pops': either a pop name or a list of pops to aggregate
    'data_pops': either a pop name or a list of pops used for data in [best, low, high] form
    'title': the title of the plot
    'ylabel': label for the y axis,
    'results_folder': where to save the plot
    'save_figs': if to save at all
    'sampled_results': a list of sampled_results to plot against each other
    't_bins': None as default plots a continuous line, if t_bins = 1 plot annual totals (weighted or summed as appropriate)
    """
    
    plot_outputs = data_pars if type(data_pars)==dict else [{ylabel: sc.promotetolist(res_pars)}]
    plot_pops = res_pops if type(res_pops)==dict else [{'Total': sc.promotetolist(res_pops)}]
    
    defcolors = sc.gridcolors(ncolors=len(results)) #, hueshift=0.3)
    colors = [(col[0],col[1],col[2],1.0) for col in defcolors] #add in a transparency channel
    lightcolors = [(col[0],col[1],col[2],0.5) for col in defcolors] #add in a transparency channel
    
    legends = ['Model (%s)'%result.name for result in results]
    
    
    fig = None
    figs = None
    handlers_to_use = []
    colors_to_use   = []
    
    if not sampled_results is None:
        assert(len(sampled_results)==len(results)), 'Results: %s, Sampled results: %s'%([res.name for res in results], ['NONE' if res is None else res.name for res in sampled_results])
        for sr, sampled_result in enumerate(sampled_results):
            rgbcol = sc.rgb2hex(defcolors[sr])
            if not sampled_result is None:  #include uncertainty as sampled
                baseline = sc.dcp(results[sr])
                mapping_function = lambda x: at.PlotData(x,outputs=plot_outputs, pops= plot_pops, pop_aggregation=pop_aggregation, t_bins=t_bins)
                ensemble = at.Ensemble(mapping_function=mapping_function)
    
                legends += ['Model (uncertainty range %s)'%baseline.name]  
                baseline.name = sampled_result[0][0].name # baseline name needs to match sampled results name
                ensemble.update(sampled_result)
                for sm in ensemble.samples:
                    sm.set_colors(rgbcol, overwrite=True)
                
                ensemble.set_baseline(baseline)
                ensemble.baseline.set_colors(rgbcol, overwrite=True)
                ensemble.name = str(res_pars)
                
    #            """the naive label for 2017 = the previous 12 months if t_bins = 1. -> manually change all the x_labels by reducing them if using t_bins"""
#                if t_bins is not None:
#                    for series in ensemble.baseline.series:
#                        series.tvec -= t_bins
#                    for pd in ensemble.samples:
#                        for series in pd.series:
#                            series.tvec -= t_bins
                fig = ensemble.plot_series(fig=fig, style='quartile', legend=False) #style='samples' is beautiful, style='quartile' for 75% interval a faster/easier option, style='ci' for 95% confidence

        if not fig is None:
            handlers_to_use = ['line' for _ in sampled_results]+['patch' for _ in sampled_results]
            colors_to_use   = colors + lightcolors
        
            if aplt.settings['legend_mode'] == 'together':
                figs = [fig]
            else:
                figs = [fig, None]
    if figs is None:  #just use the best result and plot as lines
        d = at.PlotData(results, outputs=plot_outputs, pops=plot_pops, pop_aggregation=pop_aggregation, t_bins=t_bins)
        """the naive label for 2017 = the previous 12 months if t_bins = 1. -> manually change all the x_labels by reducing them if using t_bins"""
        d.set_colors(colors=colors, results=d.results, overwrite=True)
#        if t_bins is not None:
#            for series in d.series:
#                series.tvec -= t_bins
        figs = at.plot_series(d, axis='results', plot_type='line')
        colors_to_use = colors
        handlers_to_use = ['line' for res in results]
    
    figs = addconfidence(P, figs, labels_to_use = legends, colors_to_use = colors_to_use, handlers_to_use = handlers_to_use,
                         pops=data_pops, datapars=data_pars, plot_against_data_only=plot_against_data_only) #may add further entries to the legend
    
    if not figs is None:
        fig = figs[0]
        fig.axes[0].set_ylim(bottom=0.)
        fig.axes[0].set_xlim(plot_years[0], plot_years[1])
        fig.axes[0].xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    #    if not ylabel is None: fig.axes[0].set_ylabel(ylabel)
        if not title is None: fig.axes[0].set_title(title)
        
#        if fig.axes[0].get_ylabel() in ['Proportion (proportion)', 'Proportion (proportion per year)']:
#            fig.axes[0].set_ylabel('Proportion')
        if ' (' in fig.axes[0].get_ylabel():
            fig.axes[0].set_ylabel(fig.axes[0].get_ylabel().split(' (')[0]) #get rid of all units
        
    #    label = str(res_pars).replace(':','-')
    #    print ('save figs: %s, saving to %s as %s'%(save_figs, results_folder, title))
        if save_figs: at.save_figs([fig for fig in figs if not fig is None], path=results_folder, prefix='', fnames=[title])

def run_plots(plots, plot_pars):
    language = plot_pars['plot_language'].lower()
    if language == 'english':
        lang_ind = 0
    elif language == 'french':
        lang_ind = 1
    else:
        raise Exception('ERROR: language %s not defined'%language.capitalize())
    plot_pars['lang_ind'] = lang_ind
    
    if plot_pars['plot_quality'] == 'preview':
        aplt.settings['legend_mode'] = 'together'
        aplt.settings['dpi'] = 72
        aplt.settings['transparent'] = False
    elif plot_pars['plot_quality'] == 'final':
        aplt.settings['legend_mode'] = 'separate'
        aplt.settings['dpi'] = 300
        aplt.settings['transparent'] = True
    else:
        raise Exception('plot_quality must be either \'preview\' (low quality, combined legends) or \'final\' (high quality, separate legends).')
    
    if 'demographics' in plots:
        plot_demographics(**plot_pars)
    
    if 'keypars' in plots:  #deaths, incidence
        plot_incidence(times = {'annual total': 1.}, **plot_pars)
        plot_incidence_monthly(**plot_pars)
        plot_deaths(times = {'annual total': 1.}, **plot_pars)
        
    if 'keypars_advanced' in plots:
        plot_incidence(times = {'continuous annualised': None}, **plot_pars)
        plot_deaths(times = {'continuous annualised': None}, **plot_pars)
        plot_prevalence(**plot_pars)
        
    if 'treatment' in plots:
        plot_testing_treatment(**plot_pars)
        
    if 'detailed' in plots:
#        plot_detailed_prevalence(**plot_pars)
        plot_deaths_advanced(**plot_pars)
        plot_diagnoses_advanced(**plot_pars)
        
    #    
    #
    if 'cascade' in plots:
        plot_mls_bars(**plot_pars)
        plot_treatmentcascade(**plot_pars)
        
    if 'cascade_advanced' in plots:
        plot_sym(**plot_pars)
        plot_poptrends(**plot_pars)        
        plot_mls_continuous(**plot_pars)
        # plot_treatmentstatus(**plot_pars) # fixme
    
    if 'misc' in plots:
        #various troubleshooting plots to make sure everything is working correctly
        plot_mosquito_misc(**plot_pars)
        plot_human_misc(**plot_pars)
        
    if 'programs' in plots:
        plot_LLINs(**plot_pars)
        plot_IRS(**plot_pars)
        
    if 'reconciliation' in plots:
        plot_pars['plot_against_data_only'] = False #important to plot all of these if plotting at all!
        plot_reconciliation(**plot_pars)        
        
    if 'spending' in plots:
        plot_budgets(**plot_pars)
        
    #default plots specified in the framework
    if 'default_plots' in plots:
        P.results[0].plot()

    plt.close('all')





def plot_poptrends(P, results, results_folder, save_figs, plot_years, **kwargs):
    allpops = getpops(P, 'hum')
    
    for result in results:
        d = at.PlotData(result, pops=allpops)
        figs = at.plot_series(d, plot_type="stacked") # This should look like the usual Optima-TB result
        if save_figs: at.save_figs(figs, path=results_folder, prefix='Population trends_')

def plot_mls_continuous(P, results, results_folder, save_figs, plot_years, **kwargs):
    allpops = getpops(P, 'hum')
    
    for result in results:
        d = at.PlotData(result, pops=allpops, outputs=['hSs','hLs','haPs','huPs','hsPs'])
        figs = at.plot_series(d, plot_type="stacked") # This should look like the usual Optima-TB result
        if save_figs: at.save_figs(figs, path=results_folder, prefix='Malaria-like symptoms_')
        
def plot_mls_bars(P, results, results_folder, save_figs, plot_years, **kwargs):
    allpops = getpops(P, 'hum')
    
    for pop in allpops:
        d = at.PlotData(results, pops=pop, outputs={'Non-malarial fever':['hSs'], 'Non-malarial fever with malaria parasites (latent or asymptomatic)':['hLs','haPs'], 'Malarial fever':['huPs','hsPs']}, t_bins=5)
        figs = at.plot_bars(d, outer='results', stack_outputs='all') # This should look like the usual Optima-TB result
        figs[0].axes[0].set_title('Malaria-like symptoms %s'%(pop))
        if save_figs: at.save_figs(figs, path=results_folder, prefix='Malaria-like symptoms_%s_'%(pop))

def plot_sym(P, results, results_folder, save_figs, plot_years, **kwargs):
    allpops = getpops(P, 'hum')
    
    for result in results:
        d = at.PlotData(result, pops=allpops, outputs=['all_um', 'all_sm', 'hE', 'haP', 'haPs', 'haPt', 'hL'])
        figs = at.plot_series(d, plot_type="stacked") # This should look like the usual Optima-TB result
        if save_figs: at.save_figs(figs, path=results_folder, prefix='Infectious malaria prevalence_')



def plot_demographics(P, **kwargs):
    #pop size (total and for each pop)
    allpops = getpops(P, 'hum')
    
    for pop in allpops:
        standard_plot(P, res_pars='alive', data_pars='alive', res_pops=pop, data_pops=pop,
                      ylabel='Population size', title='Population size - %s'%(pop), **kwargs)
    standard_plot(P, res_pars='alive', data_pars=None, res_pops=allpops, data_pops=None,
                      ylabel='Population size', title='Population size - total', **kwargs)
    
    
def plot_mosquito_misc(P, lang_ind, **kwargs):
    allpops = getpops(P, 'mos')
    
    pars = ['m_incub', 'm_foi', 'm_mort', 'm_prev', 'm_exp', 'm_sus', 'm_pop_seasonal', 'm_bite_factor']
    
    for par in pars:
        parlabel = P.framework.get_label(par).split('\n')[lang_ind]
        for pop in allpops:
            standard_plot(P, res_pars=par, data_pars=None, res_pops=pop, data_pops=pop,
                          ylabel=parlabel, title='%s - %s'%(parlabel, pop), **kwargs)
        
def plot_human_misc(P, lang_ind, **kwargs):
    allpops = getpops(P, 'hum')
    
    pars = ['foi', 'h_bite_factor', 'mls_seasonal', 'IStx']
    
    for par in pars:
        parlabel = P.framework.get_label(par).split('\n')[lang_ind]
        for pop in allpops:
            standard_plot(P, res_pars=par, data_pars=None, res_pops=pop, data_pops=pop,
                          ylabel=parlabel, title='%s - %s'%(parlabel, pop), **kwargs)
            
#    pair_plot(P, **kwargs)
        

def plot_incidence_monthly(P, lang_ind, **kwargs):
    allpops = getpops(P, 'hum')
    popmapping = sc.odict([(pop, pop) for pop in allpops]+[('Total', allpops)])
    
    #seasonal incidence
    times = {'monthly':1./12.} #'monthly':1./12.,
    pars = ['inc_per1000', 'h_inci', 'h_inci_diag']
    for time_label, t_bins in times.items():
        for par in pars:
            parlabel = 'Malaria cases per 1,000' #P.framework.get_label(par).split('\n')[lang_ind]
            for pop in popmapping.keys():
                data_pars = [None,'h_diag_monthly_low','h_diag_monthly_peak'] if time_label=='monthly' else None
                standard_plot(P, res_pars=par, data_pars=data_pars, t_bins=t_bins,
                              res_pops=popmapping[pop], data_pops=pop, pop_aggregation='weighted', 
                              ylabel='Number of malaria cases per 1,000', title='%s (%s) - %s'%(parlabel, time_label, pop), **kwargs)
 
def plot_incidence(P, lang_ind, times = {'continuous annualised': None}, **kwargs):
    allpops = getpops(P, 'hum')
    popmapping = sc.odict([(pop, pop) for pop in allpops]+[('Total', allpops)])
    
    def time_string(tl): return '' if tl=='annual total' else ' (%s)'%(tl)
    
    #incidence numbers over different time periods
#    times = {'continuous annualised':None, 'annual total':1.} #'monthly':1./12.,
#    for time_label, t_bins in times.items():
#        pars = ['h_u_cases', 'h_s_cases', 'h_cases']
#        for par in pars:
#            parlabel = P.framework.get_label(par).split('\n')[lang_ind]
#            for pop in allpops:
#                    standard_plot(P, res_pars=par, data_pars=None, res_pops=popmapping[pop], data_pops=pop,
#                                  ylabel='Number of malaria cases', title='%s (%s) - %s'%(parlabel, time_label, pop), t_bins=t_bins, **kwargs)
            
    #incidence per 1000
    pars = ['inc_per1000', 'h_inci']
    for time_label, t_bins in times.items():
        for par in pars:
            parlabel = 'Malaria cases per 1,000' #P.framework.get_label(par).split('\n')[lang_ind]
            for pop in popmapping.keys():
                    standard_plot(P, res_pars=par, data_pars=par, res_pops=popmapping[pop], data_pops=pop, pop_aggregation='weighted',
                                  ylabel='Number of malaria cases per 1,000', title='%s%s - %s'%(parlabel, time_string(time_label), pop), t_bins=t_bins, **kwargs)

    #incidence (annual cases)
    pars = ['h_cases']
    for time_label, t_bins in times.items():
        for par in pars:
            parlabel = 'Incident malaria cases' #P.framework.get_label(par).split('\n')[lang_ind]
            for pop in popmapping.keys():
                standard_plot(P, res_pars=par, data_pars=['h_cases', 'h_cases_low', 'h_cases_high'], res_pops=popmapping[pop],
                              data_pops=pop, pop_aggregation='sum',
                              ylabel='Number of malaria cases', title='%s%s - %s'%(parlabel, time_string(time_label), pop), t_bins=t_bins, **kwargs)



def plot_testing_treatment(P, lang_ind, **kwargs):
    allpops = getpops(P, 'hum')
    popmapping = sc.odict([(pop, pop) for pop in allpops]+[('Total', allpops)])
    
    aggpars = {'sum': ['test_annual', 'tx','smtx']} #Prevalence total numbers    
    for agg, pars in aggpars.items():
        for par in pars:
            parlabel = 'Number of people with malaria-like symptoms screened for treatment' if par == 'test_annual' else P.framework.get_label(par).split('\n')[lang_ind]
            for pop in popmapping.keys():
                    standard_plot(P, res_pars=par, data_pars=par, res_pops=popmapping[pop], data_pops=pop,
                                  pop_aggregation=agg, t_bins=1.,
                                  ylabel='Number of tests', title='%s - %s'%(parlabel, pop), **kwargs)
    
    aggpars = {'weighted': ['test_positivity']} #Prevalence total numbers    
    for agg, pars in aggpars.items():
        for par in pars:
            parlabel = P.framework.get_label(par).split('\n')[lang_ind]
#            for pop in allpops:
#                    standard_plot(P, res_pars=par, data_pars=par, res_pops=pop, data_pops=pop, pop_aggregation=agg, t_bins=None,
#                                  ylabel=parlabel, title='%s - %s'%(parlabel, pop), **kwargs)
            standard_plot(P, res_pars=par, data_pars=par, res_pops=allpops, data_pops='Total', pop_aggregation=agg, t_bins=None,
                          ylabel='Proportion of tests positive', title='%s - total'%(parlabel), **kwargs)



def plot_detailed_prevalence(P, lang_ind, **kwargs):
    allpops = getpops(P, 'hum')
    
    aggpars = {'sum': ['all_um', 'all_sm', 'all_ac', 'all_inf', 'all_par', 'all_treat']} #Prevalence total numbers    
    for agg, pars in aggpars.items():
        for par in pars:
            parlabel = P.framework.get_label(par).split('\n')[lang_ind]
            for pop in allpops:
                    standard_plot(P, res_pars=par, data_pars=par, res_pops=pop, data_pops=pop, pop_aggregation=agg,
                                  ylabel=parlabel, title='%s - %s'%(parlabel, pop), **kwargs)
            standard_plot(P, res_pars=par, data_pars=None, res_pops=allpops, data_pops=None, pop_aggregation=agg,
                          ylabel=parlabel, title='%s - total'%(parlabel), **kwargs)

def plot_prevalence(P, lang_ind, **kwargs):
    allpops = getpops(P, 'hum')
    popmapping = sc.odict([(pop, pop) for pop in allpops]+[('Total', allpops)])
    
    aggpars = {'weighted': ['inf_prev', 'ac_prev', 'sym_prev', 'lt_prev']} #Prevalence (proportional)
    for agg, pars in aggpars.items():
        for par in pars:
            parlabel = P.framework.get_label(par).split('\n')[lang_ind]
            for pop in popmapping.keys():
                    standard_plot(P, res_pars=par, data_pars=par, res_pops=popmapping[pop], data_pops=pop, pop_aggregation=agg,
                                  ylabel=parlabel, title='%s - %s'%(parlabel, pop), **kwargs)        

def plot_deaths(P, lang_ind, times={'continuous annualised':None}, **kwargs):
    allpops = getpops(P, 'hum')
    popmapping = sc.odict([(pop, pop) for pop in allpops]+[('Total', allpops)])
    
#    times = {'continuous annualised':None, 'annual total':1.} #'monthly':1./12.
    def time_string(tl): return '' if tl=='annual total' else ' (%s)'%(tl)
    
    pars = ['mal_deaths'] #deaths undiagnosed, diagnosed
    for time_label, t_bins in times.items():
        for par in pars:
#            parlabel = P.framework.get_label(par).split('\n')[lang_ind]
#            if par == 'mal_deaths':
            parlabel = 'Malaria-related deaths'
            for pop in popmapping.keys():
                    standard_plot(P, res_pars=par, data_pars=['mal_deaths', 'mal_deaths_low', 'mal_deaths_high'], res_pops=popmapping[pop], data_pops=pop, t_bins=t_bins,
                                  ylabel='Number of malaria-related deaths', title='%s%s - %s'%(parlabel, time_string(time_label), pop), **kwargs)

    pars = ['mal_deaths_per100K'] #deaths undiagnosed, diagnosed
    for time_label, t_bins in times.items():
        for par in pars:
            parlabel = 'Malaria-related deaths per 100K' #P.framework.get_label(par).split('\n')[lang_ind]
            for pop in popmapping.keys():
                    standard_plot(P, res_pars=par, data_pars=['mal_deaths_per100K', None, None], 
                                  res_pops=popmapping[pop], data_pops=pop, t_bins = t_bins, pop_aggregation= 'weighted', 
                                  ylabel='Number of malaria-related deaths per 100K', title='%s%s - %s'%(parlabel, time_string(time_label), pop), **kwargs)


    pars = ['mal_deaths_treat']
    for time_label, t_bins in times.items():
        for par in pars:
            parlabel = 'Confirmed malaria-related deaths during treatment' #P.framework.get_label(par).split('\n')[lang_ind]
            for pop in popmapping.keys():
                    standard_plot(P, res_pars=par, data_pars='mal_deaths_treat', res_pops=popmapping[pop], data_pops=pop, t_bins=t_bins,
                                  ylabel='Number of malaria-related deaths', title='%s%s - %s'%(parlabel, time_string(time_label), pop), **kwargs)


def plot_treatmentcascade(P, results, results_folder, save_figs, plot_years, **kwargs):
    allpops = getpops(P, 'hum')

    #cascade (current situation)
    fig, table = at.plot_cascade(results=results, cascade='main', show_table=False)
#    if isinstance(figs, tuple): figs = [figs[0]] #it might be (fig, table)?
    fig.set_label('current_status')
    if save_figs:
        at.save_figs(fig, path=results_folder, prefix='cascade_current_status')

    fig, table = at.plot_cascade(results=results, cascade='treatment', show_table=False)
#    if isinstance(figs, tuple): figs = [figs[0]] #it might be (fig, table)?
    fig.set_label('treatment_status')
    if save_figs:
        at.save_figs(fig, path=results_folder, prefix='cascade_treatment_status')

    for result in results:
        year = 2018
        ss_mapping = {'u':'Uncomplicated', 's':'Severe'}
        
        for ss in ['u', 's']:
            d = at.PlotData(result, pops=[{'Total':allpops}], outputs=[
                    {'New active cases':'h_'+ss+'_cases', 'Diagnosed':'casc_'+ss+'_diag', 'Undiagnosed recovery':'casc_'+ss+'_rec', 'Undiagnosed death':'casc_'+ss+'_mort',
                     'Treatment success':'casc_'+ss+'d_succ', 'Treatment incomplete/latent':'casc_'+ss+'d_lat', 'Treatment fail/LTFU':'casc_'+ss+'d_fail', 'Treatment death':'casc_'+ss+'d_mort'
                     }], t_bins=[year, year+1])
            for pop in d.pops:
                d.set_colors('#00267a', pops=pop, outputs=['New active cases', 'Diagnosed', 'Treatment success'])
                d.set_colors('#aaaaaa', pops=pop, outputs=['Undiagnosed recovery', 'Treatment incomplete/latent', 'Treatment fail/LTFU'])
                d.set_colors('#aa2626', pops=pop, outputs=['Undiagnosed death', 'Treatment death'])
            figs = at.plot_bars(d, outer='results', stack_outputs=[['New active cases'],
                                                                   ['Diagnosed', 'Undiagnosed recovery', 'Undiagnosed death'],
                                                                   ['Treatment success', 'Treatment incomplete/latent', 'Treatment fail/LTFU', 'Treatment death']])
            for fig in figs:
                fig.axes[0].set_title(ss_mapping[ss]+' treatment probabilistic outcomes')
                fig.axes[0].set_ylabel('New active malaria cases '+str(year))
                fig.axes[0].set_xticklabels(['Active malaria', 'Diagnosed', 'Treated', 'Success'])
                
                if fig.axes[0].get_legend():
                    entries = sc.odict()
                    entries['Malaria-related death'] = '#aa2626'
                    entries['Undiagnosed recovery or loss to follow-up'] = '#aaaaaa'
                    at.plot_legend(entries)
                    legendsettings = {'loc': 'center left', 'bbox_to_anchor': (1.05, 0.5), 'ncol': 1, 'framealpha':0, 'borderpad':0, }
#                    entries = sc.odict(zip(labels_to_use, colors_to_use))
                    fig = at.plot_legend(entries=entries, plot_type='patch', fig=fig, legendsettings=legendsettings)
    
            if save_figs:
                at.save_figs(figs, path=results_folder, prefix='cascade_probable_'+ss+'_'+result.name+'_')
    
    
def plot_LLINs(P, lang_ind, **kwargs):
    allpops = getpops(P, 'hum')
    pars = ['h_nets_coverage', 'h_nets_available', 'h_nets_sleeps', 'prop_nets_delivered', 'h_incoming_nets', 'h_outgoing_nets'] #deaths undiagnosed, diagnosed
    
    for par in pars:
        parlabel = P.framework.get_label(par).split('\n')[lang_ind]
        ylabel = 'Number of nets' if par in ['h_nets_available'] else 'Proportion'
        for pop in allpops:
                standard_plot(P, res_pars=par, data_pars=par, res_pops=pop, data_pops=pop,
                              ylabel=ylabel, title='%s - %s'%(parlabel, pop), **kwargs)
        agg = 'sum' if par in ['h_nets_available'] else 'weighted'
        standard_plot(P, res_pars=par, data_pars=par, res_pops=allpops, data_pops=pop, pop_aggregation=agg, 
                      ylabel=ylabel, title='%s - total'%(parlabel), **kwargs)
    
    mospops = getpops(P, 'mos')
    pars = ['m_impacted_nets'] #deaths undiagnosed, diagnosed
    for par in pars:
        parlabel = P.framework.get_label(par).split('\n')[lang_ind]
        standard_plot(P, res_pars=par, data_pars=None, res_pops=mospops, data_pops=None,
                      ylabel=parlabel, title='%s - total'%(parlabel), **kwargs)
    
    
def plot_IRS(P, lang_ind, **kwargs):
    allpops = getpops(P, 'hum')
    pars = ['h_irs_covered'] #deaths undiagnosed, diagnosed
    
    for par in pars:
        parlabel = 'IRS coverage' #P.framework.get_label(par).split('\n')[lang_ind]
#        for pop in allpops:
#                standard_plot(P, res_pars=par, data_pars=par, res_pops=pop, data_pops=pop,
#                              ylabel=parlabel, title='%s - %s'%(parlabel, pop), **kwargs)
        agg = 'weighted'
        standard_plot(P, res_pars=par, data_pars=par, res_pops=allpops, data_pops=[None, 'National - MAP', 'National - WHO'], pop_aggregation=agg, 
                      ylabel='Proportion', title='%s - total'%(parlabel), **kwargs)
    
    mospops = getpops(P, 'mos')
    pars = ['m_impacted_irs'] #deaths undiagnosed, diagnosed
    for par in pars:
        parlabel = P.framework.get_label(par).split('\n')[lang_ind]
        standard_plot(P, res_pars=par, data_pars=None, res_pops=mospops, data_pops=None,
                      ylabel=parlabel, title='%s - total'%(parlabel), **kwargs)
    


def plot_deaths_advanced(P, results, results_folder, save_figs, plot_years, **kwargs):
    allpops = getpops(P, 'hum')
#    blhpops  = natpops(P)

    times = {'continuous annualised':None, 'five-year total':5.} #'monthly':1./12.,
    for time_label, t_bins in times.items():
        #deaths by source
        d = at.PlotData(results, pops=[{'Total':allpops}], outputs=[
                {'Untreated uncomplicated malaria': ['uPdeath:flow'],
                 'Untreated severe malaria': ['sPdeath:flow'], 
                 'Uncomplicated malaria during treatment': ['uPtdeath:flow'],
                 'Severe malaria during treatment': ['sPtdeath:flow']
                 }], pop_aggregation='sum', t_bins=t_bins)
#        if t_bins is not None:
#            for series in d.series:
#                series.tvec -= t_bins
        for pop in d.pops:
            d.set_colors('Greys', pops=pop, outputs=d.outputs)
        if t_bins is None:
            figs = at.plot_series(d, plot_type='stacked', axis='outputs',data=None)
            for fig in figs:
        #        fig.axes[0].legend(['Model mortality'], **legendsettings)
                fig.axes[0].set_title(fig.axes[0].get_title().replace('parset_default-','Malaria-related deaths by source - '))
                fig.axes[0].set_xlim(plot_years[0], plot_years[1])
        else:
            figs = at.plot_bars(d, outer='results', stack_outputs='all', stack_pops='all')
            figs[0].axes[0].set_title('Malaria-related deaths by source')
        if save_figs: at.save_figs(figs, path=results_folder, prefix='mort_by_source_%s_'%(time_label))
    
        #deaths by populaion
        d = at.PlotData(results, pops=allpops, outputs=[
                {'Model deaths': ['mal_deaths']
                 }], pop_aggregation='sum', t_bins=t_bins)
#        if t_bins is not None:
#            for series in d.series:
#                series.tvec -= t_bins
        for out in d.outputs:
            d.set_colors('Reds', pops=d.pops, outputs=out)
        if t_bins is None:
            figs = at.plot_series(d, plot_type='stacked', axis='pops',data=None)
            for fig in figs:
        #        fig.axes[0].legend(['Model mortality'], **legendsettings)
                fig.axes[0].set_title(fig.axes[0].get_title().replace('parset_default','Malaria-related deaths by population'))
                fig.axes[0].set_xlim(plot_years[0], plot_years[1])
        else:
            figs = at.plot_bars(d, outer='results', stack_outputs='all', stack_pops='all')
            figs[0].axes[0].set_title('Malaria-related deaths by population')   
        if save_figs: at.save_figs(figs, path=results_folder, prefix='mort_by_pops_%s_'%(time_label))


def plot_diagnoses_advanced(P, results, results_folder, save_figs, plot_years, **kwargs):
    allpops = getpops(P, 'hum')
    
    for pop in allpops:
        d = at.PlotData(results, pops=pop, outputs={'False positive (susceptible with MLS)':['S_pos:flow'],
                                                        'Latent diagnosis (due to MLS)':['L_pos:flow'],
                                                        'Asymptomatic diagnosis (due to MLS)':['aP_pos:flow'],
                                                        'Uncomplicated malaria true positive':['uP_pos:flow'],
                                                        'Severe malaria true positive':['sP_pos:flow'],
                                                        }, t_bins=5)
        figs = at.plot_bars(d, outer='results', stack_outputs='all') # This should look like the usual Optima-TB result
        figs[0].axes[0].set_title('Diagnoses %s'%(pop))
        if save_figs: at.save_figs(figs, path=results_folder, prefix='Diagnoses %s_'%(pop))

    


def plot_reconciliation(P, results, progset_name, results_folder, calib_name = 'Calibrated', recon_name = 'Current spending', lang_ind = 'English', **kwargs):
    
    resnames = [res.name for res in results]
    if calib_name in resnames and recon_name in resnames:
        calibrated = [result for result in results if result.name == calib_name][0]
        reconciled = [result for result in results if result.name == recon_name][0]
        
        prog_headers = []
        for i in range(4):
            prog_headers.append('Program %s'%i)
            prog_headers.append('Program impact %s'%i)
            prog_headers.append('Program coverage %s'%i)
        
        recon_suggestions = [['WARNING: reconciliation suggestions do not take into account interactions with multiple programs so will need further review'],
                             ['WARNING: baseline recommendations will likely be inconsistent for different programs where multiple programs are defined'],
                             ['WARNING: the best solution is probably a combination of all options'],
                             ['For a parameter where "low is bad" like treatment probability, increasing the baseline or decreasing program effects will shrink the range and make the program less impactful'],
                             ['For a parameter where "low is bad" like treatment probability, decreasing the baseline or increasing program effects will increase the range and make the program more impactful'],
                             ['WARNING: any changes made as a result of using the suggestions will influence optimization results, so take care!'],
                             ['Parameter (full)', 'Parameter (short)', 'Population', 'Year', 'Calibrated value', 'Reconciled value', 'Mismatch', 'Suggestion (parset)', 'Suggestion (program impacts)', 'Suggestion (program coverage)', 'Suggestion (baseline)', 'Baseline']+prog_headers
                            ]

        for covout in P.progsets[progset_name].covouts.values():
            par = covout.par
            label = P.framework.get_label(par)
            pop = covout.pop
            
            calib_par = calibrated.get_variable(par, pop)[0]
            recon_par = reconciled.get_variable(par, pop)[0]
            
            assert((calib_par.t == recon_par.t).all()) #make sure they're matched for the times...
            
            #get the first year there's a mismatch between parameter and progset
            try:
                dif_ind  = list(calib_par.vals == recon_par.vals).index(False)
                dif_year = calib_par.t[dif_ind]
                
                calib_val = calib_par.vals[dif_ind]
                recon_val = recon_par.vals[dif_ind]
                diff = calib_val - recon_val
                
                if abs(diff/calib_val) > 0.5 and abs(diff) >0.01:
                    mismatch = 'SEVERE'
                elif abs(diff/calib_val) > 0.1 and abs(diff) >0.0001:
                    mismatch = 'MODERATE'
                elif abs(diff/calib_val) > 0.01 and abs(diff) >0.00001:
                    mismatch = 'LOW'
                else:
                    mismatch = ''
                
                suggest_pars     = ''
                suggest_progs    = ''
                suggest_covs     = ''
                suggest_baseline = ''

                covout = P.progsets[progset_name].covouts[(par, pop)]
                baseline = covout.baseline
                
                if mismatch != '': #make a suggestion on the parset
                    verb = 'Increase' if diff<0 else 'Decrease'
                    suggest_pars+= '%s the most recent databook value for %s to be %s\n'%(verb, label, recon_val)

                prog_covs = []
                for prog_name in covout.progs.keys():
                    
                    if P.framework.get_variable(par)[0].format=='number':
                        spending = reconciled.get_alloc()[prog_name][dif_ind]
                        prog_coverage = P.progsets[progset_name].programs[prog_name].get_capacity(tvec=dif_year, spending=spending, dt=P.settings.sim_dt)[0]
                    else:
                        prog_coverage = reconciled.get_coverage()[prog_name][dif_ind]
                    prog_effect   = covout.progs[prog_name]
                    
                    prog_covs.append(prog_name)
                    prog_covs.append(prog_effect)
                    prog_covs.append(prog_coverage)
                    
                    if mismatch != '':                        
                        if prog_coverage > 0: #if the program coverage is zero, changes to the program impact or baseline won't help at all, dont' suggest them
                            needed_new_effect = baseline + (calib_val - baseline)/prog_coverage
                            verb = 'Increase' if needed_new_effect>prog_effect else 'Decrease'
                            suggest_progs+= '%s effect for %s to be %s\n'%(verb, prog_name, abs(needed_new_effect))

                            needed_new_baseline = (calib_val - prog_coverage*prog_effect)/(1. - prog_coverage)
                            verb = 'Increase' if needed_new_baseline>baseline else 'Decrease'
                            suggest_baseline+= '%s baseline for %s to be %s\n'%(verb, prog_name, needed_new_baseline)
                        
                        needed_new_coverage = (calib_val - baseline)/(prog_effect - baseline)
                        verb = 'Increase' if needed_new_coverage>prog_coverage else 'Decrease'
                        suggest_covs+= '%s coverage for %s to be %s\n'%(verb, prog_name, needed_new_coverage)

            except:
                continue
                
            
            recon_suggestions.append([label, par, pop, dif_year, calib_val, recon_val, mismatch, suggest_pars, suggest_progs, suggest_covs, suggest_baseline, baseline]+prog_covs)
    
            sc.savespreadsheet(filename=results_folder + 'reconciliation_guide.xlsx', data=recon_suggestions)
    
    else:
        print('Results list does not contain clear parset ("%s") and progset ("%s") results for reconciliation or they had mismatched times, skipping output of suggestions.'%(calib_name, recon_name))
        
    
    for covout in P.progsets[progset_name].covouts.values():
        par = covout.par
        label = P.framework.get_label(par).split('\n')[lang_ind]
        pop = covout.pop
        standard_plot(P=P, results=results, results_folder=results_folder, res_pars=par, data_pars=None, res_pops=pop, data_pops=None,
                      ylabel=label, title='%s (%s) - %s'%(label, par, pop), **kwargs)


def plot_budgets(P, results, results_folder, save_figs, plot_years, key_year, unfunded_progs=[], time_adjustments=None, currency='$', **kwargs):
    results = [res for res in results if not res.get_equivalent_alloc() is None]
    if results == []:
        print('Tried to plot budgets but none of the results have a budget!')
        return False
    funded_progs = [prog for prog in results[0].get_equivalent_alloc().keys() if not prog in unfunded_progs]
    
#    if time_adjustments is not None:
#        num_years = np.lcm.reduce([int(np.ceil(ta['timespan'])) for ta in time_adjustments.values()])
#    else:
#        num_years = None


    #Work out an appropriate duration and number of years to plot
    durations = [x['timespan'] for x in time_adjustments.values()]
    lcm_duration = np.product(durations) #TODO actual LCM calculation
    average_years = [key_year, key_year + lcm_duration]
    indicative_years = list(np.arange(key_year, key_year + lcm_duration + 2, 1.0))
    future_years  = [2035 - lcm_duration, 2035]
    time_comparisons = {'next_funding_cycle': average_years, 'funding_cycle': indicative_years, 'future_funding_cycle': future_years}
            
    for tcn in time_comparisons.keys():
        #spending comparison of CUMULATIVE actual allocated budget - this will be incorrect for coverage scenarios
        budget_results = [res for res in results if P.scens[res.name].get_instructions(progset=res.model.progset, project=P) != sc.odict()]
        d = at.PlotData.programs(budget_results, outputs=funded_progs, quantity='spending', t_bins=time_comparisons[tcn])
    #    d.interpolate(key_year)
        figs = at.plot_bars(d, stack_outputs='all')
        figs[0].axes[0].set_ylabel(figs[0].axes[0].get_ylabel().replace('$', currency))
        figs[0].axes[0].set_title(figs[0].axes[0].get_title()+'Budgeted')
        if save_figs: 
            at.save_figs(figs, path=results_folder, prefix='spending_comparison_%s'%(tcn))
            
        #spending comparison of CUMULATIVE actual allocated budget - this will be incorrect for coverage scenarios
        budget_results = [res for res in results if P.scens[res.name].get_instructions(progset=res.model.progset, project=P) != sc.odict()]
        d = at.PlotData.programs(budget_results, outputs=funded_progs, quantity='equivalent_spending', t_bins=time_comparisons[tcn])
    #    d.interpolate(key_year)
        figs = at.plot_bars(d, stack_outputs='all')
        figs[0].axes[0].set_ylabel(figs[0].axes[0].get_ylabel().replace('$', currency))
        figs[0].axes[0].set_title(figs[0].axes[0].get_title()+'Budgeted')
        if save_figs: 
            at.save_figs(figs, path=results_folder, prefix='spending_comparison_equivalent_%s'%(tcn))
        


#    #if any of the results have time-varying budgets, do some additional plots
    if not np.array([np.array([allequal(prog_budget) for prog_budget in res.get_equivalent_alloc()[:]]).all() for res in results]).all():
        d = at.PlotData.programs(results, outputs=funded_progs, quantity='equivalent_spending')
        figs = at.plot_series(d, plot_type='stacked', axis='outputs')
        for fig in figs:
            fig.axes[0].set_xlim(plot_years[0], plot_years[1])
        if save_figs: 
            at.save_figs(figs, path=results_folder, prefix='spending_over_time_')
            
        figs = at.plot_series(d, axis='results')
        for fig in figs:
            fig.axes[0].set_xlim(plot_years[0], plot_years[1])
            fig.axes[0].set_ylabel(fig.axes[0].get_ylabel().replace('$', currency))
        if save_figs: 
            at.save_figs(figs, path=results_folder, prefix='spending_over_time_comparison_')
