# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 14:51:36 2019

@author: rowan.martin-hughes
"""
from malaria_utils import *

plot_categories = sc.odict()
plot_categories['demographics'] = ['alive']
plot_categories['keypars'] = ['mal_deaths'] #['Human incident cases', 'Incident cases (per person)', 'Human incidence (diagnosed cases per person) (monthly)',
                               #'Malaria-related deaths - best (annual total)', 'Malaria-related deaths (confirmed during treatment) (continuous annualised)',
                               #]


def _get_img_path(folder, keys, unkeys=[], case_sensitive=False, verbose=False, folder_depth=1):
    """
    Helper function to get paths
    :param results_folder: top level folder that all results should be inside
    :param keys: anything returned must contain each of the keys in the list
    :param unkeys: anything returned must NOT contain any of the unkeys in the list
    """
    keys = sc.promotetolist(keys)
    possibles = get_paths(folder, inclusions=keys, extensions=['.png'], exclusions=unkeys,
             version='all', verbose=verbose, case_sensitive=case_sensitive, folder_depth=folder_depth)
    if possibles is None:
        if verbose: print('WARNING: Could not find any files matching %s in %s.'%(keys, folder))
        return possibles
    if len(possibles)>1 and not 'legend' in keys: #exempt legend from printing as normally we want one of any legends from a folder!
        if verbose: print('WARNING: %s files in %s matching all of %s exist, returning the last result: %s'%(len(possibles), folder, keys, possibles[-1]))

    return possibles[-1]
            
def _get_slide_set(results_folder, including, distinctions, not_including=[], total_slide = True,
                  slide_title = None, title_prefix='', **kwargs):
    """Usage of this function is to generate 1-3 slides covering a specific parameter, where the plot outputs for that parameter would be specified by:
    :param including: all files must include these strings, e.g. '['incidence', 'calibration']
    :param not_including: all files must NOT include these strings
    together, including and not_including should provide enough information to precisely identify the desired output
    :param distinctions: the comparative slide will search for files matching each of these individually to put on a slide together, e.g. populations ['0-4', '5-14', 'PLHIV']
    """
    slides = []
    viablepaths = []    
    
    if slide_title is None:
        slide_title = including[0]
    slide_title = title_prefix + slide_title
    #put the total on a slide by itself
    this_path = _get_img_path(results_folder, including + ['total'], not_including + ['legend'], verbose=False)

    
    if total_slide and this_path:
        slides.append({'style': 'One',
                   'title': slide_title,
                   'img1': this_path,
                   'legend': _get_img_path(results_folder, including + ['total', 'legend'], not_including, verbose=False),
                   })
        if kwargs.get('results') and kwargs.get('parname'):
            slides[-1]['para1'] = _get_impact_string(**kwargs, as_para=True)
    

    for distinct in distinctions:
        this_path = _get_img_path(results_folder, including + [distinct], not_including + ['legend'], verbose=False)
        if this_path:
            viablepaths+=[this_path]

    mapping = {1: 'One', 2: 'Two', 3: 'Three', 4: 'Four', 5: 'Five', 6: 'Six'} #named slides in the template
    numimgs = len(viablepaths)
    if numimgs>0:
        slides.append({'style': mapping[min(numimgs, 6)],
                   'title': slide_title,
                   'legend': _get_img_path(results_folder, including + ['legend'], not_including, verbose=False)})
        for pn, path in enumerate(viablepaths):
            slides[-1]['img%s'%(pn)] = path
        if kwargs.get('results') and kwargs.get('parname') and kwargs.get('P'):
            impact_para = []
            for pop in [key for key, details in kwargs['P'].data.pops.items() if details['type']=='ind']:
                if np.array([pop in vp.split(sep)[-1] for vp in viablepaths]).any():
                    impact_para += [(0, _get_popstring(kwargs['P'], [pop]))]
                    impact_para += _get_impact_string(**kwargs, pops=[pop], as_para=True)
            if len(impact_para)<=5: #might fit onto the slide!
                slides[-1]['para1'] = impact_para
            else:
                slides.append({'style': 'Text',
                   'title': slide_title,
                   'para1': impact_para})
                
        
    return slides



def _get_category_slides(f, P, results, allpops, differentiation, variables_to_plot, uncertainty='uncertainty', title_prefix='', **kwargs):
    slides = []
    for par in variables_to_plot:
        next_slides = _get_slide_set(f, [P.framework.get_label(par), differentiation, uncertainty], allpops, total_slide = True)
        if next_slides:
            if kwargs.get('results'):
                slides[-1]['para1'] = _get_impact_string(**kwargs, parname = par, as_para=True)
            slides+=next_slides

    return slides

def _get_keypar_slides(f, allpops, differentiation, uncertainty='uncertainty', title_prefix='', **kwargs):
    slides = []
#    prop = propconf[key]
    
    slides+= _get_slide_set(f, ['Human incident cases', 'annual total', differentiation, uncertainty], allpops,
                            total_slide = False, title_prefix=title_prefix)
    
    slides+= _get_slide_set(f, ['Human incidence', 'cases per person', 'annual total', differentiation, uncertainty],
                            allpops, total_slide = False, title_prefix=title_prefix)

    slides+= _get_slide_set(f, ['Human incidence', 'monthly', differentiation, uncertainty], allpops,
                            total_slide = True, title_prefix=title_prefix)
    
    slides+= _get_slide_set(f, ['Malaria-related deaths', differentiation, uncertainty], distinctions = allpops,
                            not_including=['confirmed'], total_slide = True, title_prefix=title_prefix)
    
    slides+= _get_slide_set(f, ['Malaria-related deaths', 'confirmed', differentiation, uncertainty], distinctions = allpops,
                            total_slide = True, title_prefix=title_prefix)
        
    return slides

def _get_programs_slides(f, allpops, differentiation, uncertainty='uncertainty', title_prefix='', **kwargs):
    slides = []
#    prop = propconf[key]
    
    slides+= _get_slide_set(f, ['LLIN coverage', differentiation, uncertainty], allpops, total_slide = False, title_prefix=title_prefix)
    
    slides+= _get_slide_set(f, ['Number of LLINs available', differentiation, uncertainty], allpops, total_slide = False, title_prefix=title_prefix)

    slides+= _get_slide_set(f, ['Proportion of the population covered by IRS', differentiation, uncertainty], allpops, total_slide = True, title_prefix=title_prefix)

    return slides

def _get_demographic_slides(f, allpops, differentiation, uncertainty='uncertainty', title_prefix='', **kwargs):

    slides = []
    
    next_slides= _get_slide_set(f, ['demographics', differentiation, uncertainty], allpops, total_slide = False, title_prefix=title_prefix)
    if next_slides:
        if uncertainty:
            next_slides[-1]['text2'] = 'Model run without uncertainty on demographic projections (birth and mortality rates)'
        slides += next_slides

#    slides+= _get_keypar_slides(f, allpops, differentiation=differentiation)
#    
#    slides+= _get_slide_set(f, ['Proportion of undiagnosed active TB cases that are smear positive', differentiation, uncertainty], allpops, total_slide = True)
#
#    slides+= _get_slide_set(f, ['Case fatality ratio', differentiation, uncertainty], allpops, total_slide = True)
#    
#    slides+= _get_slide_set(f, ['Pulmonary TB case notification rate', differentiation, uncertainty], allpops, total_slide = True)

    return slides

def _get_cascade_slides(f, allpops, differentiation, uncertainty='uncertainty', title_prefix='', key_year=2019, **kwargs):
    slides = []
    differentiation = sc.promotetolist(differentiation)
    
    next_slides =  _get_slide_set(f, including=['cascade_']+differentiation, distinctions=['current_status'], not_including=['probable'], total_slide = False, slide_title='Current situation TB care cascade', title_prefix=title_prefix)
    if next_slides: next_slides[-1]['para1'] = [(0, 'This cascade view gives a “snapshot” picture, representing the status of each person in the model with active malaria at a given point in time, in this case 1 January %i. This is useful to get a picture of where new interventions could be targeted and also to understand the impact that delayed diagnosis and delayed treatment initiation has on the epidemic.'%(key_year)),
                           (0, 'The loss from active malaria to diagnosed malaria is largely due to the time it takes people to seek medical care.'),
                          ]
    slides += next_slides

    next_slides = _get_slide_set(f, ['cascade_probable']+differentiation, ['_u_bars', '_s_bars'], total_slide = False, slide_title='Probability-based malaria care cascade', title_prefix=title_prefix)
    if next_slides: slides[-1]['para1'] = [(0, 'This alternative cascade examines the probability of final outcomes from each stage in the cascade, e.g. for people with active malaria, how many of those will be diagnosed prior to natural recovery or death. This cascade is based on modeled incidence of malaria for %i and rates of flow through the stages of the cascade.'%(key_year)),
                          ]
    slides += next_slides
    
#    slides+= _get_slide_set(f, ['cascade_probable']+differentiation, ['_u', '_s'], total_slide = False, slide_title='Smear negative TB care cascade', title_prefix=title_prefix)
   
#    slides+= _get_slide_set(f, ['cascade_probable']+differentiation, ['_pd', '_pm', '_px'], total_slide = False, slide_title='Smear positive TB care cascade', title_prefix=title_prefix)

    return slides

def _get_scenario_description_slides(P, scen_name, plot_set, progset_name, allpops):
    slides = []
    relevant_scens = [P.scens[resname] for resname in plot_set['resnames'] if resname in P.scens] #actual scenario objects!
    data = None
    
    #%%Parameter scenarios
    scens = [scen for scen in relevant_scens if 
             (isinstance(scen, at.ParameterScenario) or isinstance(scen, at.CombinedScenario))
             and scen.scenario_values is not None]   
    if scens != []: #at least one parameter scenario        
        all_years = [x[n]['t'] for scen in scens for x in scen.scenario_values.values() for n in x]
#        print ('Scenario years pre-flattening: ',  all_years)
        years = sorted(set([item for pn in all_years for item in pn]))
        params_needed = sorted(set([item for scen in scens for item in scen.scenario_values.keys()]))
        #TODO remove any params where the value is the same for each population and each year under all scens
        
        data = [['Parameter', 'Populations']+['%s\n%i'%(scen.name, year) for year in years for scen in scens]]
        for param in params_needed:
            all_applied_pops = set([item for scen in scens for item in (scen.scenario_values[param].keys() if param in scen.scenario_values else [])])
            pops = [pop for pop in allpops if pop in all_applied_pops] #make sure the pops are in the same order as they were defined in the project if that matters
            applies_pops = _get_popstring(P, pops)
            
            param_str, invert = _get_inversion(P, param)
            newrow = [param_str, applies_pops]
            keeprow = False
            for year in years:
                for sn, scen in enumerate(scens):
                    effects = []
                    for pop in pops:
                        try:
                            yr_ind = scen.scenario_values[param][pop]['t'].index(year)
                            effects.append(scen.scenario_values[param][pop]['y'][yr_ind])
                        except:
                            effects.append('')
                    effect_str =  _get_summary_effect_str(P, pops, effects, invert)
                    if sn>0 and effect_str != newrow[-1]: #if there is at least one year where different scenarios have different values
                        keeprow = True 
                    newrow.append(effect_str)
            if keeprow:
                data.append(newrow)
                    
        slides.append({'style': 'Table', 'title': scen_name, 'table1': data})
        
    #%%Budget scenarios
    scens = [scen for scen in relevant_scens if 
             (isinstance(scen, at.BudgetScenario) or isinstance(scen, at.CombinedScenario))
             and scen.get_instructions(progset=P.progsets[progset_name], project=P).alloc != sc.odict()]
    if scens != []: #at least one budget scenario
        all_years = [x.t for scen in scens for x in scen.get_instructions(progset=P.progsets[progset_name], project=P).alloc.values()] #list of lists
        years = sorted(set([item for pn in all_years for item in pn])) #flattened
        if len(years)>2:
            years = [years[0], years[-1]]
        programs_needed = sorted(set([item for scen in scens for item in scen.get_instructions(progset=P.progsets[progset_name], project=P).alloc.keys()]))
        
        data = [['Program']+['%s\n%i'%(scen.name, year) for year in years for scen in scens]]
        for prog in programs_needed:
                                    
            prog_str = P.progsets[progset_name].programs[prog].label
            newrow = [prog_str]
            
            keeprow = False
            for year in years:
                for sn, scen in enumerate(scens):
                    try:
                        yr_ind = scen.get_instructions(progset=P.progsets[progset_name], project=P).alloc[prog].t.index(year)
                        effect_str = sigfigs(scen.get_instructions(progset=P.progsets[progset_name], project=P).alloc[prog].vals[yr_ind])
                    except:
                        effect_str = ''
                    if sn>0 and effect_str != newrow[-1]: #if there is at least one year where different scenarios have different values
                        keeprow = True 
                    newrow.append(effect_str)
            if keeprow:
                data.append(newrow)
#        print (data)
                    
        slides.append({'style': 'Table', 'title': scen_name, 'table1': data})
        
    #%%Coverage scenarios
    scens = [scen for scen in relevant_scens if
             (isinstance(scen, at.CoverageScenario) or isinstance(scen, at.CombinedScenario))
             and scen.get_instructions(progset=P.progsets[progset_name], project=P).coverage != sc.odict()]
    if scens != []: #at least one coverage scenario    
        all_years = [x.t for scen in scens for x in scen.get_instructions(progset=P.progsets[progset_name], project=P).coverage.values()] #list of lists
        years = sorted(set([item for pn in all_years for item in pn])) #flattened
        if len(years)>2:
            years = [years[0], years[-1]]
        programs_needed = sorted(set([item for scen in scens for item in scen.get_instructions(progset=P.progsets[progset_name], project=P).coverage.keys()]))
        
        data = [['Program']+['%s\n%i'%(scen.name, year) for year in years for scen in scens]]
        for prog in programs_needed:
                                    
            prog_str = P.progsets[progset_name].programs[prog].label
            newrow = [prog_str]
            
            keeprow = False
            for year in years:
                for sn, scen in enumerate(scens):
                    try:
                        yr_ind = scen.get_instructions(progset=P.progsets[progset_name], project=P).coverage[prog].t.index(year)
                        effect_str = sigfigs(scen.get_instructions(progset=P.progsets[progset_name], project=P).coverage[prog].vals[yr_ind])
                    except:
                        effect_str = ''
                    if sn>0 and effect_str != newrow[-1]: #if there is at least one year where different scenarios have different values
                        keeprow = True 
                    newrow.append(effect_str)
            if keeprow:
                data.append(newrow)
                    
        slides.append({'style': 'Table', 'title': scen_name, 'table1': data})
        

    return slides

def _get_popstring(P, pops, case='capitalize'):
    allpops = getpops(P)
    if pops==allpops:
        applies_pops = 'All populations'
    else:
        pop_labels = [P.data.pops[pop]['label'] for pop in pops]
        applies_pops  = str.join(', ', pop_labels)
    if case=='lower':
        applies_pops = applies_pops[0].lower()+applies_pops[1:]
    return applies_pops

def _get_compstring(P, comps, case='capitalize'):
    """Specific compartment categories can be defined"""
    comp_labels = [P.framework.get_label(comp) for comp in comps]
    
    allcats = ['Susceptible',  'Susceptible, malaria-like symptoms',  'Exposed',  'Latent',  'Latent, malaria-like symptoms',
               'Uncomplicated malaria (symptomatic)', 'Severe malaria',  'Asymptomatic malaria / natural resistance',
               'Asymptomatic malaria / natural resistance with malaria-like symptoms', 'Susceptible eligible for treatment',
               'Latent eligible for treatment', 'Uncomplicated malaria eligible for treatment', 'Severe malaria eligible for treatment',
               'Asymptomatic malaria eligible for treatment', 'Post-treatment recovery']
    mls_cats = ['Susceptible, malaria-like symptoms', 'Latent, malaria-like symptoms', 'Uncomplicated malaria (symptomatic)',
                'Severe malaria', 'Asymptomatic malaria / natural resistance with malaria-like symptoms', ]
    et_cats = ['Susceptible eligible for treatment', 'Latent eligible for treatment',
               'Uncomplicated malaria eligible for treatment', 'Severe malaria eligible for treatment',
               'Asymptomatic malaria eligible for treatment']
    ns_cats = ['Susceptible', 'Exposed',  'Latent',  'Asymptomatic malaria / natural resistance', 'Post-treatment recovery']

    if set(comp_labels)==set(allcats):
        applies_comps = 'All compartments'
    elif set(comp_labels)==set(mls_cats):
        applies_comps = 'All with malaria-like symptoms'
    elif set(comp_labels)==set(et_cats):
        applies_comps = 'All who are diagnosed and eligible for treatment'
    elif set(comp_labels)==set(ns_cats):
        applies_comps = 'All who have no current symptoms of malaria'
    else:
        applies_comps = str.join(', ', comp_labels)

    if case=='lower':
        applies_comps = applies_comps[0].lower()+applies_comps[1:]
    return applies_comps

def  _get_summary_effect_str(P, pops, effects, invert=False, case='capitalize'):
    if effects == []:
        return ''
    assert len(pops)==len(effects), '%s needs to be the same length as %s, each population must correspond to one effect'%(pops, effects)
    effects = np.array(sc.dcp(effects))
    mins = np.where(effects==min(effects))
    maxes = np.where(effects==max(effects))
    min_pops = list(np.array(pops)[mins])
    max_pops = list(np.array(pops)[maxes])
    str_min_pops = _get_popstring(P, min_pops) #, case='lower')
    str_max_pops = _get_popstring(P, max_pops) #, case='lower')
    if invert:
        effects = ['Never' if x==0. else x if isinstance(x, str) else 1./x for x in effects]
        temp = sc.dcp(str_max_pops)
        str_max_pops = sc.dcp(str_min_pops)
        str_min_pops = temp
    if min(effects) == max(effects):
        effect_string = sigfigs(effects[0])
    else:
        effect_string = '%srom %s (%s) to %s (%s)'%('F' if case=='capitalize' else 'f', sigfigs(min(effects)), str_min_pops, sigfigs(max(effects)), str_max_pops)

    return effect_string

def _get_inversion(P, par):
    par_name = P.framework.get_label(par)
    
    return par_name, False

def _get_plotdata(P, results, pops, output, tvals=[2018, 2035]):
    """ Copied and modified from _output_to_df in results.py"""
    pop_labels = {x: y for x, y in zip(results[0].pop_names, results[0].pop_labels) if x in pops}
    data = dict()
    output_name = P.framework.get_label(output)
    tvals = np.array(tvals)

    popdata = at.PlotData(results, pops=pops, outputs=output)
    assert len(popdata.outputs) == 1, 'Framework plot specification should evaluate to exactly one output series - there were %d' % (len(popdata.outputs))
    popdata.interpolate(tvals)
    for result in popdata.results:
        for pop_name in popdata.pops:
            data[(output_name, popdata.results[result], pop_labels[pop_name])] = popdata[result, pop_name, popdata.outputs[0]].vals

    # Now do a population total. Need to check the units after any aggregations
    # Check results[0].model.pops[0].comps[0].units just in case someone changes it later on
    if popdata.series[0].units in {at.FS.QUANTITY_TYPE_NUMBER, results[0].model.pops[0].comps[0].units}:
        # Number units, can use summation
        popdata = at.PlotData(results, outputs=output, pops={'total': pops}, pop_aggregation='sum')
        popdata.interpolate(tvals)
        for result in popdata.results:
            data[(output_name, popdata.results[result], 'Total (sum)')] = popdata[result, popdata.pops[0], popdata.outputs[0]].vals
    elif popdata.series[0].units in {at.FS.QUANTITY_TYPE_FRACTION, at.FS.QUANTITY_TYPE_PROPORTION,
                       at.FS.QUANTITY_TYPE_PROBABILITY} or 'prev' in output: #warning: adding 'prev' here is a bit dangerous as an assumption, but useful for characteristics
        popdata = at.PlotData(results, outputs=output, pops={'total': pops}, pop_aggregation='weighted')
        popdata.interpolate(tvals)
        for result in popdata.results:
            data[(output_name, popdata.results[result], 'Total (weighted average)')] = popdata[result, popdata.pops[0], popdata.outputs[0]].vals
    else:
        for result in popdata.results:
            data[(output_name, popdata.results[result], 'Total (unknown units)')] = np.full(tvals.shape, np.nan)
    
    return popdata

def _get_impact_string(P, results, parname, tvals=[2020, 2035], pops = None, as_para = True):
    """Return a sentence describing the difference in impact on a given parameter between multiple results"""
    if pops is None:
        pops = getpops(P)
    plotdata = _get_plotdata(P, results=results, pops=pops, output=parname, tvals=tvals)
    par_str = P.framework.get_label(parname)
    res_names = [res.name for res in results]
    plausible_baseline_resnames = [rn for rn in res_names if ('current' in rn.lower() or 'baseline' in rn.lower() or 'calibrated' in rn.lower())]
    if len(plausible_baseline_resnames)>0:
        base_resname = plausible_baseline_resnames[0]
        if len(plausible_baseline_resnames)>1:
            print('WARNING: could not determine the baseline scenario among viable results %s (assuming %s).'%(plausible_baseline_resnames, base_resname))
    else:
        base_resname = res_names[0]
        print('WARNING: could not determine the baseline scenario among results %s as none contain \'current\', \'baseline\', or \'calibrates\' (assuming the first result is the baseline: %s).'%(res_names, base_resname))
    base_ind = res_names.index(base_resname)
    
    impact_str = [(0, '%s:'%(par_str))]
    #plotdata[result, pop, output] to get something out!
    for yr, year in enumerate(tvals):
        res_outcomes = [plotdata[res_name, 'total', parname].vals[yr] for res_name in res_names]
                
        minval  = min(res_outcomes)
        maxval  = max(res_outcomes)
        baseval = res_outcomes[base_ind]
        
        impact_str.append((1, 'In %i with %s: %s'%(year, base_resname, sigfigs(baseval))))
        if minval<baseval:
            min_resname = str.join(', ', [res_name for rn, res_name in enumerate(res_names) if res_outcomes[rn]==minval])
            diff_str = '%s%%'%(sigfigs(100*(baseval-minval)/baseval))
            impact_str.append((2, '%s: %s, a reduction of %s.'%(min_resname, sigfigs(minval), diff_str)))
        if maxval>baseval:
            max_resname = str.join(', ', [res_name for rn, res_name in enumerate(res_names) if res_outcomes[rn]==maxval])
            diff_str = '%s%%'%(sigfigs(100*(maxval-baseval)/baseval))
            impact_str.append((2, '%s: %s, an increase of %s.'%(max_resname, sigfigs(maxval), diff_str)))
    
    if as_para:
        return impact_str
    else:
        return str.join('\n', ['  '*ims[0] + ims[1] for ims in impact_str])


#%%SLIDES
def make_pptx(filename=None, template_path=None, results_folder=None, sections=None, P=None, results=None,
              allpops=None, currency='$', progset_name='default', key_year = 2018, uncertainty = True,
              confidence_interval = 'ci', num_samples = 100, plot_sets={}, project_name='COUNTRY', verbose=True,):
#    key = pres_pars['book_key'].lower()
#    keydisp = keyarray[pres_pars['book_key']]
    if results_folder is None:  results_folder = ''
    if filename is None: filename = 'presentation.pptx'
    if template_path is None:
        import os
        template_path = os.path.dirname(os.path.realpath(__file__)) + sep + 'PPT template.pptx'
    file_path   = results_folder + filename 
    f = results_folder
    
    if verbose: print('Creating a slide set')
    
    if allpops is None: allpops = getpops(P)
    #NOTE: may wish to manually specify the order of the pops to get a preferred arrangement on slides 
    
    if sections is None:
        sections = ['Introduction',]
#                    'Scope of work',
#                    'Methodology',]
        sections+=list(plot_sets.keys())   #SCENARIOS ARE ALL AUTOFILLED HERE
        if progset_name: sections.append('Review of programs and costing')
        sections.append('Acknowledgements')
        
    
    slides = []
    slides.append({'style': 'Text',
                   'title': 'Investment analysis\n\n%s\n\nDraft %s'%(project_name, datetime.now().strftime('%B %Y')),            
            })
    
    slides.append({'style': 'Text',
                   'title': 'Content overview',
                   'para1': [(0, sec) for sec in sections]})
    
    if 'Introduction' in sections:
        slides.append({'style': 'Section',
                       'title': 'Introduction'
                       })
        
    

#%%Scope of work
    if 'Scope of work' in sections:
        slides.append({'style': 'Section',
                       'title': 'Scope of work'
                       })

        slides.append({'style': 'Text',
                       'title': 'Country context and rationale for study',
                       'para1': [(0, 'Rationale'),
                                ]})
        
        
        slides.append({'style': 'Text',
                       'title': 'Scope of work',
                       'para1': [(0, 'Scope of work'),
                                ]})
        

#%%Slide sections
    for section in plot_sets.keys():
        if verbose: print('Adding section %s  (%i of %i)'%(section, list(plot_sets.keys()).index(section)+1, len(plot_sets.keys())))
        plot_set = plot_sets[section]
        section_results = [res for res in results if res.name in plot_set['resnames']]
        uncertainty = 'uncertainty' if plot_set['plot_uncertainty'] else 'best'
        slides.append({'style': 'Section',
                       'title': section
                       })
        
        if len(plot_set['resnames'])>1:
            slides+= _get_scenario_description_slides(P=P, scen_name=section, plot_set=plot_set, progset_name=progset_name, allpops=allpops)
        
        for plot_category in plot_set['plots']:
            if plot_category == 'cascade':
                slides+= _get_cascade_slides(f=f, P=P, results=section_results, allpops=allpops, differentiation=section, uncertainty=uncertainty, 
                                          title_prefix = section+': ')
            elif plot_category == 'demographics':
                slides+= _get_demographic_slides(f=f, P=P, results=section_results, allpops=allpops, differentiation=section, uncertainty=uncertainty, 
                                          title_prefix = section+': ')
            elif plot_category == 'keypars':
                slides+= _get_keypar_slides(f=f, P=P, results=section_results, allpops=allpops, differentiation=section, uncertainty=uncertainty, 
                                          title_prefix = section+': ')
            elif plot_category == 'programs':
                slides+= _get_programs_slides(f=f, P=P, results=section_results, allpops=allpops, differentiation=section, uncertainty=uncertainty, 
                                          title_prefix = section+': ')
            
            elif plot_category == 'spending':
                programs = P.progsets[progset_name].programs.keys()
                slides+= _get_slide_set(f, ['spending_', section], ['comparison_equivalent_spendingbars', 'comparison_spendingbars'], not_including=programs, total_slide = False, slide_title='%s: spending comparison'%(section))
                slides+= _get_slide_set(f, ['spending_over_time', section], ['comparison_equivalent_spendingbars', 'comparison_spendingbars'], not_including=programs, total_slide = False, slide_title='%s: spending comparison'%(section))
                resnames = [res.name for res in section_results]
                slides+= _get_slide_set(f, ['spending_over_time', 'N.A.'], resnames, not_including=programs, total_slide = False, slide_title='%s: result comparison'%(section))

#                slides+= _get_slide_set(f, ['spending boxplot', section], [''], total_slide = False, slide_title='%s: spending uncertainty range'%(section))
                slides[-1]['text1'] = 'Taller bars show the optimized allocation for each scenario, with all other optimized results within 10% of the best score shown, to demonstrate that within the uncertainty of the model, alternative allocations may be feasible depending on new evidence or changing costs.'

            
            
            elif plot_category in plot_categories:
                variables_to_plot = plot_categories[plot_category]
                slides+= _get_category_slides(f, P, results, allpops, differentiation=section, uncertainty=uncertainty, 
                                              variables_to_plot = variables_to_plot, title_prefix = section+': ')

          


    
#%%Review of programs and costing   
    if 'Review of programs and costing' in sections:
        if verbose: print('Adding slides detailing programs and costings')
        slides.append({'style': 'Section',
                       'title': 'Technical annex:\nReview of programs and costing'
                       })
        
        progsplit = round(len(P.progsets[progset_name].programs)/2)
        first_half = [program.label for program in P.progsets[progset_name].programs.values()[:progsplit]]
        second_half  = [program.label for program in P.progsets[progset_name].programs.values()[progsplit:]]
        if len(second_half)<len(first_half): second_half+=['']
        total_budget = P.progsets[progset_name].get_alloc(tvec=key_year)[:].sum()
        slides.append({'style': 'Table',
                       'title': 'Program overview',
                       'table1': [['Program', '']] + [[first_half[i], second_half[i]] for i in range(progsplit)],
                       'para1': [(0, 'The most recently reported budget available for %i was estimated to be %s %s'%(key_year, currency, sc.sigfig(X=total_budget, sigfigs=2, SI=False, sep=True, keepints=True))),
                                 (0, 'This budget does not include spending on non-targeted programs including management, infrastructure, and monitoring'),
                                 ]
#                       'text2': ''
                       })
        
        #Make one slide for each program with the key details
        for prog in P.progsets[progset_name].programs.keys():
            program = P.progsets[progset_name].programs[prog]
            
            applies_pops = _get_popstring(P, program.target_pops)
            applies_comps = _get_compstring(P, program.target_comps)
            
            data = [['Parameter', 'Value', 'Comments']]
            data.append(['Annual spending %i'%(program.spend_data.t[-1]), '%s %s'%(sc.sigfig(program.spend_data.vals[-1], sep=True, keepints=True), program.spend_data.units.replace('$', currency)), ''])
            data.append(['Unit cost %i'%(program.unit_cost.t[-1]), '%s %s'%(sc.sigfig(program.unit_cost.vals[-1], sep=True, keepints=True), program.unit_cost.units.replace('$', currency)), ''])
            data.append(['Coverage %i'%(program.coverage.t[-1]), '%s %s'%(sc.sigfig(program.coverage.vals[-1], sep=True, keepints=True), program.coverage.units.replace('$', currency)), ''])
            if len(program.saturation.vals)>0:
                data.append(['Saturation %i'%(program.saturation.t[-1]), '%s%%'%(sigfigs(100*program.saturation.vals[-1])), 'Maximum percentage of the target population(s) that could be reached  (demand constraint).'])
            if len(program.capacity_constraint.vals)>0:
                data.append(['Capacity constraint %i'%(program.capacity_constraint.t[-1]), '%s'%(sigfigs(program.capacity_constraint.vals[-1])), 'Maximum number of people that could be covered per year (supply constraint).'])
            
            slides.append({'style': 'Table',
                       'title': program.label,
                       'table1': data,
                       'para1': [(0, 'Target populations: %s'%applies_pops),
                                 (0, 'Target disease status: %s'%applies_comps),
                                 (0, 'Program effects: %s'%applies_comps),
                                 ],
#                       'text2': ''
                       })
            affected_covouts = [covout for covout in P.progsets[progset_name].covouts.values() if prog in covout.progs.keys()]
            affected_pars = sorted(list(set([covout.par for covout in affected_covouts]))) #remove duplicates
            for par in affected_pars:
                #Daily treatment/diagnosis probabilities aren't 'human scale' so replace them
                affected_parameter, invert = _get_inversion(P, par)
                
                affected_pops = [cov.pop for cov in affected_covouts if cov.par == par]
                applies_pops = _get_popstring(P, program.target_pops, case = 'lower')
                

                cov_baselines = [cov.baseline for cov in affected_covouts if cov.par == par]
                baseline_str = 'Baseline: '+  _get_summary_effect_str(P, affected_pops, cov_baselines, invert)
                
                effect_levels = [cov.progs[prog] for cov in affected_covouts if cov.par == par]
                effect_str = 'Program effect: '+  _get_summary_effect_str(P, affected_pops, effect_levels, invert)
                
                if baseline_str.replace('Baseline value', 'Program effect') != effect_str:
                    slides[-1]['para1'].append((1, '%s applied to %s'%(affected_parameter, applies_pops)))
                    slides[-1]['para1'].append((2, baseline_str))
                    slides[-1]['para1'].append((2, effect_str))
                    
                if invert: 
                    if 'diagnosis' in affected_parameter: slides[-1]['text2'] = 'Average times until diagnosis/treatment include the possibility that people do not receive treatment at all.'
                    elif 'treatment' in affected_parameter: slides[-1]['text2'] = 'People not receiving official treatment may be receiving private non-standard care.'
                


#%%Acknowledgements
    if 'Acknowledgements' in sections:
        slides.append({'style': 'Section',
                       'title': 'Acknowledgements'
                       })


    sc.savepptx(filename=file_path, template=template_path, slides=slides, verbose=True)
    
