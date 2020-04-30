# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 12:46:27 2019

@author: rowan.martin-hughes
"""
import atomica as at
import sciris as sc
import numpy as np
from os import sep

from malaria_utils import try_loading, getpops, get_paths

def get_time_varying_ts(ts: at.TimeSeries = None, t: np.array = None, val: float = None, specification: dict = None) -> at.TimeSeries:
    """

    :param t: Array of time values to output spending at
    :param spend: Scalar annual spend
    :param specification: Additional options for conversion e.g. `{'start_time': 2009., 'end_time': 2035., 'cover_time': 0.25, 'timespan': 3.0, 'adjust_type': 'flat'}`
    :return: A TimeSeries suitable for inclusion as a program allocation
    """
    def _convert_ts(ts, dt=5/365, start_time = None, end_time = None, cover_time = 0.5, timespan = 3., adjust_type='flat'):
#        if start_time is None: start_time = min(ts.t)
#        if end_time is None: end_time = 2035.
        
        if adjust_type == 'flat':
            new_ts         = sc.dcp(ts)
            new_times      = list(np.arange(max([min(ts.t), start_time]), max([max(ts.t),end_time])+dt, dt)) #TODO consider whether to force new time scale or just use as supplied?
#            new_times      = [np.round(time) if time%1<0.00000001 else time for time in new_times] #floating point errors can mess things up
            new_vals_unmod = new_ts.interpolate(new_times, method='previous')
            multi          = timespan/float(cover_time)
            new_vals       = [nvu * multi if ((new_times[ind] - start_time)%timespan < cover_time) else 0. for ind, nvu in enumerate(new_vals_unmod)]
            new_ts.t       = new_times
            new_ts.vals    = new_vals
        else:
            raise Exception('No other adjust types implemented')

        return new_ts
    
    if ts is None:
        assert (t is not None and val is not None), 'Must specify either a time series or a time and spending value'
        ts = at.TimeSeries(t=list(t), vals=[val for _ in t])
    
    if specification:
        ts = _convert_ts(ts, **specification)
        
    return ts

def get_scenario_instructions(tvec, alloc: dict = None, coverage: dict = None, time_adjustments: dict = None):
    #TODO tidy up the duplication and inefficiency in this function
    tvec = sc.promotetoarray(tvec)

    instructions = at.ProgramInstructions(start_year = min(tvec))
    if alloc:
        for prog_name in alloc.keys():
            """Bunch of wrangling that could be tidied up later to take an original tvec alloc with perhaps a few years, and then 
            turn it into a bunch of time varying allocations to cover all the times in tvec
            """
            if isinstance(alloc[prog_name], at.TimeSeries):
                orig_vals = alloc[prog_name].vals
                orig_ts   = alloc[prog_name].t
                
                if time_adjustments and prog_name in time_adjustments.keys():
                    adj_ts = []
                    adj_vals = []
                    #figure out which tvecs should be adjusted for
                    for time, val in zip(orig_ts, orig_vals):
                        first_ind = np.where(tvec>=time)[0][0] #first time in the tvec that's equal to or greater than the time
                        if first_ind and time!=max(orig_ts):
                            next_time = orig_ts[np.where(np.array(orig_ts)>time)[0][0]] #first time that's larget than the time we're looking at
                            last_ind = np.where(tvec<next_time)[0][-1]
                            new_tvec = tvec[first_ind:last_ind]
                        else:
                            new_tvec = tvec[first_ind:]
                        new_ts = get_time_varying_ts(t=new_tvec, val=val, specification=time_adjustments[prog_name])
                        adj_ts += new_ts.t
                        adj_vals += new_ts.vals
                        #todo tidy this up to append things into a TimeSeries or something??
                    instructions.alloc[prog_name] = at.TimeSeries(adj_ts, adj_vals)
                else:
                    instructions.alloc[prog_name] = sc.dcp(alloc[prog_name])
            else:
                if time_adjustments and prog_name in time_adjustments.keys():
                    print ('Adjusting %s for time varying allocation but did not get a timeseries input, so assuming'%(prog_name))
                    #TODO why use [1:] in the tvec? Good question! First time point seems to be an error
                    times = tvec[1:] if len(tvec)>1 else tvec
                    instructions.alloc[prog_name] = get_time_varying_ts(t=times, val=alloc[prog_name][-1], specification=time_adjustments[prog_name])
                else:
                    #did not get a timeseries, there are no time adjustments to make...
                    instructions.alloc[prog_name] = at.TimeSeries(tvec[0], alloc[prog_name][-1])
    if coverage:
        for prog_name in coverage.keys():
            """Bunch of wrangling that could be tidied up later to take an original tvec coverage with perhaps a few years, and then 
            turn it into a bunch of time varying coverages to cover all the times in tvec
            """
            if isinstance(coverage[prog_name], at.TimeSeries):
                orig_vals = coverage[prog_name].vals
                orig_ts   = coverage[prog_name].t
                
                if time_adjustments and prog_name in time_adjustments.keys():
                    adj_ts = []
                    adj_vals = []
                    #figure out which tvecs should be adjusted for
                    for time, val in zip(orig_ts, orig_vals):
                        first_ind = np.where(tvec>=time)[0][0] #first time in the tvec that's equal to or greater than the time
                        if first_ind and time!=max(orig_ts):
                            next_time = orig_ts[np.where(np.array(orig_ts)>time)[0][0]] #first time that's larget than the time we're looking at
                            last_ind = np.where(tvec<next_time)[0][-1]
                            new_tvec = tvec[first_ind:last_ind]
                        else:
                            new_tvec = tvec[first_ind:]
                        new_ts = get_time_varying_ts(t=new_tvec, val=val, specification=time_adjustments[prog_name])
                        adj_ts += new_ts.t
                        adj_vals += new_ts.vals
                        #todo tidy this up to append things into a TimeSeries or something??
                    instructions.coverage[prog_name] = at.TimeSeries(adj_ts, adj_vals)
                else:
                    instructions.coverage[prog_name] = sc.dcp(coverage[prog_name])
            else:
                if time_adjustments and prog_name in time_adjustments.keys():
                    print ('Adjusting %s for time varying coverage but did not get a timeseries input, so assuming'%(prog_name))
                    #TODO why use [1:] in the tvec? Good question!
                    times = tvec[1:] if len(tvec)>1 else tvec
                    instructions.coverage[prog_name] = get_time_varying_ts(t=times, val=coverage[prog_name][-1], specification=time_adjustments[prog_name])
                else:
                    #did not get a timeseries, there are no time adjustments to make...
                    instructions.coverage[prog_name] = at.TimeSeries(tvec[0], coverage[prog_name][-1])
                
    return instructions



class SeasonalAdjustment(at.Adjustment):
    def __init__(self, prog_name, tvec, specification, initial, minimum_spend, maximum_spend):
        at.Adjustment.__init__(self, name=prog_name)
        self.prog_name = prog_name
        self.specification = specification # e.g. {'start_time': 2009., 'end_time': 2035., 'cover_time': 0.25, 'timespan': 3.0, 'adjust_type': 'flat'}
        self.t = tvec  # Vector of time values instructions are updated at
        self.adjustables = [at.Adjustable('spend', initial_value=initial, lower_bound=minimum_spend, upper_bound=maximum_spend)]

    def update_instructions(self, adjustable_values, instructions: at.ProgramInstructions):
        # Could make a new copy of `self.specification` and make changes if optimizing something else e.g. start_time, cover_time, timespan etc.
        instructions.alloc[self.prog_name] = get_time_varying_ts(t=self.t, val = adjustable_values[0], specification=self.specification)

class MinimumConstraint(at.Constraint):
    
    def __init__(self, limit, progs, t):
        self.limit = limit
        self.progs = progs
        self.t = t
        
    def constrain_instructions(self, instructions, hard_constraints=None):
        total = 0.0
        for prog in self.progs:
            total += instructions.alloc[prog].get(self.t)
        if total < self.limit:
            raise at.FailedConstraint
        else:
            return 0.0   

def _evaluate_model(P, optimization, parset_name='default', progset_name='default', progset_instructions=None):
    if progset_instructions is None:
        progset_instructions = at.ProgramInstructions(start_year=2019)
    result = P.run_sim(parset=parset_name, progset=progset_name, progset_instructions=progset_instructions)

    evaluation = optimization.compute_objective(result.model, [None]*len(optimization.measurables))

    return evaluation

def _safe_optimize(optim_args):
    """Try optimizing once, return None if it fails"""
    try:
        parset_name = optim_args['parset'].name
        progset_name = optim_args['progset'].name
        optimization = optim_args['optimization']
        opt_progset_instructions = at.optimize(**optim_args)
        
        #run once more to get an evaluation
        result = optim_args['project'].run_sim(parset=parset_name, progset=progset_name, progset_instructions=opt_progset_instructions)
    
        evaluation = optimization.compute_objective(result.model, [None]*len(optimization.measurables))
        
        optimized_result = {'evaluation': evaluation, 'progset_instructions': opt_progset_instructions}
        
        return optimized_result
    except:
        return None

def _repeat_optimization(optim_args, n_repeats=1, parallel=False, varying_inits=True, opt_results_folder=None, load_objects=True):
    if opt_results_folder is None: opt_results_folder='temp'+sep
    """Run asd multiple times and pick the best results"""
    seeds = list(range(n_repeats))
    opt_argset = []
    for s in seeds:
        oa = sc.dcp(optim_args)
        if varying_inits and s>0:
            """Modify the init_alloc for each run (except the first!) - ensure each program receives at least the lower limit then randomly distribute the rest"""
            alloc = sc.odict([(adjust.prog_name,adjust.adjustables[0].initial_value) for adjust in oa['optimization'].adjustments])
            budget_sum = sum(alloc[:])
            new_alloc = sc.odict([(adjust.prog_name,adjust.adjustables[0].lower_bound) for adjust in oa['optimization'].adjustments])
            constrained_sum = sum(new_alloc[:])
            additional_prop = sc.odict([(prog_name, np.random.random()) for prog_name in new_alloc.keys()])
            additional_prop[:]*=(budget_sum - constrained_sum)/sum(additional_prop[:]) #normalize so that the sum is 1, and multiply by the additional budget
            new_alloc[:]+=additional_prop[:]
            
            for a,_ in enumerate(oa['optimization'].adjustments):
                oa['optimization'].adjustments[a].adjustables[0].initial_value = sc.dcp(new_alloc[oa['optimization'].adjustments[a].prog_name])
        opt_argset.append(oa)
        
    if n_repeats==1 or not parallel:
        print ('Not running in parallel: ', n_repeats, parallel)
        for s, oa in enumerate(opt_argset):
            print('Run %s'%(s))
            #            print('using %s as initial'% sc.odict([(adjust.prog_name,adjust.adjustables[0].initial_value) for adjust in oa['optimization'].adjustments]))
            try_loading(function = _safe_optimize, fnargs={'optim_args': oa}, obj_filename='Optimization run %s'%(s),
                        obj_folder=opt_results_folder, load_objects=load_objects, run_if_unfound=True)
    else:
        print('Running %s optimizations in parallel (WARNING: not tested)'%n_repeats)
        parallel_args = [{'function': _safe_optimize, 'fnargs':{'optim_args': oa}, 'obj_filename':'Optimization run %s'%(s),
                        'obj_folder':opt_results_folder, 'load_objects':load_objects, 'run_if_unfound':True} for s, oa in enumerate(opt_argset)]
        optimized_results = sc.parallelize(func=try_loading, iterkwargs=parallel_args, maxload=0.5)
    
    return get_optimized_results(opt_results_folder, what='best_instructions')


#def verify_staircase_optimization(P, optimized_results, optimizations):
#    """Take a list of optimized budget allocations, and a list of optimization constraints
#       - test each optimized result to see if one of the adjacent optimized allocations just scaled would lead to an improved outcome
#       - if so, then rerun that optimization with an initialization of the rescaled budget"""
#    assert len(optimizations) == len(optimized_results), 'Must be one optimization object for each optimized_result'
#       
#    #TODO re-order to make sure they are in ascending budget order
#    verified = [False for res in optimized_results]
#    
#    while not np.array(verified).all():
#        unverified_ind = verified.index(False)
#        current_eval = optimized_results['evaluation']
#        for ind in [unverified_ind-1, unverified_ind +1]:
#            if ind >=0 and ind<len(verified):
#                best_known_allocation = optimized_results[unverified_ind]['progset_instructions']
#                unconstrained_bka = sc.dcp(optimized_results[unverified_ind]['progset_instructions'])
#                for adjust in optimizations:
#    
#                
#                adjacent_allocation = optimized_results[ind]['progset_instructions'].alloc
#                unconstrained_adjacent_allocation = 
#                alternative_allocation = 
#        
#        
#    
#    
#    return optimized_results


def optimize_standard(P, optim_name='Optimization', progset_name='default', parset_name='default', start_year=2019., end_year=2035., max_time = None,
                      targets = None, budget_factor=1.0, method='asd', n_repeats = 1, just_optimization_object = False,
                      unfunded_progs=[], constrained_progs=[], parallel=False, opt_results_folder = None, load_objects=True,
                      standard_lower_limit = True, standard_upper_limit = False, time_adjustments = None, prog_categories = {}
                      ):
    """
    unfunded_progs are programs that should not be included in optimizations at all:
        e.g. ART might have a budget and coverage but is not included as part of the TB budget
    constrained_progs are programs where the current budget should not be reduced
        e.g. GeneXpert testing should not be reduced.
    """
    
    if targets is None:
        targets = {'ds_prev': 1., 'mdr_prev':1000., 'xdr_prev':1000000.}
    
    def_budget = P.progsets[progset_name].get_alloc(tvec=start_year)
    tot_funded = sum([def_budget[prog][0] for prog in def_budget.keys() if not prog in unfunded_progs]) * budget_factor
    
    adjustments = []
    alloc = sc.odict()
    measurables = []
    
    all_progs = list(def_budget.keys())
    remaining_opt_progs = [prog for prog in all_progs if not prog in unfunded_progs]
    
    if not prog_categories is None:
        yrs = sc.promotetolist(start_year)
        for cat, prog_cat in prog_categories.items():
            if isinstance(prog_cat, dict):
                prog_names = list(prog_cat.keys())
                min_props       = [bounds[0] for bounds in prog_cat.values()]
                max_props       = [bounds[1] for bounds in prog_cat.values()]
            else:
                prog_names = prog_cat
                min_props       = [0. for _ in prog_names]
                max_props       = [1. for _ in prog_names]
            remaining_opt_progs = [prog for prog in remaining_opt_progs if not prog in prog_names]
            
            initial_spends  = [def_budget[prog_name][-1]*budget_factor for prog_name in prog_names] #TODO should this allow the array without the [-1]?
            tot_initial_spend = sum(initial_spends)
            #if ANY of the programs is constrained then don't allow the group funding to be reduced
            min_total_spend = tot_initial_spend if np.array([prog_name in constrained_progs for prog_name in prog_names]).any() else 0.
            max_total_spend = max(2.*tot_initial_spend, tot_funded/10.) if standard_upper_limit else tot_funded
            
            adjustments.append(at.SpendingPackageAdjustment(package_name=cat, t=yrs, prog_names=prog_names, initial_spends = initial_spends,
                                                   min_props = min_props, max_props = max_props,
                                                   min_total_spend = min_total_spend, max_total_spend = max_total_spend))
                
    #Include all programs in the spending allocation
    for prog in all_progs:
        prog_scaling = 1. if prog in unfunded_progs else budget_factor
        alloc[prog] = def_budget[prog][0] * prog_scaling
        #Add a regular spending adjustment for any other programs not specified!
        if prog in remaining_opt_progs:
            min_budget_constrained   = min(alloc[prog], def_budget[prog][0])
            min_budget_unconstrained = 0. if standard_lower_limit else min_budget_constrained/2.
            lower = min_budget_constrained if prog in constrained_progs else min_budget_unconstrained #set a minimum of 50% of the current budget on all programs
#            upper = np.inf
            upper = max(2.*def_budget[prog][0], tot_funded/10.) if standard_upper_limit else tot_funded #set a maximum of a doubling of spending or 10% of the total budget
            adjustments.append(at.SpendingAdjustment(prog_name = prog, t=start_year, initial=alloc[prog], limit_type='abs', lower=lower, upper=upper))
    
    for par in targets.keys():
        measurables.append(at.Measurable(par, t=[start_year, end_year], weight=targets[par]))

    instructions = at.ProgramInstructions(start_year=start_year, alloc=alloc)
    constraints = [at.TotalSpendConstraint(total_spend=tot_funded, t=start_year, budget_factor=1.0)]
    for constrained in constrained_progs:
        if isinstance(constrained, list) or (isinstance(constrained, str) and not constrained in remaining_opt_progs): 
#        if isinstance(constrained, list): #it's a list of programs where the sum of those should not be reduced.
            init_sum = sum([def_budget[pn][0] * prog_scaling for pn in def_budget.keys() if pn in constrained])
            constraints.append(MinimumConstraint(limit=init_sum, progs=sc.promotetolist(constrained), t=start_year))
    
    optimization = at.Optimization(name='optimized', adjustments=adjustments, measurables=measurables,constraints=constraints,
                                   method=method, maxtime=max_time, maxiters = max_time) # Evaluate from 2020 to end of simulation
    
    if just_optimization_object:
        return optimization
    
    optim_args = {'project': P, 'optimization': optimization, 'parset':P.parsets[parset_name],
                  'progset':P.progsets[progset_name],'instructions':instructions}
    optimized_instructions = _repeat_optimization(optim_args, n_repeats = n_repeats, parallel=parallel,
                                                  opt_results_folder = opt_results_folder, load_objects=load_objects)
    
    
    return optimized_instructions

#%%HANDLE RESULTS AFTER THEY HAVE BEEN CACHED

    
def verify_optimization(P, opt_results_folder, parset_name, progset_name, optimizations=None,
                        unfunded_progs = [], max_loops = 5, verbose=True, response = 'adjust'):
    """Take a list of budget allocations, and a list of optimization constraints
       - test each optimized result to see if one of the alternative optimized allocations just scaled would lead to an improved outcome
       - if so, then rerun that optimization with an initialization of the rescaled budget
       
       TODO: This assumes that the only thing we're modifying is the LAST value in the alloc, will not be applicable for time varying optimizations!
       """
    if optimizations is None:
           optimizations = P.optims
       
    optimized_results_dict = sc.odict()
    for opt_name in P.optims.keys():
        opt_result_folder = opt_results_folder%(opt_name)
        optimized_results_dict[opt_name] = get_optimized_results(opt_results_folder=opt_result_folder,what='best_object')
    
    #remove any optimization/optimized_result pairs where no valid results exist at all
    optimized_results = list(optimized_results_dict.values())
#    optimizations     = [opt for opt_ind, opt in enumerate(optimizations) if not optimized_results[opt_ind] is None]
    optimized_results = [opt_res for opt_res in optimized_results if not opt_res is None]
    
    assert len(optimizations) == len(optimized_results), 'Must be one optimization object for each optimized_result'
    
    verified = False
    tot_funded_spending = [0. for _ in optimized_results]
    for on, opt_result in enumerate(optimized_results):
        tot_funded_spending[on] = sum([spend.vals[-1] for prog, spend in opt_result['progset_instructions'].alloc.items() if not prog in unfunded_progs])
    
    loops = 0
    while not np.array(verified).all() and loops < max_loops:
        print ('Comparing optimized allocations, step %s'%(loops))
        verified = True

        #2. for each optimized_result
        for on, opt_result in enumerate(optimized_results):
            #3.for each OTHER optimized_result
            new_proposed_instructions = None
            new_proposed_evaluation = opt_result['evaluation']
            
            for on_alt, opt_result_alt in enumerate(optimized_results):
                if on_alt != on and opt_result['progset_instructions'].alloc.keys() == opt_result_alt['progset_instructions'].alloc.keys(): #don't even think about it if programs aren't the same...
                    on_name = optimized_results_dict.keys()[on]
                    on_name_alt = optimized_results_dict.keys()[on_alt]
                    
                    #4. if the proportional allocation from the OTHER optimized_result does not violate any constraints
                    scale_factor = tot_funded_spending[on] / tot_funded_spending[on_alt]
                    optim = sc.dcp(optimizations[on])
                    proposed_instructions = opt_result_alt['progset_instructions'].scale_alloc(scale_factor=scale_factor) #does return a new instruction object
                    for prog in unfunded_progs:  #TODO these shouldn't be scaled given the way these progs are set up elsewhere,
                        #but this adjustment might not be appropriate more broadly, dangerous assumption
                        proposed_instructions.alloc[prog].vals = [val / scale_factor for val in proposed_instructions.alloc[prog].vals]
#                    proposed_instructions = sc.dcp(opt_result['progset_instructions'])
#                    for prog in proposed_instructions.alloc.keys():
#                        if not prog in unfunded_progs:
#                            proposed_instructions.alloc[prog].vals[-1] *= scale_factor
                    
                    try:
                        x0,_,_ = optim.get_initialization(progset=P.progsets[progset_name], instructions=proposed_instructions)
                        hard_constraints = optim.get_hard_constraints(x0=x0, instructions=proposed_instructions)
                        constrained_penalty = optim.constrain_instructions(proposed_instructions, hard_constraints)
                    except:
                        constrained_penalty = np.inf #don't even want to know
                    
                    if constrained_penalty == 0.:
                        #5.run a model with the alternative budget
                        proposed_evaluation = _evaluate_model(P, optimizations[on], parset_name=parset_name,progset_name=progset_name,
                                                              progset_instructions=proposed_instructions)
                        if proposed_evaluation < opt_result['evaluation'] - 1e-15:
        #6.               if the result is better, then...  rerun optimization starting from alternative budget proportions but given total budget
                            print ('+++ Alloc for %s scaled using prop of %s results in a better outcome (%.2f vs %.2f)!'%(
                                   on_name, on_name_alt, opt_result['evaluation'], proposed_evaluation))
                            if proposed_evaluation < new_proposed_evaluation:
                                new_proposed_instructions = sc.dcp(proposed_instructions)
                        else:
                            if verbose: print ('--- Alloc for %s scaled using prop of %s results in a worse outcome (%.2f vs %.2f)!'%(
                                   on_name, on_name_alt, opt_result['evaluation'], proposed_evaluation))
                    else:
                        if verbose: print ('--- Alloc for %s scaled using prop of %s incomparable (would violate constraint) (penalty %s)!'%(
                                   on_name, on_name_alt, constrained_penalty))            #having made the comparisons, is anything better?
            if not new_proposed_instructions is None:
                if response is 'none':
                    print ('Optimization result (instructions) %s is not optimal, taking no action'%(on_name))
                elif response is 'replace':
                    opt_result_folder = opt_results_folder%(on_name)
                    obj_filename = 'Optimization run verification' #
                    print ('Optimization result (instructions) %s is not optimal, adding the new instructions to %s'%(on_name, opt_result_folder+obj_filename))
                    obj = {'evaluation': new_proposed_evaluation, 'progset_instructions': new_proposed_instructions}
                    sc.saveobj(obj=obj, filename=obj_filename, folder=opt_result_folder)
                    
#                    P.scens[on_name] = at.BudgetScenario(name=on_name, parsetname=parset_name, progsetname=progset_name,
#                                alloc = new_proposed_instructions.alloc)
            elif response=='reoptimize':
                raise Exception('Not implemented yet')
                verified = False #would require a further round of verification after this if making changes, in case starting from the new point finds an optimized allocation that would then be better at a different budget level
                    
    return True #TOODO this should probably return something useful?



def get_optimized_results(opt_results_folder=None, what='all'):
    """Return either all the optimized results from a folder, or just the instruction set of the best one"""
    paths = get_paths(folder=opt_results_folder, inclusions=['Optimization', 'run'], extensions='', version='all', verbose=False)
    if paths is None: 
        return None

    optimized_results = []
    for path in paths:
        optimized_results.append(sc.loadobj(path))
    if what=='all':
        return optimized_results
    else:
        evaluations  = [opt_result['evaluation'] for opt_result in optimized_results if not opt_result is None]
        instructions = [opt_result['progset_instructions'] for opt_result in optimized_results if not opt_result is None]
        best_evaluation = min(evaluations) if evaluations != [] else None
        if best_evaluation is None:
            print ('No valid evaluations found!')
            return None
        else:
#            best_instructions = instructions[evaluations.index(best_evaluation)] #change to be consistent whatever you ask for...
            best_instructions = [opt_result for opt_result in optimized_results if (opt_result is not None) and opt_result['evaluation']==best_evaluation][-1]['progset_instructions']
            print('Compared results with objective functions of %s, picked %s'%(evaluations, min(evaluations)))
            
        if what=='best_instructions':
            return best_instructions
        if what=='best_object':
            return [opt_result for opt_result in optimized_results if (opt_result is not None) and opt_result['evaluation']==best_evaluation][-1]
        else:
            raise Exception('Unknown return request from optimizations folder %s.'%(what))


def make_optim_scens(P, result, opt_results_folder, optim_names, parset_name,
                     progset_name, coverage_progs, time_adjustments, start_year=2019):
    """Return a dict of scenarios based on the optimized allocations for each optimization name"""
    scens = sc.odict()
    
    for opt_name in optim_names:
#        print ('getting results for %s from %s'%(opt_name, opt_results_folder%(opt_name)))
        opt_result_folder = opt_results_folder%(opt_name)
        
        opt_insts = get_optimized_results(opt_results_folder=opt_result_folder, what='best_instructions')
        
        if not opt_insts is None:
            #Make a scenario, but include historical spending as best available in this scenario even if it won't be used
            alloc_with_historical = opt_insts.alloc
            
            for prog in alloc_with_historical.keys():
                prog_spend = P.progsets[progset_name].programs[prog].spend_data
                hist_vals = [(yr, val) for yr, val in list(zip(prog_spend.t, prog_spend.vals)) if yr < min(opt_insts.alloc[prog].t)]
                for yr, val in hist_vals:
                    alloc_with_historical[prog].insert(t=yr, v=val)
#            print ('... alloc for IRS is %s'%(alloc_with_historical['IRS'].vals))
            scens[opt_name] = make_combined_scenario(P=P, result=result, scen_name=opt_name,
                                                       parset_name=parset_name, progset_name=progset_name, alloc = alloc_with_historical,
                                                       coverage_progs = coverage_progs, time_adjustments=time_adjustments,
                                                       start_year = start_year)
            
#            at.BudgetScenario(name=opt_name, parsetname=parset_name, progsetname=progset_name, alloc = alloc_with_historical)
#            P.make_scenario(which='budget', name=opt_key, parsetname=parset_name, progsetname=progset_name,
#                            alloc = alloc_with_historical)
        
    return scens

#%% STANDARD BASELINE SCENARIO DEFINITIONS
def scenfn_calibration(P, scen_name='Calibration', parset_name='default', **kwargs):
    scen_values = sc.odict() #change nothing in this scenario!
    scen = at.ParameterScenario(name=scen_name, parsetname=parset_name, scenario_values=scen_values)
    return scen


def _get_annual_average_coverage(result, start_year, progs=None):
#    tinds = np.where(np.logical_and(result.t>=start_year, result.t<start_year+1.))[0]
    coverage = sc.odict()
    annual_frac = result.get_coverage(quantity='annual_fraction', year=[start_year])
    if progs is None:
        progs = annual_frac.keys()
    for prog_name in progs:
        #capacity averaged over the year as values are annualized
#        d_cap    = at.PlotData.programs(result, outputs=prog_name, quantity='coverage_capacity')
#        cap_annual  = d_cap[result.name, 'N.A.', prog_name].vals[tinds].mean()
#        #eligibility summed over the year if coverage is one-off, averaged if it's continuous
#        d_elig   = at.PlotData.programs(result, outputs=prog_name, quantity='coverage_eligible')
#        fn = np.mean if '/year' in result.model.progset.programs[prog_name].unit_cost.units else np.sum
#        elig_annual = fn(d_elig[result.name, 'N.A.', prog_name].vals[tinds])
#        prog_annual_cover = cap_annual/elig_annual  #TODO should scale this somehow with capacity constraints etc
#        coverage[prog_name] = np.array([prog_annual_cover])
        coverage[prog_name] = annual_frac[prog_name]
    
    return coverage


def scenfn_current_spending(P, scen_name='Current spending', result=None, parset_name='default', progset_name='default',
                            start_year=2019, scale_year=2020, time_adjustments: dict =None, coverage_progs: list = []):
    alloc = sc.odict([(prog_name, at.TimeSeries(t=start_year, vals=P.progsets[progset_name].get_alloc(start_year)[prog_name])) for prog_name, vals in P.progsets[progset_name].get_alloc(start_year).items() if not prog_name in coverage_progs])
    
    if result is None:
        result = P.run_sim(parset_name=parset_name, progset=progset_name, progset_instructions=at.ProgramInstructions(start_year=start_year))
#    coverage = sc.odict([(prog_name, vals) for prog_name, vals in result.get_coverage(year=[start_year]).items() if prog_name in coverage_progs]) 
    coverage = _get_annual_average_coverage(result, start_year, progs=coverage_progs)
    tvec = [t for t in P.settings.tvec if t>=start_year]
    
    instructions = get_scenario_instructions(tvec = tvec, coverage = coverage, alloc = alloc,
                                             time_adjustments = time_adjustments) 
    
    scen = at.CombinedScenario(name=scen_name, parsetname=parset_name, progsetname=progset_name, instructions=instructions)
    
    return scen

def scenfn_current_coverage(P, scen_name = 'Current coverage', result=None, parset_name='default', progset_name='default',
                            start_year=2019, time_adjustments=None, **kwargs):
    if result is None:
        result = P.run_sim(progset=progset_name, progset_instructions=at.ProgramInstructions(start_year=start_year))
#    coverage = result.get_coverage(year=[start_year]) #okay for continuous programs
    #Get the target pop size over 1 year
    coverage = _get_annual_average_coverage(result, start_year)
            
    tvec = [t for t in P.settings.tvec if t>=start_year]
    instructions = get_scenario_instructions(tvec = tvec, coverage = coverage, alloc = None,
                                             time_adjustments = time_adjustments) 
    
    scen = at.CombinedScenario(name=scen_name, parsetname=parset_name, progsetname=progset_name, instructions=instructions)
    
    return scen


def scenfn_current_conditions(P, scen_name='Current conditions', result=None, parset_name='default', progset_name=None,
                              start_year=2018, scale_year=2020, time_adjustments={}):
#    allpops = getpops(P)
    
    if result is None:
        result = P.run_sim(parset=parset_name)
    
#    frs = P.framework.pars['scenario'] #keep some parameters constant
#    probability_pars = list(frs.keys()[np.where(frs=='TB treatment initiation')]) + list(frs.keys()[np.where(frs=='TB diagnosis')])
    
    scen_values = sc.odict()
    
#    for par in probability_pars:
#        scen_values[par] = sc.odict()
#        for pop in allpops:
#            popind = result.pop_names.index(pop)
#            start_ind = np.where(result.get_variable(par)[popind].t == start_year)[0][0]
##            indices = np.array([np.where(result.get_variable(par)[popind].t == start_year)[0][0],
##                       np.where(result.get_variable(par)[popind].t == scale_year)[0][0]])
#            scen_values[par][pop] = sc.odict()
#            scen_values[par][pop]['t'] = [start_year] #, scale_year]
#            scen_values[par][pop]['y'] = [result.get_variable(par)[popind].vals[start_ind]]
            
    
    scen = at.ParameterScenario(name=scen_name, parsetname=parset_name, scenario_values=scen_values)
    return scen

#%%PROGRAM SCENARIO PACKAGE
def scenfn_testprograms(P, scen_name='TP', result=None, parset_name='default', progset_name='default',
                        prog_names=None, start_year=2019, scale_year=2020, time_adjustments: dict = None, coverage_progs=[]):
    scens = sc.odict()
    def_alloc = P.progsets[progset_name].get_alloc(start_year)
    ts_alloc = sc.odict([(prog, at.TimeSeries(t=[start_year], vals=def_alloc[prog])) for prog in def_alloc.keys() if not prog in coverage_progs])
    prop_budgets = {'unfunded': 0., 'at saturation': 9e99}#zero and infinite budgets
    
    if prog_names is None:
        #Don't include any coverage programs as these wouldn't make sense
        prog_names = [prog_name for prog_name in P.progsets[progset_name].programs.keys() if not prog_name in coverage_progs]
    
#    coverage = _get_annual_average_coverage(result, start_year, progs=coverage_progs)
    tvec = [t for t in P.settings.tvec if t>=start_year]
    
    for pb in prop_budgets.keys(): 
        #For each individual program
        for prog in prog_names:
            if not prog in coverage_progs and ts_alloc[prog].vals[-1]!=prop_budgets[pb]: #likely only relevant for the zero budget - if the budget is already zero no point in comparing with a zero budget!
                alloc = sc.dcp(ts_alloc)
#                alloc[prog].t.append(scale_year)
                alloc[prog].vals = [prop_budgets[pb]]
                sname = '%s %s %s'%(scen_name, prog, pb)
                scens[sname] = make_combined_scenario(P, scen_name = sname, result=result, parset_name=parset_name,
                                                      progset_name = progset_name, alloc = alloc, coverage_progs = coverage_progs,
                                                      time_adjustments = time_adjustments, start_year = start_year)
                
#                instructions = get_scenario_instructions(tvec = tvec, coverage = coverage, alloc = alloc,
#                                             time_adjustments = time_adjustments) 
#                scens[sname] = at.CombinedScenario(name=sname, parsetname=parset_name, progsetname=progset_name, instructions=instructions)
    
        #Overall e.g. zero or infinite for ALL programs
        alloc = sc.dcp(ts_alloc)
        for prog in P.progsets[progset_name].programs.keys():
            if not prog in coverage_progs:
#                alloc[prog].t = [start_year]
                alloc[prog].vals = [prop_budgets[pb]]
                sname = '%s %s'%('All programs', pb)
                scens[sname] = make_combined_scenario(P, scen_name = sname, result=result, parset_name=parset_name,
                                      progset_name = progset_name, alloc = alloc, coverage_progs = coverage_progs,
                                      time_adjustments = time_adjustments, start_year = start_year)

                
#                instructions = get_scenario_instructions(tvec = tvec, coverage = coverage, alloc = alloc,
#                                             time_adjustments = time_adjustments) 
#                scens[sname] = at.CombinedScenario(name=sname, parsetname=parset_name, progsetname=progset_name, instructions=instructions)
    
    return scens

#%%Optimization scenarios
def make_combined_scenario(P, scen_name, result, parset_name, progset_name, alloc, coverage_progs,
                           time_adjustments, start_year = None):
    """Create a scenario where a budget alloc has been pre-determined"""
    if start_year is None:
        start_year = min([min(ap.t) for ap in alloc.values()]) #if not specified start at the first year there's a bundget

    #remove any spending allocs for programs that are defined by coverage instead
    alloc = sc.dcp(alloc)
    for prog in coverage_progs:
        if prog in alloc.keys():
            del alloc[prog]
    coverage = _get_annual_average_coverage(result, start_year, progs=coverage_progs)
    tvec = [t for t in P.settings.tvec if t>=start_year]
#    print ('%s IRS alloc %s'%(scen_name, alloc['IRS'].vals))
    instructions = get_scenario_instructions(tvec = tvec, coverage = coverage, alloc = alloc,
                                             time_adjustments = time_adjustments)
#    print ('%s IRS INSTRUCTIONS alloc %s'%(scen_name, instructions.alloc['IRS'].vals))
    
    scen = at.CombinedScenario(name=scen_name, parsetname=parset_name, progsetname=progset_name,instructions=instructions)
    scen.instructions.start_year = start_year

    return scen