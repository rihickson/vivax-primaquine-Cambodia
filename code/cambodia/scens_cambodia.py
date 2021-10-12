# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 12:46:27 2019

@author: rowan.martin-hughes
"""

from malaria_utils import *

def scenfn_primaquine(P, result = None, start_year = 2020.83, scale_year = 2020.84, target_pops = None, target_str = None, **kwargs):   # scale year when change in param occurs
    """

    :param P:
    :param result:
    :param start_year:
    :param scale_year:
    :param target_pops: populations to be given Primaquine radical cure
    :param kwargs: Makes sure this includes coverage [worst-case, expected-case, best-case],
    sensitivity_G6PDd rapid test [worst-case, expected-case, best-case], and
    efficacy of the primaquine [worst-case, expected-case, best-case]
    :return:
    """
    scens = sc.odict()  # dict of all coverage, G6PDd sensitivity, and primaquine efficacy combinations
    max_len = max(len(kwargs['coverage']), len(kwargs['sensitivity_G6PDd']), len(kwargs['efficacy']))

    for ind in range(max_len):
        if len(kwargs['coverage']) >= ind:
            cov = kwargs['coverage'][ind]
        else:
            cov = kwargs['coverage'][0]
        if len(kwargs['sensitivity_G6PDd']) >= ind:
            sens = kwargs['sensitivity_G6PDd'][ind]
        else:
            sens = kwargs['sensitivity_G6PDd'][0]
        if len(kwargs['efficacy']) >= ind:
            eff = kwargs['efficacy'][ind]
        else:
            eff = kwargs['efficacy'][0]

        scen = scenfn_calibration(P, result=result, start_year=start_year, scale_year=scale_year, **kwargs)
        scen.scenario_values = sc.odict()

        for par in ['IStx']:  # this is the parameter being changed, proportion
            scen.scenario_values[par] = sc.odict()
            for pop in target_pops:  # only males affected
                popind = result.pop_names.index(pop)
                start_ind = np.where(result.get_variable(par)[popind].t >= start_year)[0][0]
                scen.scenario_values[par][pop] = sc.odict()
                scen.scenario_values[par][pop]['t'] = [start_year] #, scale_year]
                scen.scenario_values[par][pop]['y'] = [result.get_variable(par)[popind].vals[start_ind]]


                scen.scenario_values[par][pop]['t']+= [scale_year]
    #                current_condition = scen.scenario_values[par][pop]['y'][0]
                scen.scenario_values[par][pop]['y']+= [scen.scenario_values[par][pop]['y'][0] * (1-cov) * (1-sens) * (1 - eff)] # this models the primaquine rollout by Cambodia

        scens['scenario_' + str(cov) + '_' + str(sens) + '_' + str(eff) + '_' + target_str] = scen

    return scens

# def scenfn_immediate_testing(P, scale_year = 2025, **kwargs):
#     scen = scenfn_rapid_testing(P, scale_year = scale_year, **kwargs)
#
#     for par in ['p_rdt', 'p_rdt_s']:
#         for pop in scen.scenario_values[par].keys():
#             scen.scenario_values[par][pop]['y'][-1] = [1.0] #100% daily testing probability instead of the "rapid" doubling
#
#     return scen


def scenfn_primaquine_males(**kwargs):
    return scenfn_primaquine(target_pops = ['M 15+'], target_str = 'males', **kwargs)

def scenfn_primaquine_all(**kwargs):
    return scenfn_primaquine(target_pops = ['M 15+', 'Gen'], target_str = 'all', **kwargs)