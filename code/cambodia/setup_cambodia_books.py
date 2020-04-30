# -*- coding: utf-8 -*-
"""
Created on Mon Aug 20 11:37:03 2018

@author: User
"""
import sys
sys.path.insert(0, '/Users/rihickson/repos/malaria-analyses/')
sys.path.insert(0, '/Users/rihickson/repos/atomica')
sys.path.insert(0, '/Users/rihickson/repos/sciris')
from malaria_utils import *

country = 'cambodia'
regions = ['Pursat']  # start with 1 region: start with the region where the primaquine trials have started
# regions = ['Banteay Meanchey', 'Battambang', 'Kampong Chhnang', 'Kampong Speu', 'Kampong Thom', 'Kampot', 'Kandal', 'Kep', 'Koh Kong', 'Kratie', 'Mondul Kiri', 'Oddar Meanchey', 'Pailin', 'Preah Vihear', 'Prey Veng', 'Pursat', 'Stung Treng', 'Svay Rieng', 'Takeo', 'Phnom Penh National Hospital', 'Ratanak Kiri Ratanakiri', 'Siem Reap Siemreap', 'Sihanoukville Preah Sihanouk', 'Kampong Cham Tbong Khmum']

from glob import glob
import socket
user = socket.gethostname()

project_folder = get_apps_folder(country=country)

todo = [
        'setup_databooks',
       # 'setup_progbooks',
       #  'setup_docs'
        ]

framework_path = [file for file in glob(project_folder+'*framework*.xlsx') if not '~' in file][-1]
F = at.ProjectFramework(framework_path)
print(framework_path)

if 'setup_databooks' in todo:

    data_years = [2011, 2019]

    pops = sc.odict()
    pops['M 15+'] = 'Males 18+ years'
    pops['Gen'] = 'Rest of the population'
    
    m_pops = sc.odict()
    m_pops['A. Funestus'] = 'Anopheles funestus'
    
    all_pops = sc.odict([(pn, {'label':pops[pn], 'type':'hum'}) for pn in pops.keys()]+[(pn, {'label':m_pops[pn], 'type':'mos'}) for pn in m_pops.keys()])

    transfers = sc.odict()
    transfers['age'] = 'Aging'
    
    data = at.ProjectData.new(F,np.arange(data_years[0],data_years[1]),pops=all_pops,transfers=transfers)
    
    for region in regions:
        str_ext = country + '_' + region
        data.save(project_folder + 'empty_databook%s.xlsx'%str_ext)
    
    
    
if 'setup_progbooks' in todo:
    program_years = [2015,2018]
    progs = sc.odict()
    progs['Test']    = 'Testing\nDépistage'
    progs['Treat']   = 'Treatment\nPrise en charge'
    progs['BCC']     = 'Behaviour and change communication (BCC)\nCCC'
    progs['LLIN_m']  = 'Long lasting insecticide-treated nets (LLINs) - mass distribution\nMILDA'
    progs['LLIN_r']  = 'Long lasting insecticide-treated nets (LLINs) - routine via ANCs\nMILDA'
    progs['LLIN_c']  = 'Long lasting insecticide-treated nets (LLINs) - continuous general\nMILDA'
    progs['IPTp']    = 'Intermittent preventive therapy for pregnant women (IPTp)\nTPI'
    progs['IRS']     = 'Indoor residual spraying (IRS)\nAID'
    
    databook_paths = [file for file in glob(project_folder+'*databook*.xlsx')  if not '~' in file] #all databooks in the folder
#    databook_path  = [path for path in databook_paths if book_key in path][-1] #just the most recent with the right key
    D = at.ProjectData.from_spreadsheet(databook_paths[-1],framework=F)
    
    pset = at.ProgramSet.new(tvec=np.arange(program_years[0],program_years[1]),progs=progs,framework=F,data=D)
    for region in regions:
        str_ext = country + '_' + region
        pset.save(project_folder + 'empty_progbook%s.xlsx'%str_ext)

if 'setup_docs' in todo:
    at.generate_framework_doc(framework = F, fname = project_folder +'malaria_documentation.markdown', databook_only=True)
