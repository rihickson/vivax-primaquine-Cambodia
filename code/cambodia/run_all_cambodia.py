# -*- coding: utf-8 -*-
"""

@author: RIH
"""

from run_cambodia import run_me

book_key_all = ['Pursat_low', 'Pursat_high', 'Mondul_Kiri_low', 'Mondul_Kiri_high', 'Kampong_Chhnang_low', 'Kampong_Chhnang_high', 'Battambang_low', 'Battambang_high', 'Takeo_low', 'Takeo_high', 'Pailin_low', 'Pailin_high']

# book_key_all = ['Pailin_low']

global book_key


# run the model for all provinces, for low and high baseline incidences
# with multiprocessing.Pool() as pool:
#     pool.map(run_me, book_key_all)
for prov in book_key_all:
    run_me(prov)


# aggregate and plot the results (shown in paper)
from result_aggregation import run_post_processing

run_post_processing(flag_all_prevalences = True, flag_resort_data = False, flag_make_plots = True)

