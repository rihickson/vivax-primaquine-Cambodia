# -*- coding: utf-8 -*-
"""

@author: RIH
"""

if __name__ == '__main__':
    import sciris as sc
    import gc
    from run_cambodia import run_me
    
    #split into two for memory constraints
    book_key_all = ['Pursat_low', 'Pursat_high', 'Mondul_Kiri_low', 'Mondul_Kiri_high', 'Kampong_Chhnang_low', 'Kampong_Chhnang_high', 'Battambang_low', 'Battambang_high', 'Takeo_low', 'Takeo_high', 'Pailin_low', 'Pailin_high']
    book_batches = 6 #adjust to run in batches due to memory constraints for long runs - 0 means don't (re-)run
    # book_key_all = ['Pailin_low', 'Pailin_high']
    
    global book_key
    
    # run the model for all provinces, for low and high baseline incidences
    # for i in range(book_batches):
    #     num_per_batch = int(len(book_key_all)/book_batches)
    #     book_keys = book_key_all[i * num_per_batch:(i+1)*num_per_batch]
    #     print (f'Running batch {i+1}: {book_keys}')
    #     sc.parallelize(run_me, book_keys)
    #     gc.collect() #clear memory
        
    for prov in book_key_all:
        run_me(prov)
        gc.collect() #clear memory
    
    # aggregate and plot the results (shown in paper)
    from result_aggregation import run_post_processing
    
    run_post_processing(flag_make_plots = True, flag_elimination_comparison = True)

