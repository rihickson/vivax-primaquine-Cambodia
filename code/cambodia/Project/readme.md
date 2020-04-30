# Cambodia vivax only malaria project to explore proposed Primaquine usage for adult males

`malaria_frame_work_20191028.xlsx`: Model framework provided by Rowan Martin-Hughes.

`malaria_framework_20191216_unused_senstive_mosquitoes`: Updated model framework with corrected mosquito lifespan (30 days as opposed to 30 years).

`Databook_cambodia_Pursat_uncalibrated.xlsx`: First attempt at data and parameter entry for Pursat, kept for some default values from Rowan.

`Databook_cambodia_<province>_low.xlsx`: Calibrated model for <province>, assuming the incidence for 2020-2025 will be more inline with the lower data values (e.g. for Pursat, those from before 2018), and (slightly) decreasing or flat. This is to provide a "best case scenario" for the baseline incidence from which to attempt to eliminate vivax malaria. 

`Databook_cambodia_<province>_high.xlsx`: Calibrated model for <province>, assuming the incidence for 2020-2025 will be more inline with the higher incidence value(s), and increasing. This is a "worst case scenario", where elimination of vivax is less likely because of a higher baseline incidence.

NB: The most likely scenario is that the true incidence starts high, stays fairly constant, then decreased when interventions were ramped up - and may now have stabilised or be slightly increasing from the higher incidence values. E.g. for Pursat. But these represent "best" (low) and "worst" (high) case scenarios from a control standpoint.

`<province>` = 	Pursat
				Mondul Kiri
				Kampong Chhnang
				Battambang
				Takeo
				Pailin

`Optima_malaria_results.pptx`: Report/results for Cambodia

`result_aggregation.py`: script to take the raw results from each of the 12 scenarios (6 provinces, high or low baseline each) and write them out to a csv with prevalence from the simulated scenarios with and without primaquine, for use by David Price in estimating probability of detecting at least one case for numbers of tests conducted. Relies on my populations models from `~/repos/comixing-malaria/data/Cambodia-Pengby/cambodia_population_estimates.xlsx`, which are summarised in `cambodia_population_estimates.csv` in this directory. Outputs a file to `../aggregated_results.csv` where the scenario key is `<population (`M15+` or `Gen`)>_<baseline incidence scenario (high or low)>_<intervention scenario (`no_prim` or `prim`)>```

`cambodia_population_estimates.csv`: summary of my populations models from `~/repos/comixing-malaria/data/Cambodia-Pengby/cambodia_population_estimates.xlsx`, used in `result_aggregation.py` to calculate prevalence.

