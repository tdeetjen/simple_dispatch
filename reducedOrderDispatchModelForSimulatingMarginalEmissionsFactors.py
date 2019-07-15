# -*- coding: utf-8 -*-
"""
Code for running the analysis and producing the figures for the journal article:

Deetjen, Thomas A.,  Azevedo, Ines L. (2019), A reduced-order dispatch model for simulating marginal emissions factors for the U.S. power sector. Environmental Science & Technology.

Requires custom modules "simple_dispatch" and "mefs_from_simple_dispatch", available at https://github.com/tdeetjen/simple_dispatch

"""


import pickle
import scipy
import os.path
import pandas
from simple_dispatch import generatorData
from simple_dispatch import bidStack
from simple_dispatch import dispatch
from  mefs_from_simple_dispatch import generateMefs
from  mefs_from_simple_dispatch import plotDispatch


if __name__ == '__main__':
    run_year = 2017
    #input variables. Right now the github only has 2017 data on it.
    #specific the location of the data directories
    # ferc 714 data from here: https://www.ferc.gov/docs-filing/forms/form-714/data.asp
    # ferc 714 ids available on the simple_dispatch github repository
    # egrid data from here: https://www.epa.gov/energy/emissions-generation-resource-integrated-database-egrid
    # eia 923 data from here: https://www.eia.gov/electricity/data/eia923/
    # cems data from here: ftp://newftp.epa.gov/DmDnLoad/emissions/hourly/monthly/
    # easiur data from here: https://barney.ce.cmu.edu/~jinhyok/easiur/online/
    # fuel_default_prices.xlsx compiled from data from https://www.eia.gov/
    ferc714_part2_schedule6_csv = 'Part 2 Schedule 6 - Balancing Authority Hourly System Lambda.csv'
    ferc714IDs_csv='Respondent IDs.csv'
    cems_folder_path ='C:\\Users\\tdeet\\Documents\\data\\raw\\epa\\CEMS'
    easiur_csv_path ='egrid_2016_plant_easiur.csv'
    fuel_commodity_prices_xlsx = 'fuel_default_prices.xlsx'
    if run_year == 2017:
        egrid_data_xlsx = 'egrid2016_data.xlsx'
        eia923_schedule5_xlsx = 'EIA923_Schedules_2_3_4_5_M_12_2017_Final_Revision.xlsx'
    if run_year == 2016:
        egrid_data_xlsx = 'egrid2016_data.xlsx'
        eia923_schedule5_xlsx = 'EIA923_Schedules_2_3_4_5_M_12_2016_Final_Revision.xlsx'
    if run_year == 2015:
        egrid_data_xlsx = 'egrid2014_data.xlsx'
        eia923_schedule5_xlsx = 'EIA923_Schedules_2_3_4_5_M_12_2015_Final_Revision.xlsx'
    if run_year == 2014:
        egrid_data_xlsx = 'egrid2014_data.xlsx'
        eia923_schedule5_xlsx = 'EIA923_Schedules_2_3_4_5_M_12_2014_Final_Revision.xlsx'   
    #loop through nerc regions
    #for nerc_region in ['TRE']:
    for nerc_region in ['TRE', 'MRO', 'WECC', 'SPP', 'SERC', 'RFC', 'FRCC', 'NPCC']:
        try:
            #if you've already run generatorData before, there will be a shortened pickled dictionary that we can just load in now. The 2017 pickled dictionaries can be downloaded from the simple_dispatch github repository. You can also download cems data and compile them using the generatorData object
            gd_short = pickle.load(open('generator_data_short_%s_%s.obj'%(nerc_region, str(run_year)), 'r'))
        except:
            #run the generator data object
            gd = generatorData(nerc_region, egrid_fname=egrid_data_xlsx, eia923_fname=eia923_schedule5_xlsx, ferc714IDs_fname=ferc714IDs_csv, ferc714_fname=ferc714_part2_schedule6_csv, cems_folder=cems_folder_path, easiur_fname=easiur_csv_path, include_easiur_damages=True, year=run_year, fuel_commodity_prices_excel_dir=fuel_commodity_prices_xlsx, hist_downtime=False, coal_min_downtime = 12, cems_validation_run=False)   
            #pickle the trimmed version of the generator data object
            gd_short = {'year': gd.year, 'nerc': gd.nerc, 'hist_dispatch': gd.hist_dispatch, 'demand_data': gd.demand_data, 'mdt_coal_events': gd.mdt_coal_events, 'df': gd.df}
            pickle.dump(gd_short, open('generator_data_short_%s_%s.obj'%(nerc_region, str(run_year)), 'w'))
        #now that we have the generator data cleaned up, we can build the merit order and run the dispatch
        #we can add a co2 price to the dispatch calculation
        for co2_dol_per_ton in [0]:
            #run the bidStack object - use information about the generators (from gd_short) to create a merit order (bid stack) of the nerc region's generators
            bs = bidStack(gd_short, co2_dol_per_kg=(co2_dol_per_ton / 907.185), time=30, dropNucHydroGeo=True, include_min_output=False, mdt_weight=0.5) #NOTE: set dropNucHydroGeo to True if working with data that only looks at fossil fuels (e.g. CEMS)
            #produce bid stack plots
            #bid_stack_cost = bs.plotBidStackMultiColor('gen_cost', plot_type='bar', fig_dim = (4,4), production_cost_only=True) #plot the merit order
            bid_stack_cost = bs.plotBidStackMultiColor('gen_cost', plot_type='bar', fig_dim = (4,4), production_cost_only=False) #plot the merit order
            bid_stack_co2 = bs.plotBidStackMultiColor('co2', plot_type='bar') #plot the merit order                 
            #run the dispatch object - use the nerc region's merit order (bs), a demand timeseries (gd.demand_data), and a time array (default is array([ 1,  2, ... , 51, 52]) for 52 weeks to run a whole year)
            #if you've already run and saved the dispatch, skip this step
            if not os.path.exists('simple_dispatch_%s_%s_%sco2price.csv'%(nerc_region, str(run_year), str(co2_dol_per_ton))):
                #run the dispatch object
                dp = dispatch(bs, gd_short.demand_data, time_array=scipy.arange(52)+1) #set up the object
                #dp = dispatch(bs, gd.demand_data, time_array=scipy.arange(3)+1) #test run          
                dp.calcDispatchAll() #function that solves the dispatch for each time period in time_array (default for each week of the year)
                #save dispatch results 
                dp.df.to_csv('simple_dispatch_%s_%s_%sco2price.csv'%(nerc_region, str(run_year), str(co2_dol_per_ton)), index=False)
    #now that the dispatch is run, we can calculate the marginal emissions factors and plot them            
    cedm_mefs_df = pandas.read_csv('mefs_by_decile_nerc.csv')[['year', 'region', 'dec', 'pollutant', 'factor']] #from CEDM: https://cedm.shinyapps.io/MarginalFactors/
    run_year = 2017
    #empty dataframe to hold error calculations
    error_main_df = pandas.DataFrame(columns=(['nerc', 'variable', 'co2_tot_hour_vs_rolling', 'co2_slope_sim', 'co2_slope_cedm', 'so2_tot_hour_vs_rolling', 'so2_slope_sim', 'so2_slope_cedm', 'nox_tot_hour_vs_rolling', 'nox_slope_sim', 'nox_slope_cedm']))
    #loop through the nerc regions and generate plots
    for nerc_region in ['TRE']:
    #for nerc_region in ['FRCC', 'TRE', 'WECC', 'SPP', 'MRO', 'SERC', 'RFC', 'NPCC']:
        for sim_co2_price in [0]:
            fig_suffix = '%s_%s_%sco2price.png'%(nerc_region, str(run_year), str(sim_co2_price))
            gd_short = pickle.load(open('generator_data_short_%s_%s.obj'%(nerc_region, str(run_year)), 'r'))
            #historical CEMS dispatch data
            dispatch_CEMS = gd_short['hist_dispatch']
            dispatch_CEMS.datetime = pandas.to_datetime(dispatch_CEMS.datetime) #put the datetime column into the correct type
            dispatch_CEMS.gen_cost_marg = scipy.minimum(150, dispatch_CEMS.gen_cost_marg) #remove prices larger than 150 $/MWh. Consider these outliers. For example, with 2014 TRE, R-squared between the actual data and simulated data is 0.16 while R-squared after capping the historical data at 150 $/MWh is 0.53. 150 $/MWh is something above the 99th percentile in most cases, so we are now saying that the simulation has a fit of 0.53 for 99+% of the historical data.
            gm_cems = generateMefs(dispatch_CEMS) 
            dispatch_CEMS_gm = gm_cems.df.copy(deep=True) #the script now calcultes hourly MEFs, which we want for some plotting below
            
            #dispatch solution
            dispatch_solution = pandas.read_csv('simple_dispatch_%s_%s_%sco2price.csv'%(nerc_region, str(run_year), str(sim_co2_price)))
            #dispatch_solution = pandas.read_csv('C:\\Users\\tdeet\\Documents\\analysis\\thirdParty\\PLEXOS\\2014_TRE_hourly_demand_and_fuelmix_PLEXOS_w_CEMS_demand.csv')
            gm_sim = generateMefs(dispatch_solution[list(dispatch_CEMS.columns)]) 
            #replace marginal unit informatoin with MEF calculations
            dispatch_solution[['co2_mef', 'so2_mef', 'nox_mef']] = gm_sim.df[['co2_marg', 'so2_marg', 'nox_marg']]
            
            #data below from CEDM at https://cedm.shinyapps.io/MarginalFactors/
            deciles = list(cedm_mefs_df[(cedm_mefs_df.region==nerc_region) & (cedm_mefs_df.year==run_year) & (cedm_mefs_df.pollutant=='co2')].dec.values)
            mefs_cedm_co2 = list(cedm_mefs_df[(cedm_mefs_df.region==nerc_region) & (cedm_mefs_df.year==run_year) & (cedm_mefs_df.pollutant=='co2')].factor.values)
            mefs_cedm_so2 = list(cedm_mefs_df[(cedm_mefs_df.region==nerc_region) & (cedm_mefs_df.year==run_year) & (cedm_mefs_df.pollutant=='so2')].factor.values)
            mefs_cedm_nox = list(cedm_mefs_df[(cedm_mefs_df.region==nerc_region) & (cedm_mefs_df.year==run_year) & (cedm_mefs_df.pollutant=='nox')].factor.values)
            
            #create the plotDispatch object
            pd = plotDispatch(nerc_region, dispatch_solution, dispatch_CEMS_gm, deciles, mefs_cedm_co2, mefs_cedm_so2, mefs_cedm_nox)
            #added error calculations to the error dataframe
            error_main_df = pandas.concat([error_main_df, pd.calcError(run_year)], axis=0)
            #create some figures
            fDemandEmissionsTotal = pd.plotDemandEmissions('total', figure_dimensions=(5,5))
            fDemandEmissionsTotal.savefig('fDemandEmissionsTotal' + fig_suffix, dpi=500, bbox_inches='tight')
            fDemandEmissionsMarginal = pd.plotDemandEmissions('marginal', figure_dimensions=(5,5))
            fDemandEmissionsMarginal.savefig('fDemandEmissionsMarginal' + fig_suffix, dpi=500, bbox_inches='tight')
            fDemandPrice = pd.plotDemandPrices(figure_dimensions=(5,5))
            
            #density function plots
            #price and emissions
            for s in ['WholeYear']:
                #set dates
                if s == 'WholeYear':
                    df_start = '2017-01-01 00:00:00'
                    df_end = '2018-01-01 00:00:00'
                if s == 'Summer':
                    df_start = '2014-03-01 00:00:00'
                    df_end = '2014-06-01 00:00:00'
                if s == 'Spring':
                    df_start = '2014-06-01 00:00:00'
                    df_end = '2014-09-01 00:00:00'
                #run plot and save scripts
                price_df = pd.plot_density_function('cumulative', 'data', 'gen_cost_marg', start_date=df_start, end_date=df_end, x_range=[0,100])
                price_df_err = pd.plot_density_function('probability', 'error', 'gen_cost_marg', start_date=df_start, end_date=df_end, bin_no=50, error_range=[-1,1])
                co2_df = pd.plot_density_function('cumulative', 'data', 'co2_tot', start_date=df_start, end_date=df_end)
                co2_df_err = pd.plot_density_function('probability', 'error', 'co2_tot', start_date=df_start, end_date=df_end, bin_no=50, error_range=[-1,1])   
                so2_df = pd.plot_density_function('cumulative', 'data', 'so2_tot', start_date=df_start, end_date=df_end, bin_no=50, x_range=[-scipy.inf, scipy.inf])
                so2_df_err = pd.plot_density_function('probability', 'error', 'so2_tot', start_date=df_start, end_date=df_end, bin_no=50, error_range=[-1,1])
                nox_df = pd.plot_density_function('cumulative', 'data', 'nox_tot', start_date=df_start, end_date=df_end)
                nox_df_err = pd.plot_density_function('probability', 'error', 'nox_tot', start_date=df_start, end_date=df_end, bin_no=50, error_range=[-1,1])    
                #mefs
                co2_mefs_df = pd.plot_density_function('cumulative', 'data', 'co2_marg', start_date=df_start, end_date=df_end, x_range=[0,1500])
                co2_mefs_df_err = pd.plot_density_function('probability', 'error', 'co2_marg', start_date=df_start, end_date=df_end, bin_no=50, x_range=[0,1500/0.454], error_range=[-1000,1000])
                so2_mefs_df = pd.plot_density_function('cumulative', 'data', 'so2_marg', start_date=df_start, end_date=df_end, x_range=[0,3/0.454])
                so2_mefs_df_err = pd.plot_density_function('probability', 'error', 'so2_marg', start_date=df_start, end_date=df_end, bin_no=50, x_range=[0,3/0.454], error_range=[-4,4])
                nox_mefs_df = pd.plot_density_function('cumulative', 'data', 'nox_marg', start_date=df_start, end_date=df_end, x_range=[0,2/0.454])
                nox_mefs_df_err = pd.plot_density_function('probability', 'error', 'nox_marg', start_date=df_start, end_date=df_end, bin_no=50, x_range=[0,2/0.454], error_range=[-3,3])




