# simple_dispatch

Simple Dispatch method for modeling the U.S. power sector and simulating marginal emissions factors (MEFs)

Download the files from github to run the 2017 FRCC simulation. Running simulations for additional regions will require downloading more CEMS data (see cems/ description below)

last edited: 2018-10-02

#Files included:

#Python files:

cems/ #folder containing CEMS data (one .csv per state per month). Currently, it only has Florida files (enough data to run the 2017 "FRCC" nerc region simulation. Download additional files for other states at: ftp://newftp.epa.gov/DmDnLoad/emissions/hourly/monthly/

simple_dispatch.py #runs the simple_dispatch program

mefs_from_dispatch.py #calculates hourly marginal emissions factors from historical data. Plots the simple_dispatch results

#Data files:

egrid2016_data.xlsx #see https://www.epa.gov/energy/emissions-generation-resource-integrated-database-egrid

EIA923_Schedules_2_3_4_5_M_12_2017_Early_Release.xlsx #see https://www.eia.gov/electricity/data/eia923/

Part 2 Schedule 6 - Balancing Authority Hourly System Lambda.csv #see https://www.ferc.gov/docs-filing/forms/form-714/data.asp

Respondent IDs.csv #used to cateogorize the respondents from Part 2 Schedule 6 by their appropriate NERC region

fuel_default_prices.xlsx #compiled by author, based on https://www.eia.gov/coal/markets/ , https://www.eia.gov/dnav/ng/hist/rngwhhdm.htm , https://www.eia.gov/dnav/pet/hist/RWTCD.htm , https://www.eia.gov/nuclear/

#Results examples:

'%s_%s_hourly_demand_and_fuelmix.csv'%(str(run_year), nerc_region)) #the historical dispatch: output of generatorData.hist_dispatch for a nerc region and year

'dispatch_output_weekly_%s_%s.csv'%(nerc_region, str(run_year)) #the simulated dispatch: output of dispatch.df for a nerc region and year
