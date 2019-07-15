# simple_dispatch

Simple Dispatch method for modeling the U.S. power sector and simulating marginal emissions factors (MEFs)

Published in Environmental Science and Technology:

last edited: 2019-07-15


#Files included:

#Python files:

simple_dispatch.py #runs the simple_dispatch program

mefs_from_dispatch.py #calculates and plots hourly marginal emissions factors from the dispatch data

reducedOrderDispatchModelForSimulatingMarginalEmissionsFactors.py #runs the code needed to generate the figures in the Environmental Science and Technology publication

#Data files:
    # ferc 714 data from here: https://www.ferc.gov/docs-filing/forms/form-714/data.asp
    # ferc 714 ids available on the simple_dispatch github repository
    # egrid data from here: https://www.epa.gov/energy/emissions-generation-resource-integrated-database-egrid
    # eia 923 data from here: https://www.eia.gov/electricity/data/eia923/
    # cems data from here: ftp://newftp.epa.gov/DmDnLoad/emissions/hourly/monthly/
    # easiur data from here: https://barney.ce.cmu.edu/~jinhyok/easiur/online/
    # fuel_default_prices.xlsx compiled from data from https://www.eia.gov/

#Python objects
    # generator_data_short_NERC_YEAR.obj objects contain cleaned generator data created by the simple_dispatch.generatorData class

