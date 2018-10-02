# simple_dispatch
# Thomas Deetjen
# v_21
# last edited: 2018-10-02
# class "generatorData" turns CEMS, eGrid, FERC, and EIA data into a cleaned up dataframe for feeding into a "bidStack" object
# class "bidStack" creates a merit order curve from the generator fleet information created by the "generatorData" class
# class "dispatch" uses the "bidStack" object to choose which power plants should be operating during each time period to meet a demand time-series input
# ---
# v21:
# set up to simulate the 2017 historical dispatch for the NERC regional level

import pandas
import matplotlib.pylab
import scipy
import scipy.interpolate
import datetime
import math





class generatorData(object):
    def __init__(self, nerc, egrid_fname, eia923_fname, ferc714_fname='C:\\Users\\tdeet\\Documents\\data\\raw\\ferc\\ferc_714\\Part 2 Schedule 6 - Balancing Authority Hourly System Lambda.csv', ferc714IDs_fname='C:\\Users\\tdeet\\Documents\\data\\raw\\ferc\\ferc_714\\Respondent IDs.csv', cems_folder='C:\\Users\\tdeet\\Documents\\data\\raw\\epa\\CEMS', year=2017, fuel_commodity_prices_excel_dir = 'C:\\Users\\tdeet\\Documents\\data\\processed\\eiaFuelPrices\\fuel_default_prices.xlsx'):
        """ 
        Translates the CEMS, eGrid, FERC, and EIA data into a dataframe for feeding into the bidStack class
        ---
        nerc : nerc region of interest (e.g. 'TRE', 'MRO', etc.)
        egrid_fname : a .xlsx file name for the eGrid generator data
        eia923_fname : filename of eia form 923
        ferc714_fname : filename of nerc form 714 hourly system lambda 
        ferc714IDs_fname : filename that matches nerc 714 respondent IDs with nerc regions
        year : year that we're looking at (e.g. 2017)
        fuel_commodity_prices_excel_dir : filename of national EIA fuel prices (in case EIA923 data is empty)
        """
        #read in the data. This is a bit slow right now because it reads in more data than needed, but it is simple and straightforward
        print 'Reading in data...'
        self.nerc = nerc
        egrid_year_str = str(math.floor((year / 2.0)) * 2)[2:4] #eGrid is only every other year so we have to use eGrid 2016 to help with a 2017 run, for example
        self.egrid_unt = pandas.read_excel(egrid_fname, 'UNT'+egrid_year_str, skiprows=[0]) 
        self.egrid_gen = pandas.read_excel(egrid_fname, 'GEN'+egrid_year_str, skiprows=[0])
        self.egrid_plnt = pandas.read_excel(egrid_fname, 'PLNT'+egrid_year_str, skiprows=[0])
        eia923 = pandas.read_excel(eia923_fname, 'Page 5 Fuel Receipts and Costs', skiprows=[0,1,2,3]) 
        eia923 = eia923.rename(columns={'Plant Id': 'orispl'})
        self.eia923 = eia923
        self.ferc714 = pandas.read_csv(ferc714_fname)
        self.ferc714_ids = pandas.read_csv(ferc714IDs_fname)
        self.cems_folder = cems_folder
        self.fuel_commodity_prices = pandas.read_excel(fuel_commodity_prices_excel_dir, str(year))
        self.year = year
        self.cleanGeneratorData()
        self.calcFuelPrices()
        self.addGenMinOut()
        self.addGenVom()
        self.calcDemandData()     
        self.addElecPriceToDemandData()
        self.demandTimeSeries()
        

    def cleanGeneratorData(self):
        """ 
        Converts the eGrid and CEMS data into a dataframe usable by the bidStack class.
        ---
        Creates
        self.df : has 1 row per generator unit or plant. columns describe emissions, heat rate, capacity, fuel, grid region, etc. This dataframe will be used to describe the generator fleet and merit order.
        self.df_cems : has 1 row per hour of the year per generator unit or plant. columns describe energy generated, emissions, and grid region. This dataframe will be used to describe the historical hourly demand, dispatch, and emissions
        """
        #copy in the egrid data and merge it together. In the next few lines we use the eGRID excel file to bring in unit level data for fuel consumption and emissions, generator level data for capacity and generation, and plant level data for fuel type and grid region. Then we compile it together to get an initial selection of data that defines each generator.
        print 'Compiling eGRID Data...'
        #unit-level data
        df = self.egrid_unt.copy(deep=True)
        #rename columns
        df = df[['PNAME', 'ORISPL', 'UNITID', 'PRMVR', 'FUELU1', 'HTIAN', 'NOXAN', 'SO2AN', 'CO2AN', 'HRSOP']]
        df.columns = ['gen', 'orispl', 'unit', 'prime_mover', 'fuel', 'mmbtu_ann', 'nox_ann', 'so2_ann', 'co2_ann', 'hours_on']
        df['orispl_unit'] = df.orispl.map(str) + '_' + df.unit.map(str) #orispl_unit is a unique tag for each generator unit
        #drop nan fuel
        df = df[~df.fuel.isna()]
        #gen-level data: contains MW capacity and MWh annual generation data
        df_gen = self.egrid_gen.copy(deep=True) 
        df_gen['orispl_unit'] = df_gen['ORISPL'].map(str) + '_' + df_gen['GENID'].map(str) #orispl_unit is a unique tag for each generator unit
        df_gen = df_gen[['NAMEPCAP', 'GENNTAN', 'GENYRONL', 'orispl_unit']]
        df_gen.columns = ['mw', 'mwh_ann', 'year_online', 'orispl_unit']
        #plant-level data: contains fuel, fuel_type, balancing authority, nerc region, and egrid subregion data
        df_plnt = self.egrid_plnt.copy(deep=True) 
        #fuel
        df_plnt_fuel = df_plnt[['PLPRMFL', 'PLFUELCT']]
        df_plnt_fuel = df_plnt_fuel.drop_duplicates('PLPRMFL')
        df_plnt_fuel.PLFUELCT = df_plnt_fuel.PLFUELCT.str.lower()
        df_plnt_fuel.columns = ['fuel', 'fuel_type']
        #geography
        df_plnt = df_plnt[['ORISPL', 'BACODE', 'NERC', 'SUBRGN']]
        df_plnt.columns = ['orispl', 'ba', 'nerc', 'egrid']
        #merge these egrid data together at the unit-level
        df = df.merge(df_gen, left_index=True, how='left', on='orispl_unit')  
        df = df.merge(df_plnt, left_index=True, how='left', on='orispl')  
        df = df.merge(df_plnt_fuel, left_index=True, how='left', on='fuel')  
        #keep only the units in the nerc region we're analyzing
        df = df[df.nerc == self.nerc]
        #calculate the emissions rates
        df['co2'] = scipy.divide(df.co2_ann,df.mwh_ann) * 907.185 #tons to kg
        df['so2'] = scipy.divide(df.so2_ann,df.mwh_ann) * 907.185 #tons to kg
        df['nox'] = scipy.divide(df.nox_ann,df.mwh_ann) * 907.185 #tons to kg
        ###
        #now sort through and compile CEMS data. The goal is to use CEMS data to characterize each generator unit. So if CEMS has enough information to describe a generator unit we will over-write the eGRID data. If not, we will use the eGRID data instead. (CEMS data is expected to be more accurate because it has actual hourly performance of the generator units that we can use to calculate their operational characteristics. eGRID is reported on an annual basis and might be averaged out in different ways than we would prefer.)
        print 'Compiling CEMS data...'
        #dictionary of which states are in which nerc region (b/c CEMS file downloads have the state in the filename)
        states = {'FRCC': ['fl'], 
                  'WECC': ['ca','or','wa','mn','id','wy','ut','co','az','nm','tx'],
                  'SPP' : ['nm','ks','tx','ok','la','ar','mo'],
                  'RFC' : ['wi','mi','il','in','oh','ky','wv','vi','md','pe','nj'],
                  'NPCC' : ['ny','ct','de','ri','ma','vt','nh','me'],
                  'SERC' : ['mo','ar','tx','la','ms','tn','ky','il','vi','al','fl','ga','sc','nc'],
                  'MRO': ['ia','il','mi','mn','mo','mt','nd','ne','sd','wi','wy'], 
                  'TRE': ['ok','tx']}
        #compile the different months of CEMS files into one dataframe, df_cems. (CEMS data is downloaded by state and by month, so compiling a year of data for ERCOT / TRE, for example, requires reading in 12 Texas .csv files and 12 Oklahoma .csv files)   
        df_cems = pandas.DataFrame()
        for s in states[self.nerc]:
            for m in ['01','02','03','04','05','06','07','08','09','10','11', '12']:
                print s + ': ' + m
                df_cems_add = pandas.read_csv(self.cems_folder + '\\%s\\%s%s%s.csv'%(str(self.year),str(self.year),s,m))
                df_cems_add = df_cems_add[['ORISPL_CODE', 'UNITID', 'OP_DATE','OP_HOUR','GLOAD (MW)', 'SO2_MASS (lbs)', 'NOX_MASS (lbs)', 'CO2_MASS (tons)', 'HEAT_INPUT (mmBtu)']].dropna()
                df_cems_add.columns=['orispl', 'unit', 'date','hour','mwh', 'so2_tot', 'nox_tot', 'co2_tot', 'mmbtu']
                df_cems = pandas.concat([df_cems, df_cems_add])
        #create the 'orispl_unit' column, which combines orispl and unit into a unique tag for each generation unit
        df_cems['orispl_unit'] = df_cems['orispl'].map(str) + '_' + df_cems['unit'].map(str)
        #bring in geography data and only keep generators within self.nerc
        df_cems = df_cems.merge(df_plnt, left_index=True, how='left', on='orispl')  
        df_cems = df_cems[df_cems['nerc']==self.nerc] 
        #convert emissions to kg
        df_cems.co2_tot = df_cems.co2_tot * 907.185 #tons to kg
        df_cems.so2_tot = df_cems.so2_tot * 0.454 #lbs to kg
        df_cems.nox_tot = df_cems.nox_tot * 0.454 #lbs to kg
        #calculate the hourly heat and emissions rates. Later we will take the medians over each week to define the generators weekly heat and emissions rates.
        df_cems['heat_rate'] = df_cems.mmbtu / df_cems.mwh
        df_cems['co2'] = df_cems.co2_tot / df_cems.mwh
        df_cems['so2'] = df_cems.so2_tot / df_cems.mwh
        df_cems['nox'] = df_cems.nox_tot / df_cems.mwh
        df_cems.replace([scipy.inf, -scipy.inf], scipy.nan, inplace=True) #don't want inf messing up median calculations
        #drop any bogus data. For example, the smallest mmbtu we would expect to see is 25MW(smallest unit) * 0.4(smallest minimum output) * 6.0 (smallest heat rate) = 60 mmbtu. Any entries with less than 60 mmbtu fuel or less than 6.0 heat rate, let's get rid of that row of data.
        df_cems = df_cems[(df_cems.heat_rate >= 6.0) & (df_cems.mmbtu >= 60)]
        #calculate emissions rates and heat rate for each week and each generator
        #rather than parsing the dates (which takes forever because this is such a big dataframe) we can create month and day columns for slicing the data based on time of year
        df_orispl_unit = df_cems.copy(deep=True)
        df_orispl_unit.date = df_orispl_unit.date.str.replace('/','-')
        temp = pandas.DataFrame(df_orispl_unit.date.str.split('-').tolist(), columns=['month', 'day', 'year'], index=df_orispl_unit.index).astype(float)
        df_orispl_unit['monthday'] = temp.year*10000 + temp.month*100 + temp.day
        ###
        #loop through the weeks, slice the data, and find the average heat rates and emissions rates
        #first, add a column 't' that says which week of the simulation we are in
        df_orispl_unit['t'] = 0
        for t in scipy.arange(52)+1:
            start = (datetime.datetime.strptime(str(self.year) + '-01-01', '%Y-%m-%d') + datetime.timedelta(days=7.05*(t-1)-1)).strftime('%Y-%m-%d') 
            end = (datetime.datetime.strptime(str(self.year) + '-01-01', '%Y-%m-%d') + datetime.timedelta(days=7.05*(t)-1)).strftime('%Y-%m-%d') 
            start_monthday = float(start[0:4])*10000 + float(start[5:7])*100 + float(start[8:])
            end_monthday = float(end[0:4])*10000 + float(end[5:7])*100 + float(end[8:])
            #slice the data for the days corresponding to the time series period, t
            df_orispl_unit.loc[(df_orispl_unit.monthday >= start_monthday) & (df_orispl_unit.monthday < end_monthday), 't'] = t
        #remove outlier emissions and heat rates. These happen at hours where a generator's output is very low (e.g. less than 10 MWh). To remove these, we will remove any datapoints where mwh < 10.0 and heat_rate < 30.0 (0.5% percentiles of the 2014 TRE data).
        df_orispl_unit = df_orispl_unit[(df_orispl_unit.mwh >= 10.0) & (df_orispl_unit.heat_rate <= 30.0)]    
        #aggregate by orispl_unit and t to get the heat rate, emissions rates, and capacity for each unit at each t
        temp_2 = df_orispl_unit.groupby(['orispl_unit', 't'], as_index=False).agg('median')[['orispl_unit', 't', 'heat_rate', 'co2', 'so2', 'nox']].copy(deep=True)
        temp_2['mw'] = df_orispl_unit.groupby(['orispl_unit', 't'], as_index=False).agg('max')['mwh'].copy(deep=True)
        #condense df_orispl_unit down to where we just have 1 row for each unique orispl_unit
        df_orispl_unit = df_orispl_unit.groupby('orispl_unit', as_index=False).agg('max')[['orispl_unit', 'orispl', 'ba', 'nerc', 'egrid', 'mwh']]
        df_orispl_unit.rename(columns={'mwh':'mw'}, inplace=True)
        for c in ['heat_rate', 'co2', 'so2', 'nox', 'mw']:
            temp_3 = temp_2.set_index(['orispl_unit', 't'])[c].unstack().reset_index()
            temp_3.columns = list(['orispl_unit']) + ([c + str(a) for a in scipy.arange(52)+1])
            #for any nan values (assuming these are offline generators without any ouptut data), fill nans with a large heat_rate that will move the generator towards the end of the merit order and large-ish emissions rate, so if they generator is dispatched in the model it will jack up prices but emissions won't be heavily affected (note, previously I just replaced all nans with 99999, but I was concerned that this might lead to a few hours of the year with extremely high emissions numbers that threw off the data)
            M = float(scipy.where(c=='heat_rate', 50.0, scipy.where(c=='co2', 1500.0, scipy.where(c=='so2', 4.0, scipy.where(c=='nox', 3.0, scipy.where(c=='mw', 0.0, 99.0))))))
            temp_3 = temp_3.fillna(M)
            df_orispl_unit = df_orispl_unit.merge(temp_3, on='orispl_unit', how='left')
        #merge df_orispl_unit into df. Now we have a dataframe with weekly heat rate and emissions rates for any plants in CEMS with that data. There will be some nan values in df for those weekly columns (e.g. 'heat_rate1', 'co223', etc. that we will want to fill with annual averages from eGrid for now
        orispl_units_egrid = df.orispl_unit.unique()
        orispl_units_cems = df_orispl_unit.orispl_unit.unique()
        df_leftovers = df[df.orispl_unit.isin(scipy.setdiff1d(orispl_units_egrid, orispl_units_cems))]
        #remove any outliers - fuel is solar, wind, waste-heat, purchased steam, or other, less than 25MW capacity, less than 88 operating hours (1% CF), mw = nan, mmbtu = nan
        df_leftovers = df_leftovers[(df_leftovers.fuel!='SUN') & (df_leftovers.fuel!='WND') & (df_leftovers.fuel!='WH') & (df_leftovers.fuel!='OTH') & (df_leftovers.fuel!='PUR') & (df_leftovers.mw >=25) & (df_leftovers.hours_on >=88) & (~df_leftovers.mw.isna()) & (~df_leftovers.mmbtu_ann.isna())]
        #remove any outliers that have 0 emissions (except for nuclear)
        df_leftovers = df_leftovers[~((df_leftovers.fuel!='NUC') & (df_leftovers.nox_ann.isna()))]
        df_leftovers['cf'] = df_leftovers.mwh_ann / (df_leftovers.mw *8760.)
        #drop anything with capacity factor less than 1%
        df_leftovers = df_leftovers[df_leftovers.cf >= 0.01]
        df_leftovers.fillna(0.0, inplace=True)
        df_leftovers['heat_rate'] = df_leftovers.mmbtu_ann / df_leftovers.mwh_ann
        #add in the weekly time columns for heat rate and emissions rates. In this case we will just apply the annual average to each column, but we still need those columns to be able to concatenate back with df_orispl_unit and have our complete set of generator data
        for e in ['heat_rate', 'co2', 'so2', 'nox', 'mw']:
            for t in scipy.arange(52)+1:
                df_leftovers[e + str(t)] = df_leftovers[e]
        df_leftovers.drop(columns = ['gen', 'unit', 'prime_mover', 'fuel', 'mmbtu_ann', 'nox_ann', 'so2_ann', 'co2_ann', 'mwh_ann', 'fuel_type', 'co2', 'so2', 'nox', 'cf', 'heat_rate'], inplace=True)   
        #concat df_leftovers and df_orispl_unit
        df_orispl_unit = pandas.concat([df_orispl_unit, df_leftovers])     
        #use df to get prime_mover, fuel, and fuel_type for each orispl_unit
        df_fuel = df[df.orispl_unit.isin(df_orispl_unit.orispl_unit.unique())][['orispl_unit', 'fuel', 'fuel_type', 'prime_mover']]
        df_fuel.fuel = df_fuel.fuel.str.lower()
        df_fuel.fuel_type = df_fuel.fuel_type.str.lower()
        df_fuel.prime_mover = df_fuel.prime_mover.str.lower()
        df_orispl_unit = df_orispl_unit.merge(df_fuel, on='orispl_unit', how='left')
        #if we are using, for example, 2017 CEMS and 2016 eGrid, there may be some powerplants without fuel, fuel_type, and prime_mover data. Lets assums 'ng', 'gas', and 'ct' for these units based on trends on what was built in 2017
        df_orispl_unit.loc[df_orispl_unit.fuel.isna(), ['fuel', 'fuel_type']] = ['ng', 'gas']
        df_orispl_unit.loc[df_orispl_unit.prime_mover.isna(), 'prime_mover'] = 'ct'
        df_orispl_unit.fillna(0.0, inplace=True)
        #add in some columns to aid in calculating the fuel mix
        for f_type in ['gas', 'coal', 'oil', 'nuclear', 'hydro', 'geothermal', 'biomass']:
            df_orispl_unit['is_'+f_type.lower()] = (df_orispl_unit.fuel_type==f_type).astype(int)
        self.df_cems = df_cems
        self.df = df_orispl_unit


    def calcFuelPrices(self):
        """ 
        Finds the weekly fuel prices for each orispl at the nerc level. If that data is not available, generator units will be assigned the weekly fuel profile of a generator in that nerc region of the same fuel type for which data is available. If no data for a specific fuel type is available for a nerc region we will use the national level data. If the national level data does not contain data for a certain fuel type, we will use EIA fuel commodity trends.
        ---
        Adds one column for each week of the year to self.df that contain fuel prices for each generation unit
        """   
        print 'Adding fuel prices...'
        #we use eia923, where generators report their fuel purchases
        df = self.eia923.copy(deep=True)
        df = df[['YEAR','MONTH','orispl','ENERGY_SOURCE','FUEL_GROUP','QUANTITY','Average Heat\nContent','FUEL_COST']]
        df.columns = ['year', 'month', 'orispl' , 'fuel', 'fuel_type', 'quantity', 'heat_content', 'fuel_price']
        df.fuel = df.fuel.str.lower()       
        #clean up prices
        df = df[df.fuel_price!='.']
        df.fuel_price = df.fuel_price.astype('float')/100.
        df = df.reset_index()    
        #find unique monthly prices per orispl and fuel type
        #create empty dataframe to hold the results
        df2 = self.df.copy(deep=True)[['fuel','orispl','orispl_unit']]
        orispl_prices = pandas.DataFrame(columns=['orispl_unit', 'orispl', 'fuel', 1,2,3,4,5,6,7,8,9,10,11,12, 'quantity'])
        orispl_prices[['orispl_unit','orispl','fuel']] = df2[['orispl_unit', 'orispl', 'fuel']]
        #populate the results by looping through the orispl_units to see if they have EIA923 fuel price data
        for o_u in orispl_prices.orispl_unit.unique():
            #grab 'fuel' and 'orispl'
            f = orispl_prices.loc[orispl_prices.orispl_unit==o_u].iloc[0]['fuel']
            o = orispl_prices.loc[orispl_prices.orispl_unit==o_u].iloc[0]['orispl']
            #find the weighted average monthly fuel price matching 'f' and 'o'
            temp = df[(df.orispl==o) & (df.fuel==f)][['month', 'quantity', 'fuel_price']]
            if len(temp) != 0:
                temp['weighted'] = scipy.multiply(temp.quantity, temp.fuel_price)
                temp = temp.groupby(['month'], as_index=False).sum()[['month', 'quantity', 'weighted']]
                temp['fuel_price'] = scipy.divide(temp.weighted, temp.quantity)
                temp_prices = pandas.DataFrame({'month': scipy.arange(12)+1})
                temp_prices = temp_prices.merge(temp[['month', 'fuel_price']], on='month', how='left')
                temp_prices.loc[temp_prices.fuel_price.isna(), 'fuel_price'] = temp_prices.fuel_price.median()
                #add the monthly fuel prices into orispl_prices
                orispl_prices.loc[orispl_prices.orispl_unit==o_u, orispl_prices.columns.difference(['orispl_unit', 'orispl', 'fuel'])] = scipy.append(scipy.array(temp_prices.fuel_price),temp.quantity.sum())
        #for any fuels that we have nerc_region data, apply those monthly fuel price profiles to other generators with the same fuel type but that do not have unique EIA923 data
        for f in orispl_prices.dropna().fuel.unique():
            orispl_prices_filled = orispl_prices[orispl_prices.fuel==f].dropna().drop_duplicates(subset='orispl', keep='first').sort_values('quantity', ascending=0)
            orispl_prices_empty = orispl_prices[(orispl_prices.fuel==f) & (orispl_prices[1].isna())]
            loop = 0
            loop_len = len(orispl_prices_filled) - 1
            #of the plants with EIA923 data that we are assigning to plants without EIA923 data, we will use the plant with the highest energy production first, assigning its fuel price profile to one of the generators that does not have EIA923 data. We will move on to plant with the next highest energy production and so on, uniformly distributing the available EIA923 fuel price profiles to generators without fuel price data
            for o in orispl_prices_empty.orispl.unique():
                orispl_prices.loc[orispl_prices.orispl==o, orispl_prices.columns.difference(['orispl_unit', 'orispl', 'fuel', 'quantity'])] = scipy.array(orispl_prices_filled[orispl_prices.columns.difference(['orispl_unit', 'orispl', 'fuel', 'quantity'])].iloc[loop])
                #keep looping through the generators with EIA923 price data until we have used all of their fuel price profiles, then start again from the beginning of the loop with the plant with the highest energy production
                if loop < loop_len:
                    loop += 1
                else:
                    loop = 0
        #and now we still have some nan values for fuel types that had no nerc_region EIA923 data. We'll start with the national median for the EIA923 data.
        f_array = scipy.intersect1d(orispl_prices[orispl_prices[1].isna()].fuel.unique(), df.fuel.unique())
        for f in f_array: 
            temp = df[df.fuel==f][['month', 'quantity', 'fuel_price']]
            temp['weighted'] = scipy.multiply(temp.quantity, temp.fuel_price)
            temp = temp.groupby(['month'], as_index=False).sum()[['month', 'quantity', 'weighted']]
            temp['fuel_price'] = scipy.divide(temp.weighted, temp.quantity)
            temp_prices = pandas.DataFrame({'month': scipy.arange(12)+1})
            temp_prices = temp_prices.merge(temp[['month', 'fuel_price']], on='month', how='left')
            temp_prices.loc[temp_prices.fuel_price.isna(), 'fuel_price'] = temp_prices.fuel_price.median()
            orispl_prices.loc[orispl_prices.fuel==f, orispl_prices.columns.difference(['orispl_unit', 'orispl', 'fuel'])] = scipy.append(scipy.array(temp_prices.fuel_price),temp.quantity.sum())
        #for any fuels that don't have EIA923 data at all (for all regions) we will use commodity price approximations from an excel file
        #first we need to change orispl_prices from months to weeks
        orispl_prices.columns = ['orispl_unit', 'orispl', 'fuel', 1, 5, 9, 14, 18, 22, 27, 31, 36, 40, 44, 48, 'quantity']
        scipy.array(orispl_prices.columns.difference(['orispl_unit', 'orispl', 'fuel', 'quantity']))
        test = orispl_prices.copy(deep=True)[['orispl_unit', 'orispl', 'fuel']]
        month_weeks = scipy.array(orispl_prices.columns.difference(['orispl_unit', 'orispl', 'fuel', 'quantity']))
        for c in scipy.arange(52)+1:
            if c in month_weeks:
                test['fuel_price'+ str(c)] = orispl_prices[c]
            else:
                test['fuel_price'+ str(c)] = test['fuel_price'+str(c-1)]
        orispl_prices = test.copy(deep=True)
        #now we add in the weekly fuel commodity prices
        prices_fuel_commodity = self.fuel_commodity_prices
        f_array = orispl_prices[orispl_prices['fuel_price1'].isna()].fuel.unique()
        for f in f_array:
            l = len(orispl_prices.loc[orispl_prices.fuel==f, orispl_prices.columns.difference(['orispl_unit', 'orispl', 'fuel'])])
            orispl_prices.loc[orispl_prices.fuel==f, orispl_prices.columns.difference(['orispl_unit', 'orispl', 'fuel'])] = scipy.tile(prices_fuel_commodity[f], (l,1))
        #now we have orispl_prices, which has a weekly fuel price for each orispl_unit based mostly on EIA923 data with some commodity, national-level data from EIA to supplement
        #now merge the fuel price columns into self.df
        orispl_prices.drop(['orispl', 'fuel'], axis=1, inplace=True)
        self.df = self.df.merge(orispl_prices, on='orispl_unit', how='left')    


    def addGenMinOut(self):
        """ 
        Adds fuel price and vom costs to the generator dataframe 'self.df'
        ---
        """
        df = self.df.copy(deep=True)
        #define min_out, based on the NREL Black & Veatch report (2012)
        min_out_coal = 0.4
        min_out_ngcc = 0.5
        min_out_ngst = min_out_coal #assume the same as coal boiler
        min_out_nggt = 0.5 
        min_out_oilst = min_out_coal #assume the same as coal boiler
        min_out_oilgt = min_out_nggt #assume the same as gas turbine
        min_out_nuc = 0.5
        min_out_bio = 0.4
        df['min_out_multiplier'] = scipy.where(df.fuel_type=='oil', scipy.where(df.prime_mover=='st', min_out_oilst, min_out_oilgt), scipy.where(df.fuel_type=='biomass',min_out_bio, scipy.where(df.fuel_type=='coal',min_out_coal, scipy.where(df.fuel_type=='nuclear',min_out_nuc, scipy.where(df.fuel_type=='gas', scipy.where(df.prime_mover=='gt', min_out_nggt, scipy.where(df.prime_mover=='st', min_out_ngst, min_out_ngcc)), 0.10)))))
        df['min_out'] = df.mw * df.min_out_multiplier
        self.df = df          
        
        
    def addGenVom(self):
        """ 
        Adds vom costs to the generator dataframe 'self.df'
        ---
        """
        df = self.df.copy(deep=True)
        #define vom, based on p.7 of EIA "Capital Cost Estimates for Utility Scale Electricity Generating Plants" November 2016
        vom_coal = 4.6
        vom_ngcc = 2.0
        vom_ngst = vom_coal #assume the same as coal boiler
        vom_nggt = 10.7 
        vom_oilst = vom_coal #assume the same as coal boiler
        vom_oilgt = vom_nggt #assume the same as gas turbine
        vom_nuc = 2.3 
        vom_bio = 4.2
        df['vom'] = scipy.where(df.fuel_type=='oil', scipy.where(df.prime_mover=='st', vom_oilst, vom_oilgt), scipy.where(df.fuel_type=='biomass',vom_bio, scipy.where(df.fuel_type=='coal',vom_coal, scipy.where(df.fuel_type=='nuclear',vom_nuc, scipy.where(df.fuel_type=='gas', scipy.where(df.prime_mover=='gt', vom_nggt, scipy.where(df.prime_mover=='st', vom_ngst, vom_ngcc)), 5.0)))))
        self.df = df


    def calcDemandData(self):
        """ 
        Uses CEMS data to calculate net demand (i.e. total fossil generation), total emissions, and each generator type's contribution to the generation mix
        ---
        Creates
        self.hist_dispatch : one row per hour of the year, columns for net demand, total emissions, operating cost of the marginal generator, and the contribution of different fuels to the total energy production
        """
        print 'Calculating demand data from CEMS...'
        #re-compile the cems data adding in fuel and fuel type
        df = self.df_cems.copy(deep=True)
        merge_orispl_unit = self.df.copy(deep=True)[['orispl_unit', 'fuel', 'fuel_type']]
        merge_orispl = self.df.copy(deep=True)[['orispl', 'fuel', 'fuel_type']].drop_duplicates('orispl')
        df = df.merge(merge_orispl_unit, left_index=True, how='left', on=['orispl_unit']) 
        df.loc[df.fuel.isna(), 'fuel'] = scipy.array(df[df.fuel.isna()].merge(merge_orispl, left_index=True, how='left', on=['orispl']).fuel_y)
        df.loc[df.fuel_type.isna(), 'fuel_type'] = scipy.array(df[df.fuel_type.isna()].merge(merge_orispl, left_index=True, how='left', on=['orispl']).fuel_type_y)
        #build the hist_dispatch dataframe
        #start with the datetime column
        start_date_str = (self.df_cems.date.min()[-4:] + '-' + self.df_cems.date.min()[:5] + ' 00:00')
        date_hour_count = len(self.df_cems.date.unique())*24#+1
        hist_dispatch = pandas.DataFrame(scipy.array([pandas.Timestamp(start_date_str) + datetime.timedelta(hours=i) for i in xrange(date_hour_count)]), columns=['datetime'])
        #add columns by aggregating df by date + hour
        hist_dispatch['demand'] = df.groupby(['date','hour'], as_index=False).sum().mwh
        hist_dispatch['co2_tot'] = df.groupby(['date','hour'], as_index=False).sum().co2_tot # * 2000
        hist_dispatch['so2_tot'] = df.groupby(['date','hour'], as_index=False).sum().so2_tot
        hist_dispatch['nox_tot'] = df.groupby(['date','hour'], as_index=False).sum().nox_tot
        hist_dispatch['coal_mix'] = df[(df.fuel_type=='coal') | (df.fuel=='SGC')].groupby(['date','hour'], as_index=False).sum().mwh
        hist_dispatch['gas_mix'] = df[df.fuel_type=='gas'].groupby(['date','hour'], as_index=False).sum().mwh
        hist_dispatch['oil_mix'] = df[df.fuel_type=='oil'].groupby(['date','hour'], as_index=False).sum().mwh
        hist_dispatch['biomass_mix'] = df[(df.fuel_type=='biomass') | (df.fuel=='obs') | (df.fuel=='wds') | (df.fuel=='blq') | (df.fuel=='msw') | (df.fuel=='lfg') | (df.fuel=='ab') | (df.fuel=='obg') | (df.fuel=='obl') | (df.fuel=='slw')].groupby(['date','hour'], as_index=False).sum().mwh
        hist_dispatch['geothermal_mix'] = df[(df.fuel_type=='geothermal') | (df.fuel=='geo')].groupby(['date','hour'], as_index=False).sum().mwh
        hist_dispatch['hydro_mix'] = df[(df.fuel_type=='hydro') | (df.fuel=='wat')].groupby(['date','hour'], as_index=False).sum().mwh
        hist_dispatch['nuclear_mix'] = df[df.fuel=='nuc'].groupby(['date','hour'], as_index=False).sum().mwh
        hist_dispatch.fillna(0, inplace=True)
        #fill in last line to equal the previous line
        #hist_dispatch.loc[(len(hist_dispatch)-1)] = hist_dispatch.loc[(len(hist_dispatch)-2)]
        hist_dispatch = hist_dispatch.fillna(0)
        self.hist_dispatch = hist_dispatch


    def addElecPriceToDemandData(self):
        """ 
        Calculates the historical electricity price for the nerc region, adding it as a new column to the demand data
        ---
        """
        print 'Adding historical electricity prices...'
        #We will use FERC 714 data, where balancing authorities and similar entities report their locational marginal prices. This script pulls in those price for every reporting entity in the nerc region and takes the max price across the BAs/entities for each hour.
        df = self.ferc714.copy(deep=True)
        df_ids = self.ferc714_ids.copy(deep=True)
        nerc_region = self.nerc
        year = self.year
        df_ids_bas = list(df_ids[df_ids.nerc == nerc_region].respondent_id.values)
        #aggregate the price data by mean price per hour for any balancing authorities within the nerc region
        df_bas = df[df.respondent_id.isin(df_ids_bas) & (df.report_yr==year)][['lambda_date', 'respondent_id', 'hour01', 'hour02', 'hour03', 'hour04', 'hour05', 'hour06', 'hour07', 'hour08', 'hour09', 'hour10', 'hour11', 'hour12', 'hour13', 'hour14', 'hour15', 'hour16', 'hour17', 'hour18', 'hour19', 'hour20', 'hour21', 'hour22', 'hour23', 'hour24']]
        df_bas.drop(['respondent_id'], axis=1, inplace=True)
        df_bas.columns = ['date',1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]
        df_bas_temp = pandas.melt(df_bas, id_vars=['date'])
        df_bas_temp.date = df_bas_temp.date.str[0:-7] + (df_bas_temp.variable - 1).astype(str) + ':00'
        df_bas_temp['time'] = df_bas_temp.variable.astype(str) + ':00'
        df_bas_temp['datetime'] = pandas.to_datetime(df_bas_temp.date)
        df_bas_temp.drop(columns=['date', 'variable', 'time'], inplace=True)
        #aggregate by datetime
        df_bas_temp = df_bas_temp.groupby('datetime', as_index=False).max()
        df_bas_temp.columns = ['datetime', 'price']        
        #add the price column to self.hist_dispatch
        self.hist_dispatch['gen_cost_marg'] = df_bas_temp.price


    def demandTimeSeries(self):
        """ 
        Re-formats and slices self.hist_dispatch to produce a demand time series to be used by the dispatch object
        ---
        Creates
        self.demand_data : row for each hour. columns for datetime and demand
        """
        print 'Creating "demand_data" time series...'
        #demand using CEMS data
        demand_data = self.hist_dispatch.copy(deep=True)
        demand_data.datetime = pandas.to_datetime(demand_data.datetime)
        self.demand_data = demand_data[['datetime', 'demand']]
    
    
    def cemsBoxPlot(self, plot_col):
        """ 
        Creates a box plot of the hourly CEMS data for each unique orispl_unit for the given column
        ---
        plot_col: 'co2', 'heat_rate', etc.
        """
        print 'Creating "demand_data" time series...'
        #copy the CEMS data
        cems_copy = gd.df_cems.copy(deep=True)
        #each uniqe unit tag
        ounique = cems_copy.orispl_unit.unique()
       #empty data frame for results
        result = pandas.DataFrame({'orispl_unit': ounique, plot_col+'_5': scipy.zeros_like(ounique), plot_col+'_25': scipy.zeros_like(ounique), plot_col+'_50': scipy.zeros_like(ounique), plot_col+'_75': scipy.zeros_like(ounique), plot_col+'_95': scipy.zeros_like(ounique), 'data_points': scipy.zeros_like(ounique)})
        #for each unique unit calculate the 5th, 25th, median, 75th, and 95th percentile data
        print 'Calculating quantiles...'
        for o in ounique:
            cems_e_test = cems_copy.loc[cems_copy.orispl_unit==o, plot_col]
            if len(cems_e_test) != 0:
                result.loc[result.orispl_unit==o, plot_col+'_5'] = cems_e_test.quantile(0.05)
                result.loc[result.orispl_unit==o, plot_col+'_25'] = cems_e_test.quantile(0.25)
                result.loc[result.orispl_unit==o, plot_col+'_50'] = cems_e_test.median()
                result.loc[result.orispl_unit==o, plot_col+'_75'] = cems_e_test.quantile(0.75)
                result.loc[result.orispl_unit==o, plot_col+'_95'] = cems_e_test.quantile(0.95)
                result.loc[result.orispl_unit==o, 'data_points'] = len(cems_e_test)
        #sort the results for plotting
        result = result.sort_values(by=[plot_col+'_50']).reset_index()
        #plot the results to be like a box plot
        x = scipy.arange(len(ounique))
        f, ax = matplotlib.pylab.subplots(1, figsize = (7,35))
        ax.scatter(result[plot_col+'_5'], x, marker='.', color='k', s=5)
        ax.scatter(result[plot_col+'_25'], x, marker='|', color='k', s=40)
        ax.scatter(result[plot_col+'_50'], x, marker='s', color='r', s=15)
        ax.scatter(result[plot_col+'_75'], x, marker='|', color='k', s=40)
        ax.scatter(result[plot_col+'_95'], x, marker='.', color='k', s=5)
        xmax = scipy.where(plot_col=='co2', 2000, scipy.where(plot_col=='so2', 5, scipy.where(plot_col=='nox', 3, 20)))
        for a in x:
            ax.plot([result.loc[a,plot_col+'_5'],result.loc[a,plot_col+'_95']], [a,a], lw=1, c='black', alpha=0.25)
            matplotlib.pylab.text(-xmax*0.2, a, result.loc[a, 'orispl_unit'])
        for v in [xmax*0.25, xmax*0.5, xmax*0.75]:
            ax.axvline(v, c='black', lw=1, alpha=0.5, ls=':')
        ax.set_xlim(0,xmax)
        ax.set_xticks([xmax*0.25, xmax*0.5, xmax*0.75, xmax])
        ax.get_yaxis().set_ticks([])
        if plot_col == 'heat_rate':
            ax.set_xlabel('Heat Rate [mmBtu/MWh]')
        else:
            ax.set_xlabel(plot_col.title() + ' Emissions Rate [kg/MWh]')
        return f

    
 


class bidStack(object):
    def __init__(self, gen_data_obj, time=1, dropNucHydroGeo=False, include_min_output=True):
        """ 
        1) Bring in the generator data created by the "generatorData" class.
        2) Calculate the generation cost for each generator and sort the generators by generation cost. Default emissions prices [$/kg] are 0.00 for all emissions.
        ---
        gen_data_obj : a generatorData object
        time : number denoting which time period we are interested in. Default is weeks, so time=15 would look at the 15th week of heat rate, emissions rates, and fuel prices
        dropNucHydroGeo : if True, nuclear, hydro, and geothermal plants will be removed from the bidstack (e.g. to match CEMS data)
        include_min_output : if True, will include a representation of generators' minimum output constraints that impacts the marginal generators in the dispatch. So, a "True" value here is closer to the real world.
        """
        self.df_0 = gen_data_obj.df
        self.df = self.df_0.copy(deep=True)
        self.year = gen_data_obj.year
        self.time = time
        self.include_min_output = include_min_output
        if dropNucHydroGeo:
            self.dropNuclearHydroGeo()
        self.calcGenCost()
        self.calcFullMeritOrder()
        self.addFuelColor()

    def dropNuclearHydroGeo(self):
        """ 
        Removes nuclear, hydro, and geothermal plants from self.df (since they don't show up in CEMS)
        ---
        """
        self.df = self.df[(self.df.fuel!='nuc') & (self.df.fuel!='wat') & (self.df.fuel!='geo')]
     
    
    def updateTime(self, t_new):
        """ Updates self.time
        ---
        """
        self.time = t_new
        self.calcGenCost()
        self.calcFullMeritOrder()
        
                
    def calcGenCost(self, co2_price=0.0, so2_price=0.0, nox_price=0.0):
        """ Calculate average costs that are function of generator data, fuel cost, and emissions prices.
        gen_cost ($/MWh) = (heat_rate * "fuel"_price) + (co2 * co2_price) + (so2 * so2_price) + (nox * nox_price) + vom 
        ---
        co2_price : [$/kg]
        so2_price : [$/kg]
        nox_price : [$/kg]
        """
        df = self.df.copy(deep=True)
        df['fuel_cost'] = df['heat_rate' + str(self.time)] * df['fuel_price' + str(self.time)] 
        df['co2_cost'] = df['co2' + str(self.time)] * co2_price 
        df['so2_cost'] = df['so2' + str(self.time)] * so2_price 
        df['nox_cost'] = df['nox' + str(self.time)] * nox_price 
        df['gen_cost'] = df.fuel_cost + df.co2_cost + df.so2_cost + df.nox_cost + df.vom
        df = df.append(df.loc[0]*0) #add a zero generator so that the bid stack goes all the way down to zero. This is important for calculting information for the marginal generator when the marginal generator is the first one in the bid stack.
        df.sort_values('gen_cost', inplace=True)
        df.reset_index(drop=True, inplace=True)
        df['demand'] = df.mw.cumsum()	
        df['f'] = df['demand']
        df['s'] = scipy.append(0, scipy.array(df.f[0:-1]))
        df['a'] = scipy.maximum(df.s - df.min_out*(1/0.10), 1.0)       
        self.df = df
    
    
    def addFuelColor(self):
        """ Assign a fuel type for each fuel and a color for each fuel type to be used in charts
        ---
        creates 'fuel_type' and 'fuel_color' columns
        """
        c = {'gas':'#999999', 'coal':'#bf5b17', 'oil':'#252525' , 'nuclear':'#984ea3', 'hydro':'#386cb0', 'biomass':'#7fc97f', 'geothermal':'#e31a1c', 'ofsl': '#c994c7'}
        self.df['fuel_color'] = '#bcbddc'
        for c_key in c.iterkeys():
            self.df.loc[self.df.fuel_type == c_key, 'fuel_color'] = c[c_key]        


    def returnMarginalGenerator(self, demand, return_type):
        """ Given demand and return_type inputs, return the gen (name or ID) or fuel (fuel type) of the marginal generator.
        ---
        demand : [MW]
        return : column header of self.df being returned (e.g. 'gen', 'fuel_type', etc.)
        """
        ind = scipy.minimum(self.df.index[self.df.demand <= demand][-1], len(self.df)-2)
        #if we have weekly data, return depending on which week (self.time) it is
        if return_type in ['heat_rate', 'co2', 'so2', 'nox', 'fuel_price', 'mw']:
            return self.df[return_type + str(self.time)][ind+1]
        else:
            return self.df[return_type][ind+1]
		
					
    def returnTotalCost(self, demand):
        """ Given demand input, return the integral of the bid stack generation cost (i.e. the total operating cost of the online power plants).
        ---
        demand : [MW]
        return : integral value of the bid stack cost = total operating costs of the online generator fleet [$].
        """
        ind = self.df.index[self.df.demand <= demand][-1]
        tmp = self.df.iloc[0:ind+1,]
        return sum( tmp.gen_cost * tmp['mw' + str(self.time)] ) + self.returnMarginalGenerator(demand, 'gen_cost') * (demand - self.df.iloc[ind].demand)
      
       
    def returnTotalEmissions(self, demand, emissions_type):
        """ Given demand and emissions_type inputs, return the integral of the bid stack emissions (i.e. the total emissions of the online power plants).
        ---
        demand : [MW]
        emissions_type : 'co2', 'so2', 'nox', etc.
        return : integral value of the bid stack emissions = total emissions of the online generator fleet [lbs].
        """
        ind = self.df.index[self.df.demand <= demand][-1]
        tmp = self.df.iloc[0:ind+1,]
        return sum( tmp[emissions_type + str(self.time)] * tmp['mw' + str(self.time)] ) + self.returnMarginalGenerator(demand, emissions_type) * (demand - self.df.iloc[ind].demand)
    
    
    def returnTotalFuelMix(self, demand, is_fuel_type):
        """ Given demand and emissions_type inputs, return the total MW of online generation of fuel_type is_fuel_type.
        ---
        demand : [MW]
        is_fuel_type : 'is_coal', etc.
        return : total amount of online generation of type is_fuel_type
        """
        ind = self.df.index[self.df.demand <= demand][-1]
        tmp = self.df.iloc[0:ind+1,]
        return sum(tmp[is_fuel_type] * tmp['mw' + str(self.time)]) + self.df.loc[ind+1, is_fuel_type] * (demand - self.df.loc[ind+1, 's'])
       
        
    def calcFullMeritOrder(self):
        """ Calculates the base_ and marg_ co2, so2, nox, and coal_mix, where "base_" represents the online "base load" that does not change with marginal changes in demand and "marg_" represents the marginal portion of the merit order that does change with marginal changes in demand. The calculation of base_ and marg_ changes depending on whether the minimum output constraint (the include_min_output variable) is being used.
        ---
        """
        df = self.df.copy(deep=True)
        #INCLUDING MIN OUTPUT
        if self.include_min_output:
            #emissions
            for e in ['co2', 'so2', 'nox']:
                df['full_' + e + '_base'] = 0.1*df.a.apply(self.returnTotalEmissions, args=(e,)) + 0.9*df.s.apply(self.returnTotalEmissions, args=(e,)) + df.s.apply(self.returnMarginalGenerator, args=(e,)) * df.s.apply(self.returnMarginalGenerator, args=('min_out',)) #calculate the base emissions 
                df['full_' + e + '_marg'] = ((df.s.apply(self.returnTotalEmissions, args=(e,)) - df.a.apply(self.returnTotalEmissions, args=(e,))) / (df.s-df.a) * (df.min_out/(df.f-df.s)) + df.s.apply(self.returnMarginalGenerator, args=(e,)) * (1 -(df.min_out/(df.f-df.s)))).fillna(0.0) #calculate the marginal emissions
            #fuel mix
            for f in ['gas', 'coal', 'oil', 'nuclear', 'hydro', 'geothermal', 'biomass']:
                df['full_' + f + '_mix_base'] = 0.1*df.a.apply(self.returnTotalFuelMix, args=(('is_'+f),)) + 0.9*df.s.apply(self.returnTotalFuelMix, args=(('is_'+f),)) + self.df['is_'+f] * df.s.apply(self.returnMarginalGenerator, args=('min_out',)) #calculate the base coal_mix
                df['full_' + f + '_mix_marg'] = ((df.s.apply(self.returnTotalFuelMix, args=(('is_'+f),)) - df.a.apply(self.returnTotalFuelMix, args=(('is_'+f),))) / (df.s-df.a) * (df.min_out/(df.f-df.s)) + df.s.apply(self.returnMarginalGenerator, args=(('is_'+f),)) * (1 -(df.min_out/(df.f-df.s)))).fillna(0.0) #calculate the marginal coal_mix   
        #EXCLUDING MIN OUTPUT
        if not self.include_min_output:
            #emissions
            for e in ['co2', 'so2', 'nox']:
                df['full_' + e + '_base'] = df.s.apply(self.returnTotalEmissions, args=(e,)) #calculate the base emissions, which is now the full load emissions of the generators in the merit order below the marginal unit
                df['full_' + e + '_marg'] = df.s.apply(self.returnMarginalGenerator, args=(e,)) #calculate the marginal emissions, which is now just the emissions rate of the marginal generator
            #fuel mix
            for f in ['gas', 'coal', 'oil', 'nuclear', 'hydro', 'geothermal', 'biomass']:
                df['full_' + f + '_mix_base'] = df.s.apply(self.returnTotalFuelMix, args=(('is_'+f),)) #calculate the base fuel_mix, which is now the full load coal mix of the generators in the merit order below the marginal unit
                df['full_' + f + '_mix_marg'] = df.s.apply(self.returnMarginalGenerator, args=(('is_'+f),)) #calculate the marginal fuel_mix, which is now just the fuel_type of the marginal generator
        #update the master dataframe df
        self.df = df
        
        
    def returnFullTotalValue(self, demand, col_type):
        """ Given demand and emissions_type inputs, return the total column of the online power plants in the Full model (the Full model includes the minimum output constraint).
        ---
        demand : [MW]
        col_type : 'co2', 'so2', 'nox', 'coal_mix', etc.
        return : total emissions = base emissions (marginal unit) + marginal emissions (marginal unit) * (D - s)
        """
        ind = self.df.index[self.df.demand <= demand][-1]+1
        return self.df.loc[ind,'full_' + col_type + '_base'] + self.df.loc[ind,'full_' + col_type + '_marg'] * (demand - self.df.loc[ind, 's'])
    
    
    def returnFullMarginalValue(self, demand, col_type):
        """ Given demand and col_type inputs, return the col_type (i.e. 'co2' for marginal co2 emissions rate or 'coal_mix' for coal share of the generation) of the marginal units in the Full model (the Full model includes the minimum output constraint).
        ---
        demand : [MW]
        col_type : 'co2', 'so2', 'nox', 'coal_mix', etc.
        return : full_"emissions_type"_marg as calculated in the Full model: calcFullMeritOrder
        """
        ind = self.df.index[self.df.demand <= demand][-1]+1
        return self.df.loc[ind,'full_' + col_type + '_marg']
    
    
    def plotBidStack(self, df_column, plot_type):
        """ Given a name for the df_column, plots a bid stack with demand on the x-axis and the df_column data on the y-axis. For example bidStack.plotBidStack('gen_cost', 'bar') would output the traditional merit order curve.
        ---
        df_column : column header from df, e.g. 'gen_cost', 'co2', 'so2', etc.
        plot_type : 'bar' or 'line'
        return : bid stack plot
        """	
        #color for any generators without a fuel_color entry
        empty_color = '#dd1c77'
        if df_column == 'gen_cost':
            y_lab = 'Generation Cost [$/MWh]'
            y_data = self.df[df_column]
        if df_column == 'co2':
            y_lab = 'CO$_2$ Emissions [kg/MWh]'
            y_data = self.df[df_column + str(self.time)]
        if df_column == 'so2':
            y_lab = 'SO$_2$ Emissions [kg/MWh]'
            y_data = self.df[df_column + str(self.time)]
        if df_column == 'nox':
            y_lab = 'NO$_x$ Emissions [kg/MWh]'
            y_data = self.df[df_column + str(self.time)]
        matplotlib.pylab.clf()
        f = matplotlib.pylab.figure(figsize=(4,4))
        ax = f.add_subplot(111)
        if plot_type == 'line':
            ax.plot( self.df.demand/1000, y_data, linewidth=2.5)
        elif plot_type == 'bar':
            ax.bar(self.df.demand/1000, height=y_data, width=-scipy.maximum(0.2, self.df['mw' + str(self.time)]/1000), color=self.df.fuel_color.replace('', empty_color), align='edge')
            ##add legend above chart
            #color_legend = []
            #for c in self.df.fuel_color.unique():
            #    color_legend.append(matplotlib.patches.Patch(color=c, label=self.df.fuel_type[self.df.fuel_color==c].iloc[0]))
            #ax.legend(handles=color_legend, bbox_to_anchor=(0.5, 1.2), loc='upper center', ncol=3, fancybox=True, shadow=True)
        else:
            print '***Error: enter valid argument for plot_type'
            pass
        matplotlib.pylab.ylim(ymax=y_data.quantile(0.98)) #take the 98th percentile for the y limits.
        #ax.set_xlim(7.5,38)
        #ax.set_xticks((10,15,20,25,30,35))
        #ax.set_ylim(0, 1250)
        #ax.set_ylim(0, 85)
        #ax.set_yticks((0, 20, 40, 60, 80))
        matplotlib.pylab.xlabel('Generation [GW]')
        matplotlib.pylab.ylabel(y_lab)
        matplotlib.pylab.tight_layout()
        matplotlib.pylab.show()
        return f
    
    
    
    

class dispatch(object):
    def __init__(self, bid_stack_object, demand_df, time_array=0):
        """ Read in bid stack object and the demand data. Solve the dispatch by projecting the bid stack onto the demand time series, updating the bid stack object regularly according to the time_array
        ---
        gen_data_object : a object defined by class generatorData
        bid_stack_object : a bid stack object defined by class bidStack
        demand_df : a dataframe with the demand data 
        time_array : a scipy array containing the time intervals that we are changing fuel price etc. for. E.g. if we are doing weeks, then time_array=scipy.arange(52) + 1 to get an array of (1, 2, 3, ..., 51, 52)
        """
        self.bs = bid_stack_object
        self.df = demand_df
        self.time_array = time_array
        self.addDFColumns()
        
               
    def addDFColumns(self):
        """ Add additional columns to self.df to hold the results of the dispatch. New cells initially filled with zeros
        ---
        """
        indx = self.df.index
        cols = scipy.array(('gen_cost_marg', 'gen_cost_tot' ,'co2_marg' ,'co2_tot' ,'so2_marg' ,'so2_tot' ,'nox_marg' ,'nox_tot' ,'biomass_mix' ,'coal_mix' ,'gas_mix' ,'geothermal_mix' ,'hydro_mix' ,'nuclear_mix' ,'oil_mix', 'marg_gen', 'coal_mix_marg', 'marg_gen_fuel_type'))
        dfExtension = pandas.DataFrame(index=indx, columns=cols).fillna(0)
        self.df = pandas.concat([self.df, dfExtension], axis=1)


    def calcDispatchSlice(self, bstack, start_date=0, end_date=0):
        """ For each datum in demand time series (e.g. each hour) between start_date and end_date calculate the dispatch
        ---
        bstack: an object created using the simple_dispatch.bidStack class
        start_datetime : string of format '2014-01-31' i.e. 'yyyy-mm-dd'. If argument == 0, uses start date of demand time series
        end_datetime : string of format '2014-01-31' i.e. 'yyyy-mm-dd'. If argument == 0, uses end date of demand time series
        """
        if start_date==0:
            start_date = self.df.datetime.min()
        else:
            start_date = pandas._libs.tslib.Timestamp(start_date)
        if end_date==0:
            end_date = self.df.datetime.max()
        else:
            end_date = pandas._libs.tslib.Timestamp(end_date)
        #slice of self.df within the desired dates    
        df_slice = self.df[(self.df.datetime >= pandas._libs.tslib.Timestamp(start_date)) & (self.df.datetime < pandas._libs.tslib.Timestamp(end_date))].copy(deep=True)
        #calculate the dispatch for the slice by applying the return###### functions of the bstack object
        df_slice['gen_cost_marg'] = df_slice.demand.apply(bstack.returnMarginalGenerator, args=('gen_cost',)) #generation cost of the marginal generator
        df_slice['gen_cost_tot'] = df_slice.demand.apply(bstack.returnTotalCost) #generation cost of the total generation fleet
        for e in ['co2', 'so2', 'nox']:
            df_slice[e + '_marg'] = df_slice.demand.apply(bstack.returnFullMarginalValue, args=(e,)) #emissions rate (kg/MWh) of marginal generators
            df_slice[e + '_tot'] = df_slice.demand.apply(bstack.returnFullTotalValue, args=(e,)) #total emissions (kg) of online generators
        for f in ['gas', 'oil', 'coal', 'nuclear', 'biomass', 'geothermal', 'hydro']:
            df_slice[f + '_mix'] = df_slice.demand.apply(bstack.returnFullTotalValue, args=(f+'_mix',))
        df_slice['coal_mix_marg'] = df_slice.demand.apply(bstack.returnFullMarginalValue, args=('coal_mix',))
        df_slice['marg_gen_fuel_type'] = df_slice.demand.apply(bstack.returnMarginalGenerator, args=('fuel_type',))
        #add the dispatch of the time slice, df_slice, to the results dataframe, self.df
        self.df[(self.df.datetime >= pandas._libs.tslib.Timestamp(start_date)) & (self.df.datetime < pandas._libs.tslib.Timestamp(end_date))] = df_slice
        
              
    def calcDispatchAll(self):
        """ Runs calcDispatchSlice for each time slice in the fuel_prices_over_time dataframe, creating a new bidstack each time. So, fuel_prices_over_time contains multipliers (e.g. 0.95 or 1.14) for each fuel type (e.g. ng, lig, nuc) for different slices of time (e.g. start_date = '2014-01-07' and end_date = '2014-01-14'). We use these multipliers to change the fuel prices seen by each generator in the bidStack object. After changing each generator's fuel prices (using bidStack.updateFuelPrices), we re-calculate the bidStack merit order (using bidStack.calcGenCost), and then calculate the dispatch for the slice of time defined by the fuel price multipliers. This way, instead of calculating the dispatch over the whole year, we can calculate it in chunks of time (e.g. weeks) where each chunk of time has different fuel prices for the generators. 
        Right now the only thing changing per chunk of time is the fuel prices based on trends in national commodity prices. Future versions might try and do regional price trends and add things like maintenance downtime or other seasonal factors.
        ---
        fills in the self.df dataframe one time slice at a time
        """
        #run the whole solution if self.fuel_prices_over_time isn't being used
        if scipy.shape(self.time_array) == (): #might be a more robust way to do this. Would like to say if ### == 0, but doing that when ### is a dataframe gives an error
            self.calcDispatchSlice(self.bs)
        #otherwise, run the dispatch in time slices, updating the bid stack each slice
        else:
            for t in self.time_array:
                print str(round(t/float(len(self.time_array)),3)*100) + '% Complete'
                #update the bidStack object to the current week
                self.bs.updateTime(t)
                #calculate the dispatch for the time slice over which the updated fuel prices are relevant
                start = (datetime.datetime.strptime(str(self.bs.year) + '-01-01', '%Y-%m-%d') + datetime.timedelta(days=7.05*(t-1)-1)).strftime('%Y-%m-%d') 
                end = (datetime.datetime.strptime(str(self.bs.year) + '-01-01', '%Y-%m-%d') + datetime.timedelta(days=7.05*(t)-1)).strftime('%Y-%m-%d') 
                #note that calcDispatchSlice updates self.df, so there is no need to do it in this calcDispatchAll function
                self.calcDispatchSlice(self.bs, start_date=start ,end_date=end)
                
                 



if __name__ == '__main__':
    
    #input variables. Right now the github only has 2017 data on it.
    run_year = 2017
    nerc_region = 'FRCC'
    historical_dispatch_save_folder = 'C:\\Users\\tdeet\\Documents\\data\\processed\\epa\\CEMS' #where to save the historical dispatch results
    simulated_dispatch_save_folder = 'C:\\Users\\tdeet\\Documents\\analysis\\modules\\python\\simple_dispatch' #where to save the simulated dispatch results
    #specific the location of the data directories
    #2017: (Note: uses 2017 CEMS with 2016 eGRID)
    egrid_data_xlsx = 'C:\\Users\\tdeet\\Documents\\data\\raw\\eGRID\\egrid2016_data.xlsx'
    eia923_schedule5_xlsx = 'C:\\Users\\tdeet\\Documents\\data\\raw\\eia\\eia9232017\\EIA923_Schedules_2_3_4_5_M_12_2017_Early_Release.xlsx'
    ferc714_part2_schedule6_csv = 'C:\\Users\\tdeet\\Documents\\data\\raw\\ferc\\ferc_714\\Part 2 Schedule 6 - Balancing Authority Hourly System Lambda.csv'
    ferc714IDs_csv='C:\\Users\\tdeet\\Documents\\data\\raw\\ferc\\ferc_714\\Respondent IDs.csv'
    cems_folder_path ='C:\\Users\\tdeet\\Documents\\data\\raw\\epa\\CEMS'
    fuel_commodity_prices_xlsx = 'C:\\Users\\tdeet\\Documents\\data\\processed\\eiaFuelPrices\\fuel_default_prices.xlsx'
    
    #run the generator data object
    gd = generatorData(nerc_region, egrid_fname=egrid_data_xlsx, eia923_fname=eia923_schedule5_xlsx, ferc714IDs_fname=ferc714IDs_csv, ferc714_fname=ferc714_part2_schedule6_csv, cems_folder=cems_folder_path, year=run_year, fuel_commodity_prices_excel_dir=fuel_commodity_prices_xlsx) 
    #save the historical dispatch  
    gd.hist_dispatch.to_csv(historical_dispatch_save_folder + '\\%s_%s_hourly_demand_and_fuelmix.csv'%(str(run_year), nerc_region))
            
    #run the bidStack object - use information about the generators (from gd) to create a merit order (bid stack) of the nerc region's generators
    bs = bidStack( gd, time=1, dropNucHydroGeo=True, include_min_output=True) #NOTE: set dropNucHydroGeo to True if working with data that only looks at fossil fuels (e.g. CEMS)
    #bid_stack_cost = bs.plotBidStack('gen_cost', plot_type='bar') #plot the merit order
    #bid_stack_cost.savefig('C:\\Users\\tdeet\\Documents\\media\\publications\\2018-03 MEFS Simple Dispatch\\LaTeX\\images_raw\\fStackCost%s_%s.png'%(nerc_region, str(run_year)), dpi=500, bbox_inches='tight')
        
    #run the dispatch object - use the nerc region's merit order (bs), a demand timeseries (gd.demand_data), and a time array (default is array([ 1,  2, ... , 51, 52]) for 52 weeks to run a whole year)
    dp = dispatch(bs, gd.demand_data, time_array=scipy.arange(52)+1) #set up the object
    dp.calcDispatchAll() #function that solves the dispatch for each time period in time_array (default for each week of the year)
    #save dispatch results 
    dp.df.to_csv(simulated_dispatch_save_folder = '\\dispatch_output_weekly_%s_%s.csv'%(nerc_region, str(run_year)), index=False )
        
