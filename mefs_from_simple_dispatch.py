# mefs_from_dispatch
# Thomas Deetjen
# v10
# last edited: 2018-06-08
# class "generate_mefs" creates a dataframe and plot of the MEFs for a given dispatch csv
# v11:
# updated the deciles plotting function so that it plots fuels in comparison with the CEMS historical marginal fuels
# updated the generate_mefs.df to exclude any rows where the "fuel"_mix_marg was greater than 1.01 or less than -0.01. These values shouldn't happen. A value >1.01 suggests that if demand increases by 1, generation increases by more than 1. A negative value suggests that if demand increases by 1, then generation decreases. These are due to aberrations. In the simulated data, the abberations are driven by abrupt, discrete changes in fuel price. In the empirical data, these abberations may be driven by a number of things independent of changing demand. The historical data had a few 100+ values that would clearly throw off the results of a decile that might only have 800 samples with values averaging 0.5. In that case, adding a 100 would raise the average from 0.5 to 0.625.
# added some plotting functionality. Especially 1) plotDispatch.plot_hist_vs_simulated that plots a moving average of historical and simulated data and 2) plotDispatch.plot_density_function that plots cdf or pdf of hourly data - either historical vs. simulated or the error of sim - hist / hist 
# v12:
# made a lot of tweaks and corrections to bugs
# updated plot_density_function so that MEFs error are absolute instead of relative. Relative error didn't work well with the MEFs (especially SO2 and NOX) since they are often very close to zero. If historical MEF is near zero and simulated MEF is not, then the relative error will always be very large despite the magnitude of the absolute error
# v20:
# updated the x limits of the plots the have Total Fossil Generation x axes so that the limits are from the 2.5% and 97.5% of the historical demand over the year
# changed the rolling calculations to include (center = True) so that the window centers on the input demand
# updated the way we calculate the simulated MEFs. Simulated MEFs are now based on the slope of the total emissions curve. 1) Take a scatter plot of total demand vs. total emissoins. 2) Calculate the rolling average of that data. This is now our model for predicting total emissions as a function of total demand. 3) Calculate the rolling slope of the rolling total emissions curve. In this case, we use the rolling window (p), and at a given index (i), then calcualte the change in emissions as (total_rolling_emissions(demand[i+p]) - total_rolling_emissions(demand[i-p])) and divide by the change in demand (demand[i+p] - demand[i-p]). 4) Take the rolling average of the rolling slope.
# added an error calculation function


import pandas
import matplotlib.pylab
import scipy
import scipy.linalg
import scipy.stats
import math
import sklearn.metrics




class generateMefs(object):
    def __init__(self, dispatch_df):
        """ 
        Uses a dispatch data frame to calculate the marginal emissions factors and fuel mix per hour (or per timestep)
        ---
        dispatch_df : a data frame of the dispatch (historical or simulated). From the simple_dispatch module, can use generatorData.hist_dispatch or dispatch.df data frames for the dispatch_df argument here
        """
        self.df = dispatch_df.copy(deep=True)
        self.calc_delta_g()
        self.calc_delta_g_fuel_mix(f_types = ['gas', 'coal', 'oil', 'nuclear', 'hydro', 'geothermal', 'biomass'])
        for e in ['co2', 'so2', 'nox']:
            self.calc_delta_e(e)
            self.calc_mefs(e)
        self.df.replace([scipy.inf, -scipy.inf], scipy.nan, inplace=True)
        self.df.replace(scipy.nan, 999999, inplace=True)
    
    
    def calc_delta_g(self):
        """ Calculate how generation is changing during each hour of the dispatch. Adds delta_g column to self.df
        """
        demand_t_prev = pandas.concat([pandas.Series([self.df.demand[0]]), self.df.demand])[0:-1].values
        self.df['delta_g'] = self.df.demand - demand_t_prev
    
       
    def calc_delta_g_fuel_mix(self, f_types = ['gas', 'coal', 'oil', 'nuclear', 'hydro', 'geothermal', 'biomass']):
        """ Calculate how generation of each fuel type is changing during each hour of the dispatch. Adds delta_g_['fuel'] column to self.df
        """
        for f in f_types:
            self.fuel_t = self.df[f+'_mix']
            self.fuel_t_prev = pandas.concat([pandas.Series([self.fuel_t[0]]), self.fuel_t])[0:-1]
            self.df[f+'_mix_marg'] = scipy.divide((self.fuel_t.values - self.fuel_t_prev.values), self.df.delta_g)
        
    
    def calc_delta_e(self, e_type):
        """ Calculate how emissions are changing during each hour of the dispatch. Adds delta_"e_name" columns to self.df
        ---
        e_type = name of the emissions type of interest (e.g. 'co2', 'so2', 'nox', 'pm25', etc.)
        """
        emissions_t_prev = pandas.concat([pandas.Series([self.df[e_type+'_tot'][0]]), self.df[e_type+'_tot']])[0:-1].values
        self.df['delta_'+e_type] = self.df[e_type+'_tot'] - emissions_t_prev
             
        
    def calc_mefs(self, e_type):
        """ Calculates the marginal emissions factor for each hour of the dispatch. Adds mef_"e_name" columns to self.df
        ---
        e_type = name of the emissions type of interest (e.g. 'co2', 'so2', 'nox', 'pm25', etc.)
        """
        self.df[e_type+'_marg'] = self.df['delta_'+e_type]/self.df['delta_g']
        #truncate unrealistic MEF values - i.e. MEFs cannot be below zero and cannot be larger than the dirtiest plant in the fleet
        mef_max = scipy.where(e_type=='co2', 1600.0, scipy.where(e_type=='so2', 5.0, scipy.where(e_type=='nox', 3.0, 1.0e10)))
        self.df.loc[self.df[e_type+'_marg'] <0, (e_type+'_marg')] = scipy.nan
        self.df.loc[self.df[e_type+'_marg'] >mef_max, (e_type+'_marg')] = scipy.nan
    




class plotDispatch(object):
    def __init__(self, nerc_region, dispatch_df, cems_df, deciles_cedm=[21300, 24073, 26245, 28132, 30081, 32297, 35414, 39355, 45364, 60650], mefs_cedm_co2 = [691, 623, 589, 562, 550, 529, 515, 485, 479, 489], mefs_cedm_so2 = [1.156, 0.951, 0.835, 0.745, 0.731, 0.586, 0.445, 0.324, 0.240, 0.003], mefs_cedm_nox = [0.292, 0.255, 0.238, 0.221, 0.218, 0.191, 0.198, 0.192, 0.232, 0.379]):
        """ 
        Contains a few functions used to plot the dispatch results
        ---
        dispatch_df : the dataframe version of a csv file containing the dispatch output of the simple_dispatch.dispatch
        cems_df : the CEMS '2014_%s_hourly_demand_and_fuelmix.csv'%s(nerc_region) csv file showing fuel_mix in the cems data
        deciles : the demand points that we want to calculate seperate MEFs for [MW]
        """
        self.nerc_region = nerc_region
        self.f_cedm = {'co2': scipy.interpolate.interp1d([0] + deciles_cedm + [1000000], [mefs_cedm_co2[0]] + mefs_cedm_co2 + [mefs_cedm_co2[-1]], bounds_error=False, fill_value='extrapolate'), 'so2': scipy.interpolate.interp1d([0] + deciles_cedm + [1000000], [mefs_cedm_so2[0]] + mefs_cedm_so2 + [mefs_cedm_so2[-1]], bounds_error=False, fill_value='extrapolate'), 'nox': scipy.interpolate.interp1d([0] + deciles_cedm + [1000000], [mefs_cedm_nox[0]] + mefs_cedm_nox + [mefs_cedm_nox[-1]], bounds_error=False, fill_value='extrapolate')}
        self.df = dispatch_df.copy(deep=True)
        self.df[self.df.gen_cost_marg > 150.].gen_cost_marg = 150.
        self.cems_df = cems_df
        self.cems_df[self.cems_df.gen_cost_marg > 150.].gen_cost_marg = 150.
        self.deciles_cedm = deciles_cedm
        self.xlim_tuple = (self.df.demand.quantile(0.05)*0.001, self.df.demand.quantile(0.95)*0.001) #5th and 95th percentiles of annual demand for bounding the x axes of charts with "Total Fossil Generation [GW]" as the x-axis
        self.mefs_cedm_co2 = mefs_cedm_co2
        self.mefs_cedm_so2 = mefs_cedm_so2
        self.mefs_cedm_nox = mefs_cedm_nox
        self.add_dispatch_columns()
        self.rollingCalculations()
        
        
    def add_dispatch_columns(self):
        """ Adds some information to df and cems_df that we want for plotting - namely for plotting marginal coal generation
        ---
        updates the self.df and self.cems_df
        """	
        self.cems_df['coal_mix_p'] = scipy.divide(self.cems_df.coal_mix, self.cems_df.demand)
        self.cems_df['dDemand'] = self.cems_df.demand - scipy.insert(self.cems_df.demand[0:-1].values,0,self.cems_df.demand[0])
        self.cems_df['dCoal_mix'] = self.cems_df.coal_mix - scipy.insert(self.cems_df.coal_mix[0:-1].values,0,self.cems_df.coal_mix[0])
        self.cems_df['dGas_mix'] = self.cems_df.gas_mix - scipy.insert(self.cems_df.gas_mix[0:-1].values,0,self.cems_df.gas_mix[0])
        self.cems_df['p_marginal_coal'] = scipy.maximum(scipy.minimum(scipy.divide(self.cems_df['dCoal_mix'], self.cems_df['dDemand']),1),0)
        self.df['coal_is_marginal'] = (self.df.marg_gen_fuel_type=='coal').astype(int)
        self.df['coal_mix_p'] = scipy.divide(self.df.coal_mix, self.df.demand)


    def rollingCalculations(self, p=438):
        """ Calculates rolling average of the total emissions and the rolling average of the slope of the total emissions for both the historical and the simulated dispatch data
        ---
        p : the number of values included in the rolling calculation window
        creates :
            self.hist_sorted = similar to self.cems_df but sorted by demand and including rolling total emissions and slope
            self.sim_sorted = similar to self.df but sorted by demand and including rolling total emissions and slope
        """	
        #below we will calculate the rolling slope of the demand vs. rolling emissions curves
        #define the dataframes for the hist and sim data, sorted by increasing demand. Calculate a rolling window of (demand_high - demand_low), where demand_high is (demand + p/2) and demand_low is (demand - p/2)
        hist_sorted = self.cems_df[['demand', 'datetime', 'coal_mix', 'coal_mix_p', 'co2_tot', 'so2_tot', 'nox_tot', 'gen_cost_marg']].sort_values(by='demand').copy()
        hist_sorted['dlow'] = list(scipy.repeat(0, p/2.)) + list(hist_sorted.demand.iloc[0:-p]) + list(scipy.repeat(0, p/2.))
        hist_sorted['dhigh'] = list(scipy.repeat(0, p/2.)) + list(hist_sorted.demand.iloc[p:]) + list(scipy.repeat(0, p/2.))
        hist_sorted['dd'] = hist_sorted.dhigh - hist_sorted.dlow
        sim_sorted = self.df[['demand', 'datetime', 'coal_mix', 'coal_mix_p', 'co2_tot', 'so2_tot', 'nox_tot', 'gen_cost_marg']].sort_values('demand').copy()
        sim_sorted[sim_sorted.gen_cost_marg > 150.].gen_cost_marg = 150.
        sim_sorted['dlow'] = list(scipy.repeat(0, p/2.)) + list(sim_sorted.demand.iloc[0:-p]) + list(scipy.repeat(0, p/2.))
        sim_sorted['dhigh'] = list(scipy.repeat(0, p/2.)) + list(sim_sorted.demand.iloc[p:]) + list(scipy.repeat(0, p/2.))
        sim_sorted['dd'] = sim_sorted.dhigh - sim_sorted.dlow
        #for each type of emissions
        for e in ['co2', 'so2', 'nox']:
            #calculate the rolling window of (emissions_high) - (emissions_low)
            #for the historical data
            hist_sorted[e+'_tot_rolling'] = hist_sorted[e+'_tot'].rolling(window=p*2, min_periods=p*2, center=True).mean()
            hist_sorted[e + 'low'] = list(scipy.repeat(0, p/2.)) + list(hist_sorted[e + '_tot_rolling'].iloc[0:-p]) + list(scipy.repeat(0, p/2.))
            hist_sorted[e + 'high'] = list(scipy.repeat(0, p/2.)) + list(hist_sorted[e + '_tot_rolling'].iloc[p:]) + list(scipy.repeat(0, p/2.))
            hist_sorted['d'+e] = hist_sorted[e + 'high'] - hist_sorted[e + 'low']
            #the slope of the emissions curve is d_emissions / d_demand
            hist_sorted[e+'slope_0'] = scipy.divide(hist_sorted['d'+e], hist_sorted.dd)
            hist_sorted[e+'_slope'] = hist_sorted[e+'slope_0'].rolling(window=p, min_periods=p, center=True).mean()
            #for the simulated data
            sim_sorted[e+'_tot_rolling'] = sim_sorted[e+'_tot'].rolling(window=p*2, min_periods=p*2, center=True).mean()
            sim_sorted[e + 'low'] = list(scipy.repeat(0, p/2.)) + list(sim_sorted[e + '_tot_rolling'].iloc[0:-p]) + list(scipy.repeat(0, p/2.))
            sim_sorted[e + 'high'] = list(scipy.repeat(0, p/2.)) + list(sim_sorted[e + '_tot_rolling'].iloc[p:]) + list(scipy.repeat(0, p/2.))
            sim_sorted['d'+e] = sim_sorted[e + 'high'] - sim_sorted[e + 'low']
            #the slope of the emissions curve is d_emissions / d_demand
            sim_sorted[e+'slope_0'] = scipy.divide(sim_sorted['d'+e], sim_sorted.dd)
            sim_sorted[e+'_slope'] = sim_sorted[e+'slope_0'].rolling(window=p, min_periods=p, center=True).mean()
        hist_sorted['gen_cost_marg_rolling'] = hist_sorted['gen_cost_marg'].rolling(window=p*2, min_periods=p*2, center=True).mean()
        sim_sorted['gen_cost_marg_rolling'] = sim_sorted['gen_cost_marg'].rolling(window=p*2, min_periods=p*2, center=True).mean()
        #set the sorted dataframes as self.objects
        self.hist_sorted = hist_sorted
        self.sim_sorted = sim_sorted
        

    def calcError(self, column_2):
        """ Calculates the mean standard error/hist_mean for the simulated and cedm data versus the historical data
        ---
        column_2 : the value for the second column
        p : the number of values included in the rolling calculation windows and slope windows
        return: error_df, a one-row dataframe with all of the error calculations
        """	
        #a dataframe for holding the results
        error_df = pandas.DataFrame(columns=(['nerc', 'variable', 'co2_tot_hour_vs_rolling', 'co2_slope_sim', 'co2_slope_cedm', 'so2_tot_hour_vs_rolling', 'so2_slope_sim', 'so2_slope_cedm', 'nox_tot_hour_vs_rolling', 'nox_slope_sim', 'nox_slope_cedm', 'price']))
        #now calculate the error
        #mean absolute error:
        error_list = []
        for e in ['co2', 'so2', 'nox']:
            #total emissions error: hist_hour vs sim_hour
            #error_list.append(round(sklearn.metrics.mean_absolute_error(self.hist_sorted[e+'_tot'], self.sim_sorted[e+'_tot']) / self.hist_sorted[e+'_tot'].mean(), 2))
            #total emissions error: hist_hour vs sim_rolling
            error_list.append(round(sklearn.metrics.mean_absolute_error(self.hist_sorted[self.hist_sorted[e+'_tot_rolling'].notna()][e+'_tot'], self.sim_sorted[self.sim_sorted[e+'_tot_rolling'].notna()][e+'_tot_rolling']) / self.hist_sorted[e+'_tot'].mean(), 2))
            #total emissions error: hist_rolling vs sim_rolling
            #error_list.append(round(sklearn.metrics.mean_absolute_error(self.hist_sorted[self.hist_sorted[e+'_tot_rolling'].notna()][e+'_tot_rolling'], self.sim_sorted[self.sim_sorted[e+'_tot_rolling'].notna()][e+'_tot_rolling']) / self.hist_sorted[e+'_tot'].mean(), 2))
            #emissions slope error
            error_list.append(round(sklearn.metrics.mean_absolute_error(self.hist_sorted[self.hist_sorted[e+'_slope'].notna()][e+'_slope'], self.sim_sorted[self.sim_sorted[e+'_slope'].notna()][e+'_slope']) / self.hist_sorted[self.hist_sorted[e+'_slope'].notna()][e+'_slope'].mean(), 2))
            #cedm mefs vs. historical emissions slope error
            error_list.append(round(sklearn.metrics.mean_absolute_error(self.hist_sorted[self.hist_sorted[e+'_slope'].notna()][e+'_slope'], self.hist_sorted[self.hist_sorted[e+'_slope'].notna()].demand.apply(self.f_cedm[e])) / self.hist_sorted[self.hist_sorted[e+'_slope'].notna()][e+'_slope'].mean(), 2))
        error_list.append(round(sklearn.metrics.mean_absolute_error(self.hist_sorted[self.hist_sorted['gen_cost_marg_rolling'].notna()]['gen_cost_marg'], self.sim_sorted[self.sim_sorted['gen_cost_marg_rolling'].notna()]['gen_cost_marg_rolling']) / self.hist_sorted['gen_cost_marg'].mean(), 2))
        #add error_list to error_df
        error_df.loc[len(error_df)] = [self.nerc_region, column_2] + error_list   
        return error_df


    def plotDemandEmissions(self, plot_type, figure_dimensions = (7,5)):
        """ Produces a plot of total fossil generation (x-axis) vs. emissions (y-axis)
        ---
        plot_type : if 'total' plots total emissions, if 'marginal' plots marginal emissions with CEDM MEFs
        returns figure
        """
        if plot_type == 'total':
            #total emissions plot
            f, ax = matplotlib.pylab.subplots(1, figsize=figure_dimensions) 
            #scatters
            #co2
            #ax.scatter(self.hist_sorted['demand']/1e3, self.hist_sorted['co2_tot']/1e6, c='#bcbddc', s=1, alpha=0.15)
            #so2           
            #ax.scatter(self.hist_sorted['demand']/1e3, self.hist_sorted['so2_tot']/1e3, c='#a1d99b', s=1, alpha=0.15)
            #nox            
            #ax.scatter(self.hist_sorted['demand']/1e3, self.hist_sorted['nox_tot']/1e3, c='#fdae6b', s=1, alpha=0.15)
            #lines
            #co2
            ax.plot(self.sim_sorted[self.hist_sorted['co2_tot_rolling'].notna()]['demand']/1e3, self.sim_sorted[self.sim_sorted['co2_tot_rolling'].notna()]['co2_tot_rolling']/1e6, c='#756bb1', lw=2)
            #so2
            ax.plot(self.sim_sorted[self.hist_sorted['so2_tot_rolling'].notna()]['demand']/1e3, self.sim_sorted[self.sim_sorted['so2_tot_rolling'].notna()]['so2_tot_rolling']/1e3, c='#31a354', lw=2)
            #nox
            ax.plot(self.sim_sorted[self.hist_sorted['nox_tot_rolling'].notna()]['demand']/1e3, self.sim_sorted[self.sim_sorted['nox_tot_rolling'].notna()]['nox_tot_rolling']/1e3, c='#e6550d', lw=2)
            #dotted lines
            #co2
            ax.plot(self.hist_sorted[self.hist_sorted['co2_tot_rolling'].notna()]['demand']/1e3, self.hist_sorted[self.hist_sorted['co2_tot_rolling'].notna()]['co2_tot_rolling']/1e6, c='#54278f', ls='--', lw=2)
            #so2
            ax.plot(self.hist_sorted[self.hist_sorted['so2_tot_rolling'].notna()]['demand']/1e3, self.hist_sorted[self.hist_sorted['so2_tot_rolling'].notna()]['so2_tot_rolling']/1e3, c='#006d2c', ls='--', lw=2)
            #nox
            ax.plot(self.hist_sorted[self.hist_sorted['nox_tot_rolling'].notna()]['demand']/1e3, self.hist_sorted[self.hist_sorted['nox_tot_rolling'].notna()]['nox_tot_rolling']/1e3, c='#a63603', ls='--', lw=2)
            #axis limits and labels
            ax.set_xlabel('Total Fossil Generation [GW]')
            ax.set_ylabel('Total CO$_2$ Emissions [million kg]\nTotal SO$_2$ & NO$_x$ Emissions [thousand kg]')                  
            ax.set_xlim(self.xlim_tuple)                        
            ax.set_ylim(0)
            return f
        if plot_type == 'marginal':
            #mef / emissions slope plot
            f, ax = matplotlib.pylab.subplots(1, figsize=figure_dimensions) 
            ax2 = ax.twinx()
            #cedm mefs
            ax.plot(self.hist_sorted[self.hist_sorted['co2_slope'].notna()]['demand']/1e3, self.hist_sorted[self.hist_sorted['co2_slope'].notna()]['demand'].apply(self.f_cedm['co2']), c='#cbc9e2', lw=6)
            ax2.plot(self.hist_sorted[self.hist_sorted['so2_slope'].notna()]['demand']/1e3, self.hist_sorted[self.hist_sorted['so2_slope'].notna()]['demand'].apply(self.f_cedm['so2']), c='#bae4b3', lw=6)
            ax2.plot(self.hist_sorted[self.hist_sorted['nox_slope'].notna()]['demand']/1e3, self.hist_sorted[self.hist_sorted['nox_slope'].notna()]['demand'].apply(self.f_cedm['nox']), c='#fdbe85', lw=6)
            #simulated slopes
            ax.plot(self.sim_sorted[self.hist_sorted['co2_slope'].notna()]['demand']/1e3, self.sim_sorted[self.sim_sorted['co2_slope'].notna()]['co2_slope'], c='#756bb1', lw=2)
            ax2.plot(self.sim_sorted[self.hist_sorted['so2_slope'].notna()]['demand']/1e3, self.sim_sorted[self.sim_sorted['so2_slope'].notna()]['so2_slope'], c='#31a354', lw=2) 
            ax2.plot(self.sim_sorted[self.hist_sorted['nox_slope'].notna()]['demand']/1e3, self.sim_sorted[self.sim_sorted['nox_slope'].notna()]['nox_slope'], c='#e6550d', lw=2)                         
            #historical slopes
            ax.plot(self.hist_sorted[self.hist_sorted['co2_slope'].notna()]['demand']/1e3, self.hist_sorted[self.hist_sorted['co2_slope'].notna()]['co2_slope'], c='#54278f', ls='--', lw=2)
            ax2.plot(self.hist_sorted[self.hist_sorted['so2_slope'].notna()]['demand']/1e3, self.hist_sorted[self.hist_sorted['so2_slope'].notna()]['so2_slope'], c='#006d2c', ls='--', lw=2)  
            ax2.plot(self.hist_sorted[self.hist_sorted['nox_slope'].notna()]['demand']/1e3, self.hist_sorted[self.hist_sorted['nox_slope'].notna()]['nox_slope'], c='#a63603', ls='--', lw=2)                                       
            #guide lines
            #horizontal lines
            ax.axhline(250, color='black', alpha=0.2, ls='dotted', linewidth=1)
            ax.axhline(500, color='black', alpha=0.2, ls='dotted', linewidth=1)
            ax.axhline(750, color='black', alpha=0.2, ls='dotted', linewidth=1)
            #axes limits and labels
            ax.set_xlabel('Total Fossil Generation [GW]')
            ax.set_ylabel('Marginal CO$_2$ Emissions [kg/MWh]')   
            ax2.set_ylabel('Marginal SO$_2$ & NO$_x$ Emissions [kg/MWh]')   
            ax.set_xlim(self.xlim_tuple)  
            ax.set_ylim(0,1000)
            ax.set_yticks([0, 250, 500, 750, 1000])
            ax2.set_ylim(0,1.6)
            ax2.set_yticks([0, 0.4, 0.8, 1.2, 1.6])
            return f
     
        
    def plotDemandPrices(self, figure_dimensions = (7,5)):
        """ Produces a plot of total fossil generation (x-axis) vs. price (y-axis)
        ---
        returns figure
        """        
        f, ax = matplotlib.pylab.subplots(1, figsize=figure_dimensions) 
        #ax2 = ax.twinx()
        #scatters
        #co2
        ax.scatter(self.hist_sorted['demand']/1e3, self.hist_sorted['gen_cost_marg'], c='#bcbddc', s=1, alpha=0.15)
        #lines
        #co2
        ax.plot(self.sim_sorted[self.hist_sorted['gen_cost_marg_rolling'].notna()]['demand']/1e3, self.sim_sorted[self.sim_sorted['gen_cost_marg_rolling'].notna()]['gen_cost_marg_rolling'], c='#756bb1', lw=2)
        ax.plot(self.hist_sorted[self.hist_sorted['gen_cost_marg_rolling'].notna()]['demand']/1e3, self.hist_sorted[self.hist_sorted['gen_cost_marg_rolling'].notna()]['gen_cost_marg_rolling'], c='#54278f', ls='--', lw=2)
        #axis limits and labels
        ax.set_xlabel('Total Fossil Generation [GW]')
        ax.set_ylabel('Total CO$_2$ Emissions [million kg]\nTotal SO$_2$ & NO$_x$ Emissions [thousand kg]')                  
        ax.set_xlim(self.xlim_tuple)                        
        ax.set_ylim(0, 80)
        return f


    def plot_x_demand(self):
        """ Produces a 2-panel plot with a common x axis (demand). The upper panel shows coal share of marginal generation (simulated rolling averages vs. CEMS rolling averages). The bottom panel shows the marginal emissions factors (simulated rolling averages vs. CEDM deciles)
        ---
        returns figure
        """
        #define the variables        
        x = self.x_sort/1000.
        x_deciles = list(scipy.array(self.deciles_cedm)/1000.)
        c = self.roll_c 
        cstd = self.roll_cstd 
        s = self.roll_s 
        sstd = self.roll_sstd 
        n = self.roll_n 
        nstd = self.roll_nstd 
        coal = self.roll_coal
        coalstd = self.roll_coalstd 
        cems_coal = self.roll_cems_coal 
        cems_coalstd = self.roll_cems_coalstd
        coal_total = self.roll_coal_mix_total
        coal_totalstd = self.roll_coal_mix_totalstd 
        cems_coal_total = self.roll_cems_coal_mix_total 
        cems_coal_totalstd = self.roll_cems_coal_mix_totalstd
        #set up the plot
        matplotlib.pylab.clf()
        #create a 2 panel plot sharing the same x axis
        f, axarr = matplotlib.pylab.subplots(3, sharex=True, figsize=(4,8))      
        #plot the total coal mix
        axarr[0].plot(x,cems_coal_total, c='grey', ls='--') 
        #axarr[0].fill_between(x,scipy.minimum((cems_coal_total + cems_coal_totalstd), 1),scipy.maximum((cems_coal_total - cems_coal_totalstd), 0), color='grey', alpha=0.15)   
        axarr[0].plot(x,coal_total, c='#fc8d59') 
        #axarr[0].fill_between(x,scipy.minimum((coal_total + coal_totalstd), 1),scipy.maximum((coal_total - coal_totalstd), 0), color='#fc8d59', alpha=0.15) 
        axarr[0].set_xlim(self.xlim_tuple)
        axarr[0].set_ylim(0,1.0)
        axarr[0].set_ylabel('Coal Share of \nTotal Generation')    
        #add the legend
        #coal_line = matplotlib.lines.Line2D([], [], color='#fc8d59', linewidth=2.5, label='Simulated')
        #coal_cems_line = matplotlib.lines.Line2D([], [], color='#91bfdb', linewidth=2.5, label='Historical')
        #stdline = matplotlib.patches.Patch(color='black', alpha=0.15, label='Std. Dev.')     
        #axarr[0].legend(handles=[coal_line, coal_cems_line, stdline], labelspacing=0.1) 
        #plot the betas on the lower plot
        #the betas plot will require 2 y axis, so make a twin of axarr[1]
        #plot the marginal coal mix
        axarr[1].plot(x,cems_coal, c='grey', ls='--') 
        #axarr[1].fill_between(x,scipy.minimum((cems_coal + cems_coalstd), 1),scipy.maximum((cems_coal - cems_coalstd), 0), color='grey', alpha=0.15)   
        axarr[1].plot(x,coal, c='#fc8d59') 
        #axarr[1].fill_between(x,scipy.minimum((coal + coalstd), 1),scipy.maximum((coal - coalstd), 0), color='#fc8d59', alpha=0.15) 
        axarr[1].set_xlim(self.xlim_tuple)
        axarr[1].set_ylim(0,1.0)
        axarr[1].set_ylabel('Coal Share of \nMarginal Generation')    
        #plot the marginal co2
        ax2 = axarr[2].twinx()
        axarr[2].plot(x_deciles, self.mefs_cedm_co2, c='#7570b3', ls='--')
        axarr[2].plot(x,c, c='#7570b3') 
        #axarr[2].fill_between(x,(c+cstd),scipy.maximum((c - cstd), 0), color='blue', alpha=0.1) 
        axarr[2].axhline(500, color='black', alpha=0.2, ls='dotted', linewidth=1)
        axarr[2].axhline(1000, color='black', alpha=0.2, ls='dotted', linewidth=1)
        axarr[2].set_xlim(self.xlim_tuple)
        axarr[2].set_ylim(0,1100)
        axarr[2].set_yticks([0, 500, 1000])
        axarr[2].set_xlabel('Total Fossil Generation [GW]')
        axarr[2].set_ylabel('Marg. CO$_2$ [kg/MWh]')
        dem_min = self.df.demand.min()
        dem_max = self.df.demand.max()
        if (dem_max - dem_min) < 30000:
            vline_step = 5000
        else: 
            vline_step = 10000
        vline_start = math.ceil(dem_min/vline_step)*vline_step/1000
        vline_end = math.floor(dem_max/vline_step)*vline_step/1000 + 1
        vline = scipy.arange(vline_start, vline_end, vline_step/1000.)
        for v in vline:   
            axarr[2].axvline(v, color='black', alpha=0.2, ls='dotted', linewidth=1)
        #plot the marginal so2 and nox
        ax2.plot(x_deciles, self.mefs_cedm_so2, c='#1b9e77', ls='--')       
        ax2.plot(x,s, c='#1b9e77') 
        #ax2.fill_between(x,(s+sstd),scipy.maximum((s - sstd), 0), color='green', alpha=0.1)          
        ax2.plot(x_deciles, self.mefs_cedm_nox, c='#d95f02', ls='--') 
        ax2.plot(x,n, c='#d95f02') 
        #ax2.fill_between(x,(n+nstd),scipy.maximum((n - nstd), 0), color='red', alpha=0.1)   
        ax2.set_xlim(self.xlim_tuple)
        ax2.set_ylim(0,2.2)    
        ax2.set_yticks([0, 1, 2])
        ax2.set_ylabel('Marg. SO$_2$ / NO$_x$ [kg/MWh]')
        #add the legend
        #co2line = matplotlib.lines.Line2D([], [], color='#7570b3', linewidth=2.5, label='CO$_2$')
        #so2line = matplotlib.lines.Line2D([], [], color='#1b9e77', linewidth=2.5, label='SO$_2$')
        #noxline = matplotlib.lines.Line2D([], [], color='#d95f02', linewidth=2.5, label='NO$_x$')
        #stdline = matplotlib.patches.Patch(color='black', alpha=0.15, label='Std. Dev.')
        #cedmline = matplotlib.lines.Line2D([], [], color='gray', linewidth=1.5, ls='--', label='Hist. Regress.')       
        #ax2.legend(handles=[co2line, so2line, noxline, stdline, cedmline], labelspacing=0.1)
        #return
        matplotlib.pylab.tight_layout()
        return f
    
    
    def plot_mefs(self, figure_dimensions=(5,4), rolling_window=500):
        """ Produces a plot with the CEDM MEFs, Rolling-Simulated MEFs, and 25th/75th percentile historical MEFs. The bottom panel shows the marginal emissions factors (simulated rolling averages vs. CEDM deciles)
        ---
        returns figure
        """    
        #define the variables      
        p = rolling_window
        x = self.x_sort/1000.
        x_deciles = list(scipy.array(self.deciles_cedm)/1000.)
        c = self.roll_c 
        s = self.roll_s 
        n = self.roll_n 
        cems_df_sorted = self.cems_df.sort_values(by='demand').copy()
        cems_df_sorted = cems_df_sorted.replace(999999, scipy.nan)
        c_cems_25 = cems_df_sorted.co2_marg.rolling(window=p, min_periods=20, center=True).quantile(0.25)
        c_cems_75 = cems_df_sorted.co2_marg.rolling(window=p, min_periods=20, center=True).quantile(0.75)
        s_cems_25 = cems_df_sorted.so2_marg.rolling(window=p, min_periods=20, center=True).quantile(0.25)
        s_cems_75 = cems_df_sorted.so2_marg.rolling(window=p, min_periods=20, center=True).quantile(0.75)
        n_cems_25 = cems_df_sorted.nox_marg.rolling(window=p, min_periods=20, center=True).quantile(0.25)
        n_cems_75 = cems_df_sorted.nox_marg.rolling(window=p, min_periods=20, center=True).quantile(0.75)
              
        #set up the plot
        matplotlib.pylab.clf()
        
        f, ax = matplotlib.pylab.subplots(1, figsize=figure_dimensions)      
        
        ax2 = ax.twinx()
        #plot the marginal co2
        ax.plot(x_deciles, self.mefs_cedm_co2, c='#7570b3', ls='--')
        ax.plot(x,c, c='#7570b3') 
        ax.fill_between(x,c_cems_25,c_cems_75, color='#7570b3', alpha=0.1) 
             
        ax.axhline(250, color='black', alpha=0.2, ls='dotted', linewidth=1)
        ax.axhline(500, color='black', alpha=0.2, ls='dotted', linewidth=1)
        ax.axhline(750, color='black', alpha=0.2, ls='dotted', linewidth=1)
        ax.set_xlim(self.xlim_tuple)
        ax.set_ylim(0,950)
        ax.set_yticks([0, 250, 500, 750])
        ax.set_xlabel('Total Fossil Generation [GW]')
        ax.set_ylabel('Marg. CO$_2$ [kg/MWh]')
        dem_min = self.df.demand.min()
        dem_max = self.df.demand.max()
        if (dem_max - dem_min) < 30000:
            vline_step = 5000
        else: 
            vline_step = 10000
        vline_start = math.ceil(dem_min/vline_step)*vline_step/1000
        vline_end = math.floor(dem_max/vline_step)*vline_step/1000 + 1
        vline = scipy.arange(vline_start, vline_end, vline_step/1000.)
        for v in vline:   
            ax.axvline(v, color='black', alpha=0.2, ls='dotted', linewidth=1)
        #plot the marginal so2 and nox
        ax2.plot(x_deciles, self.mefs_cedm_so2, c='#1b9e77', ls='--')       
        ax2.plot(x,s, c='#1b9e77') 
        ax2.fill_between(x,s_cems_25,s_cems_75, color='#1b9e77', alpha=0.1)       
        ax2.plot(x_deciles, self.mefs_cedm_nox, c='#d95f02', ls='--') 
        ax2.plot(x,n, c='#d95f02') 
        ax2.fill_between(x,n_cems_25,n_cems_75, color='#d95f02', alpha=0.1)  
        ax2.set_xlim(self.xlim_tuple)
        ax2.set_ylim(0,1.9)    
        ax2.set_yticks([0, 0.5, 1, 1.5])
        ax2.set_ylabel('Marg. SO$_2$ / NO$_x$ [kg/MWh]')
        #add the legend
        #co2line = matplotlib.lines.Line2D([], [], color='#7570b3', linewidth=2.5, label='CO$_2$')
        #so2line = matplotlib.lines.Line2D([], [], color='#1b9e77', linewidth=2.5, label='SO$_2$')
        #noxline = matplotlib.lines.Line2D([], [], color='#d95f02', linewidth=2.5, label='NO$_x$')
        #stdline = matplotlib.patches.Patch(color='black', alpha=0.15, label='Std. Dev.')
        #cedmline = matplotlib.lines.Line2D([], [], color='gray', linewidth=1.5, ls='--', label='Hist. Regress.')       
        #ax2.legend(handles=[co2line, so2line, noxline, stdline, cedmline], labelspacing=0.1)
        #return
        matplotlib.pylab.tight_layout()
        return f
   

    def plot_hist_vs_simulated(self, x_property='demand_smooth', y_property='co2_tot', demand_smooth_step=500, start_date='2017-01-01 00:00:00', end_date='2018-01-01 00:00:00'):
        """ Produces a plot that compares historical with simulated results. Arguments allow different results to be plotted over different time windows with different x axis values. Line shows the mean. Shading shows the +- standard deviation from the mean. Default inputs plots total co2 emissions vs total demand for the whole year.
        ---
        x_property : 'demand_smooth', 'hour', 'week', etc.
        y_property : shared columns between self.df and self.cems_df (e.g. 'co2_tot', 'coal_mix')
        w0 : first week in observation
        wf : final week in observation
        demand_smooth_step : the demand_smooth column will equal demand rounded down to the nearest demand_smooth_step
        returns figure
        """
        simu = self.df.copy(deep=True)
        hist = self.cems_df.copy(deep=True)
        #only select the data between the start_date and end_date
        simu = simu[(simu.datetime>=start_date) & (simu.datetime<=end_date)]
        hist = hist[(hist.datetime>=start_date) & (hist.datetime<=end_date)]
        #add hour, week, and demand_smooth columns, where demand_smooth is the demand rounded down to the nearest demand_smooth_step
        simu.datetime = pandas.to_datetime(simu.datetime)
        simu['hour'] = simu.datetime.dt.hour
        simu['week'] = simu.datetime.dt.week
        simu['demand_smooth'] = simu.demand // demand_smooth_step * demand_smooth_step #round down to nearest step
        hist.datetime = pandas.to_datetime(hist.datetime)
        hist['hour'] = hist.datetime.dt.hour
        hist['week'] = hist.datetime.dt.week
        hist['demand_smooth'] = hist.demand // demand_smooth_step * demand_smooth_step #round down to nearest step
        #only select the data between the start_date and end_date
        #simu = simu[(simu.datetime>=start_date) & (simu.datetime<=end_date)]
        #hist = hist[(hist.datetime>=start_date) & (hist.datetime<=end_date)]
        #aggregate to get the mean and standard deviation for the y values
        tmean = simu.groupby([x_property], as_index=False).mean()
        tstd = simu.groupby([x_property], as_index=False).std()     
        hmean = hist.groupby([x_property], as_index=False).mean()
        hstd = hist.groupby([x_property], as_index=False).std()
        #mess with the limits, ticks, labels, legend
        xmult = 1.0
        ymult = 1.0
        #y_axis values
        if y_property == 'gen_cost_marg':
            ylim_arr = (0,120)
            ytitle = 'Electricity Price [$/MWh]'
        if y_property == 'co2_tot':
            ymult = 0.000001
            ylim_arr = (0)
            ytitle = 'CO$_2$ Emissions [million kg/hr]'  
        if y_property == 'so2_tot':
            ymult = 0.001
            ylim_arr = (0)
            ytitle = 'SO$_2$ Emissions [thousand kg/hr]'
        if y_property == 'nox_tot':
            ymult = 0.001
            ylim_arr = (0)
            ytitle = 'NO$_X$ Emissions [thousand kg/hr]'
        if y_property == 'co2_marg':
            ymult = 1.0
            ylim_arr = (0)
            ytitle = 'Marginal CO$_2$ Emissions [kg/MWh]'  
        if y_property == 'so2_marg':
            ymult = 1.0
            ylim_arr = (0)
            ytitle = 'arginal SO$_2$ Emissions [kg/MWh]'
        if y_property == 'nox_marg':
            ymult = 1.0
            ylim_arr = (0)
            ytitle = 'arginal NO$_X$ Emissions [kg/MWh]'
        #x_axis values
        if x_property == 'demand_smooth':
            xmult = 0.001
            xtitle = 'Fossil Fuel Generation [GW]'
        #setup x and y for the plot
        x = tmean[x_property]*xmult
        ymean = tmean[y_property]*ymult  
        ystd = tstd[y_property]*ymult 
        xh = hmean[x_property]*xmult
        ymeanh = hmean[y_property]*ymult  
        ystdh = hstd[y_property]*ymult   
        #calculate pearson correlation between ymean and ymeanh
        #slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
        self.simu = simu
        self.hist = hist
        rs2 = scipy.stats.linregress(simu[y_property],hist[y_property]).rvalue.round(3)
        #create the plot           
        f = matplotlib.pylab.figure(figsize=(5,4))
        ax = f.add_subplot(111)
        #ax2 = ax.twinx()
        #plot the historical data first (underneath the simulated results)
        ax.plot(xh,ymeanh, c='black', ls='--')
        ax.fill_between(xh,(ymeanh + ystdh),(ymeanh - ystdh), color='gray', alpha=0.25) 
        #plot the simulated results
        ax.plot(x,ymean, 'k-')
        ax.fill_between(x,(ymean + ystd),(ymean - ystd), color='red', alpha=0.25) 
        #set limits and titles
        ax.set_xlim(self.xlim_tuple)
        ax.set_ylim(ylim_arr)
        ax.set_ylabel(ytitle)
        if x_property == 'hour':
            ax.set_xticks((0, 4, 8, 12, 16, 20, 24))
            xtitle = 'Hour of Day'
        ax.set_xlabel(xtitle)
        #add the pearson correlation text
        matplotlib.pylab.text(0.49, 0.05, 'R$^2$ Fit of Hourly Data = %.2f'%(rs2), transform=ax.transAxes, bbox=dict(facecolor='white', alpha=0.5))
        #matplotlib.pylab.text(0.54, 0.05, 'R$^2$ of Hourly Data = %.2f\nR$^2$ of Mean Lines  = %.2f'%(rs2,rs1), transform=ax.transAxes, bbox=dict(facecolor='white', alpha=0.5)) #including the R-squared of the mean lines. Not sure if that is relevant information...
        #add the legend
        #histline = matplotlib.lines.Line2D([], [], color='black', ls='--', linewidth=2.5, label='Historical Mean')
        #histstd = matplotlib.patches.Patch(color='black', alpha=0.25, label='Historical Std. Dev.')
        #simuline = matplotlib.lines.Line2D([], [], color='black', linewidth=2.5, label='Simulated Mean')
        #simustd = matplotlib.patches.Patch(color='red', alpha=0.25, label='Historical Std. Dev.')
        #ax.legend(handles=[histline, histstd, simuline, simustd], labelspacing=0.1)
        #return
        matplotlib.pylab.tight_layout()
        return f


    def plot_density_function(self, den_fun='cumulative', series='data', x_property='co2_tot', start_date='2014-01-01 00:00:00', end_date='2015-01-01 00:00:00', bin_no=100, x_range=[-scipy.inf, scipy.inf], error_range=[-1,2]):
            """ Produces a cdf for the simulated and the historical data
            ---
            den_fun : 'cumulative' for a cdf or 'probability' for pdf
            series : 'data' = plot the historical and simulated data. 'error' = plot the error between the simulated and historical data
            x_property : shared columns between self.df and self.cems_df (e.g. 'co2_tot', 'coal_mix')
            bin_no : the number of bins used in the histogram calculation
            x_range : the limits of the x axis
            error_range : the limits of the error plot's x axis (in units of %*100)
            returns figure
            """
            #copy in data
            simu = self.df.copy(deep=True)
            hist = self.cems_df.copy(deep=True)
            #slice according to dates
            simu = simu[(simu.datetime>=start_date) & (simu.datetime<=end_date)][['datetime', x_property]]
            hist = hist[(hist.datetime>=start_date) & (hist.datetime<=end_date)][x_property]
            #delete all of the rows with any 999999 values. These are from places in the historical data where MEFs were calculated to be negative or where they were calculated to be impossibly large. So they are outliers for the sake of hourly comparison
            simu['hist'] = hist
            simu = simu[simu['hist'] != 999999]
            hist = simu['hist']
            simu = simu[x_property]
            #in case we get have really high number that we want to bin at the plot tail (i.e. avoid having really long tails in the plot)
            simu[simu > x_range[1]] = x_range[1]
            hist[hist > x_range[1]] = x_range[1]
            simu[simu < x_range[0]] = x_range[0]
            hist[hist < x_range[0]] = x_range[0]        
            #mess with the limits, ticks, labels, legend
            xmult = 1.0
            #x_axis values
            if series =='data':
                xtitle = scipy.where(x_property=='gen_cost_marg', 'Hourly Electricity Price [$/MWh]', scipy.where(x_property=='co2_tot','Hourly CO$_2$ Emissions [million kg/hr]' , scipy.where(x_property=='so2_tot','Hourly SO$_2$ Emissions [thousand kg/hr]', scipy.where(x_property=='nox_tot','Hourly NO$_X$ Emissions [thousand kg/hr]', scipy.where(x_property=='co2_marg','Hourly Marginal CO$_2$ Emissions [kg/MWh]', scipy.where(x_property=='so2_marg','Hourly Marginal SO$_2$ Emissions [kg/MWh]', scipy.where(x_property=='nox_marg','Hourly Marginal NO$_X$ Emissions [kg/MWh]','X Label Error')))))))
                xmult = scipy.where(x_property=='gen_cost_marg', 1.0, scipy.where(x_property=='co2_tot',0.000001 , scipy.where(x_property=='so2_tot',0.001, scipy.where(x_property=='nox_tot',0.001,1.0))))
            if series =='error':
                xtitle = scipy.where(x_property=='gen_cost_marg', 'Hourly Electricity Price Error [%]', scipy.where(x_property=='co2_tot','Hourly CO$_2$ Emissions Error [%]' , scipy.where(x_property=='so2_tot','Hourly SO$_2$ Emissions Error [%]', scipy.where(x_property=='nox_tot','Hourly NO$_X$ Emissions Error [%]', scipy.where(x_property=='co2_marg','Hourly Marginal CO$_2$ Emissions Error [kg/MWh]', scipy.where(x_property=='so2_marg','Hourly Marginal SO$_2$ Emissions Error [kg/MWh]', scipy.where(x_property=='nox_marg','Hourly Marginal NO$_X$ Emissions Error [kg/MWh]','X Label Error')))))))
                if x_property[4:]!='marg':
                    xmult = 100.0 #converts fraction to percentage
            #error calculation
            if x_property[4:]=='marg':
                err = simu - hist #absolute error. Since marginal emissions can approach zero, relative error can be extremely large
            else:
                err = scipy.divide((simu-hist), hist) #relative error. Since total emissions do not approach zero, relative error is a good metric.
            #data for the plot
            histogram, bins = scipy.where(den_fun=='probability',scipy.histogram(simu, bins=bin_no, density=True) ,scipy.histogram(simu, bins=bin_no))
            histogram2, bins2 = scipy.where(den_fun=='probability',scipy.histogram(hist, bins=bin_no, density=True) ,scipy.histogram(hist, bins=bin_no))           
            
            #in case scipy.divide produces any -inf or inf values, or we just have really high number that we want to bin at the plot tail
            err[err > error_range[1]] = error_range[1]
            err[err < error_range[0]] = error_range[0] 
            err.dropna(inplace=True)
            histogram_error, bins_error = scipy.where(den_fun=='probability',scipy.histogram(err, bins=bin_no, density=True) ,scipy.histogram(err, bins=bin_no))
            #calculate the x and y data   
            if series == 'data':
                #cdf or pdf
                if den_fun == 'cumulative':
                    cum = scipy.cumsum(histogram)/float(len(simu))
                    cum2 = scipy.cumsum(histogram2)/float(len(simu))
                if den_fun == 'probability':
                    cum = histogram
                    cum2 = histogram2
                center = (bins[:-1] + bins[1:]) / 2 * xmult
                center2 = (bins2[:-1] + bins2[1:]) / 2 * xmult
            if series =='error':
                #cdf or pdf
                if den_fun == 'cumulative':
                    cum = scipy.cumsum(histogram_error)/float(len(simu))
                if den_fun == 'probability':
                    cum = histogram_error
                center = (bins_error[:-1] + bins_error[1:]) / 2 * xmult   
            #create the plot           
            f = matplotlib.pylab.figure(figsize=(5,4))
            ax = f.add_subplot(111)
            if series == 'data':
                #plot the historical data first (underneath the simulated results)
                ax.plot(center2,cum2, c='black', ls='--')
                #ax.set_xlim(x_range[0], x_range[1])
            #plot the simulated results
            ax.plot(center,cum, 'k-')
            #set limits and titles
            ax.set_xlabel(xtitle)
            ax.set_ylim(0)
            #set xlim for error plots
            if series == 'error':
                ax.set_xlim(error_range[0]*xmult, error_range[1]*xmult)
                if x_property[4:]=='marg':
                    if x_property[:3]!='co2':
                        matplotlib.pylab.text(0.025, 0.73, 'Error [kg/MWh]:\n5$^{th}$ Pctl.    %.2f\nMedian      %.2f\n95$^{th}$ Pctl.   %.2f'%(err.quantile(0.05)*xmult, err.median()*xmult,  err.quantile(0.95)*xmult), transform=ax.transAxes, bbox=dict(facecolor='white', alpha=0.5))
                    else:
                        matplotlib.pylab.text(0.025, 0.73, 'Error [kg/MWh]:\n5$^{th}$ Pctl.    %.0f\nMedian      %.0f\n95$^{th}$ Pctl.   %.0f'%(err.quantile(0.05)*xmult, err.median()*xmult,  err.quantile(0.95)*xmult), transform=ax.transAxes, bbox=dict(facecolor='white', alpha=0.5))
                else:
                    matplotlib.pylab.text(0.025, 0.73, 'Relative Error [%]' + '\n5$^{th}$ Pctl.    %.0f\nMedian      %.0f\n95$^{th}$ Pctl.   %.0f'%(err.quantile(0.05)*xmult, err.median()*xmult,  err.quantile(0.95)*xmult), transform=ax.transAxes, bbox=dict(facecolor='white', alpha=0.5))
            #return
            matplotlib.pylab.tight_layout()
            return f





class plotDispatchMultiple(object):
    def __init__(self, pd_list):
        """ 
        Similar plotting functions as plotDispatch but handles multiple plotDispatchObjects
        ---
        pd_list : a list of plotDispatch objects
        """
        self.pd_list = pd_list
          
        
    def plot_density_function(self, den_fun='cumulative', series='data', x_property='co2_tot', start_date='2014-01-01 00:00:00', end_date='2015-01-01 00:00:00', bin_no=100, x_range=[-scipy.inf, scipy.inf], error_range=[-1,2]):
            """ Produces a cdf for the simulated and the historical data
            ---
            den_fun : 'cumulative' for a cdf or 'probability' for pdf
            series : 'data' = plot the historical and simulated data. 'error' = plot the error between the simulated and historical data
            x_property : shared columns between self.df and self.cems_df (e.g. 'co2_tot', 'coal_mix')
            bin_no : the number of bins used in the histogram calculation
            x_range : the limits of the x axis
            error_range : the limits of the error plot's x axis (in units of %*100)
            returns figure
            """
            #create dictionary to store plot data
            colors = {0: '#e66101',1: '#fdb863', 2:'#b2abd2', 3: '#5e3c99'}
            simus = {}
            hists = {}
            errs = {}
            cums = {}
            cum2s = {}
            centers = {}
            center2s = {}
            #loop through the different plodDispatch objects
            pd_count = 0
            for pd_local in self.pd_list: 
                
                #copy in data
                simu = pd_local.df.copy(deep=True)
                hist = pd_local.cems_df.copy(deep=True)
                #slice according to dates
                simu = simu[(simu.datetime>=start_date) & (simu.datetime<=end_date)][['datetime', x_property]]
                hist = hist[(hist.datetime>=start_date) & (hist.datetime<=end_date)][x_property]
                #delete all of the rows with any 999999 values. These are from places in the historical data where MEFs were calculated to be negative or where they were calculated to be impossibly large. So they are outliers for the sake of hourly comparison
                simu['hist'] = hist
                simu = simu[simu['hist'] != 999999]
                hist = simu['hist']
                simu = simu[x_property]
                #in case we get have really high number that we want to bin at the plot tail (i.e. avoid having really long tails in the plot)
                simu[simu > x_range[1]] = x_range[1]
                hist[hist > x_range[1]] = x_range[1]
                simu[simu < x_range[0]] = x_range[0]
                hist[hist < x_range[0]] = x_range[0]          
                #error calculation
                if x_property[4:]=='marg':
                    err = simu - hist #absolute error. Since marginal emissions can approach zero, relative error can be extremely large
                else:
                    err = scipy.divide((simu-hist), hist) #relative error. Since total emissions do not approach zero, relative error is a good metric.
                #data for the plot
                histogram, bins = scipy.where(den_fun=='probability',scipy.histogram(simu, bins=bin_no, density=True) ,scipy.histogram(simu, bins=bin_no))
                histogram2, bins2 = scipy.where(den_fun=='probability',scipy.histogram(hist, bins=bin_no, density=True) ,scipy.histogram(hist, bins=bin_no))           
                
                #in case scipy.divide produces any -inf or inf values, or we just have really high number that we want to bin at the plot tail
                err[err > error_range[1]] = error_range[1]
                err[err < error_range[0]] = error_range[0] 
                err.dropna(inplace=True)
                histogram_error, bins_error = scipy.where(den_fun=='probability',scipy.histogram(err, bins=bin_no, density=True) ,scipy.histogram(err, bins=bin_no))
                
                #mess with the limits, ticks, labels, legend
                xmult = 1.0
                #x_axis values
                if series =='data':
                    xtitle = scipy.where(x_property=='gen_cost_marg', 'Hourly Electricity Price [$/MWh]', scipy.where(x_property=='co2_tot','Hourly CO$_2$ Emissions [million kg/hr]' , scipy.where(x_property=='so2_tot','Hourly SO$_2$ Emissions [thousand kg/hr]', scipy.where(x_property=='nox_tot','Hourly NO$_X$ Emissions [thousand kg/hr]', scipy.where(x_property=='co2_marg','Hourly Marginal CO$_2$ Emissions [kg/MWh]', scipy.where(x_property=='so2_marg','Hourly Marginal SO$_2$ Emissions [kg/MWh]', scipy.where(x_property=='nox_marg','Hourly Marginal NO$_X$ Emissions [kg/MWh]','X Label Error')))))))
                    xmult = scipy.where(x_property=='gen_cost_marg', 1.0, scipy.where(x_property=='co2_tot',0.000001 , scipy.where(x_property=='so2_tot',0.001, scipy.where(x_property=='nox_tot',0.001,1.0))))
                if series =='error':
                    xtitle = scipy.where(x_property=='gen_cost_marg', 'Hourly Electricity Price Error [%]', scipy.where(x_property=='co2_tot','Hourly CO$_2$ Emissions Error [%]' , scipy.where(x_property=='so2_tot','Hourly SO$_2$ Emissions Error [%]', scipy.where(x_property=='nox_tot','Hourly NO$_X$ Emissions Error [%]', scipy.where(x_property=='co2_marg','Hourly Marginal CO$_2$ Emissions Error [kg/MWh]', scipy.where(x_property=='so2_marg','Hourly Marginal SO$_2$ Emissions Error [kg/MWh]', scipy.where(x_property=='nox_marg','Hourly Marginal NO$_X$ Emissions Error [kg/MWh]','X Label Error')))))))
                    if x_property[4:]!='marg':
                        xmult = 100.0 #converts fraction to percentage
                
                #calculate the x and y data   
                if series == 'data':
                    #cdf or pdf
                    if den_fun == 'cumulative':
                        cum = scipy.cumsum(histogram)/float(len(simu))
                        cum2 = scipy.cumsum(histogram2)/float(len(simu))
                    if den_fun == 'probability':
                        cum = histogram
                        cum2 = histogram2
                    center = (bins[:-1] + bins[1:]) / 2 * xmult
                    center2 = (bins2[:-1] + bins2[1:]) / 2 * xmult
                if series =='error':
                    #cdf or pdf
                    if den_fun == 'cumulative':
                        cum = scipy.cumsum(histogram_error)/float(len(simu))
                    if den_fun == 'probability':
                        cum = histogram_error
                    center = (bins_error[:-1] + bins_error[1:]) / 2 * xmult   
                    cum2 = 0
                    center2 = 0
                #store in dictionaries  
                simus[pd_count] = simu
                hists[pd_count] = hist
                errs[pd_count] = err
                cums[pd_count] = cum 
                cum2s[pd_count] = cum2
                centers[pd_count] = center
                center2s[pd_count] = center2
                #loop counter update 
                pd_count += 1
            
            #create the plot           
            f = matplotlib.pylab.figure(figsize=(5,4))
            ax = f.add_subplot(111)
            
            for n in [0,1,2,3]:
                if series == 'data':
                    #plot the historical data first (underneath the simulated results)
                    ax.plot(center2s[n],cum2s[n], c=colors[n], ls=(0, (1, 4, 6, 4)))
                    #ax.set_xlim(x_range[0], x_range[1])
                #plot the simulated results
                ax.plot(centers[n],cums[n], c=colors[n], ls='-')
                #set limits and titles
                ax.set_xlabel(xtitle)
                ax.set_ylim(0)
                #set xlim for error plots
                if series == 'error':
                    ax.set_xlim(error_range[0]*xmult, error_range[1]*xmult)
            matplotlib.pylab.tight_layout()
            return f


    
    