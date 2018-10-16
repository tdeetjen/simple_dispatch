# mefs_from_dispatch
# Thomas Deetjen
# v12
# last edited: 2018-10-16
# class "generate_mefs" creates a dataframe and plot of the MEFs for a given dispatch csv


import pandas
import matplotlib.pylab
import scipy
import scipy.linalg
import scipy.stats
import math
import os





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
        mef_max = scipy.where(e_type=='co2', 1500.0/0.454, scipy.where(e_type=='so2', 3.0/0.454, scipy.where(e_type=='nox', 2.0/0.454, 1.0e10)))
        self.df.loc[self.df[e_type+'_marg'] <0, (e_type+'_marg')] = scipy.nan
        self.df.loc[self.df[e_type+'_marg'] >mef_max, (e_type+'_marg')] = scipy.nan
    




class plotDispatch(object):
    def __init__(self, dispatch_df, cems_df, deciles_cedm=[21300, 24073, 26245, 28132, 30081, 32297, 35414, 39355, 45364, 60650], mefs_cedm_co2 = [691, 623, 589, 562, 550, 529, 515, 485, 479, 489], mefs_cedm_so2 = [1.156, 0.951, 0.835, 0.745, 0.731, 0.586, 0.445, 0.324, 0.240, 0.003], mefs_cedm_nox = [0.292, 0.255, 0.238, 0.221, 0.218, 0.191, 0.198, 0.192, 0.232, 0.379]):
        """ 
        Contains a few functions used to plot the dispatch results
        ---
        dispatch_df : the dataframe version of a csv file containing the dispatch output of the simple_dispatch.dispatch
        cems_df : the CEMS '2014_%s_hourly_demand_and_fuelmix.csv'%s(nerc_region) csv file showing fuel_mix in the cems data
        deciles : the demand points that we want to calculate seperate MEFs for [MW]
        """
        self.df = dispatch_df.copy(deep=True)
        self.cems_df = cems_df
        self.deciles_cedm = deciles_cedm
        self.mefs_cedm_co2 = mefs_cedm_co2
        self.mefs_cedm_so2 = mefs_cedm_so2
        self.mefs_cedm_nox = mefs_cedm_nox
        self.add_dispatch_columns()
        self.rolling_calculations()
          
        
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


    def rolling_calculations(self, p=500):
        """ Calculates rolling average and standard deviations for plotting
        ---
        p : the number of values included in the rolling calculation window
        """	
        #df_sort = self.df[['demand','co2_marg', 'so2_marg', 'nox_marg', 'coal_is_marginal']].sort_values('demand') 
        df_sort = self.df[['demand','co2_marg', 'so2_marg', 'nox_marg', 'coal_mix_marg', 'coal_mix_p']].sort_values('demand') 
        cems_df_sort = self.cems_df[['demand', 'p_marginal_coal', 'coal_mix_p']].sort_values('demand') 
        self.x_sort = df_sort.demand
        self.roll_c = df_sort.co2_marg.rolling(window=p, min_periods=20).mean()
        self.roll_cstd = df_sort.co2_marg.rolling(window=p, min_periods=20).std()
        self.roll_s = df_sort.so2_marg.rolling(window=p, min_periods=20).mean()
        self.roll_sstd = df_sort.so2_marg.rolling(window=p, min_periods=20).std()
        self.roll_n = df_sort.nox_marg.rolling(window=p, min_periods=20).mean()
        self.roll_nstd = df_sort.nox_marg.rolling(window=p, min_periods=20).std()
        #self.roll_coal = df_sort.coal_is_marginal.rolling(window=p, min_periods=20).mean()
        #self.roll_coalstd = df_sort.coal_is_marginal.rolling(window=p, min_periods=20).std()
        self.roll_coal = df_sort.coal_mix_marg.rolling(window=p, min_periods=20).mean()
        self.roll_coalstd = df_sort.coal_mix_marg.rolling(window=p, min_periods=20).std()   
        self.roll_coal_mix_total = df_sort.coal_mix_p.rolling(window=p, min_periods=20).mean()
        self.roll_coal_mix_totalstd = df_sort.coal_mix_p.rolling(window=p, min_periods=20).std()    
        #cems
        self.roll_cems_coal = cems_df_sort.p_marginal_coal.rolling(window=p, min_periods=20).mean()
        self.roll_cems_coalstd = cems_df_sort.p_marginal_coal.rolling(window=p, min_periods=20).std()
        self.roll_cems_coal_mix_total = cems_df_sort.coal_mix_p.rolling(window=p, min_periods=20).mean()
        self.roll_cems_coal_mix_totalstd = cems_df_sort.coal_mix_p.rolling(window=p, min_periods=20).std()


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
        axarr[0].fill_between(x,scipy.minimum((cems_coal_total + cems_coal_totalstd), 1),scipy.maximum((cems_coal_total - cems_coal_totalstd), 0), color='grey', alpha=0.15)   
        axarr[0].plot(x,coal_total, c='#fc8d59') 
        axarr[0].fill_between(x,scipy.minimum((coal_total + coal_totalstd), 1),scipy.maximum((coal_total - coal_totalstd), 0), color='#fc8d59', alpha=0.15) 
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
        axarr[1].fill_between(x,scipy.minimum((cems_coal + cems_coalstd), 1),scipy.maximum((cems_coal - cems_coalstd), 0), color='grey', alpha=0.15)   
        axarr[1].plot(x,coal, c='#fc8d59') 
        axarr[1].fill_between(x,scipy.minimum((coal + coalstd), 1),scipy.maximum((coal - coalstd), 0), color='#fc8d59', alpha=0.15) 
        axarr[1].set_ylim(0,1.0)
        axarr[1].set_ylabel('Coal Share of \nMarginal Generation')    
        #plot the marginal co2
        ax2 = axarr[2].twinx()
        axarr[2].plot(x_deciles, self.mefs_cedm_co2, c='#7570b3', ls='--')
        axarr[2].plot(x,c, c='#7570b3') 
        axarr[2].fill_between(x,(c+cstd),scipy.maximum((c - cstd), 0), color='blue', alpha=0.1) 
        axarr[2].axhline(500, color='black', alpha=0.2, ls='dotted', linewidth=1)
        axarr[2].axhline(1000, color='black', alpha=0.2, ls='dotted', linewidth=1)
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
        ax2.fill_between(x,(s+sstd),scipy.maximum((s - sstd), 0), color='green', alpha=0.1)          
        ax2.plot(x_deciles, self.mefs_cedm_nox, c='#d95f02', ls='--') 
        ax2.plot(x,n, c='#d95f02') 
        ax2.fill_between(x,(n+nstd),scipy.maximum((n - nstd), 0), color='red', alpha=0.1)      
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
        rs1 = scipy.stats.linregress(ymean,ymeanh).rvalue.round(3)
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
            for pd in self.pd_list: 
                
                #copy in data
                simu = pd.df.copy(deep=True)
                hist = pd.cems_df.copy(deep=True)
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
# =============================================================================
#                     if x_property[4:]=='marg':
#                         if x_property[:3]!='co2':
#                             matplotlib.pylab.text(0.025, 0.73, 'Error [kg/MWh]:\n5$^{th}$ Pctl.    %.2f\nMedian      %.2f\n95$^{th}$ Pctl.   %.2f'%(errs[n].quantile(0.05)*xmult, errs[n].median()*xmult,  errs[n].quantile(0.95)*xmult), transform=ax.transAxes, bbox=dict(facecolor='white', alpha=0.5))
#                         else:
#                             matplotlib.pylab.text(0.025, 0.73, 'Error [kg/MWh]:\n5$^{th}$ Pctl.    %.0f\nMedian      %.0f\n95$^{th}$ Pctl.   %.0f'%(errs[n].quantile(0.05)*xmult, errs[n].median()*xmult,  errs[n].quantile(0.95)*xmult), transform=ax.transAxes, bbox=dict(facecolor='white', alpha=0.5))
#                     else:
#                         matplotlib.pylab.text(0.025, 0.73, 'Relative Error [%]' + '\n5$^{th}$ Pctl.    %.0f\nMedian      %.0f\n95$^{th}$ Pctl.   %.0f'%(errs[n].quantile(0.05)*xmult, errs[n].median()*xmult,  errs[n].quantile(0.95)*xmult), transform=ax.transAxes, bbox=dict(facecolor='white', alpha=0.5))
# =============================================================================
            #return
            matplotlib.pylab.tight_layout()
            return f


        

        
if __name__ == '__main__':
    
    run_year = 2017
    co2_dol_per_ton = 0 #specify co2 price in $/ton that corresponds to an existing dispatch solution produced by simple_dispatch.py. This is used for naming conventions to access the correct dispatch solution .csv file
    input_file_directory = os.getcwd()
    output_file_directory = os.getcwd()
    for nr in ['FRCC']:
        nerc_region = nr
        #dispatch solution
        dispatch_solution = pandas.read_csv(input_file_directory + '\\dispatch_output_weekly_%s_%s_%sco2.csv'%(nerc_region, str(run_year), str(co2_dol_per_ton)))
        #historical CEMS dispatch data
        dispatch_CEMS = pandas.read_csv(input_file_directory + '\\%s_%s_hourly_demand_and_fuelmix.csv'%(str(run_year), nerc_region),index_col='Unnamed: 0')
        dispatch_CEMS.datetime = pandas.to_datetime(dispatch_CEMS.datetime) #put the datetime column into the correct type
        dispatch_CEMS.gen_cost_marg = scipy.minimum(150, dispatch_CEMS.gen_cost_marg) #remove prices larger than 150 $/MWh. Consider these outliers. For example, with 2014 TRE, R-squared between the actual data and simulated data is 0.16 while R-squared after capping the historical data at 150 $/MWh is 0.53. 150 $/MWh is something above the 99th percentile in most cases, so we are now saying that the simulation has a fit of 0.53 for 99+% of the historical data.
        gm_cems = generateMefs(dispatch_CEMS) 
        dispatch_CEMS = gm_cems.df.copy(deep=True) #the script now calcultes hourly MEFs, which we want for some plotting below
        #data below from CEDM at https://cedm.shinyapps.io/MarginalFactors/
        if run_year == 2014:
            if nerc_region == 'TRE':
                deciles=[21300, 24073, 26245, 28132, 30081, 32297, 35414, 39355, 45364, 60650] 
                mefs_cedm_co2 = [691, 623, 589, 562, 550, 529, 515, 485, 479, 489]
                mefs_cedm_so2 = [1.156, 0.951, 0.835, 0.745, 0.731, 0.586, 0.445, 0.324, 0.240, 0.003]
                mefs_cedm_nox = [0.292, 0.255, 0.238, 0.221, 0.218, 0.191, 0.198, 0.192, 0.232, 0.379]
            if nerc_region == 'FRCC':
                deciles=[13247,15136,16813,18164,19421,20871,22490,24845,28172,36011] 
                mefs_cedm_co2 = [549,535,517,511,509,497,482,471,467,457]
                mefs_cedm_so2 = [0.47,0.50,0.48,0.50,0.51,0.50,0.48,0.48,0.39,0.39]
                mefs_cedm_nox = [0.22,0.25,0.24,0.26,0.28,0.30,0.33,0.38,0.46,0.61]
            if nerc_region == 'MRO':
                deciles=[12466, 13996, 15339, 16650, 17732, 18820, 20006, 21099, 22820, 28718] 
                mefs_cedm_co2 = [866, 835, 834, 823, 799, 792, 758, 727, 676, 643]
                mefs_cedm_so2 = [1.51, 1.29, 1.26, 1.33, 1.27, 1.25, 1.20, 1.06, 0.83, 0.60]
                mefs_cedm_nox = [0.87, 0.78, 0.80, 0.80, 0.76, 0.79, 0.74, 0.70, 0.64, 0.64]
        if run_year == 2017:
            if nerc_region == 'TRE':
                deciles=[18599,22053,24559,26787,29133,31471,34151,38720,46113,61452] 
                mefs_cedm_co2 = [687,639,626,608,604,587,578,556,538,518]
                mefs_cedm_so2 = [0.85,0.85,0.83,0.78,0.81,0.77,0.72,0.64,0.50,0.18]
                mefs_cedm_nox = [0.19,0.21,0.22,0.24,0.24,0.23,0.25,0.25,0.26,0.45]
            if nerc_region == 'FRCC':
                deciles=[13222,15447,17089,18459,19779,21273,23064,25557,28568,36771] 
                mefs_cedm_co2 = [504, 497, 493, 501, 504, 489, 479, 469, 457, 451]
                mefs_cedm_so2 = [0.21,0.22,0.20,0.19,0.21,0.20,0.18,0.18,0.17,0.35]
                mefs_cedm_nox = [0.18,0.17,0.18,0.18,0.21,0.21,0.23,0.27,0.34,0.51]
            if nerc_region == 'MRO':
                deciles=[10607,12321,13770,15001,16211,17429,18683,20173,22157,31205] 
                mefs_cedm_co2 = [852,834,832,822,817,813,782,732,694,628]
                mefs_cedm_so2 = [0.87,0.87,0.82,0.80,0.76,0.72,0.59,0.46,0.36,0.22]
                mefs_cedm_nox = [0.62,0.60,0.60,0.60,0.60,0.59,0.55,0.48,0.47,0.51]
            if nerc_region == 'WECC':
                deciles=[25174,29168,33728,37481,40859,43801,46845,50705,56957,83349] 
                mefs_cedm_co2 = [590,579,564,568,558,553,551,539,538,533]
                mefs_cedm_so2 = [0.21,0.21,0.19,0.19,0.18,0.17,0.16,0.14,0.13,0.18]
                mefs_cedm_nox = [0.44,0.40,0.36,0.39,0.37,0.38,0.37,0.34,0.32,0.26]
        #create the plotDispatch object
        pd = plotDispatch(dispatch_solution, dispatch_CEMS, deciles, mefs_cedm_co2, mefs_cedm_so2, mefs_cedm_nox)
        pd.rolling_calculations(p=500)
        #fuel betas plot
        fuel_betas_plot = pd.plot_x_demand()
        #save
        fuel_betas_plot.savefig(output_file_directory + '\\fDecilesFuelBeta%s_%s.png'%(nerc_region, str(run_year)), dpi=500, bbox_inches='tight')
        #rolling averages plots
        ra_start = '2017-01-01 00:00:00'
        ra_end = '2018-01-01 00:00:00'
        prices_year = pd.plot_hist_vs_simulated(x_property='demand_smooth', y_property='gen_cost_marg', demand_smooth_step=500, start_date=ra_start, end_date=ra_end)        
        co2_year = pd.plot_hist_vs_simulated(x_property='demand_smooth', y_property='co2_tot', demand_smooth_step=500, start_date=ra_start, end_date=ra_end)
        so2_year = pd.plot_hist_vs_simulated(x_property='demand_smooth', y_property='so2_tot', demand_smooth_step=500, start_date=ra_start, end_date=ra_end)
        nox_year = pd.plot_hist_vs_simulated(x_property='demand_smooth', y_property='nox_tot', demand_smooth_step=500, start_date=ra_start, end_date=ra_end)
        #save
        prices_year.savefig(output_file_directory + '\\fPricesYear%s_%s.png'%(nerc_region, str(run_year)), dpi=500)
        co2_year.savefig(output_file_directory + '\\fCo2Year%s_%s.png'%(nerc_region, str(run_year)), dpi=500, bbox_inches='tight')
        so2_year.savefig(output_file_directory + '\\fSo2Year%s_%s.png'%(nerc_region, str(run_year)), dpi=500, bbox_inches='tight')
        nox_year.savefig(output_file_directory + '\\fNoxYear%s_%s.png'%(nerc_region, str(run_year)), dpi=500, bbox_inches='tight')
        #hourly
        #ra_start = '2014-08-01 00:00:00'
        #ra_end = '2014-08-08 00:00:00'
        #prices_year = pd.plot_hist_vs_simulated(x_property='hour', y_property='gen_cost_marg', demand_smooth_step=500, start_date=ra_start, end_date=ra_end)        
        #co2_year = pd.plot_hist_vs_simulated(x_property='hour', y_property='co2_marg', demand_smooth_step=500, start_date=ra_start, end_date=ra_end)
        #so2_year = pd.plot_hist_vs_simulated(x_property='hour', y_property='so2_marg', demand_smooth_step=500, start_date=ra_start, end_date=ra_end)
        #nox_year = pd.plot_hist_vs_simulated(x_property='hour', y_property='nox_marg', demand_smooth_step=500, start_date=ra_start, end_date=ra_end)
        #density function plots
        #price and emissions
        for s in ['WholeYear']: #can plot subsets of the annual data covering only part of the time series
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
            #save       
            price_df.savefig(output_file_directory + '\\fDf%sPrice%s_%s.png'%(s,nerc_region, str(run_year)), dpi=500, bbox_inches='tight')
            price_df_err.savefig(output_file_directory + '\\fErrorDf%sPrice%s_%s.png'%(s,nerc_region, str(run_year)), dpi=500, bbox_inches='tight')
            co2_df.savefig(output_file_directory + '\\fDf%sCo2%s_%s.png'%(s,nerc_region, str(run_year)), dpi=500, bbox_inches='tight')
            co2_df_err.savefig(output_file_directory + '\\fErrorDf%sCo2%s_%s.png'%(s,nerc_region, str(run_year)), dpi=500, bbox_inches='tight')
            so2_df.savefig(output_file_directory + '\\fDf%sSo2%s_%s.png'%(s,nerc_region, str(run_year)), dpi=500, bbox_inches='tight')
            so2_df_err.savefig(output_file_directory + '\\fErrorDf%sSo2%s_%s.png'%(s,nerc_region, str(run_year)), dpi=500, bbox_inches='tight')
            nox_df.savefig(output_file_directory + '\\fDf%sNox%s_%s.png'%(s,nerc_region, str(run_year)), dpi=500, bbox_inches='tight')
            nox_df_err.savefig(output_file_directory + '\\fErrorDf%sNox%s_%s.png'%(s,nerc_region, str(run_year)), dpi=500, bbox_inches='tight')
            #mefs
            co2_mefs_df = pd.plot_density_function('cumulative', 'data', 'co2_marg', start_date=df_start, end_date=df_end, x_range=[0,1500])
            co2_mefs_df_err = pd.plot_density_function('probability', 'error', 'co2_marg', start_date=df_start, end_date=df_end, bin_no=50, x_range=[0,1500/0.454], error_range=[-1000,1000])
            so2_mefs_df = pd.plot_density_function('cumulative', 'data', 'so2_marg', start_date=df_start, end_date=df_end, x_range=[0,3/0.454])
            so2_mefs_df_err = pd.plot_density_function('probability', 'error', 'so2_marg', start_date=df_start, end_date=df_end, bin_no=50, x_range=[0,3/0.454], error_range=[-4,4])
            nox_mefs_df = pd.plot_density_function('cumulative', 'data', 'nox_marg', start_date=df_start, end_date=df_end, x_range=[0,2/0.454])
            nox_mefs_df_err = pd.plot_density_function('probability', 'error', 'nox_marg', start_date=df_start, end_date=df_end, bin_no=50, x_range=[0,2/0.454], error_range=[-3,3])
            #save
            co2_mefs_df.savefig(output_file_directory + '\\fDf%sMefCo2%s_%s.png'%(s,nerc_region, str(run_year)), dpi=500, bbox_inches='tight')
            co2_mefs_df_err.savefig(output_file_directory + '\\fErrorDf%sMefCo2%s_%s.png'%(s,nerc_region, str(run_year)), dpi=500, bbox_inches='tight')
            so2_mefs_df.savefig(output_file_directory + '\\fDf%sMefSo2%s_%s.png'%(s,nerc_region, str(run_year)), dpi=500, bbox_inches='tight')
            so2_mefs_df_err.savefig(output_file_directory + '\\fErrorDf%sMefSo2%s_%s.png'%(s,nerc_region, str(run_year)), dpi=500, bbox_inches='tight')
            nox_mefs_df.savefig(output_file_directory + '\\fDf%sMefNox%s_%s.png'%(s,nerc_region, str(run_year)), dpi=500, bbox_inches='tight')
            nox_mefs_df_err.savefig(output_file_directory + '\\fErrorDf%sMefNox%s_%s.png'%(s,nerc_region, str(run_year)), dpi=500, bbox_inches='tight')
            
# =============================================================================
#             #all four regions together
#             pdmult = plotDispatchMultiple([pd_MRO, pd_TRE, pd_FRCC, pd_WECC])
#             price_df_mult = pdmult.plot_density_function('cumulative', 'data', 'gen_cost_marg', start_date=df_start, end_date=df_end, x_range=[0,60])
#             price_df_err_mult = pdmult.plot_density_function('probability', 'error', 'gen_cost_marg', start_date=df_start, end_date=df_end, bin_no=50, error_range=[-1,1])
#             co2_df_mult = pdmult.plot_density_function('cumulative', 'data', 'co2_tot', start_date=df_start, end_date=df_end)
#             co2_df_err_mult = pdmult.plot_density_function('probability', 'error', 'co2_tot', start_date=df_start, end_date=df_end, bin_no=50, error_range=[-1,1])
#             so2_df_mult = pdmult.plot_density_function('cumulative', 'data', 'so2_tot', start_date=df_start, end_date=df_end, bin_no=50, x_range=[-scipy.inf, scipy.inf])
#             so2_df_err_mult = pdmult.plot_density_function('probability', 'error', 'so2_tot', start_date=df_start, end_date=df_end, bin_no=50, error_range=[-1,1])
#             nox_df_mult = pdmult.plot_density_function('cumulative', 'data', 'nox_tot', start_date=df_start, end_date=df_end)
#             nox_df_err_mult = pdmult.plot_density_function('probability', 'error', 'nox_tot', start_date=df_start, end_date=df_end, bin_no=50, error_range=[-1,1]) 
#             #mefs
#             co2_mefs_df_mult = pdmult.plot_density_function('cumulative', 'data', 'co2_marg', start_date=df_start, end_date=df_end, x_range=[0,1100])
#             co2_mefs_df_err_mult = pdmult.plot_density_function('probability', 'error', 'co2_marg', start_date=df_start, end_date=df_end, bin_no=50, x_range=[0,2000], error_range=[-1100,1100])
#             so2_mefs_df_mult = pdmult.plot_density_function('cumulative', 'data', 'so2_marg', start_date=df_start, end_date=df_end, x_range=[0,4])
#             so2_mefs_df_err_mult = pdmult.plot_density_function('probability', 'error', 'so2_marg', start_date=df_start, end_date=df_end, bin_no=50, x_range=[0,6], error_range=[-4,4])
#             nox_mefs_df_mult = pdmult.plot_density_function('cumulative', 'data', 'nox_marg', start_date=df_start, end_date=df_end, x_range=[0,2])
#             nox_mefs_df_err_mult = pdmult.plot_density_function('probability', 'error', 'nox_marg', start_date=df_start, end_date=df_end, bin_no=50, x_range=[0,3], error_range=[-2,2])
#             #save
#             price_df_mult.savefig(output_file_directory + '\\fDfPriceMult_%s.png'%(str(run_year)), dpi=500, bbox_inches='tight')
#             price_df_err_mult.savefig(output_file_directory + '\\fErrorDfPriceMult_%s.png'%(str(run_year)), dpi=500, bbox_inches='tight')
#             co2_df_mult.savefig(output_file_directory + '\\fDfCo2Mult_%s.png'%(str(run_year)), dpi=500, bbox_inches='tight')
#             co2_df_err_mult.savefig(output_file_directory + '\\fErrorDfCo2Mult_%s.png'%(str(run_year)), dpi=500, bbox_inches='tight')
#             so2_df_mult.savefig(output_file_directory + '\\fDfSo2Mult_%s.png'%(str(run_year)), dpi=500, bbox_inches='tight')
#             so2_df_err_mult.savefig(output_file_directory + '\\fErrorDfSo2Mult_%s.png'%(str(run_year)), dpi=500, bbox_inches='tight')
#             nox_df_mult.savefig(output_file_directory + '\\fDfNoxMult_%s.png'%(str(run_year)), dpi=500, bbox_inches='tight')
#             nox_df_err_mult.savefig(output_file_directory + '\\fErrorDfNoxMult_%s.png'%(str(run_year)), dpi=500, bbox_inches='tight')
#             co2_mefs_df_mult.savefig(output_file_directory + '\\fDfMefCo2Mult_%s.png'%(str(run_year)), dpi=500, bbox_inches='tight')
#             co2_mefs_df_err_mult.savefig(output_file_directory + '\\fErrorDfMefCo2Mult_%s.png'%(str(run_year)), dpi=500, bbox_inches='tight')
#             so2_mefs_df_mult.savefig(output_file_directory + '\\fDfMefSo2Mult_%s.png'%(str(run_year)), dpi=500, bbox_inches='tight')
#             so2_mefs_df_err_mult.savefig(output_file_directory + '\\fErrorDfMefSo2Mult_%s.png'%(str(run_year)), dpi=500, bbox_inches='tight')
#             nox_mefs_df_mult.savefig(output_file_directory + '\\fDfMefNoxMult_%s.png'%(str(run_year)), dpi=500, bbox_inches='tight')
#             nox_mefs_df_err_mult.savefig(output_file_directory + '\\fErrorDfMefNoxMult_%s.png'%(str(run_year)), dpi=500, bbox_inches='tight')
# =============================================================================


            
