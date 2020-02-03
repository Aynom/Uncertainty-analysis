
'''
This python algorithm was implemented to run a Monte Carlo simulation using the PT_GS_K
model and the sampled parameter sets. The simulation results, i.e. streamflow and fractional
snow cover area are exported to netcdf files for use in further analysis.

For details on the conceptual background, implementation, and sample outputs from this algorithm, 
the reader is referred to the following paper: 

Teweldebrhan, A. T., Burkhart, J. F., and Schuler, T. V.: Parameter uncertainty analysis for an 
operational hydrological model using residual-based and limits of acceptability approaches, 
Hydrology and Earth System Sciences, 22, 5021-5039, 2018.
'''

# import relevant modules
import os
from os import path

import sys
from matplotlib import pyplot as plt
import numpy as np

from netCDF4  import Dataset

import datetime
from matplotlib.font_manager import FontProperties
fontP = FontProperties()
fontP.set_size('9')

import datetime


from netCDF4 import num2date, date2num
from matplotlib.dates import strpdate2num
# from datetime import datetime, timedelta


def netCDF_writer(MC_sim_Q, MC_sim_scf, MC_sim_swe, sim_dates, sim_year, itern_i):
    
    MC_sim_Q_shape = np.array(MC_sim_Q).shape
    itern = MC_sim_Q_shape[0]
    daily_Q = MC_sim_Q_shape[1]

    ncfile_name = 'MC_sim_Q_era_'+np.str(sim_year)+'_'+np.str(itern_i)+'.nc'
    ncfile = Dataset(ncfile_name, 'w')
    ncfile.createDimension('itern', itern)
    ncfile.createDimension('daily_Q', daily_Q)

    sim_Q = ncfile.createVariable('sim_Q', np.float64, ('itern', 'daily_Q'))
    sim_Q[:] = MC_sim_Q #np.array(MC_sim_Q).transpose()

    ncfile.close()

    #----------------------------------------------------------------------------
    MC_sim_scf_shape = np.array(MC_sim_scf).shape
    ncfile_name = 'MC_sim_scf_era_'+np.str(sim_year)+'_'+np.str(itern_i)+'.nc'
    ncfile = Dataset(ncfile_name, 'w')
    ncfile.createDimension('itern_no', MC_sim_scf_shape[0])
    ncfile.createDimension('sim_day', MC_sim_scf_shape[1])
    ncfile.createDimension('grid_cell_id', MC_sim_scf_shape[2])
    sim_scf = ncfile.createVariable('sim_scf', np.float64, ('itern_no', 'sim_day','grid_cell_id'))
    sim_scf[:] = MC_sim_scf
    ncfile.close()

    # ----------------------------------------------------------------------------
    MC_sim_swe_shape = np.array(MC_sim_swe).shape
    ncfile_name = 'MC_sim_swe_era_'+np.str(sim_year)+'_'+np.str(itern_i)+'.nc'
    ncfile = Dataset(ncfile_name, 'w')
    ncfile.createDimension('itern_no', MC_sim_swe_shape[0])
    ncfile.createDimension('sim_day', MC_sim_swe_shape[1])
    ncfile.createDimension('grid_cell_id', MC_sim_swe_shape[2])
    sim_swe = ncfile.createVariable('sim_swe', np.float64, ('itern_no', 'sim_day', 'grid_cell_id'))
    sim_swe[:] = MC_sim_swe
    ncfile.close()


#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
if __name__=="__main__":


    sc_sim_year = 2014

    #--------------------------------------------------------------------------------------

    sys.path.insert(0,'C:\projects4\shyft_UA\shyft')
    # sys.path.insert(0, '/home/aynomtt/hycamp/software/shyft_workspace/shyft')
    # sys.path.insert(0, r'\uio\lagringshotell\hycamp\software\shyft_workspace\shyft')

    started_time = datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S')

    from shyft.repository.default_state_repository import DefaultStateRepository
    from shyft.orchestration.configuration import yaml_configs
    from shyft.orchestration.simulators.config_simulator import ConfigSimulator
    from shyft import api

    config_file_path = os.path.abspath("./nea-config/initnea/gs_yaml_era/nea_simulation_2014.yaml")
    cfg = yaml_configs.YAMLSimConfig(config_file_path, "neanidelva")

    simulator = ConfigSimulator(cfg)
    n_cells = simulator.region_model.size()
    state_repos = DefaultStateRepository(simulator.region_model.__class__, n_cells)

    simulator.run(cfg.time_axis, state_repos.get_state(0))

    discharge_ts = []
    catchment_id_map = simulator.region_model.catchment_id_map 
    for catch in range(len(catchment_id_map)):
        disch = simulator.region_model.statistics.discharge([catch])
        ts=[disch.value(t) for t in range(disch.size())]
        discharge_ts.append(ts)
    sim_times = [disch.time(i) for i in range(len(ts))] 
    sim_dates = [datetime.datetime.utcfromtimestamp(t) for t in sim_times] 
    sim_Q_org = sum([arr for arr in np.array(discharge_ts)[:]])

    MC_sim_Q=[];    MC_sim_scf=[];     MC_sim_swe=[]

    # read MC data from netcdf file
    with Dataset('Mc_calib_param_new_cv_100k.nc') as dset:
        calib_val = dset.variables['calib_MC'][:]

    # read MC data for current iteration
    mc_shape = calib_val.shape
    no_param = mc_shape[0]
    itern_len = mc_shape[1]

    print(no_param, itern_len)

    region_parameter = simulator.region_model.get_region_parameter()
    cells = simulator.region_model.get_cells()

    # get days with MODIS SC over 50%
    SC_obs_date_file = 'MODIS_SC_dates_over_50_perc_data_'+np.str(sc_sim_year)+'.txt'

    SC_obs_date = []
    data = np.genfromtxt(SC_obs_date_file, delimiter=' ', skip_header=0, dtype=None)
    for i in range(len(data)):
        line = data[i]
        b0 = line[0]
        s0 = b0.decode(encoding='UTF-8')
        sc_date = s0[0:10]
        SC_obs_date.append(sc_date)

    # sim_scf_database = [[[-9999.0 for cells in range(812)] for t in range(len(SC_obs_date))] for i in range(itern_len)]
    # sim_swe_database = [[[-9999.0 for cells in range(812)] for t in range(len(SC_obs_date))] for i in range(itern_len)]


    for i in range (itern_len):
        # i=i+63
        print('iteration: ',i)

        calib_param_i = calib_val[:, i]

        c1 =region_parameter.kirchner.c1                                        =calib_param_i[0]
        c2 = region_parameter.kirchner.c2                                       =calib_param_i[1]
        c3 =region_parameter.kirchner.c3                                        =calib_param_i[2]

        ae_scale_factor = region_parameter.ae.ae_scale_factor                   =calib_param_i[3]
        tx = region_parameter.gs.tx                                             =calib_param_i[4]
        wind_scale = region_parameter.gs.wind_scale                             =calib_param_i[5]
        max_water = region_parameter.gs.max_water                               =calib_param_i[6]
        wind_const = region_parameter.gs.wind_const                             =calib_param_i[7]
        fast_albedo_decay_rate =region_parameter.gs.fast_albedo_decay_rate      =calib_param_i[8]
        slow_albedo_decay_rate = region_parameter.gs.slow_albedo_decay_rate     =calib_param_i[9]
        surface_magnitude = region_parameter.gs.surface_magnitude               =calib_param_i[10]
        max_albedo = region_parameter.gs.max_albedo                             =calib_param_i[11]
        min_albedo = region_parameter.gs.min_albedo                             =calib_param_i[12]
        snowfall_reset_depth = region_parameter.gs.snowfall_reset_depth         =calib_param_i[13]
        snow_cv = region_parameter.gs.snow_cv                                   =calib_param_i[14]
        snow_cv_forest_factor = region_parameter.gs.snow_cv_forest_factor       =calib_param_i[15]
        snow_cv_altitude_factor = region_parameter.gs.snow_cv_altitude_factor   =calib_param_i[16]
        glacier_albedo = region_parameter.gs.glacier_albedo                     =calib_param_i[17]

        scale_factor = region_parameter.p_corr.scale_factor                     =calib_param_i[18]
        albedo = region_parameter.pt.albedo                                     =calib_param_i[19]
        alpha = region_parameter.pt.alpha                                       =calib_param_i[20]

        simulator.run(cfg.time_axis, state_repos.get_state(0))

        discharge_ts_2 = []
        for catch_2 in range(len(catchment_id_map)):
            disch_2 = simulator.region_model.statistics.discharge([catch_2])
            ts_2 = [disch_2.value(t) for t in range(disch_2.size())]
            discharge_ts_2.append(ts_2)

        sim_Q = sum([arr for arr in np.array(discharge_ts_2)[:]])
        MC_sim_Q.append(sim_Q)

        MC_sim_scf_itn_i = [];MC_sim_swe_itn_i = []

        for sc_date in SC_obs_date:
            sc_date = datetime.datetime.strptime(sc_date, '%d.%m.%Y')
            # print('sc_date: ',sc_date)
            try:
                min_sim_dates = np.min(sim_dates)
                max_sim_dates = np.max(sim_dates)
                idx = sim_dates.index(sc_date)  
                # sim_scf = [cell.rc.snow_sca.v[:] for cell in cells]
                sim_scf = [cell.rc.snow_sca.v[idx] for cell in cells]
                MC_sim_scf_itn_i.append(sim_scf)
                # sim_swe = [cell.rc.snow_swe.v[:] for cell in cells]
                sim_swe = [cell.rc.snow_swe.v[idx] for cell in cells]
                MC_sim_swe_itn_i.append(sim_swe)

            except:
                print("Error: Date index out of range")
                # idx = 0

        MC_sim_scf.append(MC_sim_scf_itn_i) # (itern, date, cells)
        MC_sim_swe.append(MC_sim_swe_itn_i)

        if i==(itern_len/4)-1 or i==(itern_len/2)-1 or i==(3*itern_len/4)-1 or i==itern_len-1:
            # check
            print(sim_Q.shape)
            print(np.array(MC_sim_Q).shape)
            print(np.array(MC_sim_scf).shape)

            netCDF_writer(MC_sim_Q, MC_sim_scf, MC_sim_swe, sim_dates,sc_sim_year, i)

            MC_sim_scf = []; MC_sim_swe = []

            print('start time :', started_time)
            print('finishing time:', datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S'))





 

