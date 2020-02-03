% Main task of this script is to conduct uncertainty analysis based on streamflow 
% as observational dataset. It generates different plots and summary statistics 
% such as mean, median, variance, and skewness of the behavioural parameter sets.
% This script was adopted from and calls to other scripts of the SAFE Toolbox (Pianosi et al, 2015)

% For details on the conceptual background, implementation, and sample outputs from this algorithm, 
% the reader is referred to the following paper: 

% Teweldebrhan, A. T., Burkhart, J. F., and Schuler, T. V.: Parameter uncertainty analysis for an 
% operational hydrological model using residual-based and limits of acceptability approaches, 
% Hydrology and Earth System Sciences, 22, 5021-5039, 2018.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all; close all
format long g
% general info
mc_calib_param_file = 'shyft/24-Mc_calib_param_new_cv_100k.nc';
LH_measure=2;% 1=NSE alone; 2=combined
threshold =0.65;% threshold value
threshold_min_NSE_for_plot = -10; %-10;% 0.0;
analysis_period_calib = 365;

% exclude the first n warm_up_days from NSE estimate
warm_up_days = 30; 
warm_up_days_2011_addnal= 30; 
%-----------------------------------------------------------------------------------------------

% SC related
LH_measure_SCA=2;% 1=CSI alone; 2=RMSE alone
mean_LH_calc_opn = 1; % 1= using weighted LH, 2=using non weighted LH
threshold_SCA = 0.17 ; % threshold value RMSE
ccf_t = 0.2; %cloud cover fraction threshold

%-----------------------------------------------------------------------------------------------
% Set up project and run simulations
my_dir = pwd ; 
cd(my_dir)

addpath('util');
addpath('visualization')
addpath('GLUE')
addpath('shyft')
addpath('x_extra')
addpath('\\lagringshotell\geofag\projects\hycamp\team\aynomtt\data\MC\snow_cv_included\')
addpath('\\lagringshotell\geofag\projects\hycamp\team\aynomtt\data\MC\nea_100k\')

N = 100000 ; % sample size parameterisations
M  = 21; % number of uncertain parameters 

% variable_calib_param_index =[1,2,3,5,6,9,10];
% S_No =   [1         2       3           4       5       6       7       8       9       10      11      12      13      14      15      16      17      18      19      20      21   ];
S_No =     {'1'       '2'       '3'       '4'      '5'    '6'     '7'     '8'     '9'      '10'   '11'    '12'    '13'    '14'    '15'    '16'    '17'    '18'    '19'    '20'    '21'   };
% X_Labels = {'c1'     'c2'       'c3'    'ae scale factor'      'tx'    'wind scale'     'max water'     'wind const'     'fast albedo decay rate'      'slow albedo decay rate'   'surface magnitude'    'max albedo'    'min albedo'    'snowfall reset depth'    'snow cv'    'snow cv forest factor'    'snow cv altitude factor'    'glacier albedo'    'scale factor'    'albedo'    'alpha'   };
X_Labels = {'c1'     'c2'       'c3'    'ae scale factor'      'tx'    'wind scale'     'max water'     'wind const'     'fast ADR'      'slow ADR'   'surface magnitude'    'max albedo'    'min albedo'    'snowfall reset depth'    'snow cv'    'snow cv forest factor'    'snow cv altitude factor'    'glacier albedo'    'scale factor'    'albedo'    'alpha'   };

% Parameter ranges:
xmin =   [-5.0,     0.0,    -0.15,      1.5,    -3.0,   1.0,    0.1,    0.1,    1.0,    20.0,   30.0,   0.9,    0.6,    5.0,    0.06,    0.0,    0.0,    0.4,    1.0,    0.2,    1.26 ] ; % minimum values
xmax =   [1.0,      1.2,    -0.05 ,     1.5,    2.0,    6.0,    0.1,    0.1,    15.0,   40.0,   30.0,   0.9,    0.6,    5.0,    0.85,    0.0,    0.0,    0.4,    1.0,    0.2,    1.26 ] ; % maximum values

%------------------------------------------------------------------------------------------------

% GLUE analysis part

variable_calib_param_index =[1,2,3,5,6,9,10,15];% parameter indices whose value to be identified
S_No =S_No(variable_calib_param_index);
X_Labels =X_Labels(variable_calib_param_index);
xmin =xmin(variable_calib_param_index);
xmax =xmax(variable_calib_param_index);

% read calibrn parameter values
mc_calib_param = ncread(mc_calib_param_file,'calib_MC');

X = mc_calib_param(:,variable_calib_param_index);

% get the monte carlo calib param value based ensembel simulations from shyft
time_obs_Q = ncread('shyft/discharge_hourly_f.nc','time');
obs_Q_data = ncread('shyft/discharge_hourly_f.nc','discharge'); 

%-----------------------------------------------------------------------------------------------
% check and fix suspected Obs_Q values
suspecteddata = [datetime('27-Mar-2011 01:00:00'),datetime('25-Mar-2012 01:00:00'),datetime('31-Mar-2013 01:00:00')]
obs_Q_data =obs_Q_data(1,:);
for i=1:length(suspecteddata)
    unix_suspecteddata = posixtime(suspecteddata(i));
    susp_index = find(time_obs_Q==unix_suspecteddata);
    susp_Q = obs_Q_data(susp_index);
    susp_Q_bef = obs_Q_data(susp_index-1);
    susp_Q_aft = obs_Q_data(susp_index+1);
    obs_Q_data(susp_index) = (susp_Q_bef+susp_Q_aft)/2;% replace by average value of the preceeding and following values
end
%-----------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------
% do analysis based on whole observed period
%-----------------------------------------------------------------------
analysis_period = analysis_period_calib; 

% convert hourly observed data to daily data
obs_Q_data = mean(reshape(obs_Q_data(:,:),24,[])); % get daily (24 hrs) mean value of discharge 
time_obs_Q = min(reshape(time_obs_Q,24,[]));

% extract obs Q and remove warm_up period values from each year
analysis_period = 365;
for yr = 1:4
    start_yr = yr-1;
    start_date = datetime(strcat('01-sep-201',num2str(start_yr)));
    unix_start_date = posixtime(start_date);
    index = find(time_obs_Q==unix_start_date);
    obs_Q =obs_Q_data(index:index+analysis_period-1);
    obs_Q_i = obs_Q';
    flow_i = obs_Q_i(warm_up_days:end);
    flow(yr,:)=flow_i;  
end
% concatenate obs Q
flow_r = cat(2,flow(1,warm_up_days_2011_addnal:end),flow(2,:),flow(3,:),flow(4,:));

% read sim Q nc
MC_sim_Q_2011_i = ncread(strcat('MC_sim_Q_era_','2011_','99999.nc'),'sim_Q'); % from hycamp\....\MC
MC_sim_Q_2012_i = ncread(strcat('MC_sim_Q_era_','2012_','99999.nc'),'sim_Q'); % from hycamp\....\MC
MC_sim_Q_2013_i = ncread(strcat('MC_sim_Q_era_','2013_','99999.nc'),'sim_Q'); % from hycamp\....\MC
MC_sim_Q_2014_i = ncread(strcat('MC_sim_Q_era_','2014_','99999.nc'),'sim_Q'); % from hycamp\....\MC

% remove warmup period
MC_sim_Q_2011 =MC_sim_Q_2011_i(warm_up_days+warm_up_days_2011_addnal-1:end,:);
MC_sim_Q_2012 =MC_sim_Q_2012_i(warm_up_days:end,:);
MC_sim_Q_2013 =MC_sim_Q_2013_i(warm_up_days:end,:);
MC_sim_Q_2014 =MC_sim_Q_2014_i(warm_up_days:end,:);

% concatenate sim_Q
MC_sim_Q = cat(1,MC_sim_Q_2011,MC_sim_Q_2012,MC_sim_Q_2013,MC_sim_Q_2014);

Qsim = MC_sim_Q';

size(obs_Q);
size(Qsim);

time_obs_Q = time_obs_Q(index:index+analysis_period-1);
norm_time_obs_Q =  datetime(time_obs_Q,'ConvertFrom','posixtime');

[nse,nse_log,sim_bias] = goodnessOfFitFuncn_shyft_May2017_for_3yrs_uncert(flow_r,Qsim);

disp(['min_and_max_nse: ',num2str(min(nse)),num2str(max(nse))]);
disp(['min_and_max_nse_log: ',num2str(min(nse_log)),num2str(max(nse_log))]);

%-----------------------------------------------------------------------
% Likelihood measure
%-------------------------------------------------------------------------

if LH_measure==1% NSE alone
    min_NSE_for_plot_idx = nse>=threshold_min_NSE_for_plot;
    Y = nse(min_NSE_for_plot_idx); %av_wted_NSE ;%nse;
    
elseif LH_measure==2% combined NSE and LnNSE    
    nse_wt = 0.54;
    nse_log_wt = 0.46;
    av_wted_NSE = nse_wt*nse+nse_log_wt*nse_log;
    min_NSE_for_plot_idx = av_wted_NSE>= threshold_min_NSE_for_plot; 
    Y = av_wted_NSE(min_NSE_for_plot_idx); 
end
%%--------------------------------------------------------------------------
Qsim=Qsim(min_NSE_for_plot_idx,:);
X = X(min_NSE_for_plot_idx,:);

% Check size of 'Y':
size(Y)
max_NSE = max(Y)%for check

[ idx, Llim, Ulim, median_predcn ] = shyft_GLUE_max(Y,threshold,Qsim,0.05) ;
[nse_median,nse_log_median,sim_bias_median] = goodnessOfFitFuncn_shyft_May2017_for_3yrs_uncert(flow_r,median_predcn');
disp(['calibrn: nse and nse_log for median crisp forecast', '  ' ,num2str(nse_median),'  ' , num2str(nse_log_median) ]);

% statistics
% calculate the containing ratio(CR)
I_Qobs=0;
for t=1:length(obs_Q)
    
    obs_Q_t = obs_Q(t);
    Llim_t = Llim(t);
    Ulim_t = Ulim(t);
    if (Llim_t<=obs_Q_t)&&(obs_Q_t<=Ulim_t)
        I_Qobs=I_Qobs+1;
    end    
end
CR_1=I_Qobs/length(obs_Q);
disp(['CR_1 = ','  ', num2str(CR_1)]);

% % Display number of behavioural simulations
disp(['Number of behavioural simulations  =  ', num2str(sum(idx))]);

%%
% -------------------------------------------------------------------------
% Statistical summary of posterior distribution (only behavioral models) for model parameters
% -------------------------------------------------------------------------

behav_MC_values = X(idx,:);
min_MC_val = round(min(X(idx,:)),2);
max_MC_val = round(max(X(idx,:)),2);
mean_MC_val = round(mean(X(idx,:)),2);
median_MC_val = round(median(X(idx,:)),2);
var_MC_val = round(var(X(idx,:)),2);
skew_MC_val = round(skewness(X(idx,:)),2);

summary_stat = [min_MC_val;max_MC_val;mean_MC_val;median_MC_val;var_MC_val;skew_MC_val];

save behav_MC_values_basedon_Q_only_4yrs_Dec2017_100k.mat behav_MC_values;

%%
% ----------------------------------------------------------------------------------
% plots

% plot NSE vs parameter and Highlight 'behavioural parameterizations' in different colour: 
shyft_scatter_plots(X,Y,[],'Likelihood',X_Labels,idx)
% disp('Press any key to continue');pause

% Parallel coordinate plots:
figure
shyft_parcoor_May_2017(X,Y,X_Labels,[],idx);
xlabel('model parameters'); ylabel('parameter range')
disp('Press any key to continue');pause

% plot 2D correlations
C = ones(sum(idx),1);% this sets the symbols to a single colour
shyft_scatter_plots_interaction_diagonalLabels_2(X(idx,:),C,Y(idx,:),12,X_Labels)

disp('Press any key to continue');pause

% plot posterior parameter distribution (histogram)
shyft_posterior_hist(X(idx,:),2,4,X_Labels,35)

%##########################################################################
% snow related
MC_scf_stat_2011 = ncread(strcat('MC_scf_stat_era_nea_100k_','2011','.nc'),'scf_stat'); 
MC_scf_stat_2012 = ncread(strcat('MC_scf_stat_era_nea_100k_','2012','.nc'),'scf_stat'); 
MC_scf_stat_2013 = ncread(strcat('MC_scf_stat_era_nea_100k_','2013','.nc'),'scf_stat'); 
MC_scf_stat_2014 = ncread(strcat('MC_scf_stat_era_nea_100k_','2014','.nc'),'scf_stat'); 

% get snow likelihood values for each model iteration and year
[mean_RMSE_2011,mean_CSI_2011] = getObsvnDayWeightedMeanRMSE_Mar2017(MC_scf_stat_2011,ccf_t,mean_LH_calc_opn);
[mean_RMSE_2012,mean_CSI_2012] = getObsvnDayWeightedMeanRMSE_Mar2017(MC_scf_stat_2012,ccf_t,mean_LH_calc_opn);
[mean_RMSE_2013,mean_CSI_2013] = getObsvnDayWeightedMeanRMSE_Mar2017(MC_scf_stat_2013,ccf_t,mean_LH_calc_opn);
[mean_RMSE_2014,mean_CSI_2014] = getObsvnDayWeightedMeanRMSE_Mar2017(MC_scf_stat_2014,ccf_t,mean_LH_calc_opn);

mean_RMSE = (mean_RMSE_2011+mean_RMSE_2012+mean_RMSE_2013+mean_RMSE_2014)/4;

Y_SCA = mean_RMSE';
Y_SCA = Y_SCA(min_NSE_for_plot_idx);% remove models with extremely poor performance in Q simuln 
idx_SCA = Y_SCA<=threshold_SCA;

behav_MC_values = X(idx_SCA,:);
save behav_MC_values_basedon_SC_only_4yrs_Dec2017_100k.mat behav_MC_values;

% SC Parallel coordinate plots:
figure
X_Labels_snow = X_Labels(4:end);
X_snow = X(:,4:end);
shyft_parcoor_snow_Feb_2017(X_snow,Y_SCA,LH_measure_SCA,X_Labels_snow,[],idx_SCA);
xlabel('model parameters'); ylabel('parameter range')

%##########################################################################

% export the figure with same format as it appears on screen
addpath('x_extra\altmany-export_fig-5be2ca4')
set(gcf, 'Color', 'w');% set figure background to white
export_fig(gcf,'param_interacn_scatter_plots.fig', '-r300');
export_fig 'param_interacn_scatter_plots' -tiff -painters -r300;

