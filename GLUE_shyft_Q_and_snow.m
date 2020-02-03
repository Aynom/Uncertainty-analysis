% This script implements residual based GLUE uncertainty analysis using a combined likelihood of 
% streamflow and fSCA as observational datasets. 
% Outputs from this script include the cross-validation plots and tabular result of the 
% evaluation metrics for streamflow and fSCA both for the calibration and validation years.
% This script was adopted from and calls to other scripts of the SAFE Toolbox (Pianosi et al, 2015

% For details on the conceptual background, implementation, and sample outputs from this algorithm, 
% the reader is referred to the following paper: 

% Teweldebrhan, A. T., Burkhart, J. F., and Schuler, T. V.: Parameter uncertainty analysis for an 
% operational hydrological model using residual-based and limits of acceptability approaches, 
% Hydrology and Earth System Sciences, 22, 5021-5039, 2018.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all
format long g

% general info on Q simuln
mc_calib_param_file = 'shyft/24-Mc_calib_param_new_cv_100k.nc';
LH_measure_Q=2;% 1=NSE alone; 2=combined
threshold_Q = 0.65 ; % threshold value
lowest_NSE_for_analysis = -10; %exclude the extremely low performing models from the analysis
analysis_period_calib = 365;
analysis_period_1 = 365;
analysis_period_2 = 365;
analysis_period_3 = 365;

% exclude the first n warm_up_days from NSE estimate
org_warm_up_days = 30; 
warm_up_days_2011 = 60; 
%-----------------------------------------------------------------------------------------------
% general info on Snow simuln
SCA_LH_wt = 0.1;
N = 100000 ; % sample size parameterisations
LH_measure_SCA=2;% 1=CSI alone; 2=RMSE alone
mean_LH_calc_opn = 1; % 1= using weighted LH, 2=using non weighted LH
threshold_1_SCA = 0.95; %0.90 ; % threshold value CSI
threshold_2_SCA = 0.17; %0.20 ; %0.15 ; % threshold value RMSE
ccf_t = 0.2;%0.2;%cloud cover fraction threshold

% **check summary_LHandCR at workspace for output

%-------------------------------------------------------------------------------------------------------------------------------------------------

%%%%%Section 1 Set up project and run simulations
my_dir = pwd ; 
cd(my_dir)

addpath('util');
addpath('visualization')
addpath('GLUE')
addpath('shyft')
addpath('x_extra')
addpath('\\lagringshotell\geofag\projects\hycamp\team\aynomtt\data\MC\snow_cv_included\')
addpath('\\lagringshotell\geofag\projects\hycamp\team\aynomtt\data\MC\nea_100k\')

M  = 21; %5 ; % number of uncertain parameters 
% S_No =   [1         2       3           4       5       6       7       8       9       10      11      12      13      14      15      16      17      18      19      20      21   ];
S_No =     {'1'       '2'       '3'       '4'      '5'    '6'     '7'     '8'     '9'      '10'   '11'    '12'    '13'    '14'    '15'    '16'    '17'    '18'    '19'    '20'    '21'   };
X_Labels = {'c1'     'c2'       'c3'    'ae scale factor'      'tx'    'wind scale'     'max water'     'wind const'     'fast ADR'      'slow ADR'   'surface magnitude'    'max albedo'    'min albedo'    'snowfall reset depth'    'snow cv'    'snow cv forest factor'    'snow cv altitude factor'    'glacier albedo'    'scale factor'    'albedo'    'alpha'   };

% Parameter ranges:
xmin =   [-5.0,     0.0,    -0.15,      1.5,    -3.0,   1.0,    0.1,    0.1,    1.0,    20.0,   30.0,   0.9,    0.6,    5.0,    0.06,    0.0,    0.0,    0.4,    1.0,    0.2,    1.26 ] ; % minimum values
xmax =   [1.0,      1.2,    -0.05 ,     1.5,    2.0,    6.0,    0.1,    0.1,    15.0,   40.0,   30.0,   0.9,    0.6,    5.0,    0.85,    0.0,    0.0,    0.4,    1.0,    0.2,    1.26 ] ; % maximum values

variable_calib_param_index =[1,2,3,5,6,9,10,15];% parameters whose value to be identified
S_No =S_No(variable_calib_param_index);
X_Labels =X_Labels(variable_calib_param_index);
xmin =xmin(variable_calib_param_index);
xmax =xmax(variable_calib_param_index);

% read calibrn parameter values
mc_calib_param = ncread(mc_calib_param_file,'calib_MC');
X = mc_calib_param(:,variable_calib_param_index);

% read sc related files    
disp(['reading SCF data for all years'])
MODIS_scf_2011 = ncread(strcat('annual_scf_for_cells_','2011','.nc'),'scf');
% MC_sim_scf_2011_1 = ncread(strcat('MC_sim_scf_era_2011_24999.nc'),'sim_scf'); % from hycamp
% MC_sim_scf_2011_2 = ncread(strcat('MC_sim_scf_era_2011_49999.nc'),'sim_scf'); % from hycamp
% MC_sim_scf_2011_3 = ncread(strcat('MC_sim_scf_era_2011_74999.nc'),'sim_scf'); % from hycamp
% MC_sim_scf_2011_4 = ncread(strcat('MC_sim_scf_era_2011_99999.nc'),'sim_scf'); % from hycamp
% MC_sim_scf_2011=cat(3,MC_sim_scf_2011_1,MC_sim_scf_2011_2,MC_sim_scf_2011_3,MC_sim_scf_2011_4);

load('MC_sim_scf_2011.mat');
MC_scf_stat_2011 = ncread(strcat('MC_scf_stat_era_nea_100k_','2011','.nc'),'scf_stat'); 
%----------------
MODIS_scf_2012 = ncread(strcat('annual_scf_for_cells_','2012','.nc'),'scf');
% MC_sim_scf_2012_1 = ncread(strcat('MC_sim_scf_era_2012_24999.nc'),'sim_scf'); % from hycamp
% MC_sim_scf_2012_2 = ncread(strcat('MC_sim_scf_era_2012_49999.nc'),'sim_scf'); % from hycamp
% MC_sim_scf_2012_3 = ncread(strcat('MC_sim_scf_era_2012_74999.nc'),'sim_scf'); % from hycamp
% MC_sim_scf_2012_4 = ncread(strcat('MC_sim_scf_era_2012_99999.nc'),'sim_scf'); % from hycamp
% MC_sim_scf_2012=cat(3,MC_sim_scf_2012_1,MC_sim_scf_2012_2,MC_sim_scf_2012_3,MC_sim_scf_2012_4);

load('MC_sim_scf_2012.mat');
MC_scf_stat_2012 = ncread(strcat('MC_scf_stat_era_nea_100k_','2012','.nc'),'scf_stat'); 
%----------------
MODIS_scf_2013 = ncread(strcat('annual_scf_for_cells_','2013','.nc'),'scf');
% MC_sim_scf_2013_1 = ncread(strcat('MC_sim_scf_era_2013_24999.nc'),'sim_scf'); % from hycamp
% MC_sim_scf_2013_2 = ncread(strcat('MC_sim_scf_era_2013_49999.nc'),'sim_scf'); % from hycamp
% MC_sim_scf_2013_3 = ncread(strcat('MC_sim_scf_era_2013_74999.nc'),'sim_scf'); % from hycamp
% MC_sim_scf_2013_4 = ncread(strcat('MC_sim_scf_era_2013_99999.nc'),'sim_scf'); % from hycamp
% MC_sim_scf_2013=cat(3,MC_sim_scf_2013_1,MC_sim_scf_2013_2,MC_sim_scf_2013_3,MC_sim_scf_2013_4);

load('MC_sim_scf_2013.mat');
MC_scf_stat_2013 = ncread(strcat('MC_scf_stat_era_nea_100k_','2013','.nc'),'scf_stat'); 

%----------------
MODIS_scf_2014 = ncread(strcat('annual_scf_for_cells_','2014','.nc'),'scf');
% MC_sim_scf_2014_1 = ncread(strcat('MC_sim_scf_era_2014_24999.nc'),'sim_scf'); % from hycamp
% MC_sim_scf_2014_2 = ncread(strcat('MC_sim_scf_era_2014_49999.nc'),'sim_scf'); % from hycamp
% MC_sim_scf_2014_3 = ncread(strcat('MC_sim_scf_era_2014_74999.nc'),'sim_scf'); % from hycamp
% MC_sim_scf_2014_4 = ncread(strcat('MC_sim_scf_era_2014_99999.nc'),'sim_scf'); % from hycamp
% MC_sim_scf_2014=cat(3,MC_sim_scf_2014_1,MC_sim_scf_2014_2,MC_sim_scf_2014_3,MC_sim_scf_2014_4);

load('MC_sim_scf_2014.mat');
MC_scf_stat_2014 = ncread(strcat('MC_scf_stat_era_nea_100k_','2014','.nc'),'scf_stat'); 

%----------------

% clearvars -except MODIS_scf_2011 MC_sim_scf_2011 MC_scf_stat_2011 MODIS_scf_2012 MC_sim_scf_2012 MC_scf_stat_2012 ...
%     MODIS_scf_2013 MC_sim_scf_2013 MC_scf_stat_2013 MODIS_scf_2014 MC_sim_scf_2014 MC_scf_stat_2014;

for m=1:4
    if m==1; start_date_calib = '01-sep-2010'; start_date_1= '01-sep-2011';start_date_2= '01-sep-2012';start_date_3= '01-sep-2013';...
            MODIS_scf_calib = MODIS_scf_2011 ;  MODIS_scf_analysis_1 = MODIS_scf_2012 ; MODIS_scf_analysis_2 = MODIS_scf_2013 ;MODIS_scf_analysis_3 = MODIS_scf_2014 ;...
            MC_sim_scf_calib = MC_sim_scf_2011 ;  MC_sim_scf_analysis_1 = MC_sim_scf_2012 ; MC_sim_scf_analysis_2 = MC_sim_scf_2013 ;MC_sim_scf_analysis_3 = MC_sim_scf_2014 ;...
            MC_scf_stat_calib = MC_scf_stat_2011 ;  MC_scf_stat_analysis_1 = MC_scf_stat_2012 ; MC_scf_stat_analysis_2 =MC_scf_stat_2013 ;MC_scf_stat_analysis_3 = MC_scf_stat_2014 ;
    end  
    if m==2; start_date_calib = '01-sep-2011'; start_date_1= '01-sep-2010';start_date_2= '01-sep-2012';start_date_3= '01-sep-2013'; ...
            MODIS_scf_calib = MODIS_scf_2012 ;  MODIS_scf_analysis_1 = MODIS_scf_2011 ; MODIS_scf_analysis_2 = MODIS_scf_2013 ;MODIS_scf_analysis_3 = MODIS_scf_2014 ;...
            MC_sim_scf_calib = MC_sim_scf_2012 ;  MC_sim_scf_analysis_1 = MC_sim_scf_2011 ; MC_sim_scf_analysis_2 = MC_sim_scf_2013 ;MC_sim_scf_analysis_3 = MC_sim_scf_2014 ;...
            MC_scf_stat_calib = MC_scf_stat_2012 ;  MC_scf_stat_analysis_1 = MC_scf_stat_2011 ; MC_scf_stat_analysis_2 = MC_scf_stat_2013 ;MC_scf_stat_analysis_3 = MC_scf_stat_2014 ;
    end
    if m==3; start_date_calib = '01-sep-2012'; start_date_1= '01-sep-2010';start_date_2= '01-sep-2011';start_date_3= '01-sep-2013'; ...
            MODIS_scf_calib = MODIS_scf_2013 ;  MODIS_scf_analysis_1 = MODIS_scf_2011 ; MODIS_scf_analysis_2 = MODIS_scf_2012 ;MODIS_scf_analysis_3 = MODIS_scf_2014 ;...
            MC_sim_scf_calib = MC_sim_scf_2013 ;  MC_sim_scf_analysis_1 = MC_sim_scf_2011 ; MC_sim_scf_analysis_2 = MC_sim_scf_2012 ;MC_sim_scf_analysis_3 = MC_sim_scf_2014 ;...
            MC_scf_stat_calib = MC_scf_stat_2013 ;  MC_scf_stat_analysis_1 = MC_scf_stat_2011 ; MC_scf_stat_analysis_2 = MC_scf_stat_2012 ;MC_scf_stat_analysis_3 = MC_scf_stat_2014 ;
    end   
    if m==4; start_date_calib = '01-sep-2013'; start_date_1= '01-sep-2010';start_date_2= '01-sep-2011';start_date_3= '01-sep-2012'; ...
            MODIS_scf_calib = MODIS_scf_2014 ;  MODIS_scf_analysis_1 = MODIS_scf_2011 ; MODIS_scf_analysis_2 = MODIS_scf_2012 ;MODIS_scf_analysis_3 = MODIS_scf_2013 ;...
            MC_sim_scf_calib = MC_sim_scf_2014 ;  MC_sim_scf_analysis_1 = MC_sim_scf_2011 ; MC_sim_scf_analysis_2 = MC_sim_scf_2012 ;MC_sim_scf_analysis_3 = MC_sim_scf_2013 ;...
            MC_scf_stat_calib = MC_scf_stat_2014 ;  MC_scf_stat_analysis_1 = MC_scf_stat_2011 ; MC_scf_stat_analysis_2 = MC_scf_stat_2012 ;MC_scf_stat_analysis_3 = MC_scf_stat_2013 ;
    end

%     summary_LHandCR = strcat('summary_LHandCR_',num2str(i));

%-----------------------------------------------------------------------------------------------
if strcmp(start_date_calib(8:end),'2010');  calib_year = '2011'; elseif strcmp(start_date_calib(8:end),'2011');  calib_year = '2012'; elseif strcmp(start_date_calib(8:end),'2012');  calib_year = '2013'; elseif strcmp(start_date_calib(8:end),'2013');  calib_year = '2014'; end
if strcmp(start_date_1(8:end),'2010');  analysis_year_1 = '2011'; elseif strcmp(start_date_1(8:end),'2011');  analysis_year_1 = '2012'; elseif strcmp(start_date_1(8:end),'2012');  analysis_year_1 = '2013'; elseif strcmp(start_date_1(8:end),'2013');  analysis_year_1 = '2014'; end
if strcmp(start_date_2(8:end),'2010');  analysis_year_2 = '2011'; elseif strcmp(start_date_2(8:end),'2011');  analysis_year_2 = '2012'; elseif strcmp(start_date_2(8:end),'2012');  analysis_year_2 = '2013'; elseif strcmp(start_date_2(8:end),'2013');  analysis_year_2 = '2014'; end
if strcmp(start_date_3(8:end),'2010');  analysis_year_3 = '2011'; elseif strcmp(start_date_3(8:end),'2011');  analysis_year_3 = '2012'; elseif strcmp(start_date_3(8:end),'2012');  analysis_year_3 = '2013'; elseif strcmp(start_date_3(8:end),'2013');  analysis_year_3 = '2014'; end


% GLUE analysis part - Q

% get obs Q
time_obs_Q = ncread('shyft/discharge_hourly_f.nc','time');
obs_Q_data = ncread('shyft/discharge_hourly_f.nc','discharge'); % (station, discharge for obs period), i.e. (6,1096)

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
    obs_Q_data(susp_index) = (susp_Q_bef+susp_Q_aft)/2;
end
%-----------------------------------------------------------------------------------------------
%%
%-----------------------------------------------------------------------
% calibration period - Q
%-----------------------------------------------------------------------
analysis_period = analysis_period_calib; 
start_date = datetime(start_date_calib);

if strcmp(start_date_calib(8:end),'2010');  warm_up_days = warm_up_days_2011; else warm_up_days = org_warm_up_days; end

% convert hourly observed data to daily time step
obs_Q_data = mean(reshape(obs_Q_data(:,:),24,[])); 
time_obs_Q = min(reshape(time_obs_Q,24,[]));

unix_start_date = posixtime(start_date);
index = find(time_obs_Q==unix_start_date);

obs_Q =obs_Q_data(index:index+analysis_period-1);

MC_sim_Q = ncread(strcat('MC_sim_Q_era_',calib_year,'_99999.nc'),'sim_Q'); % from hycamp\....\MC
Qsim = MC_sim_Q';

time_obs_Q = time_obs_Q(index:index+analysis_period-1);
norm_time_obs_Q =  datetime(time_obs_Q,'ConvertFrom','posixtime');

flow = obs_Q';

[nse,nse_log,sim_bias] = goodnessOfFitFuncn_shyft_Feb_23(obs_Q,Qsim, warm_up_days);

disp(['min_and_max_nse: ',num2str(min(nse)),'  ',num2str(max(nse))]);
disp(['min_and_max_nse_log: ',num2str(min(nse_log)),'  ',num2str(max(nse_log))]);

%-----------------------------------------------------------------------
% Likelihood measure
%-------------------------------------------------------------------------
min_NSE_for_plot = nse>=lowest_NSE_for_analysis;% exclude the extremely low performing models
if LH_measure_Q==1% NSE alone
    Y_Q = nse(min_NSE_for_plot); %av_wted_NSE ;%nse;
    
elseif LH_measure_Q==2% combined NSE and LnNSE    
    % specify the weight for NSE and LnNSE
    nse_wt = 0.54;
    nse_log_wt = 0.46;
    % calculate average weighted NSE
    av_wted_NSE = nse_wt*nse+nse_log_wt*nse_log;
    Y_Q = av_wted_NSE(min_NSE_for_plot); %av_wted_NSE ;%nse;
end
%%--------------------------------------------------------------------------
Qsim=Qsim(min_NSE_for_plot,:);

%-----------------------------------------------------------------------------------------------
% calibration period - SCF
%-----------------------------------------------------------------------
% get snow likelihood values for each model iteration
[mean_RMSE,mean_CSI] = getObsvnDayWeightedMeanRMSE_Mar2017(MC_scf_stat_calib,ccf_t,mean_LH_calc_opn);

disp(['min_and_max_mean_CSI: ',num2str(min(mean_CSI)),'  ',num2str(max(mean_CSI))]);
disp(['min_and_max_mean_RMSE: ',num2str(min(mean_RMSE)),'  ',num2str(max(mean_RMSE))]);

% Likelihood measure
%----------------------------
if LH_measure_SCA==1% CSI alone
    Y_SCA = mean_CSI';  
    threshold_SCA= threshold_1_SCA;
    Y_SCA = Y_SCA(min_NSE_for_plot); % remove models with extremely poor performance in Q simuln 
    
elseif LH_measure_SCA==2% RMSE alone    
    Y_SCA = mean_RMSE';
    threshold_SCA= threshold_2_SCA;
    Y_SCA = Y_SCA(min_NSE_for_plot);% remove models with extremely poor performance in Q simuln 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the index for the behavioral models 
% the indices from Q and SCA GLUE maximization should be the same
[ idx_Q, Llim, Ulim, median_predcn ] = shyft_GLUE_max(Y_Q,threshold_Q,Qsim,0.05) ;
no_of_behav_models_Q =  sum(idx_Q); % for summary output

idx_SCA = Y_SCA<=threshold_SCA;
no_of_behav_models_SCA =  sum(idx_SCA); % for summary output

idx_Q_SCA = idx_Q + idx_SCA;
% idx_Q_SCA_c = idx_Q_SCA > 0.9; % union, i.e. behav models on either Q or SC
idx_Q_SCA_c = idx_Q_SCA > 1.9;% intersection, i.e. behav models on both Q and SC

Y_SCA_max = 1./Y_SCA; % inverse RMSE

% normaize the weights of Q and SC between 0 and 1
Y_Q_2 = Y_Q;
Y_Q_2(~idx_Q)=0;
Y_SCA_max_2 = Y_SCA_max;
Y_SCA_max_2(~idx_SCA)=0;
Y_Q_norm_2 = (1/sum(Y_Q_2))*Y_Q_2; 
Y_SCA_norm_2 = (1/sum(Y_SCA_max_2))*Y_SCA_max_2; 

% combine the normalized weights of Q and SC
Y_Q_SCA_wted = Y_SCA_norm_2.*Y_Q_norm_2 ;

%----------------------------------------------------------------------------------
% get the index for the behavioral models as well as Llim, Ulim and median_predcn
% the indices from Q and SCA GLUE maximization should be the same
[ idx, Llim, Ulim, median_predcn ] = shyft_GLUE_max_May2017(Y_Q_SCA_wted,idx_Q_SCA_c,Qsim,0.05) ;
no_of_behav_models =  sum(idx); % for summary output
% % Display number of behavioural simulations
disp(['Number of behavioural simulations  =  ', num2str(sum(idx))]);
[nse_median,nse_log_median,sim_bias_median] = goodnessOfFitFuncn_shyft_Feb_23(obs_Q,median_predcn', warm_up_days);
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
CR=I_Qobs/length(obs_Q);
disp(['CR = ','  ', num2str(CR)]);

% %%----------------------------------------------------------------------------
[CSI_median,sim_bias_median,RMSE_median,SCF_CR]= evalMedianSCA_May2017(Y_Q_SCA_wted,idx_Q_SCA_c,MC_scf_stat_calib,MODIS_scf_calib,MC_sim_scf_calib,ccf_t,min_NSE_for_plot);
 
disp(['calibration_wted seasonal av: CSI_median, sim_bias_median and RMSE_median', '  ' ,num2str(CSI_median),'  ' , num2str(sim_bias_median),'  ' , num2str(RMSE_median) ]);

%%
% -------------------------------------------------------------------------
% Evaluation - analysis_period_1
% -------------------------------------------------------------------------
% Q related

analysis_period = analysis_period_1;
start_date = start_date_1;
if strcmp(start_date_1(8:end),'2010');  warm_up_days = warm_up_days_2011; else warm_up_days = org_warm_up_days; end

time_obs_Q = ncread('shyft/discharge_hourly.nc','time');
MC_sim_Q = ncread(strcat('MC_sim_Q_era_',analysis_year_1,'_99999.nc'),'sim_Q'); % from hycamp\....\MC

[nse_median_eval_1,nse_log_median_eval_1,sim_bias_median_eval_1,CR_1] = evalMedian_Q_May2017...
    (start_date,analysis_period,Y_Q_SCA_wted,idx_Q_SCA_c,time_obs_Q,obs_Q_data,MC_sim_Q,warm_up_days,min_NSE_for_plot);

% -------------------------------------------------------------------------
% snow related
analysis_year = analysis_year_1;

[CSI_median_1,sim_bias_median_1, RMSE_median_1,SCF_CR_1]= evalMedianSCA_May2017(Y_Q_SCA_wted,idx_Q_SCA_c,MC_scf_stat_analysis_1,MODIS_scf_analysis_1,MC_sim_scf_analysis_1,ccf_t,min_NSE_for_plot);
 
disp(['Analysis 1: wted seasonal av: CSI_median, sim_bias_median and RMSE_median', '  ' ,num2str(CSI_median_1),'  ' , num2str(sim_bias_median_1),'  ' , num2str(RMSE_median_1) ]);

%%
% -------------------------------------------------------------------------
% Evaluation - analysis_period_2
% -------------------------------------------------------------------------
% Q related

analysis_period = analysis_period_2;
start_date = start_date_2;
if strcmp(start_date_2(8:end),'2010');  warm_up_days = warm_up_days_2011; else warm_up_days = org_warm_up_days; end
time_obs_Q = ncread('shyft/discharge_hourly.nc','time');
MC_sim_Q = ncread(strcat('MC_sim_Q_era_',analysis_year_2,'_99999.nc'),'sim_Q'); % from hycamp\....\MC
[nse_median_eval_2,nse_log_median_eval_2,sim_bias_median_eval_2,CR_2] = evalMedian_Q_May2017...
    (start_date,analysis_period,Y_Q_SCA_wted,idx_Q_SCA_c,time_obs_Q,obs_Q_data,MC_sim_Q,warm_up_days,min_NSE_for_plot);

% -------------------------------------------------------------------------
% snow related
analysis_year = analysis_year_2;
[CSI_median_2,sim_bias_median_2, RMSE_median_2,SCF_CR_2]= evalMedianSCA_May2017(Y_Q_SCA_wted,idx_Q_SCA_c,MC_scf_stat_analysis_2,MODIS_scf_analysis_2,MC_sim_scf_analysis_2,ccf_t,min_NSE_for_plot);
 disp(['Analysis 2: wted seasonal av: CSI_median, sim_bias_median and RMSE_median', '  ' ,num2str(CSI_median_2),'  ' , num2str(sim_bias_median_2),'  ' , num2str(RMSE_median_2) ]);

%%
% -------------------------------------------------------------------------
% Evaluation - analysis_period_3
% -------------------------------------------------------------------------
% Q related

analysis_period = analysis_period_3;
start_date = start_date_3;
if strcmp(start_date_3(8:end),'2010');  warm_up_days = warm_up_days_2011; else warm_up_days = org_warm_up_days; end

time_obs_Q = ncread('shyft/discharge_hourly.nc','time');
MC_sim_Q = ncread(strcat('MC_sim_Q_era_',analysis_year_3,'_99999.nc'),'sim_Q'); % from hycamp\....\MC

[nse_median_eval_3,nse_log_median_eval_3,sim_bias_median_eval_3,CR_3] = evalMedian_Q_May2017...
    (start_date,analysis_period,Y_Q_SCA_wted,idx_Q_SCA_c,time_obs_Q,obs_Q_data,MC_sim_Q,warm_up_days,min_NSE_for_plot);

% -------------------------------------------------------------------------
% snow related
analysis_year = analysis_year_3;
[CSI_median_3,sim_bias_median_3, RMSE_median_3,SCF_CR_3]= evalMedianSCA_May2017(Y_Q_SCA_wted,idx_Q_SCA_c,MC_scf_stat_analysis_3,MODIS_scf_analysis_3,MC_sim_scf_analysis_3,ccf_t,min_NSE_for_plot);
disp(['Analysis 2: wted seasonal av: CSI_median, sim_bias_median and RMSE_median', '  ' ,num2str(CSI_median_2),'  ' , num2str(sim_bias_median_2),'  ' , num2str(RMSE_median_2) ]);

%%
% -------------------------------------------------------------------------
% Output Summary
% -------------------------------------------------------------------------

if m==1
summary_LHandCR_colmn_format_i = [nse_median,nse_log_median,CR,RMSE_median,CSI_median,SCF_CR,...
    nse_median_eval_1,nse_log_median_eval_1,CR_1,RMSE_median_1,CSI_median_1,SCF_CR_1,...
    nse_median_eval_2,nse_log_median_eval_2,CR_2,RMSE_median_2,CSI_median_2,SCF_CR_2,...
    nse_median_eval_3,nse_log_median_eval_3,CR_3,RMSE_median_3,CSI_median_3,SCF_CR_3,...
    no_of_behav_models];
elseif m==2
summary_LHandCR_colmn_format_i = [
    nse_median_eval_1,nse_log_median_eval_1,CR_1,RMSE_median_1,CSI_median_1,SCF_CR_1,...
    nse_median,nse_log_median,CR,RMSE_median,CSI_median,SCF_CR,...
    nse_median_eval_2,nse_log_median_eval_2,CR_2,RMSE_median_2,CSI_median_2,SCF_CR_2,...
    nse_median_eval_3,nse_log_median_eval_3,CR_3,RMSE_median_3,CSI_median_3,SCF_CR_3,...
    no_of_behav_models];
elseif m==3
    summary_LHandCR_colmn_format_i = [
    nse_median_eval_1,nse_log_median_eval_1,CR_1,RMSE_median_1,CSI_median_1,SCF_CR_1,...
    nse_median_eval_2,nse_log_median_eval_2,CR_2,RMSE_median_2,CSI_median_2,SCF_CR_2,...
    nse_median,nse_log_median,CR,RMSE_median,CSI_median,SCF_CR,...
    nse_median_eval_3,nse_log_median_eval_3,CR_3,RMSE_median_3,CSI_median_3,SCF_CR_3,...
    no_of_behav_models];
elseif m==4
    summary_LHandCR_colmn_format_i = [
    nse_median_eval_1,nse_log_median_eval_1,CR_1,RMSE_median_1,CSI_median_1,SCF_CR_1,...
    nse_median_eval_2,nse_log_median_eval_2,CR_2,RMSE_median_2,CSI_median_2,SCF_CR_2,...
    nse_median_eval_3,nse_log_median_eval_3,CR_3,RMSE_median_3,CSI_median_3,SCF_CR_3,...
    nse_median,nse_log_median,CR,RMSE_median,CSI_median,SCF_CR,...
    no_of_behav_models];
end
    
summary_LHandCR_colmn_format(m,:)=summary_LHandCR_colmn_format_i;

end