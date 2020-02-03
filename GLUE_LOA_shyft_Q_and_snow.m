% This algorithm is used to conduct uncertainty analysis based on the time relaxed limits
% of acceptability (GLUE pLoA) concept and a combined likelihood of streamflow and fSCA. This
% algorithm modifies the original formulation of GLUE LoA. Here model realizations that
% satisfy the LoA criteria above certain percentage of the observation time steps (pLoA) are
% considered behavioural. Outputs from this script include the cross-validation plots and tabular
% result of the evaluation metrics for streamflow and fSCA both for the calibration and
% validation years.
% This script was adopted from and calls to other scripts of the SAFE Toolbox (Pianosi et al, 2015)

% For details on the conceptual background, implementation, and sample outputs from this algorithm, 
% the reader is referred to the following paper: 

% Teweldebrhan, A. T., Burkhart, J. F., and Schuler, T. V.: Parameter uncertainty analysis for an 
% operational hydrological model using residual-based and limits of acceptability approaches, 
% Hydrology and Earth System Sciences, 22, 5021-5039, 2018.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all
format long g

%-----------------------------------------------------------------------------------------------
% general info
mc_calib_param_file = 'shyft/24-Mc_calib_param_new_cv_100k.nc';
analysis_period_calib = 365;
analysis_period_1 = 365;
analysis_period_2 = 365;
analysis_period_3 = 365;

org_warm_up_days = 30; 
warm_up_days_2011 = 60; 
analysis_end_indx = analysis_period_calib ;
TS_excl = [] 

% set the lower and upper acceptable limit multipliers
LA_lim = 0.75 ;
UA_lim = 1.25 ;
plots = 0;

% Specify likelhood-weighted percentiles to be returned from GLUE-LOA
Pc = [0.05 0.5 0.95];

% ** check summary_LHandCR at workspace for output
%--------------------------------------------------------------------------------------------------
% SC related
ccf_t = 0.2;% max acceptable SCA uncertainty
SCF_LA_lim = 0.75; 
SCF_UA_lim = 1.50; 

%--------------------------------------------------------------------------------------------------
%%%%%Section 1 Set up project and run simulations
my_dir = pwd ; 
cd(my_dir)

addpath('util');
addpath('visualization')
addpath('GLUE')
addpath('shyft')
addpath('x_extra')
addpath('\\lagringshotell\geofag\projects\hycamp\team\aynomtt\data\MC\snow_cv_included\')

N = 100000 ; % sample size parameterisations

M  = 21; %5 ; % number of uncertain parameters 

% S_No =   [1         2       3           4       5       6       7       8       9       10      11      12      13      14      15      16      17      18      19      20      21   ];
S_No =     {'1'       '2'       '3'       '4'      '5'    '6'     '7'     '8'     '9'      '10'   '11'    '12'    '13'    '14'    '15'    '16'    '17'    '18'    '19'    '20'    '21'   };
% X_Labels = {'c1'     'c2'       'c3'    'ae scale factor'      'tx'    'wind scale'     'max water'     'wind const'     'fast albedo decay rate'      'slow albedo decay rate'   'surface magnitude'    'max albedo'    'min albedo'    'snowfall reset depth'    'snow cv'    'snow cv forest factor'    'snow cv altitude factor'    'glacier albedo'    'scale factor'    'albedo'    'alpha'   };
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
Xrec(:,:) = X; %record of parameters

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

for m=1:4 %m=1:4

    disp(['current calibrn year is:', '  ' ,num2str(m)]);
    
    if m==1; start_date_calib = '01-sep-2010'; start_date_1= '01-sep-2011';start_date_2= '01-sep-2012';start_date_3= '01-sep-2013';...
            MODIS_scf_calib = MODIS_scf_2011 ;  MODIS_scf_analysis_1 = MODIS_scf_2012 ; MODIS_scf_analysis_2 = MODIS_scf_2013 ;MODIS_scf_analysis_3 = MODIS_scf_2014 ;...
            MC_sim_scf_calib = MC_sim_scf_2011 ;  MC_sim_scf_analysis_1 = MC_sim_scf_2012 ; MC_sim_scf_analysis_2 = MC_sim_scf_2013 ;MC_sim_scf_analysis_3 = MC_sim_scf_2014 ;...
            MC_scf_stat_calib = MC_scf_stat_2011 ;  MC_scf_stat_analysis_1 = MC_scf_stat_2012 ; MC_scf_stat_analysis_2 =MC_scf_stat_2013 ;MC_scf_stat_analysis_3 = MC_scf_stat_2014 ;
            p_LOA = 37; %25; %37; % using residual based GLUE 5-95% as ref = 37
    end  
    if m==2; start_date_calib = '01-sep-2011'; start_date_1= '01-sep-2010';start_date_2= '01-sep-2012';start_date_3= '01-sep-2013'; ...
            MODIS_scf_calib = MODIS_scf_2012 ;  MODIS_scf_analysis_1 = MODIS_scf_2011 ; MODIS_scf_analysis_2 = MODIS_scf_2013 ;MODIS_scf_analysis_3 = MODIS_scf_2014 ;...
            MC_sim_scf_calib = MC_sim_scf_2012 ;  MC_sim_scf_analysis_1 = MC_sim_scf_2011 ; MC_sim_scf_analysis_2 = MC_sim_scf_2013 ;MC_sim_scf_analysis_3 = MC_sim_scf_2014 ;...
            MC_scf_stat_calib = MC_scf_stat_2012 ;  MC_scf_stat_analysis_1 = MC_scf_stat_2011 ; MC_scf_stat_analysis_2 = MC_scf_stat_2013 ;MC_scf_stat_analysis_3 = MC_scf_stat_2014 ;
            p_LOA = 40; %30; %40;% using residual based GLUE 5-95% as ref = 40
    end
    if m==3; start_date_calib = '01-sep-2012'; start_date_1= '01-sep-2010';start_date_2= '01-sep-2011';start_date_3= '01-sep-2013'; ...
            MODIS_scf_calib = MODIS_scf_2013 ;  MODIS_scf_analysis_1 = MODIS_scf_2011 ; MODIS_scf_analysis_2 = MODIS_scf_2012 ;MODIS_scf_analysis_3 = MODIS_scf_2014 ;...
            MC_sim_scf_calib = MC_sim_scf_2013 ;  MC_sim_scf_analysis_1 = MC_sim_scf_2011 ; MC_sim_scf_analysis_2 = MC_sim_scf_2012 ;MC_sim_scf_analysis_3 = MC_sim_scf_2014 ;...
            MC_scf_stat_calib = MC_scf_stat_2013 ;  MC_scf_stat_analysis_1 = MC_scf_stat_2011 ; MC_scf_stat_analysis_2 = MC_scf_stat_2012 ;MC_scf_stat_analysis_3 = MC_scf_stat_2014 ;
            p_LOA = 30; %28; %30; % using residual based GLUE 5-95% as ref = 30
    end   
    if m==4; start_date_calib = '01-sep-2013'; start_date_1= '01-sep-2010';start_date_2= '01-sep-2011';start_date_3= '01-sep-2012'; ...
            MODIS_scf_calib = MODIS_scf_2014 ;  MODIS_scf_analysis_1 = MODIS_scf_2011 ; MODIS_scf_analysis_2 = MODIS_scf_2012 ;MODIS_scf_analysis_3 = MODIS_scf_2013 ;...
            MC_sim_scf_calib = MC_sim_scf_2014 ;  MC_sim_scf_analysis_1 = MC_sim_scf_2011 ; MC_sim_scf_analysis_2 = MC_sim_scf_2012 ;MC_sim_scf_analysis_3 = MC_sim_scf_2013 ;...
            MC_scf_stat_calib = MC_scf_stat_2014 ;  MC_scf_stat_analysis_1 = MC_scf_stat_2011 ; MC_scf_stat_analysis_2 = MC_scf_stat_2012 ;MC_scf_stat_analysis_3 = MC_scf_stat_2013 ;
            p_LOA = 30; %25; %30; % using residual based GLUE 5-95% as ref = 30
    end

    %     summary_LHandCR = strcat('summary_LHandCR_',num2str(i));
    
    %-----------------------------------------------------------------------------------------------
    % the simulated value nc files are named as 2011, 2012, 2013 and 2014
    % instead of starting from 2010. Thus following lines are mainly meant to conform
    % the nc file names with the actural simulation start dates specified above. 
    if strcmp(start_date_calib(8:end),'2010');  calib_year = '2011'; elseif strcmp(start_date_calib(8:end),'2011');  calib_year = '2012'; elseif strcmp(start_date_calib(8:end),'2012');  calib_year = '2013'; elseif strcmp(start_date_calib(8:end),'2013');  calib_year = '2014'; end
    if strcmp(start_date_1(8:end),'2010');  analysis_year_1 = '2011'; elseif strcmp(start_date_1(8:end),'2011');  analysis_year_1 = '2012'; elseif strcmp(start_date_1(8:end),'2012');  analysis_year_1 = '2013'; elseif strcmp(start_date_1(8:end),'2013');  analysis_year_1 = '2014'; end
    if strcmp(start_date_2(8:end),'2010');  analysis_year_2 = '2011'; elseif strcmp(start_date_2(8:end),'2011');  analysis_year_2 = '2012'; elseif strcmp(start_date_2(8:end),'2012');  analysis_year_2 = '2013'; elseif strcmp(start_date_2(8:end),'2013');  analysis_year_2 = '2014'; end
    if strcmp(start_date_3(8:end),'2010');  analysis_year_3 = '2011'; elseif strcmp(start_date_3(8:end),'2011');  analysis_year_3 = '2012'; elseif strcmp(start_date_3(8:end),'2012');  analysis_year_3 = '2013'; elseif strcmp(start_date_3(8:end),'2013');  analysis_year_3 = '2014'; end

    % GLUE analysis part

    % get the monte carlo calib param value based ensembel simulations from shyft
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

    %-----------------------------------------------------------------------
    % calibration period - Q
    %-----------------------------------------------------------------------
    analysis_period = analysis_period_calib; 
    start_date = datetime(start_date_calib)
    if strcmp(start_date_calib(8:end),'2010');  warm_up_days = warm_up_days_2011; else warm_up_days = org_warm_up_days; end

    % convert hourly observed data to daily time step
    obs_Q_data = mean(reshape(obs_Q_data(:,:),24,[]));  
    time_obs_Q = min(reshape(time_obs_Q,24,[]));

    unix_start_date = posixtime(start_date);
    index = find(time_obs_Q==unix_start_date);

    obs_Q =obs_Q_data(index:index+analysis_period-1);

    MC_sim_Q = ncread(strcat('MC_sim_Q_era_',calib_year,'_99999.nc'),'sim_Q'); 
    Qsim = MC_sim_Q';

    size(obs_Q);
    size(Qsim);

    time_obs_Q = time_obs_Q(index:index+analysis_period-1);
    norm_time_obs_Q =  datetime(time_obs_Q,'ConvertFrom','posixtime');

    flow = obs_Q';

    Q_LOA(:,1) = flow*LA_lim;
    Q_LOA(:,2) = flow*UA_lim;
    
    TS = [warm_up_days:analysis_end_indx] ;
    TS=TS(~ismember(TS,TS_excl));

    % Compute weigting for each time series:
    wtg = zeros(length(flow),N); 
    for ii = 1: length(flow)
        wtg(ii,:) = trapmf(Qsim(:,ii),[Q_LOA(ii,1) flow(ii)  flow(ii) Q_LOA(ii,2)] );
    end
    LMwt_Q = sum(wtg(TS,:)); % Sum generalised likelihoods (weightings).

    [idx_Q_beforeRelaxing, id_LOA, Pcnts_Q_beforeRelaxing] = GLUE_LOA(N,LMwt_Q,Q_LOA(TS,:),Qsim,TS,Pc);

    % following lines are intended for use during the iterations to relax in time when all the observations fail to
    % fall with in the prediction bounds
    no_LOA_t = p_LOA*length(TS)/100; 
    idx_relaxed_LOA = sum(id_LOA(:,:,1))>=no_LOA_t; 
    idx_Q_afterRelaxing = idx_relaxed_LOA'; 
    % just to check how many behavioral models are found satisfying the stated criteria
    no_idx_complying = sum(idx_Q_afterRelaxing)

    %-----------------------------------------------------------------------
    % calibration period - SCF related
    %-----------------------------------------------------------------------
    % get behavioral models and their weight based on their SC performance
    [idx_SC,LMwt_SC]=evalMedianSCA_LOA_Aug2017_2(MC_scf_stat_calib,MODIS_scf_calib,MC_sim_scf_calib,SCF_LA_lim,SCF_UA_lim,N,Pc,p_LOA,ccf_t,m);

    % get index of behavioral models based on both Q and SC observations
    idx = (idx_Q_afterRelaxing+idx_SC) > 1.9; %== 2.0; % get model realizations behavioral both in Q and SC simulations 
    sum_idx_test = sum(idx_Q_afterRelaxing); 
    sum_idx_SC_test = sum(idx_SC)
    sum_idx=sum(idx)

    %--------------------------------------------------------------------------------------------------------------------------

    % normaize the weights of Q and SC such that the sum of their
    % respective weights would add up to unity
    LMwt_Q_norm = (1/sum(LMwt_Q))*LMwt_Q; % normalize weight of Q 
    LMwt_SC_norm = (1/sum(LMwt_SC))*LMwt_SC; % normalize weight of SC
    LMwt = LMwt_Q_norm + LMwt_SC_norm; % combine LH measure by addition
    
    % get percentiles using the acceptable models after relaxing in time and
    % combining behavioral models in both Q and SC
    
    Pcnts(:,1:length(Pc)) = DeterminePercents(LMwt((idx(:,1)==1)),...
        Pc, Qsim((idx(:,1)==1),:));

    I_Qobs=0;
    for t=1:length(obs_Q)

        obs_Q_t = obs_Q(t);
        Llim_t = Pcnts(t,1); %Q_LOA(t,1);
        Ulim_t = Pcnts(t,3); %Q_LOA(t,2);
        if (Llim_t<=obs_Q_t)&&(obs_Q_t<=Ulim_t)
            I_Qobs=I_Qobs+1;
        end    
    end
    CR=I_Qobs/length(obs_Q);
    disp(['CR = ','  ', num2str(CR)]);

    median_predcn = Pcnts(:,2);

    [nse_median,nse_log_median,sim_bias_median] = goodnessOfFitFuncn_shyft_Feb_23(obs_Q,median_predcn', warm_up_days);
    disp(['calibrn: nse and nse_log for median crisp forecast', '  ' ,num2str(nse_median),'  ' , num2str(nse_log_median) ]);

    %--------------------------------------------------------------------------------------------------------------------------
    % SC related - calibration
    [CSI_median,sim_bias_median,RMSE_median,SCF_CR]= ...
        evalMedianSCA_LOA_Aug2017(MC_scf_stat_calib,MODIS_scf_calib,MC_sim_scf_calib,SCF_LA_lim,SCF_UA_lim,LMwt,idx,Pc,ccf_t,calib_year,m);

    disp(['calibration_wted seasonal av: CSI_median, sim_bias_median and RMSE_median', '  ' ,num2str(CSI_median),'  ' , num2str(sim_bias_median),'  ' , num2str(RMSE_median) ]);

    %--------------------------------------------------------------------------------------------------------------------------
    %%%% VISUALIZATION TOOLS
    if m==2% if calib year = 2011

    %     Plot the number of acceptable simulations based upon the LOA for the
    %     timesteps over the evaluation period
        figure
        sz = size(LMwt');
        subplot(2,1,1);plot(flow,'Color',rgb('MidnightBlue'),'LineWidth',2); tx = axis; tx(4) = 100;
        xlim([TS(1) TS(end)]);
        set(gca,'XTick',[])
        ylabel('streamflow (m^3 s^{-1})')
        subplot(2,1,2);bar(TS,sum(id_LOA')/sz(1)*100); axis(tx);
        ylabel('Acceptable LOA (%)')
        ylim([0 70]);
        xlim([TS(1) TS(end)]);
        
        startDate = datenum('01-Oct-2011');
        dateaxis('x',12,startDate)

        % suptitle('Percentage of simulations accepted at each timestep')
        disp('Press any key to continue');pause
        
        % Save the number of acceptable models in each day of current calib year
        acceptableloa = [TS;flow(TS)';sum(id_LOA')/sz(1)*100];       
        save acceptableloa_2013.mat acceptableloa;

    %     Plot uncertainty percentiles chosen
        figure(1)
        color1 = rgb('SkyBlue');
        color2 = rgb('DodgerBlue');
        [ha1 hb1 hc1]= shade_area_between_curves_LOA(Pcnts(:,3), Pcnts(:,1),color1,color1);
        hold on
        [ha2 hb2 hc2] = shade_area_between_curves_LOA(Q_LOA(:,2), Q_LOA(:,1),color2,color2);
        hold on
        h3= plot(Pcnts(:,2),'-k','LineWidth',1.5,'userdata', '50%');
        hold on
        h4 = plot(flow,'ob', 'userdata', 'flow (m3/s)','MarkerFaceColor','b','MarkerSize',3);
        hold on
        ylabel('streamflow (m^3 s^{-1})')
        xlim([TS(1) TS(end)]);% set the x range to be displayed in x axis
        legend([ha1(2) ha2(2) h3 h4],{'5-95% prediction bound','acceptable flow bound', 'median prediction','observed streamflow'});

        xlim([30 365])
         
        startDate = datenum('01-Oct-2011');
        dateaxis('x',12,startDate)
             
        %---------------------------------------------------------------------
        % trial to correct the date spacing
        figure %(1)
        
        color1 = rgb('SkyBlue');
        color2 = rgb('DodgerBlue');
        [ha1 hb1 hc1]= shade_area_between_curves_LOA(Pcnts(:,3), Pcnts(:,1),color1,color1);
        hold on
        [ha2 hb2 hc2] = shade_area_between_curves_LOA(Q_LOA(:,2), Q_LOA(:,1),color2,color2);
        hold on
        h3= plot(Pcnts(:,2),'-k','LineWidth',1.5,'userdata', '50%');
        hold on
        h4 = plot(flow,'ob', 'userdata', 'flow (m3/s)','MarkerFaceColor','b','MarkerSize',3);
        hold on
        ylabel('streamflow (m^3 s^{-1})')
        xlim([TS(1) TS(end)]);% set the x range to be displayed in x axis
        legend([ha1(2) ha2(2) h3 h4],{'5-95% prediction bound','acceptable flow bound', 'median prediction','observed streamflow'});
        startDate = datenum(start_date+warm_up_days);
        dateaxis('x',12,startDate)
        %-----------------------------
        startDate = datenum(start_date+warm_up_days); 
        finalDate = datenum(startDate+length(TS)-1);
        xData = linspace(startDate,finalDate,length(TS));
        ticks = linspace(startDate, finalDate, 6); 
       
        labels = datestr(ticks, 12); 
        set(gca, 'xtick', ticks, 'xtickLabel', labels);

    % -------------------------------------------------------------------------
    % Evaluation-1
    % -------------------------------------------------------------------------
    %% Q related
    analysis_period = analysis_period_1;
    time_obs_Q = ncread('shyft/discharge_hourly.nc','time');
    MC_sim_Q = ncread(strcat('MC_sim_Q_era_',analysis_year_1,'_99999.nc'),'sim_Q'); % from hycamp\....\MC
    start_date = start_date_1;
    if strcmp(start_date_1(8:end),'2010');  warm_up_days = warm_up_days_2011; else warm_up_days = org_warm_up_days; end

    [nse_median_eval_1,nse_log_median_eval_1,sim_bias_median_eval_1,CR_1] = ...
        evalMedian_Q_LOA_May2017(start_date,analysis_period,LMwt,idx,time_obs_Q,obs_Q_data, MC_sim_Q,warm_up_days,LA_lim,UA_lim,Pc,TS,m);

    % snow related
    [CSI_median_1,sim_bias_median_1,RMSE_median_1,SCF_CR_1]= ...
        evalMedianSCA_LOA_Aug2017(MC_scf_stat_analysis_1,MODIS_scf_analysis_1,MC_sim_scf_analysis_1,SCF_LA_lim,SCF_UA_lim,LMwt,idx,Pc,ccf_t,analysis_year_1,m); 
    disp(['Analysis 1: wted seasonal av: CSI_median, sim_bias_median and RMSE_median', '  ' ,num2str(CSI_median_1),'  ' , num2str(sim_bias_median_1),'  ' , num2str(RMSE_median_1) ]);

    % -------------------------------------------------------------------------
    % Evaluation-2
    % -------------------------------------------------------------------------
    %% Q related
    analysis_period = analysis_period_2;
    time_obs_Q = ncread('shyft/discharge_hourly.nc','time');
    MC_sim_Q = ncread(strcat('MC_sim_Q_era_',analysis_year_2,'_99999.nc'),'sim_Q'); % from hycamp\....\MC
    start_date = start_date_2;
    if strcmp(start_date_2(8:end),'2010');  warm_up_days = warm_up_days_2011; else warm_up_days = org_warm_up_days; end
    figure(2)
    [nse_median_eval_2,nse_log_median_eval_2,sim_bias_median_eval_2,CR_2] = ...
        evalMedian_Q_LOA_May2017(start_date,analysis_period,LMwt,idx,time_obs_Q,obs_Q_data, MC_sim_Q,warm_up_days,LA_lim,UA_lim,Pc,TS,m);
    
    startDate = datenum('01-Oct-2012');
    dateaxis('x',12,startDate)

    % snow related
    [CSI_median_2,sim_bias_median_2,RMSE_median_2,SCF_CR_2]= ...
        evalMedianSCA_LOA_Aug2017(MC_scf_stat_analysis_2,MODIS_scf_analysis_2,MC_sim_scf_analysis_2,SCF_LA_lim,SCF_UA_lim,LMwt,idx,Pc,ccf_t,analysis_year_2,m); 
    disp(['Analysis 2: wted seasonal av: CSI_median, sim_bias_median and RMSE_median', '  ' ,num2str(CSI_median_2),'  ' , num2str(sim_bias_median_2),'  ' , num2str(RMSE_median_2) ]);


    % -------------------------------------------------------------------------
    % Evaluation-3
    % -------------------------------------------------------------------------
    analysis_period = analysis_period_3;
    time_obs_Q = ncread('shyft/discharge_hourly.nc','time');
    MC_sim_Q = ncread(strcat('MC_sim_Q_era_',analysis_year_3,'_99999.nc'),'sim_Q'); % from hycamp\....\MC
    start_date = start_date_3;
    if strcmp(start_date_3(8:end),'2010');  warm_up_days = warm_up_days_2011; else warm_up_days = org_warm_up_days; end

    [nse_median_eval_3,nse_log_median_eval_3,sim_bias_median_eval_3,CR_3] = ...
        evalMedian_Q_LOA_May2017(start_date,analysis_period,LMwt,idx,time_obs_Q,obs_Q_data, MC_sim_Q,warm_up_days,LA_lim,UA_lim,Pc,TS,m);

    % snow related
    [CSI_median_3,sim_bias_median_3,RMSE_median_3,SCF_CR_3]= ...
        evalMedianSCA_LOA_Aug2017(MC_scf_stat_analysis_3,MODIS_scf_analysis_3,MC_sim_scf_analysis_3,SCF_LA_lim,SCF_UA_lim,LMwt,idx,Pc,ccf_t,analysis_year_3,m); 
    disp(['Analysis 2: wted seasonal av: CSI_median, sim_bias_median and RMSE_median', '  ' ,num2str(CSI_median_2),'  ' , num2str(sim_bias_median_2),'  ' , num2str(RMSE_median_2) ]);

    % -------------------------------------------------------------------------

    no_of_behav_models = sum(idx);
    % summary_LHandCR = [nse_median,nse_log_median,CR_1,nse_median_eval_1,nse_log_median_eval_1,CR_2,nse_median_eval_2,nse_log_median_eval_2,CR_3,no_of_behav_models_comb];
    %% save to summary table
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

    summary_LHandCR_colmn_format_LoA(m,:)=summary_LHandCR_colmn_format_i;

end

% put figures 1 and 2 as subplots of a single figure (fig. 3)
% N.B. the figures need be renamed to 'Figure 1' and 'Figure 2'
figure(3)% need be same with (note the 3)
ax = zeros(2,1);
for i = 1:2
    ax(i)=subplot(2,1,i);
end
for i = 1:2
    figure(i)
    h = get(gcf,'Children');
    newh = copyobj(h,3);% need be same with (note the 3 in (h,3))
    for j = length(newh)
        posnewh = get(newh(j),'Position');
        possub  = get(ax(i),'Position');
        set(newh(j),'Position',[posnewh(1) possub(2) posnewh(3) possub(4)]);
    end
delete(ax(i));
end
figure(3)% need be same with(note the 3)

% export the figure in same format as it appears on screen
addpath('x_extra\altmany-export_fig-5be2ca4')
set(gcf, 'Color', 'w');% set figure background to white
% export_fig 'Q_uncertainty_bounds_LOA' -tiff -painters -r300;
% export_fig(gcf,'Q_uncertainty_bounds_LOA.fig', '-r300');
export_fig(gcf,'Q_uncertainty_bounds_LOA_xLabelCorrected.tiff', '-r300');

