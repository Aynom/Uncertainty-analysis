% This script was implemented to assess the effect of threshold pLOA on CR and NSE
% of median prediction when using the time relaxed GLUE LoA approach(GLUE pLoA)
% Main outputs of this script are plots displaying CR and NSE against pLoA values.

% For details on the conceptual background, implementation, and sample outputs from this algorithm, 
% the reader is referred to the following paper: 

% Teweldebrhan, A. T., Burkhart, J. F., and Schuler, T. V.: Parameter uncertainty analysis for an 
% operational hydrological model using residual-based and limits of acceptability approaches, 
% Hydrology and Earth System Sciences, 22, 5021-5039, 2018.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% clear all; close all
format long g

%-----------------------------------------------------------------------------------------------
% general info
mc_calib_param_file = 'shyft/24-Mc_calib_param_new_cv_100k.nc';
analysis_period_calib = 365;
analysis_period_1 = 365;
analysis_period_2 = 365;
analysis_period_3 = 365;

% exclude the first n warm_up_days from NSE estimate
org_warm_up_days = 30; 
warm_up_days_2011 = 60; 
analysis_end_indx = analysis_period_calib ;%-50;
TS_excl = [] %[60:260]; % addtional elements to exclude with in the start and end of analysis period

% set the lower and upper acceptable limit multipliers
LA_lim = 0.75 ;%0.75; %0.8 ;%0.4; % lower LA
UA_lim = 1.25 ;%1.25;   %1.2 ;%4; % upper LA
plots = 0;

% Specify likelhood-weighted percentiles to be returned from GLUE-LOA
Pc = [0.05 0.5 0.95];

% check summary_LHandCR at workspace for output
%--------------------------------------------------------------------------------------------------
% SC related
ccf_t = 0.2;% max acceptable cloud cover fraction
SCF_LA_lim = 0.75; 
SCF_UA_lim = 1.50; 
% SCA_LH_wt = 0.1;
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

counter=1;

for m=1:4 %m=1:4
    
    for p_LOA=80:-5:5

        disp(['current calibrn year is:', '  ' ,num2str(m)]);

        if m==1; start_date_calib = '01-sep-2010'; start_date_1= '01-sep-2011';start_date_2= '01-sep-2012';start_date_3= '01-sep-2013';...
                MODIS_scf_calib = MODIS_scf_2011 ;  MODIS_scf_analysis_1 = MODIS_scf_2012 ; MODIS_scf_analysis_2 = MODIS_scf_2013 ;MODIS_scf_analysis_3 = MODIS_scf_2014 ;...
                MC_sim_scf_calib = MC_sim_scf_2011 ;  MC_sim_scf_analysis_1 = MC_sim_scf_2012 ; MC_sim_scf_analysis_2 = MC_sim_scf_2013 ;MC_sim_scf_analysis_3 = MC_sim_scf_2014 ;...
                MC_scf_stat_calib = MC_scf_stat_2011 ;  MC_scf_stat_analysis_1 = MC_scf_stat_2012 ; MC_scf_stat_analysis_2 =MC_scf_stat_2013 ;MC_scf_stat_analysis_3 = MC_scf_stat_2014 ;
%                 p_LOA = 37; %25; %37; % using residual based GLUE 5-95% as ref = 37
        end  
        if m==2; start_date_calib = '01-sep-2011'; start_date_1= '01-sep-2010';start_date_2= '01-sep-2012';start_date_3= '01-sep-2013'; ...
                MODIS_scf_calib = MODIS_scf_2012 ;  MODIS_scf_analysis_1 = MODIS_scf_2011 ; MODIS_scf_analysis_2 = MODIS_scf_2013 ;MODIS_scf_analysis_3 = MODIS_scf_2014 ;...
                MC_sim_scf_calib = MC_sim_scf_2012 ;  MC_sim_scf_analysis_1 = MC_sim_scf_2011 ; MC_sim_scf_analysis_2 = MC_sim_scf_2013 ;MC_sim_scf_analysis_3 = MC_sim_scf_2014 ;...
                MC_scf_stat_calib = MC_scf_stat_2012 ;  MC_scf_stat_analysis_1 = MC_scf_stat_2011 ; MC_scf_stat_analysis_2 = MC_scf_stat_2013 ;MC_scf_stat_analysis_3 = MC_scf_stat_2014 ;
%                 p_LOA = 40; %30; %40;% using residual based GLUE 5-95% as ref = 40
        end
        if m==3; start_date_calib = '01-sep-2012'; start_date_1= '01-sep-2010';start_date_2= '01-sep-2011';start_date_3= '01-sep-2013'; ...
                MODIS_scf_calib = MODIS_scf_2013 ;  MODIS_scf_analysis_1 = MODIS_scf_2011 ; MODIS_scf_analysis_2 = MODIS_scf_2012 ;MODIS_scf_analysis_3 = MODIS_scf_2014 ;...
                MC_sim_scf_calib = MC_sim_scf_2013 ;  MC_sim_scf_analysis_1 = MC_sim_scf_2011 ; MC_sim_scf_analysis_2 = MC_sim_scf_2012 ;MC_sim_scf_analysis_3 = MC_sim_scf_2014 ;...
                MC_scf_stat_calib = MC_scf_stat_2013 ;  MC_scf_stat_analysis_1 = MC_scf_stat_2011 ; MC_scf_stat_analysis_2 = MC_scf_stat_2012 ;MC_scf_stat_analysis_3 = MC_scf_stat_2014 ;
%                 p_LOA = 30; %28; %30; % using residual based GLUE 5-95% as ref = 30
        end   
        if m==4; start_date_calib = '01-sep-2013'; start_date_1= '01-sep-2010';start_date_2= '01-sep-2011';start_date_3= '01-sep-2012'; ...
                MODIS_scf_calib = MODIS_scf_2014 ;  MODIS_scf_analysis_1 = MODIS_scf_2011 ; MODIS_scf_analysis_2 = MODIS_scf_2012 ;MODIS_scf_analysis_3 = MODIS_scf_2013 ;...
                MC_sim_scf_calib = MC_sim_scf_2014 ;  MC_sim_scf_analysis_1 = MC_sim_scf_2011 ; MC_sim_scf_analysis_2 = MC_sim_scf_2012 ;MC_sim_scf_analysis_3 = MC_sim_scf_2013 ;...
                MC_scf_stat_calib = MC_scf_stat_2014 ;  MC_scf_stat_analysis_1 = MC_scf_stat_2011 ; MC_scf_stat_analysis_2 = MC_scf_stat_2012 ;MC_scf_stat_analysis_3 = MC_scf_stat_2013 ;
%                 p_LOA = 30; %25; %30; % using residual based GLUE 5-95% as ref = 30
        end

        %     summary_LHandCR = strcat('summary_LHandCR_',num2str(i));

        %-----------------------------------------------------------------------------------------------
        % the simulated value nc files are named as 2011, 2012, 2013 and 2014
        % instead of starting from 2010. Thus following lines are meant to conform
        % the nc file names with the actural simulation start dates specified above. 
        if strcmp(start_date_calib(8:end),'2010');  calib_year = '2011'; elseif strcmp(start_date_calib(8:end),'2011');  calib_year = '2012'; elseif strcmp(start_date_calib(8:end),'2012');  calib_year = '2013'; elseif strcmp(start_date_calib(8:end),'2013');  calib_year = '2014'; end
        if strcmp(start_date_1(8:end),'2010');  analysis_year_1 = '2011'; elseif strcmp(start_date_1(8:end),'2011');  analysis_year_1 = '2012'; elseif strcmp(start_date_1(8:end),'2012');  analysis_year_1 = '2013'; elseif strcmp(start_date_1(8:end),'2013');  analysis_year_1 = '2014'; end
        if strcmp(start_date_2(8:end),'2010');  analysis_year_2 = '2011'; elseif strcmp(start_date_2(8:end),'2011');  analysis_year_2 = '2012'; elseif strcmp(start_date_2(8:end),'2012');  analysis_year_2 = '2013'; elseif strcmp(start_date_2(8:end),'2013');  analysis_year_2 = '2014'; end
        if strcmp(start_date_3(8:end),'2010');  analysis_year_3 = '2011'; elseif strcmp(start_date_3(8:end),'2011');  analysis_year_3 = '2012'; elseif strcmp(start_date_3(8:end),'2012');  analysis_year_3 = '2013'; elseif strcmp(start_date_3(8:end),'2013');  analysis_year_3 = '2014'; end

        % GLUE analysis part
        time_obs_Q = ncread('shyft/discharge_hourly_f.nc','time');
        obs_Q_data = ncread('shyft/discharge_hourly_f.nc','discharge'); 

        %-----------------------------------------------------------------------------------------------
        % check and fix Obs_Q values suddenly dropping to zero
        suspecteddata = [datetime('27-Mar-2011 01:00:00'),datetime('25-Mar-2012 01:00:00'),datetime('31-Mar-2013 01:00:00')]
        obs_Q_data =obs_Q_data(1,:);%obs_Q_data(1,:)= nea, obs_Q_data(2,:)= tya
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

        % convert hourly observed data to daily data
        obs_Q_data = mean(reshape(obs_Q_data(:,:),24,[]));  
        time_obs_Q = min(reshape(time_obs_Q,24,[]));

        unix_start_date = posixtime(start_date);
        index = find(time_obs_Q==unix_start_date);

        obs_Q =obs_Q_data(index:index+analysis_period-1);

        MC_sim_Q = ncread(strcat('MC_sim_Q_era_',calib_year,'_99999.nc'),'sim_Q'); % from hycamp\....\MC
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
        LMwt_Q = sum(wtg(TS,:)); 

        % Call GLUE-LOA function
        [idx_Q_beforeRelaxing, id_LOA, Pcnts_Q_beforeRelaxing] = GLUE_LOA(N,LMwt_Q,Q_LOA(TS,:),Qsim,TS,Pc);

        no_LOA_t = p_LOA*length(TS)/100; 
        idx_relaxed_LOA = sum(id_LOA(:,:,1))>=no_LOA_t; 
        idx_Q_afterRelaxing = idx_relaxed_LOA'; 
        % just to check how many behavioral models are found satisfying the stated criteria
        no_idx_complying = sum(idx_Q_afterRelaxing)

        idx = idx_Q_afterRelaxing;%temporary
        
        
        if no_idx_complying>5
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
        else
            CR=nan; %-9999;
            nse_median=nan; %-9999;
        end

        p_LOA_CR_behModelNO_NSE(m,counter,1)=p_LOA;
        p_LOA_CR_behModelNO_NSE(m,counter,2)=CR;
        p_LOA_CR_behModelNO_NSE(m,counter,3)=no_idx_complying;
        p_LOA_CR_behModelNO_NSE(m,counter,4)=nse_median;
        
        counter=counter+1;

        
    end
    
    
end   
    
%--------------------------------------------------------------------------------------------------------------------------

% replace -9999 with Nan or 0 as specified
for yr=1:16      
    for pp=1:4
        for dd=1:64
            dd_v=p_LOA_CR_behModelNO_NSE(yr,dd,pp);
            if dd_v==-9999
                p_LOA_CR_behModelNO_NSE(yr,dd,pp)=nan; % 0; %nan;
            end
        end
    
    end
end

% extract values for each year of simulation
c_m_1=p_LOA_CR_behModelNO_NSE(1,1:16,:);
c_m_2=p_LOA_CR_behModelNO_NSE(2,17:32,:);
c_m_3=p_LOA_CR_behModelNO_NSE(3,33:48,:);
c_m_4=p_LOA_CR_behModelNO_NSE(4,49:64,:);



%----------------------------------------------------------------------------------------------------
% all four years

 figure              
 [ax, p11, p12]=plotyy(c_m_1(:,:,1),c_m_1(:,:,2),c_m_1(:,:,1),c_m_1(:,:,4));
 p11.Marker = '+'; p11.LineStyle = '--';p11.Color = rgb('MidnightBlue'); p11.LineWidth = 1.5;
 p12.Marker = '+';p12.LineStyle = '-'; p12.Color = rgb('Magenta'); p12.LineWidth = 1.5;
 hold(ax(1));
 p12=plot(ax(1),c_m_2(:,:,1),c_m_2(:,:,2),'Marker','o','LineStyle','--','Color',rgb('MidnightBlue'),'LineWidth',1.5);
 p13=plot(ax(1),c_m_3(:,:,1),c_m_3(:,:,2),'Marker','x','LineStyle','--','Color',rgb('MidnightBlue'),'LineWidth',1.5);
 p14=plot(ax(1),c_m_4(:,:,1),c_m_4(:,:,2),'LineStyle','--','Color',rgb('MidnightBlue'),'LineWidth',1.5);
 hold(ax(2));
 p22=plot(ax(2),c_m_2(:,:,1),c_m_2(:,:,4),'Marker','o','LineStyle','-','Color',rgb('Magenta'),'LineWidth',1.5);
 p23=plot(ax(2),c_m_3(:,:,1),c_m_3(:,:,4),'Marker','x','LineStyle','-','Color',rgb('Magenta'),'LineWidth',1.5);
 p24=plot(ax(2),c_m_4(:,:,1),c_m_4(:,:,4),'LineStyle','-','Color',rgb('Magenta'),'LineWidth',1.5);
 
 legend([p11;p12;p13;p14],'2011','2012','2013','2014')
 set(ax,{'ycolor'},{rgb('MidnightBlue');rgb('Magenta')}) 
set(ax(1),'ytick',0:0.1:1);
set(ax(2),'ytick',0:0.1:1);
 ylabel(ax(1),'CR')
 ylabel(ax(2),'NSE')
 xlabel('pLOA (%)')
 
 %----------------------------------------------------------------------------------------------------
% only years 2011 and 2012 (hydro year)
 figure              
 [ax, p11, p12]=plotyy(c_m_1(:,:,1),c_m_1(:,:,2),c_m_1(:,:,1),c_m_1(:,:,4));
 p11.Marker = '+'; p11.LineStyle = '--';p11.Color = rgb('MidnightBlue'); p11.LineWidth = 1.5;
 p12.Marker = '+';p12.LineStyle = '-'; p12.Color = rgb('Magenta'); p12.LineWidth = 1.5;
 hold(ax(1));
 p21=plot(ax(1),c_m_2(:,:,1),c_m_2(:,:,2),'Marker','o','LineStyle','--','Color',rgb('MidnightBlue'),'LineWidth',1.5);
 hold(ax(2));
 p22=plot(ax(2),c_m_2(:,:,1),c_m_2(:,:,4),'Marker','o','LineStyle','-','Color',rgb('Magenta'),'LineWidth',1.5);
 
 legend([p11;p12;p21;p22],'CR 2011','NSE 2011','CR 2012','NSE 2012')
 set(ax,{'ycolor'},{rgb('MidnightBlue');rgb('Magenta')}) 
 set(ax(1),'ytick',0:0.5:1);
 set(ax(2),'ytick',0:0.5:1);
 xlim(ax(1),[10 70]);% exclude the first pLOA (5%) result 
 xlim(ax(2),[10 70]);% exclude the first pLOA (5%) result 
 ylim(ax(1),[0.0 1.0]);%
 ylim(ax(2),[0.0 1.0]);%

 ylabel(ax(1),'CR')
 ylabel(ax(2),'NSE')
 xlabel('pLOA (%)')
 
 legend boxoff
 
 set(gcf, 'Color', 'w');% set figure background to white               
 export_fig 'CR_and_NSE_vs_pLOA_3' -tiff -painters -r300;
 export_fig(gcf,'CR_and_NSE_vs_pLOA_3.fig', '-r300');
