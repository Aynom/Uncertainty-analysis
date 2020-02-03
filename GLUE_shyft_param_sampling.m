% this algorithm is used to coordinate parameter sampling.
% Details about setup of the sampling technique such as the number of samples, parameters that
% are allowed to vary, and definition of their minimum and maximum range is specified in this
% script. It prepares the inputs and calls an existing algorithm that samples parameter values
% using the All-At-a-Time (AAT) sampling technique. The implemented algorithm retrieves the
% sampled parameter values and exports the result as a netcdf file.
% This script was adopted from and calls to other scripts of the SAFE Toolbox (Pianosi et al, 2015)

% For details on the conceptual background, implementation, and sample outputs from this algorithm, 
% the reader is referred to the following paper: 

% Teweldebrhan, A. T., Burkhart, J. F., and Schuler, T. V.: Parameter uncertainty analysis for an 
% operational hydrological model using residual-based and limits of acceptability approaches, 
% Hydrology and Earth System Sciences, 22, 5021-5039, 2018.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all
format long g

my_dir = pwd ; 

addpath('util');
addpath('visualization')
addpath('GLUE')

addpath('shyft')
 
N = 100000; % sample size parameterisations
MC_f_name = 'Mc_calib_param_garb_100k.nc';
SampStrategy = 'lhs' ; %'dds' ; % sampling strategy ('lhs' or 'dds')

M  = 21 % number of uncertain parameters 
DistrFun  = 'unif'  ; % Parameter distribution
%DistrFun  = 'norm'  ; % Parameter distribution
% Parameter ranges:
% S_No =   [1         2       3           4       5       6       7       8       9       10      11      12      13      14      15      16      17      18      19      20      21   ];
S_No =   {'1'       '2'       '3'       '4'      '5'    '6'     '7'     '8'     '9'      '10'   '11'    '12'    '13'    '14'    '15'    '16'    '17'    '18'    '19'    '20'    '21'   };

% % ff: used in running MC for nea catchment (only selected parameters)
% % S_No =   [1         2       3           4       5       6       7       8       9       10      11      12      13      14      15      16      17      18      19      20      21   ];
% xmin =   [-5.0,     0.0,    -0.15,      1.5,    -3.0,   1.0,    0.1,    0.1,    1.0,    20.0,   30.0,   0.9,    0.6,    5.0,    0.06,    0.0,    0.0,    0.4,    1.0,    0.2,    1.26 ] ; % minimum values
% xmax =   [1.0,      1.2,    -0.05 ,     1.5,    2.0,    6.0,    0.1,    0.1,    15.0,   40.0,   30.0,   0.9,    0.6,    5.0,    0.85,    0.0,    0.0,    0.4,    1.0,    0.2,    1.26 ] ; % maximum values

% ff: used in running MC for garberg (only selected parameters)
% S_No =   [1         2       3           4       5       6       7       8       9       10      11      12      13      14      15      16      17      18      19      20      21   ];
xmin =   [-5.0,    -0.1,    -0.15,      1.5,    -3.0,   1.0,    0.1,    1.0,    5.0,    20.0,   30.0,   0.9,    0.6,    5.0,    0.06,    0.0,    0.0,    0.4,    1.0,    0.2,    1.26 ] ; % minimum values
xmax =   [1.0,      1.2,    -0.05 ,     1.5,    2.0,    6.0,    0.1,    1.0,    15.0,   40.0,   30.0,   0.9,    0.6,    5.0,    0.85,    0.0,    0.0,    0.4,    1.0,    0.2,    1.26 ] ; % maximum values

variable_calib_param_index =[1,2,3,5,6,9,10,15];% parameters whose value to be identified

%%%% SAMPLING INPUT SPACE
% Sample the input space using the 'AAT_sampling' function
addpath('sampling')
% Create data structure for parameter ranges
% as required by AAT_sampling
for i=1:M; DistrPar{i} = [xmin(i) xmax(i)] ; end 
% Perform sampling:
X = AAT_sampling(SampStrategy,M,DistrFun,DistrPar,N);

% replace nan values (parameters with constant value) with xmin value
for i=1:21   
    if isnan(X(:,i))
        X(:,i)=xmin(i);
    end
end 
%write monte carlo calib parameter values to netcdf file 
my_data = X; %[xmin; xmax];
[row col]= size(X)

ncid = netcdf.create(MC_f_name,'NC_WRITE');
dimidrow = netcdf.defDim(ncid,'rows',row);
dimidcol = netcdf.defDim(ncid,'length',col);
varid = netcdf.defVar(ncid,'calib_MC','NC_DOUBLE',[dimidrow dimidcol]);
netcdf.endDef(ncid);
netcdf.putVar(ncid,varid,my_data);
netcdf.close(ncid);
ncid2 = netcdf.open(MC_f_name,'NC_NOWRITE');
data_copy = netcdf.getVar(ncid2,0);
if isequal(my_data,data_copy)
      disp('Data match');
else
      disp('Data mis-match');
end




