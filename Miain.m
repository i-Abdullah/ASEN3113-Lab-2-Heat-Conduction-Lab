%% housekeeping

clear
clc
close all


%% read data

% add data files path
addpath('./Data');

% read data
Steel = importdata('Steel_20V_185mA.txt'); % the steel sepcimen data
Brass = importdata('Brass_30V_285mA.txt'); % the brass sepcimen data
Alum = importdata('Aluminum_25V_240mA.txt'); % the aluminum specimen data

% voltage
V_Steel = 20;
V_Brass = 30;
V_Alum = 25;

% current
A_Steel = 185;
A_Brass = 285;
A_Alum = 240;


% parse the data:

time_Steel = Steel.data(:,1); % time
TC1_Steel = Steel.data(:,2); % temp reading from first thermo couple
TC2_Steel = Steel.data(:,3); % temp reading from 2nd thermo couple
TC3_Steel = Steel.data(:,4); % temp reading from 3rd thermo couple
TC4_Steel = Steel.data(:,5); 
TC5_Steel = Steel.data(:,6); 
TC6_Steel = Steel.data(:,7); 
TC7_Steel = Steel.data(:,8); 
TC8_Steel = Steel.data(:,9);



time_Brass = Brass.data(:,1); % time
TC1_Brass = Brass.data(:,2); % temp reading from first thermo couple
TC2_Brass = Brass.data(:,3); % temp reading from 2nd thermo couple
TC3_Brass = Brass.data(:,4); % temp reading from 3rd thermo couple
TC4_Brass = Brass.data(:,5); 
TC5_Brass = Brass.data(:,6); 
TC6_Brass = Brass.data(:,7); 
TC7_Brass = Brass.data(:,8); 
TC8_Brass = Brass.data(:,9);


time_Alum = Alum.data(:,1); % time
TC1_Alum = Alum.data(:,2); % temp reading from first thermo couple
TC2_Alum = Alum.data(:,3); % temp reading from 2nd thermo couple
TC3_Alum = Alum.data(:,4); % temp reading from 3rd thermo couple
TC4_Alum = Alum.data(:,5); 
TC5_Alum = Alum.data(:,6); 
TC6_Alum = Alum.data(:,7); 
TC7_Alum = Alum.data(:,8); 
TC8_Alum = Alum.data(:,9);

%%
