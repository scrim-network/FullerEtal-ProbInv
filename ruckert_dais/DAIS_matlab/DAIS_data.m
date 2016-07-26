%------------------------------------------------
% -file = DAIS_data.m
% -Antarctic Ice Sheet (AIS) model
%------------------------------------------------
% -Model can be found in Shaffer_GMDD_2014
% -Matlab codes and forcings can be found at:
% -www.dcess.dk under "DCESS models"
% -Sea-level values are reconstructed from:
% -Waelbroeck et al 2002(21-240kyr BP), Clark et al 2012(7000-21000BP),
% -Lambeck et al 2010(6000BP-1869AD), and Church & White 2011(1870-2010AD)
% -Future Sea level values and and rates are calculated using the Rahmstorf (2007) model, DEoptim, and RCP8.5
% -Future air and ocean temperatures were calculated by Robert Nicholas and Varada Vaidya using the CNRM-CMIP5 model
% -Author: Kelsey Ruckert (klr324@psu.edu)
%------------------------------------------------
% -June 17, 2014 #Updates June 10 2015 # Coded from R to Matlab March 7th 2016
%------------------------------------------------

clear all
close all

% read in forcing data
date = -239999:300;  %240 Kyr BP to 2300AD at one year intervals of the forcings
load future_GSL.txt; %Time rate of change of sea-level
load future_SL.txt;  %Antarctic temperature reduced to sea-level
load future_TA.txt;  %High latitude subsurface ocean temperature
load future_TO.txt;  %Reconstructed sea-level

GSL = future_GSL;
SL = future_SL;
Ta = future_TA;
Toc = future_TO;

clear future_GSL future_SL future_TA future_TO;

%Set model parameters at their standard values
Tf = -1.8;             %Freezing temperature of sea water
rho_w = 1030;              %Density of sea water [g/cm^3]
rho_i = 917;            %Density of ice water [g/cm^3]
rho_m = 4000;             %Density of rock [g/cm^3]
b0 = 775;                %Height of bed at the center of the continent [m]
slope = 0.0006;              %Slope of the bed
f0 = 1.2;                %Constant of proportionality for ice speed [m/yr]
h0 = 1471;               %Initial value for runoff line calculation [m]
c = 95;                 %Second value for runoff line calculation [m]
Rad0 = 1.8636e6;          %Steady state AIS radius for present day Ta and SL [m]
Volo = 2.4789e16;        %Steady state AIS volume for present day Ta and SL [m^3]
Toc_0 = 0.72;              %Present day high latitude ocean subsurface temperature [K]

del = rho_w/rho_i;
eps1 = rho_i/(rho_m - rho_i);
eps2 = rho_w/(rho_m - rho_i);

%Initial condition for integration
R = Rad0;

%Setup AIS melt ranges and specific dates so there is no use of magic numbers
%in the code
last_interglacial = [date(110000),date(120000),date(130000)];
last_glacialmax = [date(219000),date(220000),date(221000)];
holocene = [date(233800),date(234000),date(234200)];
SL_1961_1990 = 239961:239990;
kyrbp_25 = 220000;
kyrbp_6 = 234000;
AD_1880 = 239880;
present = 240002;
endsat = 240000;
enddate = 240010;
