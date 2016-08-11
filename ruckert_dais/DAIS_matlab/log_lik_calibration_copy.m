%------------------------------------------------
% -file = DAIS_IceFlux_model.R
% -Antarctic Ice Sheet (AIS) model
%------------------------------------------------
% -Function of the Shaffer DAIS model 2014
% -This function/DAIS model estimates the sea-level equivalence of Antarctic ice sheet melt
% -Model can be found in Shaffer_GMDD_2014
% -Matlab codes and forcings can be found at:
% -www.dcess.dk under "DCESS models"
%
% -Author: Kelsey Ruckert (klr324@psu.edu)
%------------------------------------------------
% -June 17, 2014 #Updates June 10 2015 # Coded from R to Matlab March 7th 2016
%------------------------------------------------

function llik = iceflux(theta)

global data;

% theta = parameters
% data = forcings to run the model

Gamma = theta(1) ;           %sensitivity of ice flow to sea level
alpha = theta(2) ;           %sensitivity of ice flow to ocean subsurface temperature
mu = theta(3) ;              %Profile parameter related to ice stress [m^(1/2)]
nu = theta(4) ;              %Relates balance gradient to precipitation [m^(-1/2) yr^(-1/2)]
P0 = theta(5);               %Precipitation at 0C [m of ice/yr]
kappa = theta(6) ;           %Relates precipitation to temperature [K^-1]
f0 = theta(7) ;              %Constant of proportionality for ice speed [m/yr]
h0 = theta(8) ;              %Initial value for runoff line calculation [m]
c = theta(9) ;              %Second value for runoff line calculation [m]
b0 =  theta(10) ;            %Height of bed at the center of the continent [m]
slope =  theta(11) ;             %Slope of the bed

var_y = theta(12);         % variance

Ta = data(:,1) ;            %Air temperature
Toc = data(:,2) ;            %Ocean temperature
GSL = data(:,3) ;           %Global sea-level rate of change
SL = data(:,4) ;            %Global sea-level rise

%standards                      %standards = [Tf, rho_w, rho_i, rho_m, Toc_0, Rad0, Volo]
Tf = -1.8;                     % Freezing temperature of sea water
rho_w = 1030;                  % Density of sea water [g/cm^3]
rho_i = 917;                   % Density of ice water [g/cm^3]
rho_m = 4000;                  % Density of rock [g/cm^3]
del  = rho_w/rho_i;            % Ratio sea water and ice density [-]
eps1 = rho_i/(rho_m - rho_i);  % Ratio ice density and density difference between rock and ice [-]
eps2 = rho_w/(rho_m - rho_i);  % Ratio sea water density and density difference between rock and ice [-]
Toc_0 = 0.72;                  % Present day high latitude ocean subsurface temperature [K]
Volo = 2.4789e16;              % Steady state AIS volume for present day Ta and SL [m^3]
Rad0 = 1.8636e6;               % Steady state AIS radius for present day Ta and SL [m]
R = Rad0;                      % Initial condition for integration

%Simple time stepping integration with one year time step; solves for ice sheet
% radius R and from R calculates ice sheet volume Vol and sea level equivalent SLE
Rad(1:length(Ta))=NaN;            %Radius of ice sheet
Vol(1:length(Ta))=NaN;            %Ice volume
SLE(1:length(Ta))=NaN;            %Sea-level equivalent [m]

% Run model
for i=1:length(Ta)
    % function used in calculating Ice speed and ice flux (modified from equation 11)
    f = f0*((1-alpha)+alpha*((Toc(i)-Tf)/(Toc_0-Tf))^2)/((slope*Rad0-b0)^(Gamma-1));
    hr= h0+c*Ta(i); % equation 5
    rc = (b0-SL(i))/slope; % application of equation 1 (paragraph after eq3)
    P = P0*exp(kappa*Ta(i)); % equation 6
    beta = nu*P^0.5; % equation 7 (corrected with respect to text)
    rR = R-(abs(hr-b0+slope*R)*(hr-b0+slope*R))/mu; % Distance from the continent center to where the runoff line intersects the ice sheet surface.
    
    if R<=rR && R<=rc
        % Total mass accumulation on ice sheet (equation 8)
        Btot = pi*P*R*R; 
        % In case there is no marine ice sheet / grounding line
        F = 0;    % no ice flux
        ISO = 0;  % (third term equation 14) NAME?
        fac = pi*(1+eps1)*(4/3*mu^0.5*R^1.5-slope*R*R);  % ratio dV/dR
    
    elseif R>rR && R<=rc
        % Total mass accumulation on ice sheet minus runoff
        Btot = pi*P*R*R-pi*beta*(hr-b0+slope*R)*(R*R-rR*rR)-...
            (4*pi*beta*mu^0.5)/5*(R-rR)^2.5+...
            (4*pi*beta*mu^0.5)/3*R*(R-rR)^1.5;
        % In case there is no marine ice sheet / grounding line
        F = 0;    % no ice flux
        ISO = 0;  % (third term equation 14) NAME?
        fac = pi*(1+eps1)*(4/3*mu^0.5*R^1.5-slope*R*R);  % ratio dV/dR
    
    elseif R<=rR && R>=rc
        % Total mass accumulation with marine ice sheet / grounding line
        Btot = pi*P*R*R;
        Hw = slope * R - b0 + SL(i);  % (equation 10)
        F = 2*pi*R*f*del*Hw^(Gamma+1); % Ice flux (equation 9)
        ISO = 2*pi*eps2*(slope*rc*rc-b0/slope*rc)*GSL(i); % (third term equation 14) NAME?
        fac = pi*(1+eps1)*(4/3*mu^0.5*R^1.5-slope*R*R)-2*pi*eps2*(slope*R*R-b0*R);
    
    else
        % Total mass accumulation minus runoff with marine ice sheet / grounding line
        Btot = pi*P*R*R-pi*beta*(hr-b0+slope*R)*(R*R-rR*rR)-...
            (4*pi*beta*mu^0.5)/5*((R-rR)^2.5)+...
            (4*pi*beta*mu^0.5)/3*(R*(R-rR)^1.5);
        Hw = slope * R - b0 + SL(i);  % (equation 10)
        F = 2*pi*R*f*del*Hw^(Gamma+1); % Ice flux (equation 9)
        ISO = 2*pi*eps2*(slope*rc*rc-b0/slope*rc)*GSL(i); % (third term equation 14) NAME?
        fac = pi*(1+eps1)*(4/3*mu^0.5*R^1.5-slope*R*R)-2*pi*eps2*(slope*R*R-b0*R);
    end
    
    dR = (Btot-F+ISO)/fac;
    R = R+dR; % Estimate new radius
    V = 8/15*pi*mu^0.5*R^2.5-1/3*pi*slope*R^3; 
    
    % Calculate sea volume
    Vsea = pi*(2/3*slope*(R^3-rc^3)-b0*(R^2-rc^2));
            
    % Calulate the volume change over time
    if R<=rc
        Volt = (1+eps1)*V;
    else
        Volt = (1+eps1)*V-eps2*Vsea;
    end
    
    % Ice sheet volume (equation 13)
    Rad(i) = R;
    Vol(i) = Volt;
    SLE(i) = 57*(1-Vol(i)/Volo); % Takes steady state present day volume to correspond to 57m SLE
end

% Estimate the residuals
resid_y(1) = 3.9 - (SLE(120000) - mean(SLE(239961:239990)));
resid_y(2) = -11.35 - (SLE(220000) - mean(SLE(239961:239990)));
resid_y(3) = -2.625 - (SLE(234000) - mean(SLE(239961:239990)));
resid_y(4) = 0.002 - (SLE(240002) - mean(SLE(239961:239990)));

sterr_y = [2.1, 4.45, 1.375, 0.0009];

% Caluculate the log likelihood assuming heteroskedastic observation errors
% and independent and identically distributed residuals.
llik = sum(log(normpdf(resid_y, repmat(0,1,4), sqrt(var_y + sterr_y.^2))));






