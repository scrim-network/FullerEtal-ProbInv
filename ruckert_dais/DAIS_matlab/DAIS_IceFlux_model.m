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

function SLE = iceflux(iceparam,forcings,standards)

%iceflux = function(iceflux, forcings, standards){ #R syntax
%  model.p = length(iceflux)
  
    Gamma = iceparam(1);           %sensitivity of ice flow to sea level
    alpha = iceparam(2);           %sensitivity of ice flow to ocean subsurface temperature
    mu = iceparam(3);              %Profile parameter related to ice stress [m^(1/2)]
    nu = iceparam(4);              %Relates balance gradient to precipitation [m^(-1/2) yr^(-1/2)]
    P0 = iceparam(5);               %Precipitation at 0C [m of ice/yr]
    kappa = iceparam(6);           %Relates precipitation to temperature [K^-1]
    f0 = iceparam(7);              %Constant of proportionality for ice speed [m/yr]
    h0 = iceparam(8);              %Initial value for runoff line calculation [m]
    c = iceparam(9);              %Second value for runoff line calculation [m]
    b0 =  iceparam(10);            %Height of bed at the center of the continent [m]
    slope =  iceparam(11);             %Slope of the bed
  
    Ta = forcings(:,1);            %Air temperature
    Toc = forcings(:,2);            %Ocean temperature
    GSL = forcings(:,3);           %Global sea-level rate of change
    SL = forcings(:,4);            %Global sea-level rise
    
    %standards                      %standards = [Tice, eps1, del, eps2, TOo, Volo, Roa, R]
    Tf = standards(1);
    rho_w = standards(2);
    rho_i = standards(3);
    rho_m = standards(4);
    Toc_0 = standards(5);
    Rad0 = standards(6);
    Volo = standards(7);
    
    % Initialize intermediate parameters.
    del  = rho_w/rho_i;            % Ratio sea water and ice density [-]
    eps1 = rho_i/(rho_m - rho_i);  % Ratio ice density and density difference between rock and ice [-]
    eps2 = rho_w/(rho_m - rho_i);  % Ratio sea water density and density difference between rock and ice [-]
    R  = Rad0;                     % gets updated at end of loop
  
    %Simple time stepping integration with one year time step; solves for ice sheet
    % radius R and from R calculates ice sheet volume Vol and sea level equivalent SLE   
    Rad(1:length(Ta))=NaN;            %Radius of ice sheet
    Vol(1:length(Ta))=NaN;            %Ice volume
    SLE(1:length(Ta))=NaN;            %Sea-level equivalent [m]
    
  % Run model
    for i=1:length(Ta)
    % function used in calculating Ice speed and Ice Flux (modified from equation 11)
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
    
end;


