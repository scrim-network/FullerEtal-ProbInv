# file = DaisDE.R
#Function of the Shaffer DAIS model 2014
# This function/DAIS model estimates the sea-level equivalence of Antarctic ice sheet melt
modelfn = function(parameters, forcings, standards, windows, obs.years){
  
  gamma = parameters[1] #sensitivity of ice flow to sea level
  alpha = parameters[2] #sensitivity of ice flow to ocean subsurface temperature
  mu = parameters[3]    #Profile parameter related to ice stress [m^(1/2)]
  eta = parameters[4]   #Relates balance gradient to precipitation [m^(-1/2) yr^(-1/2)]
  Po = parameters[5]    #Precipitation at 0C [m of ice/yr]
  kappa = parameters[6] #Relates precipitation to temperature [K^-1]
  fo = parameters[7]    #Constant of proportionality for ice speed [m/yr]
  ho = parameters[8]    #Initial value for runoff line calculation [m]
  co = parameters[9]    #Second value for runoff line calculation [m]
  bo =  parameters[10]  #Height of bed at the center of the continent [m]
  s =  parameters[11]   #Slope of the bed

  
  TA = forcings[,1]
  TO = forcings[,2]
  GSL = forcings[,3]
  SL = forcings[,4]
  
  Rad = rep(NA,end)
  Vol = rep(NA,end)
  SLE = rep(NA,end)
  
  #Setup the model
  for(i in 1:enddate){ 
    #  print(i)
    #Calculate total ice flux across the grounding line
    f = fo*((1-alpha)+alpha*((TO[i]-Tice)/(TOo-Tice))^2)/((s*Roa-bo)^(gamma-1))
    hr= ho+co*TA[i] #Calculate ice sheet surface height
    #Calculate the distance from the continent center to where the ice sheet enters the sea
    rcon = (bo-SL[i])/s
    P = Po*exp(kappa*TA[i]) #Calculate precipitation
    beta = eta*P^1/2 #calculate the mass balance gradient
    #Calculate the distance from the continent center to where the runoff line intersects the ice sheet surface
    rR = R-(abs(hr-bo+s*R)*(hr-bo+s*R))/mu
    
    #Apply whether or not the ice sheet has components of a marine ice sheet
    if (R<=rR && R<=rcon) {
      Btot = pi*P*R*R #Total mass accumulation on ice sheet
      F = 0
      ISO = 0
      fac = pi*(1+eps1)*(4/3*mu^0.5*R^1.5-s*R*R)
    } else if (R>rR && R<=rcon) {
      Btot = pi*P*R*R-pi*beta*(hr-bo+s*R)*
        (R*R-rR*rR)-(4*pi*beta*mu^0.5)/5*(R-rR)^2.5+
        (4*pi*beta*mu^0.5)/3*R*(R-rR)^1.5
      F = 0
      ISO = 0
      fac = pi*(1+eps1)*(4/3*mu^0.5*R^1.5-s*R*R)
    } else if (R<=rR && R>=rcon) {
      Btot = pi*P*R*R
      F = 2*pi*R*f*del*(s*R-bo+SL[i])^(gamma+1)
      ISO = 2*pi*eps2*(s*rcon*rcon-bo/s*rcon)*GSL[i]
      fac = pi*(1+eps1)*(4/3*mu^0.5*R^1.5-s*R*R)-2*pi*eps2*(s*R*R-bo*R)
    } else {
      Btot = pi*P*R*R-pi*beta*(hr-bo+s*R)*
        (R*R-rR*rR)-(4*pi*beta*mu^0.5)/5*
        ((R-rR)^2.5)+(4*pi*beta*mu^0.5)/3*(R*(R-rR)^1.5)
      F = 2*pi*R*f*del*(s*R-bo+SL[i])^(gamma+1)
      ISO = 2*pi*eps2*(s*rcon*rcon-bo/s*rcon)*GSL[i]
      fac = pi*(1+eps1)*(4/3*mu^0.5*R^1.5-s*R*R)-2*pi*eps2*(s*R*R-bo*R)
    }
    dR = (Btot-F+ISO)/fac
    R = R+dR #calculate radius
    V = 8/15*pi*mu^0.5*R^2.5-1/3*pi*s*R^3 #Calculate volume
    #calculate sea volume
    Vsea = pi*(2/3*s*(R^3-rcon^3)-bo*(R^2-rcon^2))
    #calulate the volume change over time
    if(R<=rcon){
      Volt = (1+eps1)*V
    }
    else{
      Volt = (1+eps1)*V-eps2*Vsea
    }
    Rad[i] = R
    Vol[i] = Volt
    SLE[i] = 57*(1-Volt/Volo) #Takes steady state present day volume to correspond to 57m SLE
  }
  #Set up equation to find the residuals  
  wresid <- rep(NA,length(obs.years)) #Create a vector of the weighted residuals
  for(i in 1:length(obs.years)){
    wresid[i] <- (median(windows[i,])-(SLE[obs.years[i]]-
                                         mean(SLE[SL.1961_1990]))) #/sd(windows[i,])
  }
  #Find the weighted root means square error
  wrmse = sqrt(mean(wresid^2))
  #Return the weighted root mean square error
  return(wrmse)
}
