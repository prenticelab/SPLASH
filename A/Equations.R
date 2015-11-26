#Net.rad is Net Radiation
#A.H.C is Air specific heat capacity
#B.L.C is Boundary layer conductance which is proportional to sqrt(u/d), with u being wind speed and d being leaf width

Temp.change <-( Net.rad - lambda.E)/(A.H.C*B.L.C)


#alpha is the Priestley-Taylor coefficient 
#s is the slope of the Clausius-Clapeyron relationship 
#gamma is the psychrometer constant
#vapour.water is vapourisation of water
lambda.E<-vapour.water*Evapotranspiration = Net.rad*alpha*(s/s+gamma) 

T.change <- (Net.rad - Net.rad*alpha*s/(s+gamma))/(A.H.C*B.L.C)
T.change <- Net.rad(1- alpha*s/(s+gamma))/A.H.C*B.L.C

water.vapourisation<-2.26
# latent heat of vapourisation

MW.Ratio <- 0.622
# MW.Ratio = molecular weight ratio of H2O/dry air

#air.specific.heat depends on atmospheric.pressure
gamma <- atmospheric.pressure*air.specific.heat/water.vapourisation*MW.Ratio


