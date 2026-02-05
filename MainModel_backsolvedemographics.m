%%% Backsolving equations 
%%% This routine identify a number of demographic parameters from targets

rho_times_muh0              = 1/D*(2+Chi)/2/(1+Chi) ;
lambda_M0_over_M            = Chi*rho_times_muh0;

s_total                     = Agg_Muni_Supply/Trade_size;
s                           = Agg_Muni_Supply/Trade_size/N_retail ;


%M1                          = BD_share*Agg_Muni_Supply/Trade_size;
M1                          = CD_turnover*Agg_Muni_Supply/Trade_size/rho_times_muh0 ; 
M1_over_M                   = 1 - M0_over_M ;
M                           = M1/M1_over_M ;
M0                          = M-M1 ;

m0                          = M0/N_retail; 
m1                          = M1/N_retail ;
m                           = M/N_retail;

lambda                      = lambda_M0_over_M/M0*M;
rho                         = 1/Customer_time_to_sell/m0 ;%1/delay/(0.5*m0+0.5*m1) ;%


muh0                        = rho_times_muh0/rho;
piH                         = s-m1+muh0*(1-m1/m0);
piL                         = 1-piH ;
muh1                        = piH - muh0;
gamma                       = rho*muh0*m1*m0/(piH*piL*m0 - muh0*(piH*m1+piL*m0)) ;
mul1                        = s - m1 - (piH-muh0);
mul0                        = piL - mul1;
trader_per_dealer           = M/dd ;