%%%%%%%%%%%%%%%%%%%%%%%%
%%% Data Moments     %%%
%%%%%%%%%%%%%%%%%%%%%%%%

YieldSpread_ABX     = 1.4; %target yield spread, from Ang Bhansali and Xing


%%%% SOME CHAIN LENGTH STATISTICS FROM LS
K_LS                        = 7 ;
SizeOfMarkup_LS             = [1.85 1.94 2.26 2.92 3.26 3.57 3.71]' ;

breakdown_markup_LS         = [1.85 0.00 0.00 0.00 0.00 0.00 0.00 ; ...
                       0.84 1.10 0.00 0.00 0.00 0.00 0.00 ; ...
                       0.66 0.52 1.08 0.00 0.00 0.00 0.00 ; ...
                       0.64 0.60 0.55 1.13 0.00 0.00 0.00 ; ...
                       0.63 0.30 0.82 0.40 1.11 0.00 0.00 ; ...
                       0.60 0.27 0.45 0.85 0.27 1.14 0.00 ; ...
                       0.62 0.23 0.43 0.53 0.46 0.29 1.15] ;
                   
DistributionOfMarkupMat_LS  = breakdown_markup_LS./(SizeOfMarkup_LS*ones(1,7)) ;

N_obs_per_chain_length      = [8808119;...
                               1511196;...
                               866459;...
                               173579;...
                               37229;
                               6866;
                               822] ;
                           
freq_chain_length_LS        = N_obs_per_chain_length'/sum(N_obs_per_chain_length) ;

vec_chain_length            = (1:7)' ;

avg_chain_length_LS         = freq_chain_length_LS*vec_chain_length ;

avg_chain_length_sq_LS      = freq_chain_length_LS*(vec_chain_length.^2) ;

var_chain_length_LS         = avg_chain_length_sq_LS - avg_chain_length_LS^2 ;

cov_markup_chain_length_LS  = freq_chain_length_LS*( SizeOfMarkup_LS.*(vec_chain_length-avg_chain_length_LS) ) ;

beta_markup_chain_length_LS = cov_markup_chain_length_LS/var_chain_length_LS ;

avg_markup_LS               = freq_chain_length_LS*SizeOfMarkup_LS ;





%%% Average Inventory Duration (year)
D                           = 3.3/250;
%%% Intermediation chain parameters (lambda M0/M)/(rho muh0)
avg_length                  = 1.3466 ;
Chi                         = fsolve( @(x) (1+1./x).*log(1+x) - avg_length,1) ;



%%% Fraction of assets in the dealer sector (M1*Q/A)
BD_share                    = 0.01 ;
%%% CD Turnover 
CD_turnover                 = 0.411 ;
%%% Muni Market Cap (dollars)
Agg_Muni_Supply             = 2308598605189;
%%% Retail investor base
N_retail                    = 54187500;

%%% Average trade size for DD transaction (dollars)
%Trade_size                  = 110000;  %avg DD, non-seasoned according to LS
%Trade_size                  = 40000;  %avg DD, non-seasoned according to LS
Trade_size                  = 206989 ; %avg DD Green Hollifield and Shueroff
%%% Average time to sell for a customer
Customer_time_to_sell       = 5/250 ;
%%% Number of dealer firms: LS supplied # on page 6
dd                          = 2257-63-117;

%%% Solve for the m0 that will make piH=s

rho_times_muh0              = (2+Chi)/2/D/(1+Chi) ;
M0_over_M_balanced          = (1+rho_times_muh0*Customer_time_to_sell)/(1+2*rho_times_muh0*Customer_time_to_sell) ;
