%%% A parameter set for a model w/o heterogeneity

theta                = 0.971; % DC bargaining power
theta1               = 0.5 ; % DD bargaining power
yL                   = 0.022846; % low utility flow for customer
yH                   = r ; % high utility flow for customer
deltaLow             = 0.11*yH ; % lower bound of utility flow for dealer
deltaHigh            = 0.8*yH; % upper bound of utility flow for dealer
beta_delta_a         = 1 ; % first parameter of the beta function for dealer utility flow
beta_delta_b         = 1 ; % second parameter of the beta function for dealer utility flow
epsiHigh             = 0 ; % upper bound of expected match specific flow for dealer
beta_epsi_a          = 1; % first parameter of the beta function for expected match specific utility flow
beta_epsi_b          = 1 ; % second parameter of the beta function for expected match specific utility flow
M0_over_M            = M0_over_M_balanced; % fraction of asset-less dealer in the dealer sector
theta0               = 1-theta1 ; % DD bargaining power
