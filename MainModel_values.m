%%% Beta inverse function for delta(x)
beta_delta_pdf     = @(delta_normalized) pdf('beta',delta_normalized,beta_delta_a,beta_delta_b) ;
beta_delta_icdf    = @(x) icdf('beta',x,beta_delta_a,beta_delta_b) ;
deltaFunc          = @(x) deltaLow + (deltaHigh-deltaLow)*beta_delta_icdf(x) ;
deltaFunc_prime    = @(x) (deltaHigh-deltaLow)./beta_delta_pdf(beta_delta_icdf(x)) ;

%%% Beta inverse function for epsi(x)
beta_epsi_pdf      = @(epsi_normalized) pdf('beta',epsi_normalized,beta_epsi_a,beta_epsi_b) ;
beta_epsi_icdf     = @(x) icdf('beta',x,beta_epsi_a,beta_epsi_b) ;
epsiFunc           = @(x) epsiLow + (epsiHigh-epsiLow)*beta_epsi_icdf(x) ;
epsiFunc_prime     = @(x) (epsiHigh-epsiLow)./beta_epsi_pdf(beta_epsi_icdf(x)) ;

%%% Linear system of equation for DWH, DWL and DVlowerbar

A                  = zeros(3,3); % Matrix of coefficients for DWH_bar, DWL, and DV(0), respectively
% Equation for DWH_average
A(1,1)             = r + gamma*piL + rho*(1-theta)*m1Eqm ;
A(1,2)             = -gamma*piL ;
A(1,3)             = - rho*(1-theta)*m1Eqm ;
% Equation for DWL
A(2,1)             = -gamma*piH ;
A(2,2)             = r + gamma*piH + rho*(1-theta)*m0Eqm ;
A(2,3)             = - rho*(1-theta)*m0Eqm ;
% Equation for DVlowerbar
A(3,1)             = - rho*theta*muh0Eqm ;
A(3,2)             =  -rho*theta*mul1Eqm ;
A(3,3)             = r + rho*theta*muh0Eqm + rho*theta*mul1Eqm ;

C                  = zeros(3,1) ; %% Vector of right-hand side coefficients

sigma              = @(x) (deltaFunc_prime(x) + rho*muh0Eqm*theta*epsiFunc_prime(x))...
    ./(r + rho*theta*muh0Eqm + rho*theta*mul1Eqm + lambda*theta1*(m0Eqm-Phi0(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm))/m + lambda*theta0*Phi1(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm)/m) ;
% Constant for DVH
C(1,1)             = yH + rho*(1-theta)*integral( @(x) (m1Eqm-Phi1(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm)).*(sigma(x)-epsiFunc_prime(x)),LowX,1,'AbsTol',TolX);
% Constant for DVL
C(2,1)             = yL + rho*(1-theta)*integral( @(x) (m0Eqm-Phi0(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm)).*sigma(x),LowX,1,'AbsTol',TolX);

% Constant for DVlow
C(3,1)             = deltaFunc(0)+ lambda/m*theta1*integral( @(x) (m0Eqm-Phi0(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm) ).*sigma(x),LowX,1,'AbsTol',TolX) ;




Delta = zeros(3,1) ; %% Solve for Delta
Delta = A\C ;

DWH   = Delta(1) ;
DWL   = Delta(2) ;

DVmat = zeros(1,length(gridsize)) ;

for n = 1:gridsize ;
    DVmat(n) = Delta(3,1) + integral( @(x) sigma(x),LowX,xMat(n),'AbsTol',TolX) ;
end;

DVlow = DVmat(1) ;
DVhigh = DVmat(end) ;

% Interpolation
DVfunc = @(x) interp1(xMat,DVmat,x,'linear','extrap') ;


