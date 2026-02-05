%%% Function to calculate m1 using Cardano's method

m1EqmFunc_balanced = @(rho,s,m,piH,piL,gamma) m/2 ;

a               = @(rho,s,m,piH,piL,gamma) (piH-s) - (m + gamma*piH./rho-gamma*piL./rho) ;
b               = @(rho,s,m,piH,piL,gamma) -gamma*piL*m./rho - (piH-s)*(m + gamma*piH./rho - gamma*piL./rho) - 2*gamma*piH*piL./rho ;
c               = @(rho,s,m,piH,piL,gamma) gamma*piL*m*s./rho ;

p               = @(rho,s,m,piH,piL,gamma)  b(rho,s,m,piH,piL,gamma) - (a(rho,s,m,piH,piL,gamma).^2)/3 ;
q               = @(rho,s,m,piH,piL,gamma)  c(rho,s,m,piH,piL,gamma) - a(rho,s,m,piH,piL,gamma).*b(rho,s,m,piH,piL,gamma)/3  + 2/27*(a(rho,s,m,piH,piL,gamma).^3);
Delta           = @(rho,s,m,piH,piL,gamma) -4*(p(rho,s,m,piH,piL,gamma).^3) - 27*(q(rho,s,m,piH,piL,gamma).^2) ;

j               = exp(2*i*pi/3) ;

u0              = @(rho,s,m,piH,piL,gamma) ( 0.5*(-q(rho,s,m,piH,piL,gamma) + (-Delta(rho,s,m,piH,piL,gamma)./27).^(0.5) ) ).^(1/3) ;
v0              = @(rho,s,m,piH,piL,gamma) ( 0.5*(-q(rho,s,m,piH,piL,gamma) - (-Delta(rho,s,m,piH,piL,gamma)./27).^(0.5) ) ).^(1/3) ;
x0              = @(rho,s,m,piH,piL,gamma) u0(rho,s,m,piH,piL,gamma)+v0(rho,s,m,piH,piL,gamma) - a(rho,s,m,piH,piL,gamma)/3;

u1              = @(rho,s,m,piH,piL,gamma) j*u0(rho,s,m,piH,piL,gamma) ;
v1              = @(rho,s,m,piH,piL,gamma) j^(-1)*v0(rho,s,m,piH,piL,gamma) ;
x1              = @(rho,s,m,piH,piL,gamma) real(u1(rho,s,m,piH,piL,gamma)+v1(rho,s,m,piH,piL,gamma) - a(rho,s,m,piH,piL,gamma)/3);

u2              = @(rho,s,m,piH,piL,gamma) j^2*u0(rho,s,m,piH,piL,gamma) ;
v2              = @(rho,s,m,piH,piL,gamma) j^(-2)*v0(rho,s,m,piH,piL,gamma) ;
x2              = @(rho,s,m,piH,piL,gamma) real(u2(rho,s,m,piH,piL,gamma)+v2(rho,s,m,piH,piL,gamma) - a(rho,s,m,piH,piL,gamma)/3);

m1EqmFunc_unbalanced ...
                = @(rho,s,m,piH,piL,gamma) (Delta(rho,s,m,piH,piL,gamma)<0).*x0(rho,s,m,piH,piL,gamma) ...
                                           + (Delta(rho,s,m,piH,piL,gamma)>0)...
                                           .*( x0(rho,s,m,piH,piL,gamma).*(x0(rho,s,m,piH,piL,gamma)<m).*(x0(rho,s,m,piH,piL,gamma)>0)  ...
                                             + x1(rho,s,m,piH,piL,gamma).*(x1(rho,s,m,piH,piL,gamma)<m).*(x1(rho,s,m,piH,piL,gamma)>0) ... 
                                             + x2(rho,s,m,piH,piL,gamma).*(x2(rho,s,m,piH,piL,gamma)<m).*(x2(rho,s,m,piH,piL,gamma)>0) ) ;

m1EqmFunc = @(rho,s,m,piH,piL,gamma) (abs(piH-s+m/2)<eps).*m1EqmFunc_balanced(rho,s,m,piH,piL,gamma) ...
                                   + (abs(piH-s+m/2)>=eps).*m1EqmFunc_unbalanced(rho,s,m,piH,piL,gamma) ;

%%% Function to calculate the distribution of assets in the dealer sector
APhi                = 1 ;
BPhi                = @(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm)  rho*m/lambda*(mul1Eqm+muh0Eqm) + m0Eqm - m*x ;
CPhi                = @(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm) -rho*m/lambda*mul1Eqm*m*x ;

dBPhi               = @(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm)  - m ;
dCPhi               = @(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm) -rho*m/lambda*mul1Eqm*m ;
DeltaPhi            = @(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm) BPhi(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm).^2 - 4*APhi.*CPhi(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm) ;
Phi1                = @(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm) 1/2/APhi*( - BPhi(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm) + DeltaPhi(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm).^0.5 ) ;
dPhi1               = @(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm) 1/2/APhi*( - dBPhi(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm) ...
                                                                         +  0.5*(DeltaPhi(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm).^-0.5)...
                                                                         .*(2*BPhi(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm)...
                                                                         .*dBPhi(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm)-4*APhi*dCPhi(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm))); 
Phi0                = @(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm) m*x - Phi1(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm) ;
dPhi0               = @(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm) m - dPhi1(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm);

%%% Functions to calculate markup statistics

%%% Calculate the expected reservation value conditional on initial type
lambda1                                         = @(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm) lambda/m*(m0Eqm - Phi0(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm)) ;
lambda0                                         = @(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm) lambda/m*Phi1(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm) ;
Lambda                                          = @(x,y,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm) log( lambda1(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm) + rho*muh0Eqm) - log( lambda1(y,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm) + rho*muh0Eqm) ;
dlambda1                                        = @(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm) -lambda/m*dPhi0(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm) ;



% The distribution of x1 conditional on chain length
ConditionalDistributionOfx1                     = @(x1,k,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm) (-dlambda1(x1,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm))./(lambda1(x1,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm)+rho*muh0Eqm)...
                                                .*k.*Lambda(x1,1,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm).^(k-1)...
                                                ./Lambda(0,1,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm).^k ;
                                            
% The joint distribution of x1 and xj conditional on chain length                                            
ConditionalDistributionOfx1xj                   = @(x1,xj,j,k,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm)...
                                                  (-dlambda1(x1,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm))./(lambda1(x1,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm)+rho*muh0Eqm)...
                                                .*Lambda(x1,xj,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm).^(j-2)./factorial(j-2)...
                                                .*(-dlambda1(xj,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm))./(lambda1(xj,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm)+rho*muh0Eqm)...
                                                .*Lambda(xj,1,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm).^(k-j)./factorial(k-j)...
                                                .*factorial(k)./Lambda(0,1,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm).^k ;
                                            

% THE PENALTY FUNCTION FOR VIOLATING HLW PATTERNS
penalty_func = @(x) 1000*(abs(max(x,0)).^2) ; % create a penalty for violating the HLW patterns 
                             

