m1Eqm               = m1EqmFunc(rho,s,m,piH,piL,gamma) ;

m0Eqm               = m-m1Eqm;
mul1Eqm             = (gamma*piL*piH*m1Eqm)/(rho*m1Eqm*m0Eqm + gamma*piH*m1Eqm + gamma*piL*m0Eqm) ;
muh0Eqm             = (gamma*piL*piH*m0Eqm)/(rho*m1Eqm*m0Eqm + gamma*piH*m1Eqm + gamma*piL*m0Eqm) ;
mul0Eqm             = piL - mul1Eqm;
muh1Eqm             = piH- muh0Eqm;

%%% Calculate DD and CD/DC volumes, and DD share of total sales

VolDD               = rho*muh0Eqm*m1Eqm*( (1+rho*mul1Eqm./(lambda*m1Eqm/m)).*log(1+(lambda.*m1Eqm/m)./(rho.*mul1Eqm))-1) ;
VolCD               = rho*mul1Eqm*m0Eqm;
DealerShare         = VolDD/(VolDD + 2*VolCD);

DD_SellingShare     = lambda*dPhi1(xMat,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm).*(m0Eqm-Phi0(xMat,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm))./VolDD ;

CD_Share            = dPhi0(xMat,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm)/m0Eqm ;
DC_Share            = dPhi1(xMat,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm)/m1Eqm ;

