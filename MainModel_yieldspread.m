% Inter dealer price and expected inter dealer price
inter_dealer        = @(x,y) (theta0*DVfunc(x) + theta1*DVfunc(y))...
                                .*(x<=y).*dPhi1(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm).*dPhi0(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm);
scal_inter_dealer   = @(x,y) (x<=y).*dPhi1(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm).*dPhi0(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm);
    
expected_inter      = integral2( @(x,y) inter_dealer(x,y),0,1,0,1,'AbsTol',TolX)...
                         /integral2( @(x,y) scal_inter_dealer(x,y),0,1,0,1,'AbsTol',TolX) ;

yield_spread     = ( yH/expected_inter - r) ;
                     