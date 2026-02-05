if piH >= s & piH*piL<s
    
    epsiHigh_W = epsiHigh*(rho*m1+gamma)/(rho*m1+gamma*piL) ;
    
    % The walrasian price, assuming that the price capitalize the match specific utility of the asset best holder    
    p_walras            = 1/r*yH ;
    % The monopoly ask and bid price
    p_monopoly_high     = 1/r*( r/(r+gamma)*(yH+(r+gamma)*epsiHigh) + gamma/(r+gamma)*(piH*yH+piL*yL)) ;
    p_monopoly_low      = 1/r*( r/(r+gamma)*(yL) + gamma/(r+gamma)*(piH*yH+piL*yL)) ;
 
    % Welfare in the Walrasian equilibrium assuming the asset is allocated to its best holder
    welfare_walras      = s*yH/r + (r+gamma)*epsiHigh_W*piH*piL/r ;
    
    welfare_walras_nomatch = s*yH/r ;
    
    
    % Welfare in autarky. The asset is randomly allocated, and as a result
    % customers derive no match specific utility
    welfare_autarky     = s*(piH*yH+piL*yL)/r ;
    
    % In the search equilibrium, dmu1Eqm is the measure of customers who
    % were matched with type-x dealers, and are still deriving match
    % specific utility (i.e. they have not been hit by a gamma shock yet)
    dmuh1EqmEpsiHigh    =  @(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm) rho*muh0Eqm/gamma*dPhi1(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm) ;
    
    % The function to integrate in order to calculate welfare via sum of
    % flow utilities
    welfare_func        = @(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm) muh1Eqm*yH +dmuh1EqmEpsiHigh(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm).*epsiFunc(x)*(r+gamma)...
                                                                + dPhi1(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm).*deltaFunc(x) ...
                                                                + mul1Eqm*yL ; 
    
    muh1EqmEpsiHigh     = rho*muh0Eqm*m1/gamma ;

    welfare_otc_nochain...
                        =  (muh1Eqm*yH + muh1EqmEpsiHigh*integral(@(x) epsiFunc(x),0,1)*(r+gamma)...
                                                                + m1*integral(@(x) deltaFunc(x),0,1) ...
                                                                + mul1Eqm*yL)/r ; 
    
    match_quality_no_chain ...
                        = muh1EqmEpsiHigh*integral(@(x) epsiFunc(x),0,1)*(r+gamma)/r;
    match_quality_with_chain ... 
                        = integral( @(x) dmuh1EqmEpsiHigh(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm).*epsiFunc(x)*(r+gamma),0,1)/r ;
    match_quality_total = muh1EqmEpsiHigh*epsiHigh_W*(r+gamma)/r;
    
    frac_match_quality_welfare = (match_quality_with_chain-match_quality_no_chain)/(match_quality_total-match_quality_no_chain) ;
    % The integration to calculate welfare
	welfare_otc         = integral(@(x) welfare_func(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm)/r,0,1,'AbsTol',TolX) ;
    
    % The function to integrate in order to calculate the value of a high valuation customer with zero assets
    WH0_func            = @(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm) (r + r*gamma*piL/(r+gamma*piH))^(-1) ... 
                                                                 *rho*(1-theta)*(DWH+epsiFunc(x)-DVfunc(x)).*dPhi1(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm) ;
                                                             
    WH0_func2            = @(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm) (r + r*gamma*piL/(r+gamma*piH))^(-1) ... 
                                                                *rho*(1-theta)...
                                                                *( m1*(DWH - DVlow) ...
                                                                   + (epsiFunc_prime(x) - sigma(x)).*(m1-Phi1(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm)) ) ;
  
    % The value of a high-valuation customer with zero asset                                                                
	WH0                 = integral( @(x) WH0_func(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm),0,1,'AbsTol',TolX) ;
    % The value of a low valuation customer with zero asset
    WL0                 = WH0*gamma*piH/(r+gamma*piH) ;
    
    % The function to integrate in order to calculate the welfare of
    % customers
    welfare_otc_cust_func ...
                        = @(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm) piH*WH0 ...
                                                                + muh1Eqm*DWH ...
                                                                + dmuh1EqmEpsiHigh(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm).*epsiFunc(x) ...
                                                                + piL*WL0 + mul1Eqm*DWL ;
    % The welfare of customers
    welfare_otc_cust    = integral( @(x) welfare_otc_cust_func(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm),0,1,'AbsTol',TolX) ;
    
    % The value of a dealer with x=0 and zero asset
    V00                 = 1/r*rho*mul1Eqm*theta*(DVlow - DWL) ;
    V0prime             = @(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm) (lambda*theta0*Phi1(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm)/m+rho*theta*mul1Eqm)...
                                                                 .*sigma(x)/r ;
    welfare_otc_dealer_func ...
                        = @(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm) m*V00  + m*V0prime(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm).*(1-x) ...
                                                               + m1*DVlow + sigma(x).*(m1-Phi1(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm)) ;
    % The welfare of dealers                                                        
    welfare_otc_dealer  =  integral(@(x) welfare_otc_dealer_func(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm),0,1,'AbsTol',TolX)  ;
    
    
    
    % Distribution of the last dealer in a chain
    distribution_last_dealer_in_chain ...
                        = @(x) -rho*muh0Eqm*dlambda1(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm) ...
                               ./(lambda1(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm)+rho*muh0Eqm).^2 ...
                               .*exp(Lambda(0,1,x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm)) ;
    % Ask and expected ask
    ask                 = @(x) (1-theta)*DVfunc(x) + theta*DWH + theta*epsiFunc(x) ; 
    expected_ask        = integral(@(x) ask(x).*dPhi1(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm)/m1,0,1,'AbsTol',TolX) ;
    
    % Bid and expected bid
    bid                 = @(x) (1-theta)*DVfunc(x) + theta*DWL ;
    expected_bid        = integral(@(x) bid(x).*dPhi0(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm)/m0,0,1,'AbsTol',TolX) ;
    
    % Inter dealer price and expected inter dealer price
    inter_dealer        = @(x,y) (theta0*DVfunc(x) + theta1*DVfunc(y))...
                                .*(x<=y).*dPhi1(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm).*dPhi0(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm);
    scal_inter_dealer   = @(x,y) (x<=y).*dPhi1(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm).*dPhi0(x,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm);
    
    expected_inter      = integral2( @(x,y) inter_dealer(x,y),0,1,0,1,'AbsTol',TolX)...
                         /integral2( @(x,y) scal_inter_dealer(x,y),0,1,0,1,'AbsTol',TolX) ;
                            
    disp(' '); 
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    disp('%% Prices and welfare relative to frictionless economy                  %%');
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    disp(' ') ;
    disp('Found that piH>s') ;

    disp(' ') ;
    disp(['The competitive price without search frictions is         = ',num2str(p_walras)]) ;
    disp(['The monopoly    price without search frictions is         = ',num2str(p_monopoly_high)]) ;
    disp(['The autarchy value of high valuation customer  is         = ',num2str(p_monopoly_high)]) ;
    disp(['The autarchy value of low  valuation customer  is         = ',num2str(p_monopoly_low)]) ;
    disp(['The autarchy value of spread btw customers     is         = ',num2str(10000*(p_monopoly_high/p_monopoly_low-1)),' bps']) ;
    disp(['The res. value of low  valuation customer      is         = ',num2str(DWL)]) ;
    disp(['The res. value of high valuation customer      is         = ',num2str(DWH)]) ;
    disp(['The res. value spread btw        customers     is         = ',num2str(10000*(DWH/DWL-1)),' bps']) ;
    disp([' ']) ;
    disp(['The average ask price with    search frictions is         = ',num2str(expected_ask)]) ;
    disp(['The average int-d price  with search frictions is         = ',num2str(expected_inter)]) ;
    disp(['The average bid price with    search frictions is         = ',num2str(expected_bid)]) ;
    %disp(['Yield spread implied by expected ask           is         = ',num2str(10000*((yH+(r+gamma)*epsiHigh)/expected_ask-r)),' bps']) ;
    
    disp([' ']) ;
    disp(['The following calculation assumes that the Walrasian price correspond to a zero yield spread, i.e., cash flows for the marginal investor equal to yH']) ;  
    disp(['Yield spread implied by expected inter-dealer         = ',num2str(10000*(yH/expected_inter-r)),' bps']) ;
    

    disp([' ']) ; 
    disp(['The market with search frictions achieves                 = ',num2str(100*(welfare_otc-welfare_autarky)/(welfare_walras-welfare_autarky)),' percent of frictionless gains from trade']);
    
    disp([' ']) 
    disp(['Welfare in the frictionless market                        = ',num2str(welfare_walras)]) ;
    disp(['Welfare in autarky                                        = ',num2str(welfare_autarky)]) ;
    disp(['Welfare of customers in the otc market                    = ',num2str(welfare_otc_cust)]) ;
    disp(['Welfare of dealers in the otc market                      = ',num2str(welfare_otc_dealer)]) ;
    disp(['Total welfare in the otc market, adding value functions   = ',num2str(welfare_otc_cust+welfare_otc_dealer)]); 
    disp(['Total welfare in the otc market, adding flow utils        = ',num2str(welfare_otc)]); 
    disp(['Numerical error                                           = ',num2str(abs(welfare_otc_cust+welfare_otc_dealer-welfare_otc)/welfare_otc*100),' percent']);
    disp(['Welfare of in the competitive market                      = ',num2str(welfare_walras)]) ;
    
    disp([' ']);
    disp(['Customers appropriate                                     = ',num2str((welfare_otc_cust-welfare_autarky)/(welfare_otc-welfare_autarky)*100),' percent of gains from trade']) ;
    disp(['Dealers   appropriate                                     = ',num2str((welfare_otc_dealer)/(welfare_otc-welfare_autarky)*100),' percent of gains from trade']) ;
else
    disp('Found that piH<s or piH*piL>s and exit')   
end ;


