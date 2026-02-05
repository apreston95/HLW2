 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISPLAY RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('%% Demographic parameters                %%');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp(' ') ;
disp(' ') ;
disp(['Supply per capita                          = ',num2str(s)]) ;
disp(['Relative size of dealer sector             = ',num2str(m)]) ;
disp(['Type switching intensity                   = ',num2str(gamma),' per year']) ;
disp(['Proba of switch to high                    = ',num2str(piH)]) ;
disp(['Intensity of customer to dealer contact    = ',num2str(rho*m),' per year']); 
disp(['Intensity of dealer to dealer contact      = ',num2str(lambda),' per year']) ;
 
disp(' '); 
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('%% Valuation parameters                  %%');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp(' ') ;
disp(' ') ;
disp(['Discount rate                               = ',num2str(r)]);
disp(' ') ; 
disp(['DC bargaining power, theta                  = ',num2str(theta)]);
disp(['DD bargaining power, dealer-owner, theta1   = ',num2str(theta1)]);
%disp(['DD bargaining power, dealer non-owner,        = ',num2str(theta0)]);
disp(' ') ;
disp(['Flow utility of customer in H state   = ',num2str(yH)]);
disp(['Flow utility of customer in L state         = ',num2str(yL/yH),' times yH']);
disp(' ') ;
disp(['Flow utility of dealer, upper bound         = ',num2str(deltaHigh/yH),' times yH']);
disp(['Flow utility of dealer, lower bound         = ',num2str(deltaLow/yH),' times yH']);
disp(['Flow utility of dealer, beta distribution a = ',num2str(beta_delta_a)]) ; 
disp(['Flow utility of dealer, beta distribution b = ',num2str(beta_delta_b)]) ; 
disp(' ') ;
disp(['Match specific value, lower bound           = ',num2str(epsiLow)]) ; 
disp(['Match specific value, upper bound           = ',num2str(epsiHigh)]) ; 
disp(['Match specific value, beta distribution a   = ',num2str(beta_epsi_a)]) ; 
disp(['Match specific value, beta distribution b   = ',num2str(beta_epsi_b)]) ; 
disp(' '); 
disp(['Fraction of asset-less dealers M0/M         = ',num2str(M0_over_M)]) ; 
%%% THE LIST OF VARIABLE
% theta                = vec_var(1) ;
% theta1               = vec_var(2) ;
% deltaH               = vec_var(3) ;
% deltaHigh            = vec_var(4) ; 
% deltaLow             = vec_var(5) ;
% beta_delta_a         = vec_var(6) ;
% beta_delta_b         = vec_var(7) ;
% epsiHigh             = vec_var(8) ;
% beta_epsi_a          = vec_var(9) ;
% beta_epsi_b          = vec_var(10); 


disp([' ']);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('%% Basic Checks                          %%');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp(' ') ;
disp(['m0 according to back solving          = ',num2str(m0)]);
disp(['m0 according to cubic equation        = ',num2str(m0Eqm)]);
disp(['Sum of muh1, mul1, and m1             = ',num2str(muh1Eqm+mul1Eqm+m1Eqm)]); 
disp(['Aggregate asset supply                = ',num2str(s)]) ;
disp([' ']);

disp([' ']);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('%% Demographics                          %%');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp(['mu_h0                                 = ',num2str(muh0Eqm)]); 
disp(['mu_l1                                 = ',num2str(mul1Eqm)]); 
disp(['m0                                    = ',num2str(m0Eqm)]); 
disp(['m1                                    = ',num2str(m1Eqm)]); 
disp([' ']);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('%% Contact times                        %%');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp(['Dealer time to contact a h0 customer  = ',num2str(1/rho/muh0Eqm*250),' days']) ;
disp(['Dealer time to contact a l1 customer  = ',num2str(1/rho/mul1Eqm*250),' days']) ;
disp(['Dealer time to contact a m0 dealer    = ',num2str(1/lambda/m0*m*250),' days']) ;
disp(['Dealer time to contact a m1 dealer    = ',num2str(1/lambda/m1*m*250),' days']) ;

disp(['Customer time to contact a m0 dealer  = ',num2str(1/rho/m0*250),' days']) ;
disp(['Customer time to contact a m1 dealer  = ',num2str(1/rho/m1*250),' days']) ;

disp([' '])
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%') ;
disp('%% The following should be in ascending order for HLW patterns to hold %%') ;
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%') ;
disp(['Reservation value of a low  customer  = ',num2str(DWL)]) ;
disp(['Reservation value of a low  dealer    = ',num2str(DVlow)]) ;
disp(['Reservation value of a high dealer    = ',num2str(DVhigh)]) ;
disp(['Reservation value of a high customer  = ',num2str(DWH)]) ;
BP = 10000*(Delta(1,1)/Delta(2,1)-1) ;
disp(['There are ',num2str(BP),' basis points between the reservation value of a low and a high type']); 

disp([' '])
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%') ;
disp('%% The Size of Markup Conditional on Chain length                      %%') ;
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%') ;
disp('In column: the  total markup in a chain of length k=1,2,..., in basis points');
disp(round(SizeOfMarkup*10000));
disp([' '])
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%') ;
disp('%% The Matrix of the Share of Markup Along Chain                       %%') ;
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%') ;
disp('In rows: the conditioning variable, chain length, k') ;
disp('In column: the share of total markup appropriated by dealer j, for j=1,2,...,k, in percent');
disp(round(DistributionOfMarkupMat*100));
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%') ;
disp('%% Two moments about markup and chain length                           %%') ;
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%') ;
disp(['The average markup in LS               = ',num2str(avg_markup_LS*100),' bps']) ;
disp(['The average markup in  model           = ',num2str(avg_markup*10000),' bps']) ;
disp(['The average chain length in LS         = ',num2str(avg_chain_length_LS)]) ;
disp(['The average chain length in model      = ',num2str(avg_chain_length)]) ;
disp(['beta of markup on  length in LS        = ',num2str(beta_markup_chain_length_LS*100),' bps']) ;
disp(['beta of markup on  length in model     = ',num2str(beta_markup_chain_length*10000),' bps']) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FIGURES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
 figure(1) ; %Asset distribution in the dealer sector
 hold on ;
 plot(xMat,Phi1(xMat,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm)/m,'b-');
 plot(xMat,Phi0(xMat,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm)/m,'r-');
 plot(xMat,ones(size(xMat)),'k--');
 grid on ;
 axis tight
 hold off
 xlabel('x');
 legend('\Phi_1(x)/m','\Phi_0(x)/m','Location','Northwest');
 legend('Boxoff'); 

% 
 figure(2) %Buying and selling intensities
 subplot(2,1,1)
 hold on;
 plot(xMat, rho*mul1Eqm + lambda0(xMat,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm),'b-') ;
 plot(xMat, rho*muh0Eqm+  lambda1(xMat,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm),'r-') ;

 xlabel('x') ;
 legend('DC + DD Buying intensity','DC + DD Selling intensity','Location','Northwest')
 legend('Boxoff');
 hold off
 subplot(2,1,2) ;
 hold on;
 plot(xMat,lambda1(xMat,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm),'b-') ;
 plot(xMat,lambda0(xMat,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm),'r-') ;
 xlabel('x') ;
 legend('DD Buying intensity','DD Selling intensity','Location','Northwest')
 legend('Boxoff');
 hold off

 figure(3)
contribution_to_volume =  2*(rho*mul1Eqm+lambda0(xMat,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm))...
          .*(rho*muh0Eqm+lambda1(xMat,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm))...
          ./( rho*mul1Eqm + lambda0(xMat,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm)...
             +rho*muh0Eqm + lambda1(xMat,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm) ) ;
plot(xMat,contribution_to_volume,'b-') ;
grid on
title('Contribution to Volume')


figure(4) %Reservation values
hold on;
plot(xMat,(DWH)*ones(1,gridsize),'b-') ;
plot(xMat,DVmat,'r-') ;
plot(xMat,DWL*ones(1,gridsize),'g-') ;
%axis([deltaLow,deltaHigh,Delta(2,1)*0.97,Delta(1,1)*1.03])
axis tight
grid on
hold off
xlabel('x')
legend('\Delta V_h','\Delta  V(x)','\Delta  V_l','Location','East'); 
legend('Boxoff')

%%
figure(5)
hold on
b = bar(1:5,[freq_chain_length_LS(1:5)',FrequencyOfChain(1:5)']*100) ;
set(b(1),'FaceColor',rgb('darkgreen')) ;
set(b(2),'FaceColor',rgb('darkorange')) ;
set(gca,'XTick',[1 2 3 4 5]) ;
title('Distribution of chain length (percent)') ;
legend('empirical distribution','model distribution') ;
legend('boxoff')

grid on
hold off
