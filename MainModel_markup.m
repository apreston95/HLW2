ExpRatioResValueToBid                           = zeros(K,K);
ExpInverseOfBid                                 = zeros(K,1) ;
ExpectedMarkupMat                               = zeros(K,K) ;  
ExpRatioExtraToBid                              = zeros(K,1) ;                                            
                                            
% The first function to integrate: 
% the ratio of reservation value of x1 to bid, multiplied by the distribution of x1 conditional on chain length
TheFirstFunctionToIntegrateForjEqual1           = @(x1,k)   DVfunc(x1)./(theta*DWL + (1-theta).*DVfunc(x1)).*ConditionalDistributionOfx1(x1,k,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm) ;
% The second function to integrate: 
% the inverse of the bid, multiplied by the distribution of x1 conditional on chain length
TheSecondFunctionToIntegrateForjEqual1           = @(x1,k)   1./(theta*DWL + (1-theta).*DVfunc(x1)).*ConditionalDistributionOfx1(x1,k,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm) ;
% The third function to integrate:
% the ratio of reservation value of xj to bid, multiplied by the joint distribution of x1 and xj conditional on chain length                                   
TheThirdFunctionToIntegrateForOtherj             = @(x1,xj,j,k) DVfunc(xj)./(theta*DWL + (1-theta).*DVfunc(x1)).*ConditionalDistributionOfx1xj(x1,xj,j,k,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm) ;
% The fourth function to intergrate
% the ratio of the extra value of trading to bid, multiplied by the joint distribution of x1 and xk conditional on chain length, if k>1
TheFourthFunctionToIntegrate                    = @(x1,xk,k) epsiFunc(xk)./(theta*DWL + (1-theta).*DVfunc(x1)).*ConditionalDistributionOfx1xj(x1,xk,k,k,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm) ;
% The fifth function to intergrate
% the ratio of the extra value of trading to bid, multiplied by the joint distribution of x1 and xk conditional on chain length, if k=1
TheFifthFunctionToIntegrateForjEqual1           = @(x1) epsiFunc(x1)./(theta*DWL + (1-theta).*DVfunc(x1)).*ConditionalDistributionOfx1(x1,1,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm) ;

% The function parameterizing the domain of integration 
% x1 varies in [deltaLow,deltaHigh], while xj varies in [x1,deltaHigh]
LowerBoundxj                                    = @(x1) x1 ;


% Calculation of the relevant expectations
for k=1:K
    % Expected ratio of reservation value of x1 relative to bid
    ExpRatioResValueToBid(k,1)        = integral( @(x1) TheFirstFunctionToIntegrateForjEqual1(x1,k),LowX,1,'AbsTol',TolX) ;
    % Expected inverse of the bid
    ExpInverseOfBid(k,1)              = integral( @(x1) TheSecondFunctionToIntegrateForjEqual1(x1,k),LowX,1,'AbsTol',TolX) ;
    
    if k==1
       ExpRatioExtraToBid(1,1)           = integral( @(x1) TheFifthFunctionToIntegrateForjEqual1(x1),LowX,1,'AbsTol',TolX) ;
    else
        % Expected ratio of extra value relative to bid
        ExpRatioExtraToBid(k,1)           = integral2( @(x1,xk) TheFourthFunctionToIntegrate(x1,xk,k),LowX,1,LowerBoundxj,1,'AbsTol',TolX) ;
    end
    for j=2:k    
    % Expected ratio of reservation value of xj relative to bid
    ExpRatioResValueToBid(k,j)        = integral2( @(x1,xj) TheThirdFunctionToIntegrateForOtherj(x1,xj,j,k),LowX,1,LowerBoundxj,1,'AbsTol',1e-15) ;
    end;
end

% Calculation of expected markups
for k=1:K 
    if k==1
        ExpectedBuyOverBid                        = theta*DWL*ExpInverseOfBid(1,1) + (1-theta)*ExpRatioResValueToBid(1,1) ;
        ExpectedSellOverBid                       = theta*DWH*ExpInverseOfBid(1,1) + theta*ExpRatioExtraToBid(1,1) + (1-theta)*ExpRatioResValueToBid(1,1) ;
        ExpectedMarkupMat(1,1)                    = ExpectedSellOverBid-ExpectedBuyOverBid ;
       
    elseif k>=2
        % Start with j=1
        ExpectedBuyOverBid                        = theta*DWL*ExpInverseOfBid(k,1) + (1-theta)*ExpRatioResValueToBid(k,1) ;
        ExpectedSellOverBid                       = theta0*ExpRatioResValueToBid(k,1)+theta1*ExpRatioResValueToBid(k,2) ;
        ExpectedMarkupMat(k,1)                    = ExpectedSellOverBid-ExpectedBuyOverBid ;
        % Continue with j=k
        ExpectedBuyOverBid                        = theta0*ExpRatioResValueToBid(k,k-1) + theta1*ExpRatioResValueToBid(k,k) ;
        ExpectedSellOverBid                       = theta*DWH*ExpInverseOfBid(k,1) + theta*ExpRatioExtraToBid(k,1) + (1-theta)*ExpRatioResValueToBid(k,k) ;
        ExpectedMarkupMat(k,k)                    = ExpectedSellOverBid-ExpectedBuyOverBid ;
        % Fill in the intermediate markups for j = 2,3,...k-1 (so for this we need k>=3)
        if k>=3
            for j=2:k-1
                ExpectedBuyOverBid                  = theta0*ExpRatioResValueToBid(k,j-1) + theta1*ExpRatioResValueToBid(k,j) ;
                ExpectedSellOverBid                 = theta0*ExpRatioResValueToBid(k,j) + theta1*ExpRatioResValueToBid(k,j+1) ;
                ExpectedMarkupMat(k,j)              = ExpectedSellOverBid-ExpectedBuyOverBid ;
            end;
        end
    end ;
end


DistributionOfMarkupMat = ExpectedMarkupMat./(sum(ExpectedMarkupMat,2)*ones(1,K)) ;
SizeOfMarkup            = sum(ExpectedMarkupMat,2) ;
FrequencyOfChain        = rho*muh0Eqm/(lambda*m0/m)*( Lambda(0,1,rho,lambda,m,m0Eqm,mul1Eqm,muh0Eqm).^(1:K) ) ./ (factorial(1:K)) ;

vec_chain_length        = (1:K)' ;

avg_chain_length        = FrequencyOfChain*vec_chain_length ;

avg_chain_length_sq     = FrequencyOfChain*(vec_chain_length.^2) ;

avg_markup              = FrequencyOfChain*SizeOfMarkup ; 
var_chain_length        = avg_chain_length_sq - avg_chain_length^2 ;

cov_markup_chain_length...
                        = FrequencyOfChain*( SizeOfMarkup.*(vec_chain_length-avg_chain_length) ) ;

beta_markup_chain_length...
                         = cov_markup_chain_length/var_chain_length ;

