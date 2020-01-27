% FI6121 - 40% ASSIGNMENT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART-2 

% Initial Price of S&P 500 Index
S0 = 3076.74;
% Strike
K = 3076.74;

% Total Trading Days in a year
trade_days = 365;
% 90 Day-up call option
no_of_time_steps = 90;
time = no_of_time_steps/trade_days;
delT = 1/365;
% US Rate Mmkt
rf = 0.01904;
% Dividend Yieldclear
div = 0.01963;
% Annualized Total Return
mu = 0.00567;
% Heston Parameters for the heston stochastic volatility model
% Initial Variance
V0 = 0.011;
% Long Run Variance
theta = 0.048;
% Correlation
rho = -0.7401;
% Volatility of Volatility
sigma = 1.30129;
% Mean Reversion
kappa = 4.2721;

No_of_sims = 10000;
NormRand1 = randn(no_of_time_steps, No_of_sims);
NormRand2 = randn(no_of_time_steps, No_of_sims);
St = [S0*ones(1,No_of_sims); zeros(no_of_time_steps-1,No_of_sims)];
V = [V0*ones(1,No_of_sims); zeros(no_of_time_steps-1,No_of_sims)];

W1 = NormRand1;
W2 = zeros(no_of_time_steps,No_of_sims);

for i = 1:1:No_of_sims
    for j = 1:1:no_of_time_steps
        W2(j,i) = rho*NormRand1(j,i)+sqrt(1-rho^2)*NormRand2(j,i);
    end
end

for i = 1:1:No_of_sims
    for j = 2:1:no_of_time_steps
        St(j,i) = St(j-1,i) + (rf-div)*St(j-1,i)*delT + St(j-1,i)*sqrt(V(j-1,i))*W1(j,i)*sqrt(delT);
        V(j,i) = V(j-1,i) + kappa*(theta - V(j-1,i))*delT + sigma*sqrt(V(j-1,i))*W2(j,i)*sqrt(delT);
        V(j,i) = abs(V(j,i));
    end
end

Barriers = [1.00:0.01:1.30].*S0;
Heston_KOBarrierCall_Premium = zeros(length(Barriers),1);
 Call_Payoff = zeros(1,No_of_sims);

for i = 1:1:length(Barriers)
    KO_Barrier = Barriers(i);
    
    %Call_Payoff = zeros(1,No_of_sims);
    for j = 1:No_of_sims
        Call_Payoff(1,j) = max(0,St(end,j)-K);
        if max(St(:,j)) >= KO_Barrier == 1
            Call_Payoff(1,j) = 0;
        else
            Call_Payoff(1,j) = exp(-rf*time)*Call_Payoff(1,j);
        end
    end
    Heston_KOBarrierCall_Premium = mean(Call_Payoff);
end

Heston_VanillaCall_Premium = exp(-rf*time)*mean(max(0,St(end,:)-K));

figure
plot(St)
hold all
Barrier1 = plot(KO_Barrier*ones(no_of_time_steps),'-','Color','blue','LineWidth',3);
Strike1 = plot(K*ones(no_of_time_steps),'-','Color','red','LineWidth',3);
title(['Heston Monte Carlo Simulations: 10000 Simulations'])
xlabel('T (Days)')
ylabel('S&P 500 Index Price')

figure
scatter(V(end,:),St(end,:)./K,10,'filled','d','MarkerFaceColor',[0 0.7 0.7],'MarkerEdgeColor',[0 0.5 0.5])
hold all
Strike = refline([0 K/S0]);
set(Strike,'Color','r','LineWidth',2)
Barrier = refline([0 KO_Barrier/S0]);
set(Barrier,'Color','b','LineWidth',2)
title(' Moneyness v/s Heston Volatilities at Maturity')
xlabel('Heston Volatility')
ylabel('Moneyness')
legend([Barrier Strike],'130% Barrier','Strike','location','NorthEast')

% Monte Carlo Simulation : Black Scholes - CS
St2 = [S0*ones(1,No_of_sims);zeros(no_of_time_steps-1,No_of_sims)];
%Bloomberg Quoted Volatility
BS = 0.13020; 

E = randn(no_of_time_steps,No_of_sims);
Call_PayoffBS = zeros(1,No_of_sims);
BSKO_Barrier = S0*1.3;


for i=1:1:No_of_sims
    for j=2:1:no_of_time_steps
        St2(j,i) = St2(j-1,i)*exp((mu-0.5*BS^2)*delT+BS*sqrt(delT)*E(j,i));
        Call_PayoffBS(1,i)=max(0,St2(end,i)-K);
        if max(St2(:,i)) >= BSKO_Barrier == 1
            Call_PayoffBS(1,i) = 0;
        else
            Call_PayoffBS(1,i) = exp(-rf*time)*Call_PayoffBS(1,i);
        end
    end
    BS_KOBarrierCall_Premium = mean(Call_PayoffBS);
end

BS_VanillaCall_Premium = exp(-rf*time)*mean(max(0,St2(end,:)-K));

returns=(St(end,:)-St(1,:))./St(1,:);
%figure
%disp('Kurtosis and Skewness for Negative Correlation Case')
%display(num2str(kurtosis(returns)))
%display(num2str(skewness(returns)))
%[f,xi]=ksdensity(returns);
%plot(xi,f,'--r');
%hleg1 = legend('\rho = 0.0', '\rho =-0.5');
%zoom on
%toc
% A check to see that the correlation between the simulated returns and
% (changes in) variance is as specified in the Cholesky Decomposition above
%format bank
%corrcoef(price2ret(St),diff(V))