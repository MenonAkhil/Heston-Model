% FI6121 - 40% ASSIGNMENT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART-1 

FilePath = 'C:\Users\menon\Documents\MATLAB\VolSurf.xlsx';
ImportOptions = detectImportOptions(FilePath, 'NumHeaderLines', 2);
T = readtable(FilePath, ImportOptions);
T.Properties.VariableNames = {'Expiry' 'Years' 'ImpFwd' 'M030' 'M033' 'M035' 'M038' 'M040' 'M043' 'M045' 'M048' 'M050' 'M053' 'M055' 'M058' 'M060' 'M063' 'M065' 'M068' 'M070' 'M073' 'M075' 'M078' 'M080' 'M083' 'M085' 'M088' 'M090' 'M093' 'M095' 'M098' 'M100' 'M103' 'M105' 'M108' 'M110' 'M113' 'M115' 'M118' 'M120' 'M123' 'M125' 'M128' 'M130' 'M133' 'M135' 'M138' 'M140' 'M143' 'M145' 'M150' 'M153' 'M155' 'M158' 'M160'};
X= T.Years(2:end);
Y= table2array( T(1, 4:end) )';
Z = table2array( T(2:end,4:end) )';

% Local Volatility Smiles For the 5 Options
figure; hold on
a1 = plot(Y,Z(1:end,1),'LineWidth',4); M1 = "08-11-2019";
a2 = plot(Y,Z(1:end,2),'LineWidth',4); M2 = "15-11-2019";
a3 = plot(Y,Z(1:end,3),'LineWidth',4); M3 = "22-11-2019";
a4 = plot(Y,Z(1:end,4),'LineWidth',4); M4 = "29-11-2019";
a5 = plot(Y,Z(1:end,5),'LineWidth',4); M5 = "06-12-2019";
title('Local Volatility Smiles','FontSize',30)
xlabel('Moneyness','FontSize',20)
ylabel('Local Volatility','FontSize',20)
legend([a1;a2;a3;a4;a5],M1,M2,M3,M4,M5,'FontSize',20);


listed_expiry_dates = [1/252;7/252;14/252;21/252;28/252];
lv_matrix = [1.2186	1.1678	1.1353	1.0883	1.0579	1.0135	0.9846	0.9420	0.9140	0.8725	0.8451	0.8042	0.7769	0.7359	0.7084	0.6677	0.6385	0.5952	0.5656	0.5197	0.4877	0.4371	0.4010	0.3416	0.2968	0.2156	0.1692	0.1344	0.0829	0.0962	0.1118	0.1258	0.1311	0.1392	0.1442	0.1512	0.1556	0.1619	0.1659	0.1715	0.1752	0.1803	0.1837	0.1884	0.1915	0.1959	0.1987	0.2055	0.2094	0.2119	0.2155	0.2178;
    0.4938	0.4763	0.4652	0.4491	0.4388	0.4239	0.4142	0.4001	0.3909	0.3774	0.3685	0.3554	0.3468	0.3340	0.3256	0.3130	0.3046	0.2921	0.2837	0.2711	0.2627	0.2499	0.2412	0.2283	0.2167	0.1885	0.1601	0.1109	0.0808	0.0834	0.1006	0.1287	0.1426	0.1582	0.1663	0.1758	0.1808	0.1863	0.1894	0.1938	0.1967	0.2009	0.2035	0.2074	0.2099	0.2135	0.2158	0.2214	0.2246	0.2267	0.2297	0.2317;
    0.5438	0.5257	0.5142	0.4977	0.4871	0.4718	0.4619	0.4475	0.4381	0.4244	0.4155	0.4023	0.3936	0.3808	0.3724	0.3599	0.3516	0.3393	0.3311	0.3188	0.3100	0.2926	0.2775	0.2484	0.2251	0.1880	0.1629	0.1229	0.0953	0.0898	0.0939	0.1063	0.1211	0.1411	0.1524	0.1669	0.1752	0.1860	0.1923	0.2009	0.2056	0.2122	0.2162	0.2211	0.2241	0.2284	0.2311	0.2378	0.2416	0.2440	0.2476	0.2499;
    0.4950	0.4798	0.4701	0.4563	0.4475	0.4347	0.4265	0.4145	0.4068	0.3955	0.3881	0.3773	0.3703	0.3599	0.3531	0.3430	0.3363	0.3265	0.3198	0.3072	0.2963	0.2755	0.2585	0.2297	0.2099	0.1795	0.1565	0.1203	0.0979	0.0884	0.0908	0.0938	0.0997	0.1146	0.1243	0.1374	0.1451	0.1553	0.1614	0.1695	0.1744	0.1810	0.1850	0.1904	0.1938	0.1983	0.2012	0.2069	0.2101	0.2122	0.2152	0.2171;
    0.5251  0.5047	0.4941	0.4788	0.4691	0.4549	0.4458	0.4326	0.4240	0.4104	0.4031	0.3911	0.3832	0.3715	0.3638	0.3525	0.3449	0.3317	0.3213	0.3028	0.2884	0.2645	0.2480	0.2244	0.2089	0.1816	0.1591	0.1250	0.1052	0.0899	0.0921	0.0958	0.1033	0.1207	0.1319	0.1469	0.1558	0.1678	0.1750	0.1848	0.1908	0.1989	0.2038	0.2106	0.2148	0.2206	0.2242	0.2325	0.2369	0.2394	0.2432	0.2446];
simulation_dates =1/252:1/252:28/252;


% Spline Interpolated Matrix (Smoothened Volatility Surface Matrix)
interpolated_matrixlv=interp1(listed_expiry_dates,lv_matrix,simulation_dates,'spline');

A = transpose(simulation_dates);
B = [0.30	0.33	0.35	0.38	0.40	0.43	0.45	0.48	0.50	0.53	0.55	0.58	0.60	0.63	0.65	0.68	0.70	0.73	0.75	0.78	0.80	0.83	0.85	0.88	0.90	0.93	0.95	0.98	1.00	1.03	1.05	1.08	1.10	1.13	1.15	1.18	1.20	1.23	1.25	1.28	1.30	1.33	1.35	1.38	1.40	1.43	1.45	1.50	1.53	1.55	1.58	1.60];
C = B';
E = transpose(interpolated_matrixlv);


% Building Local Volatility Surface : Cubic Spline
figure
surf(A,C,E)
title('Cubic Spline Interpolation of Local Volatility : Fine Mesh ','FontSize',30)
xlabel('T (days) ','FontSize',20)
ylabel('Moneyness ','FontSize',20)
zlabel('Local Volatility ','FontSize',20)
axis tight
rotate3d on


%Implied forward prices for 5 options
Imp = [3075.60;3074.88;3074.83;3074.39;3073.89];

% Interpolating Forward Prices for 28 days on both linear and spline
lv_1 = interp1(listed_expiry_dates,Imp,simulation_dates,'spline');
ImpFwd_spline = transpose(lv_1);



% MONTE CARLO SIMULATION - CUBIC SPLINE
St_price = 3076.74;
dt=1/252;
N = 10000;
e = randn(28,N);
mu =0.00567;
simulations = ones(28,N);
simulations(1,:) = simulations(1,:)*St_price;
forwardPrice_spline = ImpFwd_spline;
LocalVolMoneyness = B;
LocalVolSurface_spline = interpolated_matrixlv;
for i=1:1:27
lv_spline = interp1(LocalVolMoneyness,LocalVolSurface_spline(i+1,:),simulations(i,:)./forwardPrice_spline(i),'spline','extrap');
simulations(i+1,:) = simulations(i,:).*exp((mu-.5*lv_spline.^2)*dt+lv_spline.*sqrt(dt).*e(i,:));
end

USt_spline = [St_price zeros(1,28-1)];
LSt_spline = [St_price zeros(1,28-1)];
ESt_spline = [St_price zeros(1,28-1)];
lv2_spline(1,1)=interp1(LocalVolMoneyness,LocalVolSurface_spline(1,:),ESt_spline(1,1)./forwardPrice_spline(1),'spline','extrap');

%95 percent confidence interval%

for i=1:1:27
    ESt_spline(i+1)=simulations(1,1)*exp(mu*i*dt);
    USt_spline(i+1)=ESt_spline(i+1)+1.96*lv2_spline(1,i)*simulations(1,1)*sqrt(i*dt);
    LSt_spline(i+1)=ESt_spline(i+1)-1.96*lv2_spline(1,i)*simulations(1,1)*sqrt(i*dt);
    lv2_spline(1,i+1)=interp1(LocalVolMoneyness,LocalVolSurface_spline(i+1,:),ESt_spline(1,i+1)./forwardPrice_spline(i+1),'spline','extrap');
end


% Plotting Monte Carlo Simulation
figure
plot(simulations);
title('Local Volatility Monte Carlo Simulations : 10000 Simulations');
ylabel('S&P 500 INDEX PRICE');
xlabel('T (days)');
hold on;
plot(ESt_spline,'--','color','black','LineWidth',6);
plot(USt_spline,'--','color','yellow','LineWidth',3);
plot(LSt_spline,'--','color','yellow','LineWidth',3);

%ret_step = zeros(28,10000);
for i = 1:1:10000
    for j = 1:1:27
ret_step(j+1,i)=log(simulations(j+1,i))-log(simulations(j,i));
    end
end

% Calculating Returns at maturity
ret = zeros(1,10000);
for i = 1:1:10000
    ret(1,i) = (simulations(28,i) - simulations(1,i))/simulations(1,i);
end


% PLOT 5 :Histogram of Simulated S&P 500 Index Prices at Maturity
fig_1 = simulations(28,1:10000);
figure; 
histfit(fig_1); 
title('Histogram of Simulated S&P 500 Index Prices at Maturity'); 
xlabel('Simulated S&P 500 Index Prices'); 
ylabel('Frequency');

% PLOT 6 :Histogram of Returns of Simulated S&P 500 Index Prices at Maturity
[g,xi] = ksdensity(ret);
kurt = kurtosis(ret);
ske = skewness(ret);
figure;
plot(xi,g,'LineWidth',2);
title('Returns Distribution of Simulated S&P 500 Index Prices');
xlabel('Returns');
ylabel('Density');
annotation('textbox', [0.7, 0.3, 0.1, 0.1], 'String', "Kurtosis = " + kurt);
annotation('textbox', [0.7, 0.2, 0.1, 0.1], 'String', "Skewness = " + ske);

% CBOE Skew Corroboration
xo = randn(N,1);
trans_ret = ret;
CurrentSkewIndex = 125.2;
bl = -3*std(trans_ret) - mean(trans_ret);
b2 = -2*std(trans_ret) - mean(trans_ret);
CBOEsigma = interp1([125,130],[0.0163,0.0192],CurrentSkewIndex,'spline');
ModelSigma = ksdensity(trans_ret,bl,'function','cdf');
MSigma = ksdensity(trans_ret,b2,'function','cdf');


skewint = interp1([0.0074,0.0104],[110,115],ModelSigma,'spline');





% Linear Values
% ImpFwd_linear - Interpolated Implied Forward Prices
% B - LocalVolMoneyness
% interpolated_lvmatrix - LocalVolSurfaces

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TO COMPLETE IN PART-1 OF PROJECT
% STEP 1. Need to find sigma(Volatility) value??
% STEP 2. Run Monte-Carlo Simulation - 10,000 times to simulate 28 trading days in the future
% STEP 3. Calculate mean and standard deviation for the Implied return (Obtained
% from MC-Simulation of 10,000 days) 
% STEP 4. to impute the tail probability for a return outlier 3 or more standard deviation below the mean
% STEP 5.Then benchmark this with the corresponding value computed for the interpolated Skew index value