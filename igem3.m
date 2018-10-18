% Fluorescence/OD Analysis at Various IPTG Concentrations

x = linspace(380,820);
XIPTG = [400:100:800];
FluorOD = [32565 31556 31470 30481 30010];
ymid = [74455.40014 75360.77596 74239.54373 74749.81694 73752.76481];
ypos = [523.23 759.055 887.198 877.331 1000.92];
yneg = [450.3 641.575 679.32 713.635 824.29];

%%
%p1*x + p2
%p1 =      -6.185  (-8.575, -3.795)
%p2 =   3.493e+04  (3.345e+04, 3.64e+04)

hold on
scatter(XIPTG, FluorOD, 'x', 'black')
plot(x, -6.185*x + 3.493e+04, 'r')
ylabel('Fluorescence/OD')
xlabel('IPTG Concentration (mM)')
xlim([351 849])
ylim([2.93*10^4 3.32*10^4])

%%
% Total Protein

hold on
ylabel('IDF')
xlabel('IPTG Concentration (mM)')
bar(XIPTG, ymid, 'white')
errorbar(XIPTG, ymid, yneg, ypos, '+')
ylim([7.11*10^4 7.8*10^4])
xlim([351 849])
