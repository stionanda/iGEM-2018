%Optical Density Analysis of E. coli at Various IPTG Concentrations

hold on
x = linspace(380,820);
XIPTG = [400:100:800];
yleft = [2.293 2.366 2.376 2.441 2.458];
ypos1 = [0.015 0.022 0.026 0.026 0.030];
yneg1 = [0.015 0.022 0.024 0.026 0.031];
%yyaxis left
%plot(x, 0.000405*x + 2.144)
%ylabel('Maximum optical cell density [1/b]')
%errorbar(XIPTG, yleft, yneg1, ypos1, 'sq')
%ylim([2.1501 2.5499])

yright = [0.4215 0.4649 0.4735 0.4773 0.4759];
ypos2 = [0.0076 0.0132 0.0147 0.0154 0.0177];
yneg2 = [0.0077 0.0132 0.0147 0.0153 0.0176];
%yyaxis right
%plot(x, 0.1606 + 0.0009303*x - 6.743e-07*x.^2)
ylabel('Population growth rate [k]')
errorbar(XIPTG, yright, yneg2, ypos2, 'sq')
ylim([0.300 0.600])

xlabel('IPTG Concentration (mM)')
xlim([351 849])