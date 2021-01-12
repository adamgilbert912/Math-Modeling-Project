function MathModelingProject

syms x

%%Ethanol Simulation Data for mixture

partialPEthanolMixture = [30.32, 46.52, 74.2, 238.48, 3111.41, 4006.98, 4773.9, 5578, 6517.6, 6776.2];

EthanolMoleculesMixture = [1.42, 2.67, 5.34, 9.82, 12.84, 13.18, 13.52, 13.81, 14.03, 14.32];

f1 = 1.629e-10*x^3 - 2.045e-06*x^2 + 0.007971*x + 3.934; %third degree polynomial for q1 using cftool in command window

figure(1)
ezplot(f1, [0.1e+04, 2e+04])
title('Ethanol Mixture Simulation Data (q1)')
xlabel('Partial Pressure Ethanol (Pa)')
ylabel('Molecules/uc')

piEthanolMixture = [];
P1O = [];

%getting different values of pi by varying P1O
for i = 30:100:2000
    P1O = [P1O, i];
    newVal = double(int(f1, 0, i)) * 298 * 8.314;
    piEthanolMixture = [piEthanolMixture, newVal];
end

assignin('base','piEthanolMixture', piEthanolMixture)
assignin('base','P1O', P1O)

fPi1 = 8.59e+05*x^3 - 3.822e+094*x^2 + 5.875e+12*x - 1.995e+15; %third degreee polynomial for pi using cftool in command window

%%Water Simulation Data for Mixture

partialPWaterMixture = [30.32, 46.52, 75.67, 238.49, 3111.41, 4006.98, 4773.9, 5470.51, 5799.55, 6517.55, 6909.34, 7324.68];
waterMoleculesMixture = [0.398, 0.568, 0.966, 1.25, 0.739, 0.682, 0.568, 0.455, 0.341, 0.227, 0.057, 0.057];

f2 = -2.47e-08*x^2 + 7.406e-05*x + 0.7792; %second degree polynomial for q2 using cftool in command window

figure(2)
ezplot(f2, [0, 7000])
title('Water Mixture Simulation Data (q2)')
xlabel('Partial Pressure Ethanol (Pa)')
ylabel('Molecules/uc')

piWaterMixture = [];
P2O = [];

%getting different values of pi by varying P2O
for i = 30:100:2000
    P2O = [P2O, i];
    newVal = double(int(f2, 0, i)) * 298 * 8.314;
    piWaterMixture = [piWaterMixture, newVal];
end

assignin('base','piWaterMixture',piWaterMixture)
assignin('base','P2O',P2O)

fPi2 = -123.4*x^3 + 2234*x^2 - 1.01e+04*x - 7.7374e-10; %third degree polynomial for pi using cftool in command window

guess = [10000, 10000, 0.002];

options = optimoptions('fsolve');
options.MaxFunctionEvaluations = 9999999;
options.MaxIterations = 999999;

P1O = [];
P2O = [];
xEthanol = [];
qTot = [];

%solving system of equations using varying values of Partial Pressure of
%Ethanol
for i = 30:1000:100000

answer = fsolve(@(X)IAST_System(X, i),guess, options);

P1O = [P1O, answer(1)];
P2O = [P2O, answer(2)];
xEthanol = [xEthanol, answer(3)];

newQTot = (1.629e-10*i^3 - 2.045e-06*i^2 + 0.007971*i + 3.934)/answer(3);

qTot = [qTot, newQTot];

end

lambdaEthanol = [];
lambdaWater = [];
P1ORAST = [];
P2ORAST = [];

%Solving a system of equations using RAST_System to solve for activity
%coefficients
for i = 1:length(partialPEthanolMixture)

IASTAnswer = fsolve(@(X)IAST_System(X, partialPEthanolMixture(i)), guess, options);

qTot = (1.629e-10*i^3 - 2.045e-06*i^2 + 0.007971*i + 3.934)/IASTAnswer(3);

guess2 = [1, 1, IASTAnswer(1), IASTAnswer(2)];

RASTAnswer = fsolve(@(X)RAST_System(X, partialPEthanolMixture(i), IASTAnswer(3), qTot), guess2, options)

lambdaEthanol = [lambdaEthanol, RASTAnswer(1)];
lambdaWater = [lambdaWater, RASTAnswer(2)];
P1ORAST = [P1ORAST, RASTAnswer(3)];
P2ORAST = [P2ORAST, RASTAnswer(4)];
end

assignin('base','lambdaEthanol',lambdaEthanol);
assignin('base','lambdaWater',lambdaWater);

figure(3)
plot([30:1000:100000], xEthanol)
title('Ethanol Mol Fraction in Adsorbed Phase for Ethanol/Water Mixture')
xlabel('Partial Pressure Ethanol (Pa)')
ylabel('Mol Fraction')

figure(4)
plot([30:1000:100000], qTot)
title('Total Molecules in Adsorbed Phase per Unit Cell for Ethanol/Water Mixture')
xlabel('Partial Pressure Ethanol (Pa)')
ylabel('Molecules/uc')

figure(5)
plot([30:1000:100000], P1O)
title('P1O for Ethanol/Water Mixture (Ethanol)')
xlabel('Partial Pressure Ethanol (Pa)')
ylabel('Ethanol Resevoir Pressure (Pa)')

figure(6)
plot([30:1000:100000], P2O)
title('P2O for Ethanol/Water Mixture (Water)')
xlabel('Partial Pressure Ethanol (Pa)')
ylabel('Water Resevoir Pressure (Pa)')

figure(7)
plot(partialPEthanolMixture, lambdaEthanol)
title('Activity Coefficient for Ethanol')
xlabel('Partial Pressure Ethanol (Pa)')
ylabel('Activity Coefficient')

figure(8)
plot(partialPEthanolMixture, lambdaWater)
title('Activity Coefficient for Water')
xlabel('Partial Pressure Ethanol (Pa)')
ylabel('Activity Coefficient')

figure(9)
plot(partialPEthanolMixture, P1ORAST)
title('P1O for Ethanol/Water Mixture using RAST (Ethanol)')
xlabel('Partial Pressure Ethanol (Pa)')
ylabel('Ethanol Resevoir Pressure (Pa)')

figure(10)
plot(partialPEthanolMixture, P2ORAST)
title('P2O for Ethanol/Water Mixture using RAST (Water)')
xlabel('Partial Pressure Ethanol (Pa)')
ylabel('Water Resevoir Pressure (Pa)')
end