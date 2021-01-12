function F = RAST_System(X, partialPEthanol, xEthanol, qTot)

lambdaEthanol = X(1); lambdaWater = X(2); PEthanol = X(3); PWater = X(4);

P = 101325;

% Pi(P1O) - Pi(P2O)
F(1) = (5.576*PEthanol^2 + 1.291e+04*PEthanol - 4.806e+05) - (0.03177*PWater^2 + 1977*PWater - 7275);

% P1 - P1O*x1
F(2) = partialPEthanol - PEthanol*xEthanol*lambdaEthanol;

% P2 - P2O*x2
F(3) = P - partialPEthanol - PWater*(1 - xEthanol)*lambdaWater;

%1/qTot - x1/q1(P1O) - x2/q2(P2O)
F(4) = (1/qTot) - (xEthanol/(1.629e-10*PEthanol^3 - 2.045e-06*PEthanol^2 + 0.007971*PEthanol + 3.934)) - ((1 - xEthanol)/(-2.47e-08*PWater^2 + 7.406e-05*PWater + 0.7792)); 

end