function F = IAST_System(X, partialPEthanol)

PEthanol = X(1); PWater = X(2); xEthanol = X(3);

P = 101325; %pressure

% Pi(P1O) - Pi(P2O)
F(1) = (5.576*PEthanol^2 + 1.291e+04*PEthanol - 4.806e+05) - (0.03177*PWater^2 + 1977*PWater - 7275);

% P1 - P1O*x1
F(2) = partialPEthanol - PEthanol*xEthanol;

% P2 - P2O*x2
F(3) = P - partialPEthanol - PWater*(1 - xEthanol);
end