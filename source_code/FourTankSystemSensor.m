function y=FourTankSystemSensor(x,p)
A = p(5:8,1); % Tank cross sectional areas [cm2]
rho = p(12,1); % Density of water [g/cm3]
y = zeros(4,1);
y = x./(rho*A);% Liquid level in each tank [cm]
end