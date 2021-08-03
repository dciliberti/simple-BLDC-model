function [Imot,Pmot,Pload_array,Qload,omega,eff] = motorCalc(V,Kv,I0,Rm,Pmax)
% MOTORCALC(V,KV,R;,PLOAD) calculates motor current, power, angular speed,
% and torque from motor characteristics and applied voltage.
% A virtual load is applied from zero to the maximum possible value.

Piron = V * I0;
Pload_max = 0.999*(V^2/(4*Rm) - Piron);

if Pload_max > Pmax
    Pload_max = Pmax;
end

steps = 100;
Pload_array = linspace(0,Pload_max,steps);

c = 0;
Imot = zeros(1,steps);
Pmot = Imot;
omega = Imot;
Qload = Imot;
eff = Imot;
for Pload = Pload_array
    c = c + 1;
    
    Imot(c) = (V - sqrt(V^2 - 4*Rm*(Piron + Pload))) / (2*Rm);
    
    Pmot(c) = V * Imot(c);
    RPM = Kv * (V - Rm*Imot(c));
    omega(c) = 2*pi/60 * RPM;
    Qload(c) = Pload/omega(c);
    eff(c) = Pload/Pmot(c);
end

end