function [Imot,Pmot,Pload_array,Qload,omega,eff] = motorCalc(V,Kv,I0,Rm,Imax)
% MOTORCALC(V,KV,I0,Rm) calculates motor current, power, angular speed,
% and torque from motor characteristics and applied voltage.
% MOTORCALC(V,KV,I0,Rm,Imax) truncate output values when I > Imax.
% A virtual load is applied from zero to the maximum possible value to
% generate an array of output.

%     Copyright (C) 2021 Danilo Ciliberti danilo.ciliberti@unina.it
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program. If not, see <https://www.gnu.org/licenses/>.

% Check if maximum current has been assigned
if nargin < 5
   Imax = 0; 
end

Piron = V * I0;
Pload_max = 0.999*(V^2/(4*Rm) - Piron);

steps = 1000;
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

% Delete motor data when motor current is above max current
if Imax > 0
    idx = find(Imot < Imax);
    Pload_array = Pload_array(idx);
    Imot = Imot(idx);
    Pmot = Pmot(idx);
    omega = omega(idx);
    Qload = Qload(idx);
    eff = eff(idx);
end

end