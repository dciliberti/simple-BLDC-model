function [Imot,Pelec,Pshaft_array,Qtorque,omega,eff] = ...
    motorCalc(Vmax,Kv,I0,Vref,Rm,phi,Imax)
% MOTORCALC(Vmax,KV,I0,Vref,Rm) calculates motor current, power, angular speed,
% and torque from motor characteristics and applied voltage.
% MOTORCALC(Vmax,KV,I0,Vref,Rm,phi) scales V with throttle value phi (0 to 1).
% MOTORCALC(Vmax,KV,I0,Rm,Vref,phi,Imax) truncates output values when I > Imax.
% A virtual load is applied from zero to the maximum possible value to
% generate an array of output.

%     Copyright (C) 2024 Danilo Ciliberti danilo.ciliberti@unina.it
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

% Check if throttle value has been assigned and it has the correct range
if nargin < 5
    phi = 1.0;
elseif phi < 0 || phi > 1
    error('Throttle value out of range (0,1)')
end

% Check if maximum current has been assigned
if nargin < 6
   Imax = 0;
end

V = phi*Vmax;           % applied voltage
I0 = I0*sqrt(V/Vref);   % scale no-load current with applied voltage
Pno_load = V * I0;      % no-load power losses
Pshaft_max = 0.999*(V^2/(4*Rm) - Pno_load); % maximum shaft power

% Initialize shaft power array
steps = 1000;
Pshaft_array = linspace(0,Pshaft_max,steps);

c = 0;
Imot = zeros(1,steps);
Pelec = zeros(1,steps);
omega = zeros(1,steps);
Qtorque = zeros(1,steps);
eff = zeros(1,steps);
for Pload = Pshaft_array
    c = c + 1;
    
    % Small root of the second order equation (physical solution)
    Imot(c) = (Vmax - sqrt(Vmax^2 - 4*Rm*(Pno_load + Pload))) / (2*Rm);
    
    Pelec(c) = Vmax * Imot(c);
    RPM = Kv * (V - Rm*Imot(c));
    omega(c) = 2*pi/60 * RPM;
    Qtorque(c) = Pload/omega(c);
    eff(c) = Pload/Pelec(c);
end

% Delete motor data when motor current is above max current
if Imax > 0
    idx = find(Imot < Imax);
    Pshaft_array = Pshaft_array(idx);
    Imot = Imot(idx);
    Pelec = Pelec(idx);
    omega = omega(idx);
    Qtorque = Qtorque(idx);
    eff = eff(idx);
end

end