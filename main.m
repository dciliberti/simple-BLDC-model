% Simple BLDC motor model
% A simple (and probably inaccurate) brushless DC motor modelling script
% See companion PDF file generated from the MATLAB live script for details
%
%     Copyright (C) 2023 Danilo Ciliberti danilo.ciliberti@unina.it
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

close all; clearvars; clc

%% Input data
motorName = 'Hacker A30-12XL V4';
V = 15.0;               % Voltage
Imax = 60;              % max current, Ampere (useful for axis limits)

% Motor constants
Kv = 149;               % Motor RPM/Volt constant
I0 = 1.6;               % Motor idle (no load) current, A
Rm = 0.020;             % Motor internal resistance, Ohm

%% Calculate and plot motor performance

[Imot, Pmot, Pload, Qload, omega, eff] = motorCalc(V,Kv,I0,Rm,Imax);
RPM = omega*60/(2*pi);

t = tiledlayout(3,1);
ax1 = nexttile;
plot(Imot,Pmot,'LineWidth',2)
ylabel('Power, W')
hold on
plot(Imot,Pload,'LineWidth',2)
hold off
legend('Electric power','Shaft power','Location','best')
grid on

ax2 = nexttile;
yyaxis left
plot(Imot,Qload,'LineWidth',2)
ylabel('Torque, Nm')
%
yyaxis right
plot(Imot,RPM,'LineWidth',2)
ylabel('RPM')
ax2.YAxis(2).Exponent = 0;
% ytickformat('%.0f')
grid on


ax3 = nexttile;
plot(Imot,eff,'LineWidth',2)
grid on
ylabel('Efficiency')

linkaxes([ax1,ax2,ax3],'x');
xticklabels([ax1,ax2],{})
xlabel(t,'Motor current, A')
title(t,[motorName, ' (',num2str(V,'%.1f'),' V)'],'interpreter','none')
t.TileSpacing = 'compact';

%% Build a motor map

volt = linspace(V/10,V,10);    % give reasonable voltage range
maxTorque = zeros(1,numel(volt));
revMaxTorque = zeros(1,numel(volt));
x = []; %#ok<*AGROW> 
y = [];
z = [];
for idx = 1:numel(volt)

    [~, ~, ~, torque, radSec, etaMot] = motorCalc(volt(idx),Kv,I0,Rm,Imax);
    revMinute = radSec*60/(2*pi);

    [maxTorque(idx), pos] = max(torque);
    revMaxTorque(idx) = revMinute(pos);

    x = [x, revMinute];
    y = [y, torque];
    z = [z, etaMot];

end

figure, hold on

xArray = [revMaxTorque, fliplr(revMinute(1:end-1))];
yArray = [maxTorque, fliplr(torque(1:end-1))];
area(xArray,yArray,'FaceAlpha',0.5)

[xNorm, Cx, Sx] = normalize(x);
[yNorm, Cy, Sy] = normalize(y);
[zNorm, Cz, Sz] = normalize(z);
[X, Y] = meshgrid(linspace(min(xNorm),max(xNorm)), linspace(min(yNorm),max(yNorm)));  % x: RPM, y: torque
Z = griddata(xNorm,yNorm,zNorm,X,Y);
levels = 0.70:0.05:0.95;
[C,H] = contour(X*Sx+Cx,Y*Sy+Cy,Z*Sz+Cz,levels,'black');
clabel(C,H,levels)

hold off, grid on
xlabel('RPM'), ylabel('Torque, Nm')
title([motorName, ' (map)'],'interpreter','none')