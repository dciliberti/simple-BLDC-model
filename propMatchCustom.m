% MOTOR-PROPELLER MATCH
% Script to match user-imported propeller data at constant flow speed (and
% variable RPM) with motor data.
%
% This script makes use of the functions set POLYFITZERO:
%
% Mark Mikofski (2022). polyfitZero 
% (https://www.mathworks.com/matlabcentral/fileexchange/35401-polyfitzero), 
% MATLAB Central File Exchange. Retrieved June 9, 2022.
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

%% Initial input data
Vmax = 15.0;            % max voltage (limited by the motor)
Kv = 700;               % motor RPM/Volt constant
I0 = 1.5;               % motor idle (no load) current, A
Imax = 60;              % max current, A
Rm = 0.034;             % motor internal resistance, Ohm

% Import propeller data from csv file. Columns order as follows:
% RPM, Power (W), Thrust (N), Torque (Nm)
filename = 'therm-prop-flowspeed-35ms.csv';
propTable = readtable(filename);
propData = propTable.Variables;

% Propeller data curves fitting with RPM as the independent variable
% (usually a quadratic fit is good, but we need to intercept zero with a
% horizontal tangent). It expects data columns in this order:
% RPM, Power (W), Thrust (N), Torque (Nm)
funPropPower = polyfitB0(propData(:,1),propData(:,2),2,0);
funPropThrust = polyfitB0(propData(:,1),propData(:,3),2,0);
funPropTorque = polyfitB0(propData(:,1),propData(:,4),2,0);

%% Cycle over several voltage (simulate motor throttle)

figure(1), hold on
xProp = linspace(0,max(propData(:,1))); % x array of propeller RPM, from 0
yProp = polyval(funPropPower,xProp);    % y array of propeller power
plot(xProp,yProp,'k','LineWidth',2,'DisplayName','Prop Power')

figure(2), hold on
yProp = polyval(funPropTorque,xProp);    % y array of propeller torque
plot(xProp,yProp,'k','LineWidth',2,'DisplayName','Prop Torque')

c = 0;  % counter
volt = linspace(4,Vmax,10); % give a reasonable lower voltage
for V = volt

    throttle = V/Vmax*100;  % assume throttle linear with voltage

    [Imot, Pelec, Pshaft, Torque, omega, ~] = motorCalc(V,Kv,I0,Rm,Imax);
    RPM = omega*60/(2*pi);

    % Motor data curves fitting with RPM as the independent variable
    % (usually a linear fit is very good)
    funMotElecPower = polyfit(RPM,Pelec,1);
    funMotShaftPower = polyfit(RPM,Pshaft,1);
    funMotTorque = polyfit(RPM,Torque,1);
    funMotCurrent = polyfit(RPM,Imot,1);

    % Search for curves intersection
    xval = linspace(min(RPM),max(RPM)); % limit the search to the available motor RPM data
    
    if any(diff(sign( polyval(funPropPower,xval) - polyval(funMotShaftPower,xval) )))
        
        xPowerMatch = fzero(@(x) polyval(funPropPower,x)-polyval(funMotShaftPower,x),...
            [min(RPM), max(RPM)]);

        c = c + 1; %#ok<*SAGROW> 
        matchingRPM(c) = xPowerMatch;
        matchingThrottle(c) = throttle;
        matchingElecPower(c) = polyval(funMotElecPower,xPowerMatch);
        matchingShaftPower(c) = polyval(funMotShaftPower,xPowerMatch);
        matchingCurrent(c) = polyval(funMotCurrent,xPowerMatch);

        figure(1)
        plot(RPM,Pshaft,'LineWidth',2, ...
            'DisplayName',[num2str(V,'%.1f'),'V, \Phi=', num2str(throttle,'%.0f'), '%'])
        plot(xPowerMatch,polyval(funPropPower,xPowerMatch),'ko','MarkerSize',6,...
            'HandleVisibility','off')


        xTorqueMatch = fzero(@(x) polyval(funPropTorque,x)-polyval(funMotTorque,x),...
            [min(RPM), max(RPM)]);

        figure(2)
        plot(RPM,Torque,'LineWidth',2, ...
            'DisplayName',[num2str(V,'%.1f'),'V, \Phi=', num2str(throttle,'%.0f'), '%'])
        plot(xTorqueMatch,polyval(funPropTorque,xTorqueMatch),'ko','MarkerSize',6,...
            'HandleVisibility','off')

    else
        disp(['No intersection found, no matching for V = ', num2str(V), ...
            ' Volt and Imax = ', num2str(Imax), ' Ampere'])
    end

end

figure(1)
grid on, hold off, legend(Location="best")
xlabel('RPM'), ylabel('Shaft Power, W')
title(filename,'Interpreter','none')

figure(2)
grid on, hold off, legend(Location="best")
xlabel('RPM'), ylabel('Torque, Nm')
title(filename,'Interpreter','none')

%% Summary figure
matchingThrust = polyval(funPropThrust,matchingRPM);

figure(3)
t = tiledlayout(4,1);
ax1 = nexttile;
plot(matchingRPM,matchingThrottle,'LineWidth',2)
ylabel('Throttle, %')
grid on

ax2 = nexttile;
plot(matchingRPM,matchingCurrent,'LineWidth',2)
ylabel('Current, A')
grid on

ax3 = nexttile; hold on
plot(matchingRPM,matchingShaftPower,'LineWidth',2)
plot(matchingRPM,matchingElecPower,'LineWidth',2)
hold off; legend('Shaft Power','Electric Power','Location','southeast')
ylabel('Power, W')
grid on

ax4 = nexttile;
plot(matchingRPM,matchingThrust,'LineWidth',2)
ylabel('Thrust, N')
grid on

linkaxes([ax1,ax2,ax3,ax4],'x');
xticklabels([ax1,ax2,ax3],{})
xlabel(t,'RPM')
title(t,filename,'interpreter','none')
t.TileSpacing = 'compact';

%% Limited power
Pmax = 1200/5;            % max electric power available per motor, Watt
yline(ax3,Pmax,'--','color',ax3.ColorOrder(2,:),'LineWidth',2,'DisplayName','Limited Power')
RPMmax = interp1(matchingElecPower,matchingRPM,Pmax);
xline(ax1,RPMmax,'--','color',ax3.ColorOrder(2,:),'LineWidth',2,'HandleVisibility','off')
xline(ax2,RPMmax,'--','color',ax3.ColorOrder(2,:),'LineWidth',2,'HandleVisibility','off')
xline(ax3,RPMmax,'--','color',ax3.ColorOrder(2,:),'LineWidth',2,'HandleVisibility','off')
xline(ax4,RPMmax,'--','color',ax3.ColorOrder(2,:),'LineWidth',2,'HandleVisibility','off')

maxUsefulThrottle = interp1(matchingRPM,matchingThrottle,RPMmax);
usefulThrottle = [matchingThrottle(matchingThrottle<maxUsefulThrottle), maxUsefulThrottle, 100];
limitedRPM = [matchingRPM(matchingThrottle<maxUsefulThrottle), RPMmax, RPMmax];
figure(4), hold on
plot(matchingThrottle,matchingRPM,'LineWidth',2)
plot(usefulThrottle,limitedRPM,'LineWidth',2)
xlabel('Throttle, %'), ylabel('RPM'), grid on
title(filename,'Interpreter','none')
legend('Unlimited Power','Limited Power','Location','northwest')