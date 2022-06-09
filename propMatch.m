% Prop Match
% Script to match user-imported propeller data with motor data
%
% This script makes use of the functions set POLYFITZERO:
% Mark Mikofski (2022). polyfitZero 
% (https://www.mathworks.com/matlabcentral/fileexchange/35401-polyfitzero), 
% MATLAB Central File Exchange. Retrieved June 9, 2022.
%
%     Copyright (C) 2022 Danilo Ciliberti danilo.ciliberti@unina.it
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
Rm = 0.034;             % motor internal resistance, Ohm

% Import propeller data from csv file. Columns order as follows:
% RPM, Power (W), Thrust (N), Torque (Nm)
propTable = readtable(['dep-prop-flowspeed-35ms.csv']);
propData = propTable.Variables;

% Propeller data curves fitting with RPM as the independent variable
% (usually a quadratic fit is good, but we need to intercept zero with a
% horizontal tangent)
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

    Imax = 1200/V/5;        % max current, Ampere (usually limited by power supply)
    % (in my case I need 5 motors on a power supply of 1200 W working at voltage V)
    % (it should be the max current sustainable by the motor, but in my
    % case constraints are different. WARNING: motor curves will be limited
    % by this Imax, and that is what I want)

    [Imot, Pmot, Pload, Qload, omega, eff] = motorCalc(V,Kv,I0,Rm,Imax);
    RPM = omega*60/(2*pi);

    % Motor data curves fitting with RPM as the independent variable
    % (usually a linear fit is very good)
    funMotPower = polyfit(RPM,Pload,1);
    funMotTorque = polyfit(RPM,Qload,1);
    funMotCurrent = polyfit(RPM,Imot,1);

    % Search for curves intersection
    xval = linspace(min(RPM),max(RPM)); % limit the search to the available motor RPM data
    
    if any(diff(sign( polyval(funPropPower,xval) - polyval(funMotPower,xval) )))
        
        xPowerMatch = fzero(@(x) polyval(funPropPower,x)-polyval(funMotPower,x),...
            [min(RPM), max(RPM)]);

        c = c + 1; %#ok<*SAGROW> 
        matchingRPM(c) = xPowerMatch;
        matchingCurrent(c) = polyval(funMotCurrent,xPowerMatch);
        matchingPower(c) = polyval(funMotPower,xPowerMatch); 
        matchingThrottle(c) = throttle;

        figure(1)
        plot(RPM,Pload,'LineWidth',2, ...
            'DisplayName',[num2str(V,'%.1f'),'V, \Phi=', num2str(throttle,'%.0f'), '%'])
        plot(xPowerMatch,polyval(funPropPower,xPowerMatch),'ko','MarkerSize',6,...
            'HandleVisibility','off')


        xTorqueMatch = fzero(@(x) polyval(funPropTorque,x)-polyval(funMotTorque,x),...
            [min(RPM), max(RPM)]);

        figure(2)
        plot(RPM,Qload,'LineWidth',2, ...
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
xlabel('RPM'), ylabel('Power, W')
title('Propeller-motor matching on shaft power')

figure(2)
grid on, hold off, legend(Location="best")
xlabel('RPM'), ylabel('Torque, Nm')
title('Propeller-motor matching on torque')

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

ax3 = nexttile;
plot(matchingRPM,matchingPower,'LineWidth',2)
ylabel('Shaft Power, W')
grid on

ax4 = nexttile;
plot(matchingRPM,matchingThrust,'LineWidth',2)
ylabel('Thrust, N')
grid on

linkaxes([ax1,ax2,ax3,ax4],'x');
xticklabels([ax1,ax2],{})
xlabel(t,'RPM')
title(t,'Propeller-motor performance')
t.TileSpacing = 'compact';