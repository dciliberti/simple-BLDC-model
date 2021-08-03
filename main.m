close all; clearvars; clc

%% Input data
V = 36.0;               % Voltage
Pmax = 300;             % max shaft power, W (use also for axis limits)

% Motor constants
Kv = 300;               % Motor RPM/Volt constant
I0 = 1.8;               % Motor idle (no load) current, A
Rm = 0.032;             % Motor internal resistance, Ohm

%% Calculate and plot motor performance

[Imot, Pmot, Pload, Qload, omega, eff] = motorCalc(V,Kv,I0,Rm,Pmax);
RPM = omega*60/(2*pi);

t = tiledlayout(3,1);
ax1 = nexttile;
yyaxis left
plot(Imot,Pmot,'LineWidth',2)
ylabel('Electric Power, W')
%
yyaxis right
plot(Imot,Pload,'LineWidth',2)
ylabel('Shaft Power, W')
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
title(t,['Motor performance at ',num2str(V,'%.1f'),' V'])
t.TileSpacing = 'compact';