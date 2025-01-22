% MOTOR-PROPELLER MATCH WITH APC PROPELLER
% Script to match the performance of an electric motor with an APC
% propeller, imported from the file released on the APC propeller website,
% provided at constant RPM and variable airspeed.
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
Vmax = 44.0;            % max voltage (limited by the motor)
Kv = 200;               % motor RPM/Volt constant
I0ref = 2.37;           % motor idle (no load) current at reference voltage, A
Vref = 20.0;            % reference voltage for the no-load current, volt
Imax = 95;              % max current, A
Rm = 0.0084;            % motor internal resistance, Ohm
velArray = 60:20:200;   % expected airspeed operating range in km/h

addpath('apc')
[filename, pathname] = uigetfile({'*.dat'},'Select APC propeller performance file');
apcPropPerfData = readAPCperf([pathname,filename]);

%% Data interpolation initialization
dataset = size(apcPropPerfData,2);

propRPM = zeros(dataset,1);
for idx = 1:dataset
    propRPM(idx) = unique(apcPropPerfData{idx}.RPM);
end

% Concatenate all datasets vertically and add RPM as final column
allDataset = [];
for idx = 1:dataset
    temp = table2array(apcPropPerfData{idx});
    allDataset = [allDataset; temp]; %#ok<*AGROW>
end

% Enable a correspondence between variable names and their positions
header = apcPropPerfData{1}.Properties.VariableNames;

% Use headFun function to link numeric array position to table variable
% names (and hopefully improve code readability and avoid errors)

%% Interpolate among propeller curves at several airspeeds
% APC propellers data are given at constant RPM and variable speed, while
% electric motors data are provided at variable RPM (they cannot directly
% sense the airspeed). Thus, a preliminary interpolation must be done on
% the propeller data at different airspeeds, sweeping the available data at
% constant RPM.

queryRPM = linspace(0,max(allDataset(:,headFun(header,'RPM')))); % desired RPM array

xVector = allDataset(:,headFun(header,'V_kph'));
yVector = allDataset(:,headFun(header,'RPM'));

[xNorm, Cx, Sx] = normalize(xVector);
[yNorm, Cy, Sy] = normalize(yVector);

funPropPower = cell(size(velArray));
funPropThrust = cell(size(velArray));
funPropTorque = cell(size(velArray));

f1 = figure(1);
f2 = figure(2);
f1.Visible = 'off';
f2.Visible = 'off';

cProp = 0; % counter
for queryVelocity = velArray % cycle over velocity array

    % Thrust interpolant
    vVector = allDataset(:,headFun(header,'T_N'));
    [vNorm, Cv, Sv] = normalize(vVector);
    fThrust = scatteredInterpolant(xNorm, yNorm, vNorm, 'linear', 'none');

    % Interpolate thrust among velocity and RPM
    thrust = fThrust(repmat((queryVelocity-Cx)/Sx,1,100),(queryRPM-Cy)/Sy);
    posDataIndex = thrust*Sv+Cv > 0; % do not show data with negative thrust
    interpThrust = [queryRPM(posDataIndex)', thrust(posDataIndex)'*Sv+Cv'];


    % Torque interpolant
    vVector = allDataset(:,headFun(header,'Q_Nm'));
    [vNorm, Cv, Sv] = normalize(vVector);
    fTorque = scatteredInterpolant(xNorm, yNorm, vNorm, 'linear', 'none');

    % Interpolate torque among velocity and RPM
    torque = fTorque(repmat((queryVelocity-Cx)/Sx,1,100),(queryRPM-Cy)/Sy);
    posDataIndex = torque*Sv+Cv > 0; % do not show data with negative thrust
    interpTorque = [queryRPM(posDataIndex)', torque(posDataIndex)'*Sv+Cv];


    % Power interpolant
    vVector = allDataset(:,headFun(header,'P_W'));
    [vNorm, Cv, Sv] = normalize(vVector);
    fPower = scatteredInterpolant(xNorm, yNorm, vNorm, 'linear', 'none');

    % Interpolate power among velocity and RPM
    power = fPower(repmat((queryVelocity-Cx)/Sx,1,100),(queryRPM-Cy)/Sy);
    posDataIndex = power*Sv+Cv > 0; % do not show data with negative thrust
    interpPower = [queryRPM(posDataIndex)', power(posDataIndex)'*Sv+Cv];


    % Propeller data curves fitting with RPM as the independent variable
    % (usually a quadratic fit is good, but we need to intercept zero with a
    % horizontal tangent).
    cProp = cProp + 1;
    %     funPropPower{cProp} = polyfitB0(interpPower(:,1),interpPower(:,2),2,0);
    %     funPropThrust{cProp} = polyfitB0(interpThrust(:,1),interpThrust(:,2),2,0);
    %     funPropTorque{cProp} = polyfitB0(interpTorque(:,1),interpTorque(:,2),2,0);
    %     funPropPower{cProp} = polyfit(interpPower(:,1),interpPower(:,2),2);
    %     funPropThrust{cProp} = polyfit(interpThrust(:,1),interpThrust(:,2),2);
    %     funPropTorque{cProp} = polyfit(interpTorque(:,1),interpTorque(:,2),2);
    funPropPower{cProp} = @(rpm) interp1(interpPower(:,1),interpPower(:,2),rpm,'linear','extrap');
    funPropThrust{cProp} = @(rpm) interp1(interpThrust(:,1),interpThrust(:,2),rpm,'linear','extrap');
    funPropTorque{cProp} = @(rpm) interp1(interpTorque(:,1),interpTorque(:,2),rpm,'linear','extrap');

    % Erase propeller data with decreasing power and torque with RPM
    yPropPower = funPropPower{cProp}(queryRPM);
    yPropTorque = funPropTorque{cProp}(queryRPM);
    posDataIndex = diff(yPropPower) > 0 & diff(yPropTorque) > 0;
    yPropPower = yPropPower(posDataIndex);
    yPropTorque = yPropTorque(posDataIndex);
    queryRPMgood = queryRPM(posDataIndex);

    set(0,'CurrentFigure',f1), hold on
%     yProp = polyval(funPropPower{cProp},queryRPM);    % y array of propeller power
    plot(queryRPMgood,yPropPower,'LineWidth',2,...
        'DisplayName',['Prop at ', num2str(queryVelocity), ' km/h'])

    set(0,'CurrentFigure',f2), hold on
%     yProp = polyval(funPropTorque{cProp},queryRPM);    % y array of propeller torque
    plot(queryRPMgood,yPropTorque,'LineWidth',2,...
        'DisplayName',['Prop at ', num2str(queryVelocity), ' km/h'])

end

%% Cycle over several throttle values
cMot = 0;  % motor-prop matching counter
for phi = 0:0.1:1

    V = Vmax*phi;
    throttle = phi*100;

    [Imot, motPelec, motPshaft, motTorque, motOmega, ~] = ...
    motorCalc(Vmax,Kv,I0ref,Vref,Rm,phi,Imax);
    motRPM = motOmega*60/(2*pi);

    % Motor data curves fitting with RPM as the independent variable
    % (usually a linear fit is very good)
    [funMotElecPower, ~, muElecPower] = polyfit(motRPM,motPelec,1);
    [funMotShaftPower, ~, muShaftPower] = polyfit(motRPM,motPshaft,1);
    [funMotTorque, ~, muTorque] = polyfit(motRPM,motTorque,1);
    [funMotCurrent, ~, muCurrent] = polyfit(motRPM,Imot,1);

    set(0,'CurrentFigure',f1)
    plot(motRPM,motPshaft,'--','LineWidth',2, ...
        'DisplayName',[num2str(V,'%.1f'),'V, \Phi=', num2str(throttle,'%.0f'), '%'])

    set(0,'CurrentFigure',f2)
    plot(motRPM,motTorque,'--','LineWidth',2, ...
        'DisplayName',[num2str(V,'%.1f'),'V, \Phi=', num2str(throttle,'%.0f'), '%'])

    % Search for curves intersection
    xval = linspace(min(motRPM),max(motRPM)); % limit the search to the available motor RPM data

%     cMot = 0; % motor-prop matching counter (reset for each voltage)
    for cProp = 1:length(velArray)
        if any(diff(sign( funPropPower{cProp}(xval) - polyval(funMotShaftPower,xval,[],muShaftPower) )))     % any(diff(sign( polyval(funPropPower{cProp},xval) - polyval(funMotShaftPower,xval) )))

            cMot = cMot + 1; %#ok<*SAGROW>

%             xPowerMatch = fzero(@(x) polyval(funPropPower{cProp},x)-polyval(funMotShaftPower,x),...
%                 [min(motRPM), max(motRPM)]);
            xPowerMatch = fzero(@(x) funPropPower{cProp}(x)-polyval(funMotShaftPower,x,[],muShaftPower),...
                [min(motRPM), max(motRPM)]);

            matchingRPM = xPowerMatch;
            matchingThrottle = throttle;
            matchingElecPower = polyval(funMotElecPower,xPowerMatch,[],muElecPower);
            matchingShaftPower = polyval(funMotShaftPower,xPowerMatch,[],muShaftPower);
            matchingCurrent = polyval(funMotCurrent,xPowerMatch,[],muCurrent);

            set(0,'CurrentFigure',f1)
%             plot(xPowerMatch,funPropPower{cProp},xPowerMatch),...
%                 'ko','MarkerSize',6,'HandleVisibility','off')
            plot(xPowerMatch,funPropPower{cProp}(xPowerMatch),...
                'ko','MarkerSize',6,'HandleVisibility','off')


%             xTorqueMatch = fzero(@(x) polyval(funPropTorque{cProp},x)-polyval(funMotTorque,x),...
%                 [min(motRPM), max(motRPM)]);
            xTorqueMatch = fzero(@(x) funPropTorque{cProp}(x)-polyval(funMotTorque,x,[],muTorque),...
                [min(motRPM), max(motRPM)]);

            set(0,'CurrentFigure',f2)
%             plot(xTorqueMatch,polyval(funPropTorque{cProp},xTorqueMatch),...
%                 'ko','MarkerSize',6,'HandleVisibility','off')
            plot(xTorqueMatch,funPropTorque{cProp}(xTorqueMatch),...
                'ko','MarkerSize',6,'HandleVisibility','off')

            % Collect throttle, airspeed, RPM, power, and thrust into an array
            xPowerMatchArray(cMot,:) = [throttle,velArray(cProp),...
                xPowerMatch,funPropPower{cProp}(xPowerMatch),funPropThrust{cProp}(xPowerMatch)];

        else
            disp(['No intersection found for motor at ', num2str(V), ...
                ' volt and propeller airspeed ', num2str(velArray(cProp)), ' m/s'])
        end % if-else

    end % for-loop sweeping airspeed array

end % for-loop sweeping voltage array

% Check if no intersection between curves has been found
if cMot == 0
    error('No intersection found between motor and propeller curves')
end

% Erase data with negative power and thrust
posDataIndex = xPowerMatchArray(:,4) > 0 & xPowerMatchArray(:,5) > 0;
xPowerMatchArray = xPowerMatchArray(posDataIndex,:);

figure(1)
ylim([0,inf]), grid on, hold off, legend(Location="best")
xlabel('RPM'), ylabel('Shaft Power, W')
title(filename,'Interpreter','none')

figure(2)
ylim([0,inf]), grid on, hold off, legend(Location="best")
xlabel('RPM'), ylabel('Torque, Nm')
title(filename,'Interpreter','none')

%% Plotting vs airspeed

% Distinguish each throttle value
for idx = 1:length(xPowerMatchArray)
    sepIdx = find(diff(xPowerMatchArray(:,1)));
end

% Make different curves for different throttle values (power vs airspeed)
figure(3), hold on
plot(xPowerMatchArray(1:sepIdx(1),2), ...
    xPowerMatchArray(1:sepIdx(1),4), ...
    'LineWidth',2,...
    'DisplayName',['\Phi=', num2str(xPowerMatchArray(sepIdx(1),1),'%.0f'), '%'])
for idx = 2:length(sepIdx)
    plot(xPowerMatchArray(sepIdx(idx-1)+1:sepIdx(idx),2), ...
        xPowerMatchArray(sepIdx(idx-1)+1:sepIdx(idx),4), ...
        'LineWidth',2, ...
        'DisplayName',['\Phi=', num2str(xPowerMatchArray(sepIdx(idx),1),'%.0f'), '%'])
end
plot(xPowerMatchArray(sepIdx(end)+1:end,2), ...
    xPowerMatchArray(sepIdx(end)+1:end,4), ...
    'LineWidth',2, ...
    'DisplayName',['\Phi=', num2str(xPowerMatchArray(end,1),'%.0f'), '%'])
hold off, grid on, legend(Location="best")
xlabel('Velocity, km/h'), ylabel('Shaft power, W')
title('Coupled propeller-motor power')


% Make different curves for different throttle values (thrust vs airspeed)
figure(4), hold on
plot(xPowerMatchArray(1:sepIdx(1),2), ...
    xPowerMatchArray(1:sepIdx(1),5), ...
    'LineWidth',2,...
    'DisplayName',['\Phi=', num2str(xPowerMatchArray(sepIdx(1),1),'%.0f'), '%'])
for idx = 2:length(sepIdx)
    plot(xPowerMatchArray(sepIdx(idx-1)+1:sepIdx(idx),2), ...
        xPowerMatchArray(sepIdx(idx-1)+1:sepIdx(idx),5), ...
        'LineWidth',2, ...
        'DisplayName',['\Phi=', num2str(xPowerMatchArray(sepIdx(idx),1),'%.0f'), '%'])
end
plot(xPowerMatchArray(sepIdx(end)+1:end,2), ...
    xPowerMatchArray(sepIdx(end)+1:end,5), ...
    'LineWidth',2, ...
    'DisplayName',['\Phi=', num2str(xPowerMatchArray(end,1),'%.0f'), '%'])
hold off, grid on, legend(Location="best")
xlabel('Velocity, km/h'), ylabel('Thrust, N')
title('Coupled propeller-motor thrust')


return
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