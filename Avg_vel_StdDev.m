%% COMPUTE STATS
% This file has been created at Cranfield University in order to
% provide supplementary materials for the Data Analysis and Uncertainty
% module assignment (N-CFD-DAPP)
%
% The script reads data obtained from direct numerical simulation of
% a low Reynolds number channel flow. It reads the excel data sheet
% named "Reynolds_stresses.xlsx" including spatio-temporal averages 
% computed from 200,000 time steps. The data sheet also contains the 
% parameters of the channel flow.
%
% Spatio-temporal averages are computed by computing the time-averaged 
% value at every cell thorughout the simulations and then computing the
% spatial averages of the time-averaged data along planes parallel to the
% wall. Spatial-averaging is denoted by <> in this script.
%
% The script also reads the 'time_samples.hdf5' file which contains time
% samples of the velocity components and the pressure along a wall-normal
% line from the bottom wall to the centre line of the channel. This data is
% processed using time-averaging symbolised by <>_t.
%
% The data is made non-dimensional based on the fluid density, the average
% streamwise velocity, and the channel half-height.
%
% For equations Kundu's book has been used (fifth edition)
%
% written by Tamás István Józsa
% 21/09/2023
clear all; close all; clc;

%% READ CHANNEL PARAMETERS AND SPATIO-TEMPORAL AVERAGES

% read parameters
% Lz, Lx, Ly, nu, Delta p
params = xlsread('Reynolds_stresses.xlsx','parameters');

Lz = params(1); Lx = params(2); Ly = params(3);
nu = params(4); % kinematic viscosity

% these values are equal to unity because they are the reference quantities
% used to make the data dimensionless
u_b = 1.0; % bulk velocity (average velocity in the entire channel)
rho = 1.0; % density
delta = Lx/2; % boundary layer thickness  =  channel half-height

% bulk Reynolds number based on channel half height and mean velocity
Re_b  =  u_b*delta/nu;

% read wall-normal coordinate and spatio-temporal averages
% x, <w>, <w'w'>, <u'u'> , <v'v'>, <u'w'>
ST_ave_dat = xlsread('Reynolds_stresses.xlsx','Reynolds_stresses');


%% READ TIME SAMPLES AT PROBES PLACED ALONG A WALL_NORMAL LINE

hinfo  =  hdf5info('time_samples.hdf5');

% sampling time
t_smpl = hdf5read(hinfo.GroupHierarchy.Datasets(1));
% wall-normal location of the samples
x_smpl = hdf5read(hinfo.GroupHierarchy.Datasets(5))+1.0;

% sampled velocity components
% each row represents a time instant as dictated by t_smpl
% each column represents a spatial location as dictated by y_smpl
w_smpl = hdf5read(hinfo.GroupHierarchy.Datasets(2));
u_smpl = hdf5read(hinfo.GroupHierarchy.Datasets(3));
v_smpl = hdf5read(hinfo.GroupHierarchy.Datasets(4));


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


% Code to get the average for one column
%format shortE; % Set the format to scientific notation
%disp([(1:size(u_smpl, 1))', u_smpl(:, 1)]);
%format;

%avg = mean(u_smpl(:, 1));
%fprintf('%.4e\n', avg);


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


cols = size(u_smpl, 2);
avg_u = zeros(1, cols); % array in u to store averages
avg_v = zeros(1, cols); % array in v to store averages
avg_w = zeros(1, cols); % array in v to store averages


stddev_u = zeros(1, cols); % array in u to store std dev
stddev_v = zeros(1, cols); % array in v to store std dev
stddev_w = zeros(1, cols); % array in v to store std dev

uplim_u = zeros(1, cols);
lwlim_u = zeros(1, cols);

uplim_v = zeros(1, cols);
lwlim_v = zeros(1, cols);

uplim_w = zeros(1, cols);
lwlim_w = zeros(1, cols);




for i = 1:cols
    avg_u(i) = mean(u_smpl(:, i));
    stddev_u(i)= std(u_smpl(:, i));
    uplim_u(i) = avg_u(i) + stddev_u(i)
    lwlim_u(i) = avg_u(i) - stddev_u(i)

    avg_v(i) = mean(v_smpl(:, i));
    stddev_v(i)= std(v_smpl(:, i));
    uplim_v(i) = avg_v(i) + stddev_v(i)
    lwlim_v(i) = avg_v(i) - stddev_v(i)

    avg_w(i) = mean(w_smpl(:, i));
    stddev_w(i)= std(w_smpl(:, i));
    uplim_w(i) = avg_w(i) + (stddev_w(i))
    lwlim_w(i) = avg_w(i) - stddev_w(i)
    
    

    
end




[max_stddev_u, idx_max_stddev_u] = max(stddev_u);
[max_stddev_v, idx_max_stddev_v] = max(stddev_v);
[max_stddev_w, idx_max_stddev_w] = max(stddev_w);


value_at_max_stddev_w = avg_w(idx_max_stddev_w);





%% instantaneous velocity plots

% % % % % % % % % % % % % % % FIGURE 1 AVG AND INSTANTANEOUS % % % % % % % % % % % % % % %




figure(1)
hold on
plot(x_smpl,u_smpl(1,:),':r',LineWidth=1)
plot(x_smpl,v_smpl(1,:),':b',LineWidth=1)
plot(x_smpl,w_smpl(1,:),':k',LineWidth=1)

plot(x_smpl,avg_u(1,:),'-r',LineWidth=2)
plot(x_smpl,avg_v(1,:),'-b',LineWidth=2)
plot(x_smpl,avg_w(1,:),'-k',LineWidth=2)


plot(x_smpl,u_smpl(100000,:),':r',LineWidth=1)
plot(x_smpl,u_smpl(200000,:),':r',LineWidth=1)

plot(x_smpl,v_smpl(100000,:),':b',LineWidth=1)
plot(x_smpl,v_smpl(200000,:),':b',LineWidth=1)

plot(x_smpl,w_smpl(100000,:),':k',LineWidth=1)
plot(x_smpl,w_smpl(200000,:),':k',LineWidth=1)


xlabel('$x/\delta$','Interpreter','latex','FontSize',12,'FontName','Times New Roman')
xlim([0,1])
legend('$u_s$','$v_s$','$w_s$', '$u_{s_{avg}}$','$v_{s_{avg}}$','$w_{s_{avg}}$',...
    'Interpreter','latex','FontSize',12,'Location','east','FontName','Times New Roman')
xticks(0:0.1:1);
yticks(-0.2:0.25:1.5);
grid on
set(gca,'TickLabelInterpreter','latex','FontSize',12)
saveas(gcf, 'Instantaneous_average.png');



% % % % % % % % % % % % % % FIGURE 2 % % % % % % % % % % % % % % % % % % % %

figure(2)
hold on
plot(x_smpl,stddev_u(1,:),'-r',LineWidth=2)
plot(x_smpl,stddev_v(1,:),'-b',LineWidth=2)
plot(x_smpl,stddev_w(1,:),'-k',LineWidth=2)

% algnu=-0.01
% algnv=-0.02
% algnw=-0.03

scatter(x_smpl(idx_max_stddev_u), max_stddev_u, 50, 'Marker', 'o', 'MarkerEdgeColor', 'r', 'LineWidth', 1);
scatter(x_smpl(idx_max_stddev_v), max_stddev_v, 50, 'Marker', 'o', 'MarkerEdgeColor', 'b', 'LineWidth', 1);
scatter(x_smpl(idx_max_stddev_w), max_stddev_w, 50, 'Marker', 'o', 'MarkerEdgeColor', 'k', 'LineWidth', 1);

line([x_smpl(idx_max_stddev_u), x_smpl(idx_max_stddev_u)], [0, max_stddev_u], 'Color', 'r', 'LineStyle', ':', 'LineWidth', 1);
line([x_smpl(idx_max_stddev_v), x_smpl(idx_max_stddev_v)], [0, max_stddev_v], 'Color', 'b', 'LineStyle', ':', 'LineWidth', 1);
line([x_smpl(idx_max_stddev_w), x_smpl(idx_max_stddev_w)], [0, max_stddev_w], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 1);



text(x_smpl(idx_max_stddev_u), max_stddev_u + 0.0025, sprintf('x/\\delta=%.2f', x_smpl(idx_max_stddev_u)), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 10, 'Color', 'k');
text(x_smpl(idx_max_stddev_v), max_stddev_v + 0.0025, sprintf('x/\\delta=%.2f', x_smpl(idx_max_stddev_v)), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 10, 'Color', 'k');
text(x_smpl(idx_max_stddev_w), max_stddev_w + 0.0025, sprintf('x/\\delta=%.2f', x_smpl(idx_max_stddev_w)), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 10, 'Color', 'k');



xlabel('$x/\delta$','Interpreter','latex','FontSize',12,'FontName','Times New Roman')
%xlim([0,1])
legend('$u_{std}$','$v_{std}$','$w_{std}$','$u_{std_{max.}}$','$v_{std_{max.}}$','$w_{std_{max.}}$',...
    'Interpreter','latex','FontSize',12,'Location','east','FontName','Times New Roman')
%xticks(0:0.1:1);
%yticks(-0.2:0.25:1.5);
grid on
set(gca,'TickLabelInterpreter','latex','FontSize',12)
saveas(gcf, 'Standard_Deviation.png');

% % % % % % % % % % % % % % FIGURE 3 % % % % % % % % % % % % % % % % % % % %





figure(4)
hold on
plot(x_smpl,avg_u(1,:),'--r',LineWidth=2)

plot(x_smpl,avg_v(1,:),'-.b',LineWidth=2)

plot(x_smpl,avg_w(1,:),'-k',LineWidth=2)



xlabel('$x/\delta$','Interpreter','latex','FontSize',12,'FontName','Times New Roman')
%xlim([0,1])
legend('$u_{s_{avg}}$','$v_{s_{avg}}$','$w_{s_{avg}}$',...
    'Interpreter','latex','FontSize',12,'Location','east','FontName','Times New Roman')
%xticks(0:0.1:1);
%yticks(-0.2:0.25:1.5);
grid on
set(gca,'TickLabelInterpreter','latex','FontSize',12)
saveas(gcf, 'Averages_Only.png');





