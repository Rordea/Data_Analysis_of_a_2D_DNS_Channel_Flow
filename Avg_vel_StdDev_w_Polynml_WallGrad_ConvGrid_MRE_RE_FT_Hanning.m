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





figure(3)
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


degree = 8; 
coefficients = polyfit(x_smpl, avg_w, degree);
disp('Coefficients of the polynomial:');
disp(coefficients);

y_fit = polyval(coefficients, x_smpl);

figure(4);
hold on;
plot(x_smpl, avg_w, '--k', 'LineWidth', 3);

plot(x_smpl, y_fit, ':','Color', [0.5, 0, 0.5], 'LineWidth', 2);
xlabel('$x/\delta$', 'Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel('$w_{{avg}}$', 'Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times New Roman');
legend('$w_{s_{avg}}$','Analytical Solution',...
    'Interpreter','latex','FontSize',12,'Location','east','FontName','Times New Roman')
legend('show', 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'east', 'FontName', 'Times New Roman');

grid on;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
saveas(gcf, 'Analytical_Solution.png');


x_arr = linspace(0, x_smpl(end), 25);

%not used
x_arr2 = linspace(0, x_smpl(end), 50);



x_arr3 = linspace(0, x_smpl(end), 150);



%not used
x_arr4 = linspace(0, x_smpl(end), 100);


x_arr5 = linspace(0, x_smpl(end), 400);

x_arr6 = linspace(0, x_smpl(end), 900);




% Add the polynomial equation as text
poly_eq = sprintf('%.4fx^8 + %.4fx^7 + %.4fx^6 + %.4fx^5 + %.4fx^4 + %.4fx^3 + %.4fx^2 + %.4fx + %.4f', ...
     coefficients(1), coefficients(2), coefficients(3), coefficients(4), ...
     coefficients(5), coefficients(6), coefficients(7), coefficients(8), coefficients(9));
text(0.5, 1.0, poly_eq, 'Interpreter', 'latex', 'FontSize', 5, 'FontName', 'Times New Roman');
% 
% ...


dx = x_smpl(end)/(length(x_arr)-1)
dx2 = x_smpl(end)/(length(x_arr2)-1)
dx3 = x_smpl(end)/(length(x_arr3)-1)
dx4 = x_smpl(end)/(length(x_arr4)-1)
dx5 = x_smpl(end)/(length(x_arr5)-1)
dx6 = x_smpl(end)/(length(x_arr6)-1)

syms x;

y=-81.8081*x^8 + 423.6784*x^7 + -925.5110*x^6 + 1109.3943*x^5 + -795.2143*x^4 + 348.2765*x^3 + -91.4918*x^2 + 13.8670*x + -0.0180;
dw_dy = gradient(y,x)

% Plot the derivative of the polynomial
figure(5);
hold on;
dy_fun = matlabFunction(dw_dy);

% Plot the derivative of the polynomial using fplot
fplot(dy_fun, [0, x_smpl(end)], 'k', 'LineWidth', 2);
title('Wall Gradient');
xlabel('$x/\delta$', 'Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel('$\frac{dy}{dx}$', 'Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times New Roman');
ylim([0, 15]);
yticks(0:2:14);
grid on;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);






x_arr = linspace(0, x_smpl(end), 25);

%not used
x_arr2 = linspace(0, x_smpl(end), 50);



x_arr3 = linspace(0, x_smpl(end), 150);



%not used
x_arr4 = linspace(0, x_smpl(end), 100);


x_arr5 = linspace(0, x_smpl(end), 400);

x_arr6 = linspace(0, x_smpl(end), 900);




% Add the polynomial equation as text
poly_eq = sprintf('%.4fx^8 + %.4fx^7 + %.4fx^6 + %.4fx^5 + %.4fx^4 + %.4fx^3 + %.4fx^2 + %.4fx + %.4f', ...
     coefficients(1), coefficients(2), coefficients(3), coefficients(4), ...
     coefficients(5), coefficients(6), coefficients(7), coefficients(8), coefficients(9));
text(0.5, 1.0, poly_eq, 'Interpreter', 'latex', 'FontSize', 5, 'FontName', 'Times New Roman');
% 
% ...


dx = x_smpl(end)/(length(x_arr)-1)
dx2 = x_smpl(end)/(length(x_arr2)-1)
dx3 = x_smpl(end)/(length(x_arr3)-1)
dx4 = x_smpl(end)/(length(x_arr4)-1)
dx5 = x_smpl(end)/(length(x_arr5)-1)
dx6 = x_smpl(end)/(length(x_arr6)-1)

syms x;

y=-81.8081*x^8 + 423.6784*x^7 + -925.5110*x^6 + 1109.3943*x^5 + -795.2143*x^4 + 348.2765*x^3 + -91.4918*x^2 + 13.8670*x + -0.0180;
dy = diff(y)

y0 = zeros(1, length(x_arr));
y2 = zeros(1, length(x_arr2));
y3 = zeros(1, length(x_arr3));
y4 = zeros(1, length(x_arr4));
y5 = zeros(1, length(x_arr5));
y6 = zeros(1, length(x_arr6));


% forward finite differences


for i = 1:length(x_arr)    
    y0(i)= subs(y, x, x_arr(i)) +  ( subs(dy, x, x_arr(i)) * dx )
end    
for i = 1:length(x_arr2)    
    y2(i)= subs(y, x, x_arr2(i)) +  ( subs(dy, x, x_arr2(i)) * dx2 )
end

for i = 1:length(x_arr3)    
    y3(i)= subs(y, x, x_arr3(i)) +  ( subs(dy, x, x_arr3(i)) * dx3 )
end

for i = 1:length(x_arr4)    
    y4(i)= subs(y, x, x_arr4(i)) +  ( subs(dy, x, x_arr4(i)) * dx4 )
end

for i = 1:length(x_arr5)    
    y5(i)= subs(y, x, x_arr5(i)) +  ( subs(dy, x, x_arr5(i)) * dx5 )
end

for i = 1:length(x_arr6)    
    y6(i)= subs(y, x, x_arr6(i)) +  ( subs(dy, x, x_arr6(i)) * dx6 )
end


disp(['Size of x_arr: ', num2str(size(x_arr))]);

disp(['Size of y1: ', num2str(size(y0))]);



figure(6);
hold on;

plot(x_smpl, y_fit, '-k', 'LineWidth', 2);
plot(x_arr, y0(1,:), ':', 'Color', [0.7, 0, 0], 'LineWidth', 2);  % Dark Red
%plot(x_arr2, y2(1,:), '--', 'Color', [0, 0.7, 0], 'LineWidth', 2);  % Dark Green
plot(x_arr3, y3(1,:), ':', 'Color', [0, 0, 0.7], 'LineWidth', 2);  % Dark Blue
%plot(x_arr4, y4(1,:), '--', 'Color', [0.7, 0.7, 0], 'LineWidth', 2);  % Dark Yellow
plot(x_arr5, y5(1,:), ':', 'Color', [0.7, 0, 0.7], 'LineWidth', 2);  % Dark Magenta
plot(x_arr6, y6(1,:), ':', 'Color', [0, 0.7, 0.7], 'LineWidth', 2);  % Dark Cyan

xlabel('$x/\delta$', 'Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel('$w_{{avg}}$', 'Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times New Roman');
legend('Analytical Solution', ...
    [ num2str(length(x_arr)), ' Points'], ...
    [ num2str(length(x_arr3)), ' Points'], ...
    [ num2str(length(x_arr5)), ' Points'], ...
    [ num2str(length(x_arr6)), ' Points'], ...
    'Interpreter', 'latex', 'FontSize', 12, 'Location', 'east', 'FontName', 'Times New Roman');
grid on;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);

grid on;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
saveas(gcf, 'Analytical_Solution_discrete.png');



yy0 = zeros(1, length(x_arr));
yy3 = zeros(1, length(x_arr3));
yy5 = zeros(1, length(x_arr5));
yy6 = zeros(1, length(x_arr6));


for i = 1:length(x_arr)    
    yy0(i)= subs(y, x, x_arr(i))
end    

for i = 1:length(x_arr3)    
    yy3(i)= subs(y, x, x_arr3(i))
end

for i = 1:length(x_arr5)    
    yy5(i)= subs(y, x, x_arr5(i)) 
end

for i = 1:length(x_arr6)    
    yy6(i)= subs(y, x, x_arr6(i)) 
end

abs_relative_errors_0 = abs((yy0 - y0) ./ y0);
abs_relative_errors_3 = abs((yy3 - y3) ./ y3);
abs_relative_errors_5 = abs((yy5 - y5) ./ y5);
abs_relative_errors_6 = abs((yy6 - y6) ./ y6);

mean_relative_error_0 = (1/length(abs_relative_errors_0)) * sum(abs_relative_errors_0) * 100;
mean_relative_error_3 = (1/length(abs_relative_errors_3)) * sum(abs_relative_errors_3) * 100;
mean_relative_error_5 = (1/length(abs_relative_errors_5)) * sum(abs_relative_errors_5) * 100;
mean_relative_error_6 = (1/length(abs_relative_errors_6)) * sum(abs_relative_errors_6) * 100;

% Create an array to store mean relative errors
all_mean_relative_errors = [];
points=[25,50,100,600]
% Append mean relative errors to the array
all_mean_relative_errors = [all_mean_relative_errors, mean_relative_error_0];
all_mean_relative_errors = [all_mean_relative_errors, mean_relative_error_3];
all_mean_relative_errors = [all_mean_relative_errors, mean_relative_error_5];
all_mean_relative_errors = [all_mean_relative_errors, mean_relative_error_6];



figure(8);
plot(points, all_mean_relative_errors, '-o', 'LineWidth', 2, 'MarkerSize', 5);

% Set labels and title
xlabel('Number of Points','Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel('MRE (%)','Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times New Roman');
title('MRE vs Grid spacing','Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times New Roman');

% Set the x-axis tick labels
%xticks(1:length(all_mean_relative_errors));
%xticklabels({[ num2str(length(x_arr))], [ num2str(length(x_arr3))], [ num2str(length(x_arr5))], [ num2str(length(x_arr6))]});

% Display the mean relative errors next to the data points
%text(1:length(all_mean_relative_errors), all_mean_relative_errors + 0.5, ...
    %sprintf('%.4f%%', all_mean_relative_errors), ...
    %'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% Show grid
grid on;

xlim([0,700])
% Set the font size
set(gca, 'FontSize', 12);

% Display the plot
saveas(gcf, 'All_Mean_Relative_Errors_Line.png');

% Calculate the refinement ratio (assuming uniform mesh refinement)
dx_coarse = x_arr(2) - x_arr(1);
dx_medium = x_arr3(2) - x_arr3(1);
dx_fine = x_arr5(2) - x_arr5(1);

r = dx_coarse / dx_medium;

% Calculate the order of convergence
P = ( log((y0(1) - y3(1)) ./ (y3(1) - y5(1))) / log(r) );
%P = ( log((y3 - y5) ./ (y0 - y3)) / log(r) );
disp(['Size of yy0: ', num2str(size(yy0))]);
disp(['Size of yy3: ', num2str(size(yy3))]);
disp(['Size of yy5: ', num2str(size(yy5))]);

fprintf('yy0(1): %.4f\n', y0(1));
fprintf('yy3(1): %.4f\n', y3(1));
fprintf('yy5(1): %.4f\n', y5(1));


% Display the result
disp(['Orders of Convergence: ', num2str(P)]);
disp(['Change of h: ', num2str(r)]);




yy_extrapolated = y5 + ( (1 ./ (r.^P - 1)) .* (y5(1) - y3(1)) );
%yy_extrapolated = y5 + ( (1 ./ (r.^P - 1)) .* (y5 - y3) );

% Plot the results
figure(9);
hold on;
plot(x_arr5, yy_extrapolated, '-x', 'LineWidth', 2, 'MarkerSize', 1);
plot(x_arr, y0, '-o', 'LineWidth', 2, 'MarkerSize', 1);
plot(x_arr3, y3, '-o', 'LineWidth', 2, 'MarkerSize', 1);
plot(x_arr5, y5, '-o', 'LineWidth', 2, 'MarkerSize', 1);
plot(x_smpl, y_fit, "-k", 'LineWidth', 2, 'MarkerSize', 1);
xlabel('$x/\delta$', 'Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel('Extrapolated Value', 'Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times New Roman');

legend('Richardson Extrapolation', ...
    [ num2str(length(x_arr)), ' Points'], ...
    [ num2str(length(x_arr3)), ' Points'], ...
    [ num2str(length(x_arr5)), ' Points'], ...
    "Analytical Solution", ...
    'Interpreter', 'latex', 'FontSize', 12, 'Location', 'east', 'FontName', 'Times New Roman');

title('Richardson Extrapolation', 'Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times New Roman');
grid on;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);


saveas(gcf, 'RichardsonExtrapolation.png');


x_target = 0.0610;
%disp(min(abs(x_smpl - x_target)))
[~, idx_x_target] = min(abs(x_smpl - x_target)); %looks for the closest value to the x number specified.

% Extract the wall gradient values at x = 0.06
w_smpl_x06 = w_smpl( : , idx_x_target ) ;

dw_dy_sampl_x06 = gradient(w_smpl_x06,x_target)








figure(10);
hold on;
plot(t_smpl, w_smpl_x06, ':', 'LineWidth', 2);

% Set labels and title
xlabel('Time sample [s]','Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel('Instantaneous w vel ','Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times New Roman');
title('Instantaneous w velocity vs Time','Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times New Roman');
% Show grid
grid on;
% Set the font size
set(gca, 'FontSize', 12);

% Display the plot
saveas(gcf, 'x0_06_sample.png');




figure(11);
hold on;
plot(t_smpl, dw_dy_sampl_x06, ':r', 'LineWidth', 2);

% Set labels and title
xlabel('Time sample [s]','Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel('Wall Gradient ','Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times New Roman');
title('Instantaneous w velocity vs Time','Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times New Roman');
% Show grid
grid on;
% Set the font size
set(gca, 'FontSize', 12);

% Display the plot
saveas(gcf, 'wallgrad_x0_06_sample.png');




num_timesteps = length(t_smpl);
averaging_interval = 20;

% Initialize variables to store averaged values and corresponding timesteps
avg_dw_dy = zeros(floor(num_timesteps/averaging_interval), 1);
avg_timesteps = zeros(floor(num_timesteps/averaging_interval), 1);

% Calculate the averaged values and corresponding timesteps
for i = 1:floor(num_timesteps/averaging_interval)
    start_idx = (i-1) * averaging_interval + 1;
    end_idx = i * averaging_interval;
    
    % Extract the wall gradient values within the specified interval
    dw_dy_interval = dw_dy_sampl_x06(start_idx:end_idx);
    
    % Calculate the average and store it in the result vector
    avg_dw_dy(i) = mean(dw_dy_interval);
    
    % Save the corresponding timestep
    avg_timesteps(i) = mean(t_smpl(start_idx:end_idx));
end






% Plot the averaged wall gradient values against corresponding timesteps
figure(12);
hold on;
plot(avg_timesteps, avg_dw_dy, ':r', 'LineWidth', 2 );

% Set labels and title
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel('Average Wall Gradient', 'Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times New Roman');
title('Average Wall Gradient vs Time', 'Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times New Roman');

% Show grid
grid on;

% Set the font size
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);

% Display the plot
saveas(gcf, 'Average_Wall_Gradient_vs_Time.png');

disp(['Size : ', num2str(size(t_smpl))]);
disp(['Size : ', num2str(size(w_smpl_x06))]);
disp(['Size : ', num2str(size(avg_timesteps))]);
disp(['Size : ', num2str(size(avg_dw_dy))]);


% Assuming w_smpl_x06 is the signal you want to analyze
N = length(w_smpl_x06);  % Length of the signal
Fs = 1 / (t_smpl(2) - t_smpl(1));  % Sampling frequency

% Perform FFT
fft_result = fft(dw_dy_sampl_x06);
frequencies = Fs*(0:(N/2))/N;  % Frequency axis

% Compute the one-sided spectrum
P1 = abs(fft_result(1:N/2+1));





% Assuming w_smpl_x06 is the signal you want to analyze
N = length(avg_dw_dy);  % Length of the signal
Fs = 1 / (t_smpl(2) - t_smpl(1));  % Sampling frequency

% Apply Hanning window
hanning_window = hanning(N);
windowed_signal = avg_dw_dy .* hanning_window';

% Perform FFT on the windowed signal
fft_result = fft(windowed_signal);
frequencies = Fs * (0:(N/2)) / N;  % Frequency axis

% Compute the one-sided spectrum
P1 = abs(fft_result(1:N/2+1));

amp=2/N * P1

% Reference amplitude (you can adjust this based on your preference)
reference_amplitude = 1.0;

% Convert amplitude to decibels
amplitude_dB = 20 * log10(amp / reference_amplitude);

% Plot the results
figure(13);
hold on;


%subplot(2,1,2);
plot(frequencies, amplitude_dB, 'LineWidth', 2);
xlabel('Frequency [Hz]', 'Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel('Amplitude Spectrum', 'Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times New Roman');
%title('Amplitude Spectrum', 'Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times New Roman');
grid on;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);

sgtitle(' Hanning Window ', 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');

% Display the plot
saveas(gcf, 'Hanning_Window_Effect_Dbs.png');



figure(14);
% subplot(2,1,1);
% plot(avg_timesteps, avg_dw_dy, 'LineWidth', 2);
% xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times New Roman');
% ylabel('Original Signal', 'Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times New Roman');
% title('Original Signal vs Time', 'Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times New Roman');
% grid on;
% set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);

%subplot(2,1,2);
plot(frequencies, amp, 'LineWidth', 2);
xlabel('Frequency [Hz]', 'Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel('Amplitude Spectrum', 'Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times New Roman');
%title('Amplitude Spectrum', 'Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times New Roman');
grid on;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);

sgtitle(' Hanning Window ', 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');

% Display the plot
saveas(gcf, 'Hanning_Window_Effect.png');





