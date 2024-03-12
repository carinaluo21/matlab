%% Plot the resolved spectra

%Note: All data is from Wednesday Group 1, and I simply uploaded the data
%as an array instead of pointing to a file path. 

figure; 
plot(table2array(resfreq),table2array(resspec),'LineWidth', 4); 
set(gca, 'XDir','reverse','FontSize', 20,'FontName','Times New Roman')
xlim([2070,2120])
xlabel('Frequency (Hz)','FontSize',40,'FontName','Times New Roman');
ylabel('Peak Intensity (arbitrary units)','Color', 'black','FontName', 'Times New Roman','FontSize',35);
title('1D Spectrum','FontSize',40,'FontName','Times New Roman','FontWeight', 'bold');


%% Plot the nutation curve
time = 2:2:40;
figure;
hold on
for i = 1:20
    subplot(1,20,i)
    plot(calispec(i,:),'LineWidth', 4)
    set(gca,'box','off','YColor', 'black','FontName', 'Times New Roman','FontSize',20)
    ylim([-4.5e6,4.5e6])
    xlim([1000,2000])
    set(gca,'XTick',[])
    set(gca, 'YColor', 'none');
    xlabel(time(i), 'HorizontalAlignment', 'left','Color', 'black','FontSize', 22,'FontName', 'Times New Roman')
    if i == 1
        set(gca, 'YColor', 'black','FontName', 'Times New Roman','FontSize',20);
        ylabel('Peak Intensity (arbitrary units)','Color', 'black','FontName', 'Times New Roman','FontSize',40)
        continue
    else
        set(gca,'YTick',[])
    end
end
hold off 
annotation('textbox', [0.3, 0.02, 0.4, 0.06], 'String', 'Time (\museconds)', 'HorizontalAlignment', 'center', 'EdgeColor', 'none','FontSize', 40,'FontName', 'Times New Roman');
annotation('textbox', [0.3, 0.95, 0.5, 0.05], 'String', 'Nutation Curve', 'HorizontalAlignment', 'center', 'EdgeColor', 'none','FontSize', 40,'FontName', 'Times New Roman');

%% Fit the T1 data to a function with an exponential offset

T1_time = transpose([0.0625,0.125,0.25,0.5,1,2,4,8,16,32]); 
T1_int_peak1 = transpose([-17.4,-16,-16.3,-14.1,-11.4,-6.6,0.76,9.53,15.9,17.7]); 
T1_int_peak2 = transpose([-50,-48,-45.2,-39.4,-30.9,-15.8,6.49,31.2,47.2,51.2]); 

figure;
set(gca,'FontSize', 30,'FontName','Times New Roman')
title('T1 Inversion Recovery','FontSize',40,'FontName','Times New Roman')
hold on
scatter(T1_time,T1_int_peak1, 'b','DisplayName', 'Peak 1','LineWidth',4)
scatter(T1_time,T1_int_peak2, 'r','DisplayName', 'Peak 2','LineWidth',4)

T1fittype = fittype(@(a,b,c,x) a*exp(-b*x)+c);
initial_guess_1 = [-35.17,0.1821,17.79]
initial_guess_2 = [-101.3, 0.2066, 51.03];
fit_T1_peak1 = fit(T1_time,T1_int_peak1,T1fittype, 'StartPoint',initial_guess_1);
fit_T1_peak2 = fit(T1_time,T1_int_peak2,T1fittype,'StartPoint', initial_guess_2);

legend('T1 Peak 1','T1 Peak 2','Location', 'best','FontSize',24,'FontName','Times New Roman');
xlabel('Time (seconds)','FontSize',40,'FontName','Times New Roman')
ylabel('Signal Intensity (arbitrary units)','FontSize',40,'FontName','Times New Roman')

%Add error bars to the plot
T1p1coeff = coeffvalues(fit_T1_peak1);
T1p2coeff = coeffvalues(fit_T1_peak2);

T1p1ci = confint(fit_T1_peak1);
T1p2ci = confint(fit_T1_peak2);

T1p1fit = T1fittype(T1p1coeff(1),T1p1coeff(2),T1p1coeff(3), T1_time);
T1p2fit = T1fittype(T1p2coeff(1),T1p2coeff(2),T1p2coeff(3), T1_time);

T1p1_error = abs(T1fittype(T1p1ci(2),T1p1ci(4),T1p1ci(6), T1_time) - T1fittype(T1p1ci(1),T1p1ci(3),T1p1ci(5), T1_time));
T1p2_error = abs(T1fittype(T1p2ci(2),T1p2ci(4),T1p2ci(6), T1_time) - T1fittype(T1p2ci(1),T1p2ci(3),T1p2ci(5), T1_time));

% Plot fit with error bars
errorbar(T1_time, T1p1fit, T1p1_error,'k--', 'LineWidth', 2,'DisplayName', 'T1 Peak 1 Fit');
errorbar(T1_time, T1p2fit, T1p2_error,'k--', 'LineWidth', 2,'DisplayName', 'T1 Peak 2 Fit');
hold off

%% Fit the T2 for both peaks
T2_time = transpose([0.025,0.05,0.1,0.2,0.4,0.8,1.6,3.2,6.4,12.8,25.6])
T2_int_peak1 = transpose([35.6,35,34.6,33.2,32.1,29.4,24.5,16.8,7.4,1.16,0.03])
T2_int_peak2 = transpose([101,101,99.8,97,93.3,85,70.6,48.4,22.3,4.02,0.05])

figure;
set(gca,'FontSize', 30,'FontName','Times New Roman')
title('T2 Spin Echo','FontSize',40,'FontName','Times New Roman')
hold on
scatter(T2_time,T2_int_peak1, 'b','DisplayName', 'Peak 1','LineWidth',4)
scatter(T2_time,T2_int_peak2, 'r','DisplayName', 'Peak 2','LineWidth',4)

T2fittype = fittype(@(a,b,c,x) a*exp(-b*x)+c);
initial_guess_3 = [35.77,0.2339,-0.3498]
initial_guess_4 = [102.8, 0.2314, -0.7598];
fit_T2_peak1 = fit(T2_time,T2_int_peak1,T2fittype, 'StartPoint', initial_guess_3);
fit_T2_peak2 = fit(T2_time,T2_int_peak2,T2fittype, 'StartPoint', initial_guess_4);

legend('T2 Peak 1','T2 Peak 2','Location', 'best','FontSize',24,'FontName','Times New Roman');
xlabel('Time (seconds)','FontSize',40,'FontName','Times New Roman')
ylabel('Signal Intensity (arbitrary units)','FontSize',40,'FontName','Times New Roman')

%Add error bars to the plot
T2p1coeff = coeffvalues(fit_T2_peak1);
T2p2coeff = coeffvalues(fit_T2_peak2);

T2p1ci = confint(fit_T2_peak1);
T2p2ci = confint(fit_T2_peak2);

T2p1fit = T2fittype(T2p1coeff(1),T2p1coeff(2),T2p1coeff(3), T2_time);
T2p2fit = T2fittype(T2p2coeff(1),T2p2coeff(2),T2p2coeff(3), T2_time);

T2p1_error = abs(T2fittype(T2p1ci(2),T2p1ci(4),T2p1ci(6), T2_time) - T2fittype(T2p1ci(1),T2p1ci(3),T2p1ci(5), T2_time));
T2p2_error = abs(T2fittype(T2p2ci(2),T2p2ci(4),T2p2ci(6), T2_time) - T2fittype(T2p2ci(1),T2p2ci(3),T2p2ci(5), T2_time));

% Plot fit with error bars
errorbar(T2_time, T2p1fit, T2p1_error,'k--', 'LineWidth', 2,'DisplayName', 'T2 Peak 1 Fit');
errorbar(T2_time, T2p2fit, T2p2_error,'k--', 'LineWidth', 2,'DisplayName', 'T2 Peak 2 Fit');
hold off