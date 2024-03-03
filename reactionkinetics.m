%% Parse the data to extract the time and intensity components

time = WednesdayGroup2Kinetics430nm{:,1};
kin370 = WednesdayGroup2Kinetics370nm{:,2};
kin430 = WednesdayGroup2Kinetics430nm{:,2};

wv = WednesdayGroup2BlankSpectrum{:,1};
specr = WednesdayGroup2BlankSpectrum{:,2};
specp = WednesdayGroup3ProductSpectrum{:,2};

%% Plot the reactant and product spectra

figure; 
hold on; 
plot(wv, specr); 
plot(wv, specp); 
xlabel('Wavelength (nm)')
ylabel('Absorbance (dimensionless)')
title('Reaction Spectra')
legend('Reactant', 'Intermediate and Product')
grid on; 
hold off; 
xlim([190, 1100]);


%% Plot the kinetic traces

figure; 
hold on; 
plot(time, kin370); 
plot(time, kin430); 
xlabel('Time (s)')
ylabel('Absorbance (dimensionless)')
title('Kinetics Spectra')
legend('Kinetics at 370 nm', 'Kinetics at 430 nm')
grid on; 
hold off; 
xlim([1.4, 2497.9]);


%% Subtract reactant component from 430 nm Data
%Figure out how intense the reactant is at 430 nm relative to 370 nm

wv430 = find(wv>=430,1); %Finds the index (row number) where the wavelength is 430 nm
wv370 = find(wv>=370,1); %Finds the index (row number) where the wavelength is 370 nm

ratio430to370 = (specr(wv430)/specr(wv370));

kin430corr = kin430-ratio430to370.*kin370; %The period means element-wise

%% Fit data to biexponentials

Roft = fittype( @(a,b,c,d,x) a.*(exp(b.*x)+c.*exp(d.*x)) );
Ioft = fittype( @(a,b,c,x) a.*(exp(b.*x)-exp(c.*x)) );
%Compare to equations 5 and 6 from the Halpern Exerpt 

%Perform a rough fit to set the initial conditions 
%This helps avoid local minima

testfit = fit(time,kin370,'exp2');

%We know lambda 2 should be larger

if testfit.b<testfit.d
    lam1 = testfit.b;
    lam2 = testfit.d;
    A = testfit.c/testfit.a;
    C = testfit.a;
else
    lam1 = testfit.d;
    lam2 = testfit.b;
    A = testfit.a/testfit.c;
    C = testfit.c;
end

%Run the fits

fit370 = fit(time,kin370,Roft,'Lower',[0,-.1,0,-.1],'Upper',[100,-0.000001,100,-.000001],'StartPoint',[C,lam1,A,lam2]);
fit430 = fit(time,kin430corr,Ioft,'Lower',[0,-.1,-.1],'Upper',[100,-0.000001,-.000001],'StartPoint',[C,lam1,lam2]);

%Plot these fits on top of the kinetic traces 
% Pay attention to error bars

figure; 
hold on; 
plot(time, kin370); 
plot(fit370,'--');
plot(time, kin430corr); 
plot(fit430,'--');
xlabel('Time (s)')
ylabel('Absorbance (dimensionless)')
title('Fitted Kinetics Spectra')
legend('Kinetics at 370 nm', '370 nm fit','Kinetics at 430 nm','430 nm fit')
grid on; 
hold off; 
xlim([1.4, 2497.9]);