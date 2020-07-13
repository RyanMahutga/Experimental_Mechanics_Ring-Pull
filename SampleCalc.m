close all
clear all
clc

%% Experimental data
posLoad = xlsread('sampleData.xlsx'); % position and load data

delta = posLoad(1:end,2); % displacement
F = posLoad(1:end,1); % force

%% Tissue measurements (in mm)
Ri = 12.65/2; % inner tissue radius (diameter/2)
t = 1.84; % tissue thickness
w = 7.224; % tissue width
pin_rad = 2.4; % pin radius
x0 = 10.25; % pin to pin centerline distance at start

%% Fit ranges
lam_lo_lims = [1.1 1.25]; % small-strain stretch range
PK1_lo_lims = []; % small-strain stress range (unused in this example)
lam_hi_lims = []; % large-strain stretch range (unused in this example)
PK1_hi_lims = [0.2 0.35]; % large strain stress range in MPa

%% Finding uniaxial stretch and transverse stretch 
[ lam, lam_t ] = lam_transverse( Ri,pin_rad, t, x0, delta);

strain = (lam-1); % engineering strain
stress = F./(2*w*t); % PK1 in MPa (force in N, measurements in mm)

%% Fitting small-strain and large strain moduli, and transition stretch
[SS_mod, lock_mod, lam_cross, PK1_cross, p_low, p_high, error] =...
    fitall(lam', stress,lam_lo_lims, PK1_lo_lims, lam_hi_lims, PK1_hi_lims);

eps_cross = lam_cross-1; % engineering transition strain

fprintf(strcat('Small-Strain Modulus: \n', num2str(SS_mod),'\n\n'))
fprintf(strcat('Large-Strain Modulus: \n', num2str(lock_mod),'\n\n'))
fprintf(strcat('Transition Stretch: \n', num2str(lam_cross),'\n\n'))

%% calculating fit
PK_lo_fit = polyval(p_low,lam);
PK_hi_fit = polyval(p_high,lam);

%% Plotting stress-strain curve
figure
p1=plot(lam, PK_lo_fit.*1000,'r--','LineWidth',2);
hold on
plot( lam, PK_hi_fit.*1000,'r--','LineWidth',2)
p3=plot(lam,stress.*1000,'k-','LineWidth',2);
legend([p1 p3],'Bilinear Fit','Experimental Data')
ylim([-50 500])
xlabel('Stretch, [mm/mm]')
ylabel('PK1 Stress, [kPa]')
set(gca,'FontSize',16)


