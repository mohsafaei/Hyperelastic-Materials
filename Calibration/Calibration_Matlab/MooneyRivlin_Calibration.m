close all;
clear; clc;

% -------------------------------------------------------------
% Mooney–Rivlin parameter identification from test data
%% -------------------------------------------------------------
% This script identifies the two Mooney–Rivlin constants (c1, c2) for an
% incompressible hyperelastic material by minimizing a least–squares
% objective function based on experimental stress–stretch data.
%
% Tests used:
%   - Uniaxial tension
%   - Pure shear
% (Equi–biaxial test is included as a function and can be activated.)
%
% The optimization is performed with fmincon, and the calibrated model is
% compared against the experimental data in combined plots.
%
%
% Author: Mohammad Ali Safaei
% Google Scholar: https://scholar.google.com/citations?user=jD_-4JcAAAAJ&hl
% email: ma.safaei@ut.ac.ir
%
%-------------------------------------------------------------
resp = input('Include equi-biaxial test? (1 = yes, 0 = no): ');

useEquiBiaxial = (resp == 1);   %boolean
% -------------------------------------------------------------
% Experimental data: uniaxial tension
% -------------------------------------------------------------
uniaxial_data = [ ...
    1.00   0.00;
    1.05   0.12;
    1.10   0.21;
    1.15   0.28;
    1.20   0.325;
    1.30   0.37;
    1.40   0.39;
    1.50   0.40;
    1.55   0.41;
    1.60   0.43;
    1.70   0.46;
    1.75   0.48;
    1.80   0.50;
    1.85   0.53;
    1.90   0.56;
    2.00   0.63;
    2.10   0.70;
    2.15   0.74;
    2.20   0.79;
    2.25   0.84;
    2.30   0.89;
    2.35   0.95;
    2.40   1.02;
    2.45   1.10;
    2.50   1.18;
    2.55   1.27;
    2.60   1.37;
    2.65   1.47;
    2.70   1.59;
    2.75   1.71;
    2.85   2.00;
    2.90   2.15];

A1  = uniaxial_data(:,1);
SiggSE_exp  = uniaxial_data(:,2);

%% -------------------------------------------------------------
%% Experimental data: pure shear
%% -------------------------------------------------------------
shear_data = [ ...
    1.000   0.000;
    1.020   0.040;
    1.030   0.070;
    1.040   0.090;
    1.050   0.110;
    1.070   0.160;
    1.080   0.180;
    1.090   0.210;
    1.100   0.230;
    1.110   0.250;
    1.125   0.270;
    1.140   0.290;
    1.155   0.320;
    1.170   0.340;
    1.185   0.360;
    1.200   0.370;
    1.215   0.390;
    1.225   0.410;
    1.240   0.420;
    1.260   0.450;
    1.275   0.460;
    1.300   0.480;
    1.320   0.500;
    1.340   0.515;
    1.360   0.530;
    1.375   0.540;
    1.400   0.550;
    1.420   0.560;
    1.440   0.575;
    1.460   0.590;
    1.480   0.600;
    1.500   0.610];

C1 = shear_data(:,1);
SiggP_exp = shear_data(:,2);


%{

%% 
%%A  is  principal stretches in simple Extention
filename1 = 'simpleExtention.xlsx';
A = xlsread(filename1,'D2:E114');
A1 = A(:,1);     %%Lambda 1 in Simple extension test
A2 = 1./(sqrt(A1));    %%Lambda 2 in Simple extension test
A3 = sqrt(1./A1);      %%Lambda 3 in Simple extension test
SiggSE_exp = (A(:,2))*10e6;
M_SiggS_exp = mean(SiggSE_exp);

%% 
%%B  is  principal stretches in Equi-Biaxial
filename2 = 'EquiBiaxial.xlsx';
B = xlsread(filename2,'A2:B87');
B1 = B(:,1);      %%Lambda 1 in Equi-biaxial test
B2 = B(:,1);      %%Lambda 1 in Equi-biaxial test
B3 = 1./(B1.^2);   %%Lambda 3 in Equi-biaxial test
SiggE_exp = (B(:,2))*10e6;
M_SiggE_exp = mean(SiggE_exp);

 %% 
%%C  is  principal stretches in Pure Shear
filename3 = 'PureShear.xlsx';
C = xlsread(filename3,'A2:B80');
C1 = C(:,1);    %%Lambda 1 in Pure Shear test
C2 = 1./(C1);   %%Lambda 2 in Pure Shear test
C3 = ones(79,1);   %%Lambda 3 in Pure Shear test
SiggP_exp = (C(:,2))*10e6;
M_SiggP_exp = mean(SiggP_exp);
%}

% -------------------------------------------------------------
% Objective function for Mooney–Rivlin calibration
% -------------------------------------------------------------
% The total cost TT(x) is defined as the sum of squared residuals between
% experimental and model stresses over the selected tests.
% Here: uniaxial + pure shear (equi-biaxial can be included if needed).
% -------------------------------------------------------------
if useEquiBiaxial
    % Uniaxial + pure shear + equi-biaxial
    TT = @(x) MR_Eb(B1,SiggE_exp,x(1),x(2)) + ...
                MR_PSh(C1,SiggP_exp,x(1),x(2)) + ...
                MR_Ten(A1,SiggSE_exp,x(1),x(2));
else
    % Uniaxial + pure shear only
    TT = @(x)  MR_PSh(C1,SiggP_exp,x(1),x(2)) + MR_Ten(A1,SiggSE_exp,x(1),x(2));

end



% Lower and upper bounds on the Mooney–Rivlin parameters [c1 c2]
lb = [0.01 0.001];
ub = [1.4*10e3 10e3];
% Run constrained nonlinear optimization (fmincon)
[xcte, y] = fmincon(TT,[0.1 0.1],[],[],[],[],lb,ub);
c1 = xcte(1);  c2 = xcte(2);

nnn = 'The constant arrays of the strain energy function of Mooney-Rivlin will be:';
disp(nnn);
disp('xcte');
disp(xcte);

bbb = 'The value of the RSS objective function will be:';
disp(bbb);
disp('y');
disp(y);


%% 
tension_f = 2*c1*((A1)-(1./((A1).^2)))+2*c2*((1-(1./(A1).^3)));     % First Piola Equation

SSE_SE  = sum((SiggSE_exp - tension_f).^2);
perc_RMSE = sqrt(SSE_SE/numel(SiggSE_exp))/mean(abs(SiggSE_exp));
str = sprintf('Relative RMSE = %.2f %%', perc_RMSE*100);
disp(str);

%EquiB_f = 2*c1*((B1)-(1./(B1).^5))+ 2*c2*(((B1).^3)-(1./(B1).^3));  % First Piola Equation
PureShear_f = 2*c1*(C1-(1./C1))+2*c2*((1./(C1))-(1./((C1).^3)));    % First Piola Equation


figure('Position',[100 100 1200 400]);   % adjust to taste [web:66]

% Create a tight 1x3 layout:
t = tiledlayout(1,3);
t.TileSpacing = 'tight';                 % or 'compact' [web:61][web:62]
t.Padding     = 'tight';                 % or 'compact' [web:61][web:65]

% --- Tile 1: Simple extension ---
ax1 = nexttile;
Pfit(1) = plot(A1,tension_f,'r','LineWidth',3); hold on; grid on;
PData(1) = plot(A1,SiggSE_exp,'Marker','o','MarkerFaceColor', 'b','Color','b');
xlabel('stretch');
ylabel('stress  (in Pa)');
title('Uniaxial');
legend([Pfit(1),PData(1)],'Approximation function','Uniaxial Data','Location','southeast');
ylim([0 1.2*max(SiggSE_exp)]);

% --- Tile 2: Equi-biaxial ---

% Tile 2: equi-biaxial (only if available)
if useEquiBiaxial
    ax2 = nexttile;
    Pfit(2) = plot(B1,EquiB_f,'r','LineWidth',3); hold on; grid on;
    PData(2) = plot(B1,SiggE_exp,'Marker','o','MarkerFaceColor', 'b','Color','b');
    xlabel('stretch');
    ylabel('stress  (in Pa)');
    title('Equi-biaxial');
    legend([Pfit(2),PData(2)],'Approximation function','Equi-biaxial Data','Location','southeast');
    ylim([0 1.2*max(SiggE_exp)]);
else
    % skip or use nexttile to keep layout consistent, e.g. empty tile
    nexttile; axis off;
end

% --- Tile 3: Pure shear ---
ax3 = nexttile;
Pfit(3) = plot(C1,PureShear_f,'r','LineWidth',3); hold on; grid on;
PData(3) = plot(C1,SiggP_exp,'Marker','o','MarkerFaceColor', 'b','Color','b');
xlabel('stretch');
ylabel('stress  (in Pa)');
title('Pure shear');
legend([Pfit(3),PData(3)],'Approximation function','Pure Shear Data','Location','southeast');
ylim([0 1.2*max(SiggP_exp)]);

%End of the script

% -------------------------------------------------------------
%   Helper functions for Mooney–Rivlin parameter identification
% -------------------------------------------------------------
% Each function below evaluates the sum of squared errors (SSE) between
% experimental nominal stresses and the Mooney–Rivlin model prediction
% for a specific loading mode. These SSE values are used by the main
% objective function TT in the fmincon optimization.
% -------------------------------------------------------------

function [SiggP_for] = MR_PSh(C1,SiggP_exp,x,y)   %Pure Shear test
z = numel(C1);
FFP = zeros(z,1);
TTP1 = zeros(z,1);
Sig_ex = SiggP_exp;
for i=1:z
   FFP(i) = (1/(C1(i)))*(2*x*((C1(i)^2)-1)+2*y*(1-((1/C1(i))^2)));
   TTP1(i) = (((Sig_ex(i))- FFP(i))^2);
end
SiggP_for = sum(TTP1(i));    %%True Stress formulation for Pure Shear in Pa unit
end

%% ----------------
% ----------------
% ----------------

function [SiggE_for] = MR_Eb(B1,SiggE_exp,x,y) %Equi-biaxial test
z = numel(B1);
FFT = zeros(z,1);
 MME1 = zeros(z,1);
Sig_ex = SiggE_exp;
for i= 1:z
    FFT(i) =(1/(B1(i)))*(2*x*((B1(i)^2)-(1/(B1(i)^4)))+(2*y*(((B1(i))^4)-(1/(B1(i))^2))));
    MME1(i) =  ((Sig_ex(i) - FFT(i))^2);
end 
SiggE_for = sum(MME1(i));   %%True Stress formulation for Equi-biaxial in Pa unit
end

%% ----------------
% ----------------
% ----------------

function [SiggS_for] = MR_Ten(A1,SiggS_exp,x,y)   %Tension test
z = numel(A1);
FFT = zeros(z,1);
MMS1 = zeros(z,1);
Sig_ex=SiggS_exp;

for i=1:z
   FFT(i) =  (1/(A1(i)))*(2*x*((A1(i)^2)-(1/A1(i)))+2*y*((A1(i)-(1/(A1(i))^2))));
    MMS1(i) =  ((Sig_ex(i) - FFT(i))^2);
end
SiggS_for = sum(MMS1);   %%True Stress formulation for Simple Extention in Pa unit
end
