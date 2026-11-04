close all;
clear; clc;

% -------------------------------------------------------------
% Yeoh parameter calibration from test data
%% -------------------------------------------------------------
% This script identifies the three Yeoh constants (c10, c20, c30) for an
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
disp('--- YEOH MODEL PARAMETER CALIBRATION ---');
disp('**Yeoh hyperelastic model** parameter identification started...\n');

resp = input('Include equi-biaxial test? (1 = yes, 0 = no): ');
useEquiBiaxial = (resp == 1);   %boolean

res2 = input('Include pure hear test? (1 = yes, 0 = no): ');
usePureShear = (res2 == 1);   %boolean
% -------------------------------------------------------------
% Experimental data: uniaxial tension
% -------------------------------------------------------------
uniaxial_data = xlsread('Benam_Uniaxial.xlsx','A1:B70');



A1  = uniaxial_data(:,1);
A2 = sqrt(1./A1);    %%Lambda 2 in Simple extension test
A3 = sqrt(1./A1);    %%Lambda 3 in Simple extension test
disp(numel(A1));
IS1 = zeros(numel(A1),1);
for i= 1:numel(A1)
    IS1(i) = ((A1(i,1))^2)+((A2(i,1))^2)+((A3(i,1))^2);
end
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



SiggP_exp = shear_data(:,2);
C1 = shear_data(:,1);
C3 = 1./(C1);      %%Lambda 2 in Pure Shear test
C2 = ones(79,1);   %%Lambda 3 in Pure Shear test
IP1 = zeros(numel(C1),1);
for i= 1:numel(C1)
    IP1(i) = ((C1(i,1))^2)+((C2(i,1))^2)+((C3(i,1))^2);
end


%{
%% 
%%A  is  pricipal stretches in simple Extention
filename1 = 'simpleExtention.xlsx';
A = xlsread(filename1,'D2:E114');
A1 = A(:,1);         %%Lambda 1 in Simple extension test
A2 = sqrt(1./A1);    %%Lambda 2 in Simple extension test
A3 = sqrt(1./A1);    %%Lambda 3 in Simple extension test
zz1 = numel(A1);
for i= 1:zz1
IS1(i) = ((A1(i,1))^2)+((A2(i,1))^2)+((A3(i,1))^2);
end
IS1;     %%the first invariant of the Simple Extension stretches
SiggSE_exp = (A(:,2))*10e6;
M_SiggSE_exp = mean(SiggSE_exp);


%%B  is  pricipal stretches in Equi-Biaxial
filename2 = 'EquiBiaxial.xlsx';
sheet = 1;
B = xlsread(filename2,sheet,'A2:B87');
B1 = B(:,1);     %%Lambda 1 in equi-biaxial test
B2 = B(:,1);     %%Lambda 2 in equi-biaxial test
B3 = 1./((B1).^2);      %%Lambda 3 in equi-biaxial test
zz2 = numel(B1);
 
for i= 1:zz2
IE1(i) = ((B1(i,1))^2)+((B2(i,1))^2)+((B3(i,1))^2);
end
IE1;       %%the first invariant of the equi-biaxial stretches
SiggE_exp = (B(:,2))*10e6;
M_SiggE_exp = mean(SiggE_exp);

 
%%C  is  pricipal stretches in Pure Shear
filename3 = 'PureShear.xlsx';
C = xlsread(filename3,'A2:B80');
C1 = C(:,1);      %%Lambda 1 in Pure Shear test
C2 = 1./(C1);      %%Lambda 2 in Pure Shear test
C3 = ones(79,1);   %%Lambda 3 in Pure Shear test
zz3 = numel(C1);
for i= 1:zz3
IP1(i) = ((C1(i,1))^2)+((C2(i,1))^2)+((C3(i,1))^2);
end
IP1;       %%the first invariant of the Pure Shear stretches
SiggP_exp = (C(:,2))*10e6;
M_SiggP_exp = mean(SiggP_exp);
%}

%%
% -------------------------------------------------------------
% Objective function for Yeoh calibration
% -------------------------------------------------------------
% The total cost TT(x) is defined as the sum of squared residuals between
% experimental and model stresses over the selected tests.
% Here: uniaxial + pure shear (equi-biaxial can be included if needed).
% -------------------------------------------------------------
if useEquiBiaxial
    % --- Outer IF: Equi-biaxial included ---
    
    if usePureShear
        % Uniaxial + pure shear + equi-biaxial
        TT = @(x) Yeoh_Uniaxial(A1,A3,SiggSE_exp,IS1,x(1),x(2),x(3))+ ...
            Yeoh_EquiBiaxial(B1,B3,IE1,SiggE_exp,x(1),x(2),x(3))+ ...
            Yeoh_PureShear(C1,C3,IP1,SiggP_exp,x(1),x(2),x(3));
    else
        % Uniaxial + pure shear only
        TT = @(x) Yeoh_Uniaxial(A1,A3,SiggSE_exp,IS1,x(1),x(2),x(3))+ ...
                  Yeoh_EquiBiaxial(B1,B3,IE1,SiggE_exp,x(1),x(2),x(3));
    end

else
    % Outer ELSE: No equi-biaxial

    if usePureShear
        % Only uniaxial + pure shear
        TT = @(x) Yeoh_Uniaxial(A1, A3, SiggSE_exp, IS1, x(1), x(2), x(3)) + ...
                  Yeoh_PureShear(C1, C3, IP1, SiggP_exp, x(1), x(2), x(3));
    else
        % Only uniaxial
        TT = @(x) Yeoh_Uniaxial(A1,A3,SiggSE_exp,IS1,x(1),x(2),x(3))
    end
end

format longE     % scientific notation with more precision


%lb= [0.1*10e5 -0.1*10e5 0.1*10e5];
%ub =[5*10e5 -1000 5*10e5];

lb = [0 -0.052 -1];
ub = [0.243 0.5 0.01];



% Calculate a midpoint for the initial guess
%x0 = (lb + ub) / 2;
x0 = [2.0e-01    -4.7e-02     8.5e-03];
%[xcte, y] = fmincon(TT,[0.1 -0.1 0.1],[-1 1 -1],0,[],[],lb,ub);
options = optimoptions('fmincon', ...
    'Algorithm','sqp', ...                       % often better for bound-only problems
    'StepTolerance', 1e-9, ...
    'OptimalityTolerance', 1e-9, ...
    'FunctionTolerance', 1e-9, ...
    'FiniteDifferenceType','central', ...        % more accurate gradients
    'FiniteDifferenceStepSize', [1e-4 1e-4 1e-4], ... % scale to your variable magnitudes
    'ScaleProblem',true, ...                     % automatic scaling
    'MaxIterations', 2000, ...
    'MaxFunctionEvaluations', 1e5, ...
    'Display','iter-detailed', ...
    'PlotFcn', {'optimplotfval','optimplotx','optimplotfirstorderopt'});

[xcte, y] = fmincon(TT,x0,[],[],[],[],lb,ub,[],options);
c10 = xcte(1);  c20 = xcte(2); c30 = xcte(3);

disp('The constant arrays of the strain energy function of Yeoh will be:');
xcte
disp('The value of the RSS objective function will be:');
y

%%
Uniaxial_fun = zeros(numel(A1),1);
for i=1:numel(A1)
 Uniaxial_fun(i) = (1/(A1(i)))*2*(((A1(i))^2)-(A3(i))^2)*(c10+(2*c20*((IS1(i))-3))+(3*c30*(IS1(i)-3)^2));
end

figure('Position',[100 100 1200 400]);   % adjust to taste [web:66]
% Create a tight 1x3 layout:
t = tiledlayout(1,3);
t.TileSpacing = 'tight';                 % or 'compact' [web:61][web:62]
t.Padding     = 'tight';                 % or 'compact' [web:61][web:65]

% --- Tile 1: Simple extension ---
ax1 = nexttile;
nnn1 = plot(A1,Uniaxial_fun,'Color','r','marker','*');
xlabel('stretch');
ylabel('stress  in Pa');
title('Simple extention fitting with RRS Optimazation');
hold on
grid on
nnn2 = plot(A1,SiggSE_exp,'Color','b');
ffff1 = legend([nnn1,nnn2],'Approximation function','Experiment data');
ylim([0 1.2*max(SiggSE_exp)]);

hold off
%
if useEquiBiaxial
    EquiB_fun = zeros(zz2,1);
    for i= 1:zz2
        EquiB_fun(i) = (1/(B1(i)))*2*((B1(i)^2)-(B3(i)^2))*((c10+(2*c20*(IE1(i)-3))+(3*c30*(IE1(i)-3)^2))) ; 
    end
    ax2 = nexttile;

    nnn3 = plot(B1,EquiB_fun,'Color','r','marker','*');
    xlabel('stretch');
    ylabel('stress  in Pa');
    title('Equi-biaxial tension fitting with RRS Optimazation');
    hold on
    grid on
    nnn4 = plot(B1,SiggE_exp,'Color','b');
    ffff2 = legend([nnn3,nnn4],'Approximation function','Experiment data');
    hold off
else
    % skip or use nexttile to keep layout consistent, e.g. empty tile
    nexttile; axis off;
end
    
    %

if usePureShear

    ax3 = nexttile;
    PureShear_fun = zeros(numel(C1),1);
    for i= 1:numel(C1)
       PureShear_fun(i) = (1/(C1(i)))*2*((C1(i)^2)-(C3(i)^2))*((c10+(2*c20*(IP1(i)-3))+(3*c30*(IP1(i)-3)^2))) ; 
    end
    nnn5 = plot(C1,PureShear_fun,'Color','r','marker','*');
    xlabel('stretch');
    ylabel('stress  in Pa');
    title('Pure Shear tension fitting with RRS Optimazation');
    hold on
    grid on
    nnn6 = plot(C1,SiggP_exp,'Color','b');
    ffff3 = legend([nnn5,nnn6],'Approximation function','Experiment data');
    hold off

end
%
%End of the script
%

% -------------------------------------------------------------
%   Helper functions for Mooney–Rivlin parameter identification
% -------------------------------------------------------------
% Each function below evaluates the sum of squared errors (SSE) between
% experimental nominal stresses and the Mooney–Rivlin model prediction
% for a specific loading mode. These SSE values are used by the main
% objective function TT in the fmincon optimization.
% -------------------------------------------------------------

function [SiggS_for] = Yeoh_Uniaxial(A1,A3,IS1,SiggSE_exp,x,y,z)
zz = numel(A1);
FFT=zeros(zz,1);
TTS1=zeros(zz,1);
Sig_ex=SiggSE_exp;
%M_Sig_ex = mean(Sig_ex);
numfor = size(A1,1);
for i = 1:numfor  
    FFT(i) =  (1/(A1(i)))*2*(((A1(i))^2)-(A3(i))^2)*(x+(2*y*((IS1(i))-3))+(3*z*(IS1(i)-3)^2));         
    TTS1(i) = (Sig_ex(i) - FFT(i))^2;
end
SiggS_for = sum(TTS1(i));   %%True Stress formulation for Simple Extention in Pa unit
end

%% ----------------
% ----------------
% ----------------
function [SiggE_for] = Yeoh_EquiBiaxial(B1,B3,IE1,SiggE_exp,x,y,z)
zz = numel(B1);
FFT = zeros(zz,1);
TTS1 = zeros(zz,1);
Sig_ex = SiggE_exp;
%M_Sig_ex = mean(Sig_ex);
for i = 1:zz 
    FFT(i) = (1/(B1(i)))*2*((B1(i)^2)-(B3(i)^2))*((x+(2*y*(IE1(i)-3))+(3*z*(IE1(i)-3)^2))) ;         
    TTS1(i) = (Sig_ex(i) - FFT(i))^2;
end
SiggE_for = sum(TTS1(i));   %%True Stress formulation for Equi-biaxial Extention in Pa unit
end

%% ----------------
% ----------------
% ----------------

function [SiggP_for] = Yeoh_PureShear(C1,C3,IP1,SiggP_exp,x,y,z)
zz = numel(C1);
FFT = zeros(zz,1);
TTS1 = zeros(zz,1);
Sig_ex = SiggP_exp;
%M_Sig_ex = mean(Sig_ex);
for i = 1:zz 
    FFT(i) = (1/(C1(i)))*2*((C1(i)^2)-(C3(i)^2))*((x+(2*y*(IP1(i)-3))+(3*z*(IP1(i)-3)^2))) ;         
    TTS1(i) = (Sig_ex(i) - FFT(i))^2;
end
SiggP_for = sum(TTS1(i));   %%True Stress formulation for Pure Shear in Pa unit
end

