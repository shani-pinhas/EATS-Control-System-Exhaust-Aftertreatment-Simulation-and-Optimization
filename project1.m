%% 
%% *1. CATALYST DYNAMICS*

clear all 
% Initialize parameters for VTC thermal model
close all;

% Load the file
load('catalystModel.mat');

% Initial variables
s = 1;  % Index for plotting
Tinit = [190, 220, 220];  % Initial temperatures
Xint = [0, 15, 0];        % Initial states

timeStep = 0.1;              % Time step
totalTime = 180;             % Total time
numSteps = totalTime / timeStep;  % Number of steps

% Constant vectors over time
emf = 0.1 * ones(1, numSteps);        % Mass flow
inNOx = 0.2 * ones(1, numSteps);      % NOx flow
inNH3_demand = 0.2 * ones(1, numSteps);  % NH3 demand

% Loop for each set of parameters
for i = 1:length(Tinit)

    % Reset the model
    CatalystModel.reset(Xint(i), Tinit(i), totalTime);

    % Create temperature vector for each time step
    T = Tinit(i) * ones(1, numSteps);

    % Time-step loop
    for j = 1:numSteps
        CatalystModel.makeStep(emf(j), T(j), inNOx(j), inNH3_demand(j));
    end

    % Plotting
    CatalystModel.plotBufferEmissions(s);
    CatalystModel.plotBufferState(s + 1);

    s = s + 2;

    % Get logs
    CatalystLog = CatalystModel.getLogs;
end

%% *1.2 MINIMIZE EMISSIONS*

clear all 
close all;

% Load the file
load('catalystModel.mat');

timeStep = 0.1;              
totalTime = 180;              
numSteps = totalTime / timeStep;  

% Constant vectors
emf = 0.1 * ones(1, numSteps);       
inNOx = 0.2 * ones(1, numSteps);     

Tincr_demand = 0 * ones(1, numSteps);
for i = 1:numSteps
    % If time < 13 seconds, set NH3 demand to 1
    if i * timeStep < 13
        inNH3_demand(i) = 1;        
    else
        inNH3_demand(i) = 0.2;      
    end
end

T = 220 * ones(1, numSteps);           

% Initialize model
CatalystModel.reset(0, 220, totalTime);

% Simulation
for j = 1:numSteps
    CatalystModel.makeStep(emf(j), T(j), inNOx(j), inNH3_demand(j));
end

% Plot results
CatalystModel.plotEmissionsReduction(1);
CatalystModel.plotBufferState(2);

%% *2.1 TERMO_1_a*

% Case 1
clear all 
close all;

% Load VTC thermal model
load('ThermalModel.mat');

% Simulation parameters
timeStep = 0.1;                  
totalTime = 15 * 60;             
numSteps = totalTime / timeStep; 

% Constant values
emf = 0.1 * ones(1, numSteps);      
Texh = 210 * ones(1, numSteps);     
Tincr_demand = 0 * ones(1, numSteps); 
Tinit = 180;                       

% Initialize model
ThermalModel.reset(Tinit, totalTime);

% Simulation loop
for i = 1:numSteps
    ThermalModel.makeStep(emf(i), Texh(i), Tincr_demand(i));  
end

% Get logs
ThermalLog = ThermalModel.getLogs;                  

% Plot temperature readings
figure
plot(ThermalLog.t, ThermalLog.Ts1, 'r:', 'DisplayName', 'Ts1 sensor')  
hold on
plot(ThermalLog.t, ThermalLog.Ts2, 'g--', 'DisplayName', 'Ts2 sensor') 
plot(ThermalLog.t, ThermalLog.Ts3, 'b-.', 'DisplayName', 'Ts3 sensor') 
hold off

ylabel('Temperature [°C]')
xlabel('Time [s]')
title('Thermal sensor readings');
legend show
grid on

% Plot catalyst temperature
figure
plot(ThermalLog.t, ThermalLog.Tcatalyst, 'm-', 'DisplayName', 'Catalyst temperature')
ylabel('Temperature [°C]')
xlabel('Time [s]')
title('Catalyst Temperature');
legend show
grid on

%% *2.1 TERMO_1_b*

close all;
load('ThermalModel.mat');

timeStep = 0.1;                  
totalTime = 15 * 60;             
numSteps = totalTime / timeStep; 

emf = 0.1 * ones(1, numSteps);      
Texh = 210 * ones(1, numSteps);     
Tincr_demand = zeros(1, numSteps);  
Tinit = 180;                        

% Adjust temperature demand
for i = 1:numSteps
    if i * timeStep < 180
        Tincr_demand(i) = 0;        
    else
        Tincr_demand(i) = 120;      
    end
end

ThermalModel.reset(Tinit, totalTime);

for i = 1:numSteps
    ThermalModel.makeStep(emf(i), Texh(i), Tincr_demand(i));  
end

ThermalLog = ThermalModel.getLogs;                  

% Plot sensor data
figure
plot(ThermalLog.t, ThermalLog.Ts1, 'r:', 'DisplayName', 'Ts1 sensor')
hold on
plot(ThermalLog.t, ThermalLog.Ts2, 'g--', 'DisplayName', 'Ts2 sensor')
plot(ThermalLog.t, ThermalLog.Ts3, 'b-.', 'DisplayName', 'Ts3 sensor')
hold off

title('Thermal sensor readings');
ylabel('Temperature [°C]')
xlabel('Time [s]')
legend show
grid on

figure
plot(ThermalLog.t, ThermalLog.Tcatalyst, 'm-', 'DisplayName', 'Catalyst temperature')
ylabel('Temperature [°C]')
xlabel('Time [s]')
title('Catalyst Temperature');
legend show
grid on

%% *2.2 TERMO_2*

close all;
load('ThermalModel.mat');

timeStep = 0.1;                  
totalTime = 15 * 60;             
numSteps = totalTime / timeStep; 

emf = 0.8 * ones(1, numSteps);      
Texh = 210 * ones(1, numSteps);     
Tincr_demand = zeros(1, numSteps);  
Tinit = 180;                        

for i = 1:numSteps
    if i * timeStep < 180
        Tincr_demand(i) = 0;        
    else
        Tincr_demand(i) = 120;      
    end
end

ThermalModel.reset(Tinit, totalTime);

for i = 1:numSteps
    ThermalModel.makeStep(emf(i), Texh(i), Tincr_demand(i));  
end

ThermalLog = ThermalModel.getLogs;                  

figure
plot(ThermalLog.t, ThermalLog.Ts1, 'r:', 'DisplayName', 'Ts1 sensor')  
hold on
plot(ThermalLog.t, ThermalLog.Ts2, 'g--', 'DisplayName', 'Ts2 sensor') 
plot(ThermalLog.t, ThermalLog.Ts3, 'b-.', 'DisplayName', 'Ts3 sensor') 
hold off

title('Thermal sensor readings');
ylabel('Temperature [°C]')
xlabel('Time [s]')
legend show
grid on

figure
plot(ThermalLog.t, ThermalLog.Tcatalyst, 'm-', 'DisplayName', 'Catalyst temperature')
ylabel('Temperature [°C]')
xlabel('Time [s]')
title('Catalyst Temperature');
legend show
grid on

%% *3.1 OVERALL EMISSION DYNAMICS*

close all;

load('CatalystModel.mat');  
load('ThermalModel.mat');     
load('InputData_Cycle1.mat');

time = InputData.time;       
emf = InputData.emf;         
Texh = InputData.Texh;       
inNOx = InputData.eoNOx;     

numSteps = length(time);  

Tincr_demand = zeros(1, numSteps);  
inNH3_demand = zeros(1, numSteps);  

CatalystModel.reset(0,25,900);
ThermalModel.reset(25,900);

T_limit = 200;
Tincr_demand = zeros(1, numSteps);  

for i = 1:numSteps-1
    if Texh(i) >= T_limit
        Tincr_demand(i) = 0;  
    else
        Tincr_demand(i) = T_limit - Texh(i);  
    end

    inNH3_demand(i) = inNOx(i);

    ThermalModel.makeStep(emf(i), Texh(i), Tincr_demand(i));  
    ThermalValues = ThermalModel.getCurrentValues;
    CatalystModel.makeStep(emf(i), ThermalValues.Tcatalyst, inNOx(i), inNH3_demand(i));
end

iFig = 14;
Tlim = 200;
iFig = ThermalModel.plotInputs(iFig,Tlim);
iFig = ThermalModel.plotSensorData(iFig,Tlim);
iFig = ThermalModel.plotCatalystTemperature(iFig,Tlim);
iFig = ThermalModel.plotPower(iFig);
iFig = CatalystModel.plotBufferInputs(iFig);
iFig = CatalystModel.plotBufferState(iFig);
iFig = CatalystModel.plotBufferEmissions(iFig);
iFig = CatalystModel.plotEmissionsReduction(iFig);

figure
plot(time, Tincr_demand, 'b-', 'DisplayName', '∆Tdem');
ylabel('∆Tdem [°C]');
xlabel('Time [s]');
legend show;
title('Control Signal: ∆Tdem');
hold off

figure
plot(time, inNH3_demand, 'g-', 'DisplayName', 'inNH3\_dem');
ylabel('inNH3\_dem [kg/s]');
xlabel('Time [s]');
legend show;
title('Control Signal: inNH3\_dem');
hold off

%% *4.1.a FINDING ALPHA*

alpha = 0.8;  
emf = 0.1;  
T_exh = 210;  
T_sim = 900;  
T0 = 180;  

s = tf('s');
stone_tf  = (1 / ((1 / (alpha * emf)) * s + 1))^15;

time = 0:0.1:900;  
input_signal = ones(size(time)) * (T_exh - T0);

[T_response, T_time] = lsim(stone_tf, input_signal, time);  
T_response = T_response + T0;

plot(T_time, T_response, 'k-', 'LineWidth', 1.5, 'DisplayName', '15-Stone Model Response');
xlabel('Time [s]');
ylabel('Temperature [°C]');
title('15-Stone Model Response');
legend('Stone model');
grid on;
hold off

%% *4.1.b FINDING ALPHA (VTC comparison)*

clear all 
close all;

load('ThermalModel.mat');

timeStep = 0.1;                  
totalTime = 15 * 60;             
numSteps = totalTime / timeStep; 

emf = 0.1 * ones(1, numSteps);      
Texh = 210 * ones(1, numSteps);     
Tincr_demand = 0 * ones(1, numSteps); 
Tinit = 180;                       

ThermalModel.reset(Tinit, totalTime);

for i = 1:numSteps
    ThermalModel.makeStep(emf(i), Texh(i), Tincr_demand(i));  
end

ThermalLog = ThermalModel.getLogs;                  
plot(ThermalLog.t, ThermalLog.Ts3, 'b-.', 'DisplayName', 'Ts3 sensor')
hold on

alpha = 0.5;  
emf = 0.1;  
T_exh = 210;  
T_sim = 900;  
T0 = 180;  

s = tf('s');
stone_tf  = (1 / ((1 / (alpha * emf)) * s + 1))^15;

time = 0:0.1:900;  
input_signal = ones(size(time)) * (T_exh - T0);

[T_response, T_time] = lsim(stone_tf, input_signal, time);  
T_response = T_response + T0;

plot(T_time, T_response, 'k-', 'LineWidth', 1.5, 'DisplayName', '15-Stone Model Response');
xlabel('Time [s]');
ylabel('Temperature [°C]');
title('15-Stone Model Response');
legend('VTC','Stone model');
grid on;

%% 4.3 FINDING A, B MATRICES

num_stones = 15;         
emf = 0.1;
alpha = 0.04;

A = sym(zeros(num_stones, num_stones)); 
B = sym(zeros(num_stones, 1));          

for i = 1:num_stones
    if i == 1
        A(i, i) = -alpha * emf;          
    else
        A(i, i) = -alpha * emf;          
        A(i, i-1) = alpha * emf;         
    end
end

B(1) = alpha * emf;  

disp('Matrix A (symbolic):');
disp(A);

disp('Matrix B (symbolic):');
disp(B);

%% 4.4 FINDING Ad, Bd (DISCRETE MATRICES)

alpha = 0.5;              
emf = 0.1;                
num_stones = 15;          
Ts = 5;                   

A = sym(zeros(num_stones, num_stones));
B = sym(zeros(num_stones, 1));

for i = 1:num_stones
    if i == 1
        A(i, i) = -alpha * emf;          
    else
        A(i, i) = -alpha * emf;          
        A(i, i-1) = alpha * emf;         
    end
end

B(1) = alpha * emf;  

A_d = expm(A * Ts);  
B_d = simplify((expm(A * Ts) - eye(num_stones)) / A * B);  

A_d_numeric = vpa(A_d, 10);  
B_d_numeric = vpa(B_d, 10);  

disp('Discrete-time matrix A_d (numeric):');
disp(A_d_numeric);

disp('Discrete-time matrix B_d (numeric):');
disp(B_d_numeric);

%% Compare continuous vs discrete responses

alpha = 0.5;            
emf = 0.1;              
num_stones = 15;        
Ts = 5;                 
total_time = 1000;      
time = 0:Ts:total_time; 

T_exh = 210;            
T0 = 180;               

input_signal = ones(size(time)) * (T_exh - T0); 

A1 = -alpha * emf * eye(num_stones);
for i = 2:num_stones
    A1(i, i-1) = alpha * emf;  
end
B1 = [alpha * emf; zeros(num_stones-1, 1)];

X0 = T0 * ones(num_stones, 1);  

sys1 = ss(A1, B1, eye(num_stones), zeros(num_stones, 1));  
sys2 = c2d(sys1, Ts, 'zoh');  

[y1, t1] = step(sys1, time);  
[y2, t2] = step(sys2, time);  

y1 = y1 * (T_exh - T0) + T0; 
y2 = y2 * (T_exh - T0) + T0; 

figure;
hold on;
stairs(t2, y2(:, 15), 'ro', 'LineWidth', 2);  
title('Step Responses of Continuous vs. Discrete Systems for T15');
plot(t1, y1(:, 15), 'b', 'LineWidth', 2);  

xlabel('Time [s]');
ylabel('Temperature of T15 [°C]');
legend('Discrete','Continuous');
grid on;
hold off;
