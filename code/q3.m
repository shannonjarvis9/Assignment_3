% Change bottleneck width and calculate the current density after 50
% itterations 

clear all 
clear

% Initialize fundamental constants 
global C
C.temp = 300;                       % Initial temperature 
C.kb = 1.3806504e-23;               % Boltzmann constant
C.m_0 = 9.10938215e-31;             % Electron mass
C.m_e = 0.26*C.m_0;                 % Effective mass 
C.q = 1.60217662e-19;               % Charge of electron


Wb = 5:5:50;             % Initalize bottleneck width
Lb = 40; % 20:5:50;       % Initialize bottleneck height 

current = zeros(size(Wb,2)*size(Lb,2),1);         % Vector of calculated currents 
calc_fc_current = zeros(size(Wb,2)*size(Lb,2),1);         % Vector of calculated currents 
bottleneck_area = zeros(size(Wb,2)*size(Lb,2),1); % Area of bottlenect
                               

%--------------------------------------------------------------------------
%Calculate current, changing the bottleneck dimensions
%--------------------------------------------------------------------------

idx = 1;                        %Index for loop itteration
for i= 1:size(Lb,2)
    for j=1:size(Wb,2)
        bottleneck_area(idx,1) = Lb(i)*Wb(j);
        [current(idx,1), calc_fc_current(idx,1)] = calc_current_density(Lb(i), Wb(j));
        idx = idx +1
    end
end



%-------------------------------------------------------------------------- 
% Plot the current and area:
%-------------------------------------------------------------------------- 
figure('Name','Question 3b: Current density versus Bottleneck Width');  %Optional
plot(Wb, current, '*');
grid;
title('Current density versus Bottleneck Width', 'FontSize',20);
xlabel('Bottleneck Width (nm)','FontSize',13);
ylabel('Current Density (A/m)','FontSize',13);


%-------------------------------------------------------------------------- 
figure('Name','Question 3b: Current versus Bottleneck Width');  %Optional
plot(Wb, calc_fc_current, '*');
grid;
title('Current versus Bottleneck Width', 'FontSize',20);
xlabel('Bottleneck Width (nm)','FontSize',13);
ylabel('Current (A)','FontSize',13);



