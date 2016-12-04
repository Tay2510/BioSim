%% Script for simulating p53 network %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @author: Chao-Ming (Jeremy) Yen
% @model: p53...MDM2...ARF...E2F1...Rb
% @algorithm: Gillespie SSA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;
clc
tic
fprintf('----------------------[Start Simulation]----------------------\n')

%% General Setting
% [ General ]
M = 1;                                                                      % number of trials
Time_End = 1000;                                                           % unit: ?   (end point of time line)  

% [ Parameters ]
% degradation
gamma_x = 0.0004;                                                           % unit: 1/s
gamma_z = 0.0004;                                                           % unit: 1/s    

% mean burst size
b_x = 20;                                                                   % unit: none                                                                  
b_z = 20;                                                                   % unit: none    

a_x = 50;
a_z = 50;

% production
kx_ON = gamma_x*a_x;                                                              
kz_ON = gamma_x*a_x;                                                              

% Basal
basal_z = 0.1;

% [ State of input X ]
px = 1; % for modifying the state, ON = 1, OFF = 0.1 etc;

% [ Hill functions ]
% Hill coefficient
h = 2;                                                                    
Kxz = 80;                                                                  
                                                            
Hxz_act = @(x) x^h/(Kxz^h + x^h);

%% Model parameters
% cell volume
volume = 100;

% Initial Condition
%DYRK2_initial = 2*volume;
p53helper_initial = 0.39*volume;
p53killer_initial = 0.39*volume;
MDM2_initial = 2.53*volume;
ARF_initial = 2.79*volume;
ARF_MDM2_initial = 1.99*volume;
E2F1_initial = 0.9*volume;
RBp_initial = 1.88*volume;
RB_E2F1_initial = 0*volume;
RB_initial = 0*volume;

% constant boundary
total_E2F1 = 1*volume;
total_RB = 2*volume;

% coefficients
k_out = 0.0003;
k_sp53 = 0.5;
kp_dp53 = 0.1;
k_pp53 = 1;
k_dpp53 = 0.5;
j_pp53 = 0.1;
j_dpp53 = 0.1;
kp_sMDM2 = 0.02;
kpp_sMDM2 = 0.02;
kp_dMDM2 = 0.1;
k_asAM = 10;
k_dsAM = 2;
kp_sARF = 0.01;
kpp_sARF = 0.3;
k_dARF = 0.1;
k_asRE = 5;
k_dsRE = 1;
k_pRB = 1;
k_dpRB = 0.5;
j_pRB = 0.1;
j_dpRB = 0.1;

kp_sE = 0.01;
kpp_sE = 0.5;
k_dE = 0.12;
k_asp21E = 1;
k_dsp21E = 10;
k_dp21 = 0.2;
kp_sp21 = 0.03;
kpp_sp21 = 0.3;
kppp_sp21 = 0.01;


% volume modification of on coefficients?
k_out = k_out*volume;
k_sp53 = k_sp53*volume;
kp_dp53 = kp_dp53*volume;
k_pp53 = k_pp53*volume;
k_dpp53 = k_dpp53*volume;
j_pp53 = j_pp53*volume;
j_dpp53 = j_dpp53*volume;
kp_sMDM2 = kp_sMDM2*volume;
kpp_sMDM2 = kpp_sMDM2*volume;
kp_dMDM2 = kp_dMDM2*volume;
k_asAM = k_asAM*volume;
k_dsAM = k_dsAM*volume;
kp_sARF = kp_sARF*volume;
kpp_sARF = kpp_sARF*volume;
k_dARF = k_dARF*volume;
k_asRE = k_asRE*volume;
k_dsRE = k_dsRE*volume;
k_pRB = k_pRB*volume;
k_dpRB = k_dpRB*volume;
j_pRB = j_pRB*volume;
j_dpRB = j_dpRB*volume;


%%  Preallocation 
[tg,ng] = deal(cell(1,M));                                           

dist_protein_x_OFF1 = zeros(1,M);                                          
dist_protein_z_OFF1 = zeros(1,M); 

dist_p53helper = zeros(1,M);
dist_p53killer = zeros(1,M);
dist_MDM2 = zeros(1,M);
dist_ARF = zeros(1,M);
dist_ARF_MDM2 = zeros(1,M);
dist_E2F1 = zeros(1,M);
dist_RBp = zeros(1,M);
dist_RB_E2F1 = zeros(1,M);
dist_RB = zeros(1,M);

n_p53helper = p53helper_initial;
n_p53killer = p53killer_initial;
n_MDM2 = MDM2_initial;
n_ARF = ARF_initial;
n_ARF_MDM2 = ARF_MDM2_initial;
n_E2F1 = E2F1_initial;
n_RBp = RBp_initial;
n_RB_E2F1 = RB_E2F1_initial;
n_RB = RB_initial;  

%% State description 
% [state vector]
%-------------------------------------------
stateVector = [
               n_p53helper;  ... 1
               n_p53killer;  ... 2
               n_MDM2;       ... 3
               n_ARF;        ... 4
               n_ARF_MDM2;   ... 5
               n_E2F1;       ... 6
               n_RBp;        ... 7
               n_CycE;       ... 8
               n_p21;        ... 9
               n_RB_E2F1;    ... 10 (modified by boundary)
               n_RB];        ... 11 (modified by boundary)                  
%-------------------------------------------
N = length(stateVector);
initialState = stateVector; % store for future initialization

% events (columns of state-change matrix)
    
% (1: p53helper) 
% + 
dp_p53helper_born = zeros(N,1); 
dp_p53helper_born(2) = 1;
dp_p53helper_pp = zeros(N,1);
dp_p53helper_pp(2) = 1;
%dp_p53helper_pp(3) = -1;  %???
    
% - 
dn_p53helper_die = zeros(N,1);
dn_p53helper_die(2) = -1;
dn_p53helper_DYRK2 = zeros(N,1); 
dn_p53helper_DYRK2(2) = -1;
    
% (2: p53killer)
% + 
dp_p53killer_DYRK2 = zeros(N,1);
dp_p53killer_DYRK2(3) = 1;
    
% - 
dn_p53killer_die = zeros(N,1);
dn_p53killer_die(3) = -1;
dn_p53killer_pp = zeros(N,1);
dn_p53killer_pp(3) = -1;
%dn_p53killer_pp(2) = 1;   %???    
    
% (3: MDM2)
% + 
dp_MDM2_born = zeros(N,1);
dp_MDM2_born(4) = 1;
dp_MDM2_p53T = zeros(N,1);
dp_MDM2_p53T(4) = 1;
dp_MDM2_disso = zeros(N,1);
dp_MDM2_disso(4) = 1;

% -
dn_MDM2_die = zeros(N,1);
dn_MDM2_die(4) = -1;
dn_MDM2_ARF = zeros(N,1);
dn_MDM2_ARF(4) = -1;

% (4: ARF)
% +
dp_ARF_born = zeros(N,1);
dp_ARF_born(5) = 1;
dp_ARF_E2F1 = zeros(N,1);
dp_ARF_E2F1(5) = 1;
dp_ARF_disso = zeros(N,1);
dp_ARF_disso(5) = 1;

% - 
dn_ARF_die = zeros(N,1);
dn_ARF_die(5) = -1;
dn_ARF_asso = zeros(N,1);
dn_ARF_asso(5) = -1;

% (5: ARF/MDM2)
% +
dp_ARF_MDM2_asso = zeros(N,1);
dp_ARF_MDM2_asso(6) = 1;

% -
dn_ARF_MDM2_disso = zeros(N,1);
dn_ARF_MDM2_disso(6) = -1;
dn_ARF_MDM2_dieMDM2 = zeros(N,1);
dn_ARF_MDM2_dieMDM2(6) = -1;
dn_ARF_MDM2_dieARF = zeros(N,1);
dn_ARF_MDM2_dieARF(6) = -1;

% (6: E2F1)
% +
dp_E2F_disso = zeros(N,1);
dp_E2F_disso(7) = 1;

% -
dn_E2F_asso = zeros(N,1);
dn_E2F_asso(7) = -1;

% (7: RBp)
% +
dp_RBp_born = zeros(N, 1);
dp_RBp_born(8) = 1;

% - 
dn_RBp_die = zeros(N,1);
dp_RBp_born(8) = -1;    

% (8: CycE)
% +
dp_CycE_born = zeros(N,1);
dp_CycE_born(9) = 1;

% - 
dn_CycE_die = zeros(N,1);
dn_CycE_die(9) = -1;
dn_CycE_asso = zeros(N,1);
dn_CycE_asso(9) = -1;

% (9: p21)
% +
% - 

% (RB/E2F1)
% (RB)
    
% State-change matrix
D = [
     dp_p53helper_born, dp_p53helper_pp, dn_p53helper_die dn_p53helper_DYRK2,...
     dp_p53killer_DYRK2, dn_p53killer_die, dn_p53killer_pp, ...
     dp_MDM2_born, dp_MDM2_p53T, dp_MDM2_disso, dn_MDM2_die, dn_MDM2_ARF, ...
     dp_ARF_born, dp_ARF_E2F1, dp_ARF_disso, dn_ARF_die, dn_ARF_asso, ...
     dp_ARF_MDM2_asso, dn_ARF_MDM2_disso, dn_ARF_MDM2_dieMDM2, dn_ARF_MDM2_dieARF, ...
     dp_E2F_disso, dn_E2F_asso, ...
     dp_RBp_born, dn_RBp_die, ...
     
    ]; 


%% Operating Center
for i = 1:M
    
    % [Initialization]
    t = 0;                                                                  % starting time     
    n_x = 0;                                                                % starting particle number
    n_z = 0;
    n_xz = [n_x ;...
            n_z];          

    stateVector = initialState;    
        
    tg{i} = [];                                                             % claim a empty [] for recording time
    ng{i} = [];                                                             % claim a empty [] for recording particle number
    
    dist_protein_x_OFF1_temp = 0;
    dist_protein_z_OFF1_temp = 0;
    
    % [Gillespie SSA]
    % ====================================================================
    while t < Time_End                                                      % [Loop] Gillespie main framework        
        % sampling
        r1 = rand;                                                          % 1st dice        
        r2 = rand;                                                          % 2nd dice    

        % update propensity matrix
        %a_matrix = [px*kx_ON, n_xz(1,1)*gamma_x, (basal_z*kz_ON + Hxz_act(n_xz(1,1))*(1-basal_z)*kz_ON), n_xz(2,1)*gamma_z];
        
        a_matrix = [px*kx_ON, n_xz(1,1)*gamma_x, (basal_z*kz_ON + Hxz_act(n_xz(1,1))*(1-basal_z)*kz_ON), n_xz(2,1)*gamma_z];
        
        
        % State changing matrix
        %D = [round(exprnd(b_x)) , (-1) , 0 , 0 ;...
             %0 , 0 , round(exprnd(b_z)) , (-1)];

        % [State recording]
        tg{i} = [tg{i} t];                                                       
        %ng{i} = [ng{i} n_xz];
        ng{i} = [ng{i} stateVector];

        % [Stochastic part]
        astr = sum(a_matrix);
        if ~astr , break ,end
        tau =  (log(1/r1))/astr;
        u = find(cumsum(a_matrix)>=astr*r2,1);
        t = t + tau;                                                        % updating time

        % [State Changing]
        n_xz = n_xz + D(:,u);
        stateVector = stateVector + D(:,u);

        % [Negative value filters]
        %{
        if n_xz(1,1) <0
            n_xz(1,1) =0;
        end

        if n_xz(2,1) <0
            n_xz(2,1) = 0;
        end
        %}
        
        for j = 1 : length(stateVector)
            if stateVector(j) < 0
                stateVector(j) = 0;
            end
        end

        % [Data recording for steady state analysis]
        dist_protein_x_OFF1_temp = n_xz(1,1);
        dist_protein_z_OFF1_temp = n_xz(2,1);
    
    end                                                                     % [End] Gillespie main framework
    % ====================================================================
    
    % [Statistic data extraction]
    % OFF_1
    %dist_protein_x_OFF1(1,i) = dist_protein_x_OFF1_temp;
    %dist_protein_z_OFF1(1,i) = dist_protein_z_OFF1_temp;

    
    % [update waitbar]
    fprintf('progress: %d %% \r', 100*i/M);

end                                                                         % [End] of one trial 


%% Plot panel
% Sample trajectory
subplot(2,3,[1 2 3])
stairs(tg{1},ng{1}(1,:));
hold all
stairs(tg{1},ng{1}(2,:),'--r');
hold off
title('one of sample trajectories')
xlabel('time(sec)')
ylabel('particle numbers')
legend('X particle number','Z particle number')

%{
% Steady state distribution
bin_num = 20;
% OFF_1
subplot(2,3,4)
[elements_x_OFF1,centers_x_OFF1] = hist(dist_protein_x_OFF1,bin_num);
plot(centers_x_OFF1,elements_x_OFF1,'b','LineWidth',2);
hold on
[elements_z_OFF1,centers_z_OFF1] = hist(dist_protein_z_OFF1,bin_num);
plot(centers_z_OFF1,elements_z_OFF1,'--r', 'LineWidth',2);
hold off
title('Steady state distribution (OFF state 1)')
xlabel('particle number')
ylabel('frequency')
legend('X','Z')

%}
%% Statistics Analysis Output
% [Steady state distribution]
% OFF_1
%mean_SSdist_X_OFF1 = mean(dist_protein_x_OFF1)
%var_SSdist_X_OFF1 = var(dist_protein_x_OFF1)
%CV2_SSdist_X_OFF1 = var_SSdist_X_OFF1 / mean_SSdist_X_OFF1^2
%Fano_SSdist_X_OFF1 = var_SSdist_X_OFF1 / mean_SSdist_X_OFF1

%mean_SSdist_Z_OFF1 = mean(dist_protein_z_OFF1)
%var_SSdist_Z_OFF1 = var(dist_protein_z_OFF1)
%CV2_SSdist_Z_OFF1 = var_SSdist_Z_OFF1 / mean_SSdist_Z_OFF1^2
%Fano_SSdist_Z_OFF1 = var_SSdist_Z_OFF1 / mean_SSdist_Z_OFF1



%% Miscellaneous
fprintf('-----------------------[End Simulation]-----------------------\n')

toc