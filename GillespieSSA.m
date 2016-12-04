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
Time_End = 100;                                                           % unit: ?   (end point of time line)  

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
V = 100;

% constant boundary
total_E2F1 = 1*V;
total_RB = 2*V;

% Initial Condition
init_p53helper = 0.39*V;
init_p53killer = 0.39*V;
init_MDM2 = 2.53*V;
init_ARF = 2.79*V;
init_ARF_MDM2 = 1.99*V;
init_CycE = 3.83*V; 
init_E2F1 = 0.9*V;
init_RBp = 1.88*V;
init_p21 = 1.32*V;
init_p21_CycE = 1.25*V;

init_RB_E2F1 = total_E2F1 - init_E2F1;
init_RB = total_RB - init_RBp - init_RB_E2F1;

% coefficients
k_out = 0.0003;
k_sp53 = 0.5;
kp_dp53 = 0.1;
k_pp53 = 1;
k_dpp53 = 0.5; % increase in virus B-pathway
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
k_dsRE = 1; % increase in virus A-pathway
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

% volume modification on coefficients?
k_out = k_out*V;
k_sp53 = k_sp53*V;
kp_dp53 = kp_dp53*V;
k_pp53 = k_pp53*V;
k_dpp53 = k_dpp53*V;
j_pp53 = j_pp53*V;
j_dpp53 = j_dpp53*V;
kp_sMDM2 = kp_sMDM2*V;
kpp_sMDM2 = kpp_sMDM2*V;
kp_dMDM2 = kp_dMDM2*V;
k_asAM = k_asAM*V;
k_dsAM = k_dsAM*V;
kp_sARF = kp_sARF*V;
kpp_sARF = kpp_sARF*V;
k_dARF = k_dARF*V;
k_asRE = k_asRE*V;
k_dsRE = k_dsRE*V;
k_pRB = k_pRB*V;
k_dpRB = k_dpRB*V;
j_pRB = j_pRB*V;
j_dpRB = j_dpRB*V;

kp_sE = kp_sE*V;
kpp_sE = kpp_sE*V;
k_dE = k_dE*V;
k_asp21E = k_asp21E*V;
k_dsp21E = k_dsp21E*V;
k_dp21 = k_dp21*V;
kp_sp21 = kp_sp21*V;
kpp_sp21 = kpp_sp21*V;
kppp_sp21 = kppp_sp21*V;


%%  Preallocation 
[tg,ng] = deal(cell(1,M));                                           

dist_protein_x_OFF1 = zeros(1,M);                                          
dist_protein_z_OFF1 = zeros(1,M); 

dist_p53helper = zeros(1,M);
dist_p53killer = zeros(1,M);

n_p53helper = init_p53helper;
n_p53killer = init_p53killer;
n_MDM2 = init_MDM2;
n_ARF = init_ARF;
n_ARF_MDM2 = init_ARF_MDM2;
n_CycE = init_CycE;
n_E2F1 = init_E2F1;
n_RBp = init_RBp;
n_p21 = init_p21;
n_p21_CycE = init_p21_CycE;

n_RB_E2F1 = init_RB_E2F1;
n_RB = init_RB;  

%% State description 
% [state vector]
%-------------------------------------------
stateVector = [n_p53helper;  ... 1
               n_p53killer;  ... 2
               n_MDM2;       ... 3
               n_ARF;        ... 4
               n_ARF_MDM2;   ... 5
               n_CycE;       ... 6
               n_E2F1;       ... 7
               n_RBp;        ... 8               
               n_p21;        ... 9
               n_p21_CycE;   ... 10
               n_RB_E2F1;    ... 11 (modified by boundary)
               n_RB];        ... 12 (modified by boundary)                  
%-------------------------------------------
N = length(stateVector);

initialState = stateVector; % store for future initialization

% Reaction channels (events) = columns of the state-change matrix    
% (1: p53helper) 
d_p53helper_born = zeros(N,1); 
d_p53helper_born(1) = 1;

d_p53helper_die = zeros(N,1);
d_p53helper_die(1) = -1;
    
% (2: p53killer)
d_p53killer_born = zeros(N,1);
d_p53killer_born(2) = 1;

d_p53killer_die = zeros(N,1);
d_p53killer_die(2) = -1;

d_p53killer_PP = zeros(N,1);
d_p53killer_PP(1) = 1;
d_p53killer_PP(2) = -1;

% (3: MDM2)
d_MDM2_born = zeros(N,1);
d_MDM2_born(3) = 1;

d_MDM2_die = zeros(N,1);
d_MDM2_die(3) = -1;

% (4: ARF)
d_ARF_born = zeros(N,1);
d_ARF_born(5) = 1;

d_ARF_die = zeros(N,1);
d_ARF_die(5) = -1;

% (5: ARF/MDM2)
d_ARF_MDM2_asso = zeros(N,1);
d_ARF_MDM2_asso(3) = -1;
d_ARF_MDM2_asso(4) = -1;
d_ARF_MDM2_asso(5) = 1;

d_ARF_MDM2_disso = zeros(N,1);
d_ARF_MDM2_disso(3) = 1;
d_ARF_MDM2_disso(4) = 1;
d_ARF_MDM2_disso(5) = -1;

d_ARF_MDM2die = zeros(N,1);
d_ARF_MDM2die(3) = -1;
d_ARF_MDM2die(4) = 1;
d_ARF_MDM2die(5) = -1;

d_ARFdie_MDM2_ARFdie = zeros(N,1);
d_ARFdie_MDM2_ARFdie(3) = 1;
d_ARFdie_MDM2_ARFdie(4) = -1;
d_ARFdie_MDM2_ARFdie(5) = -1;

% (6: CycE)
d_CycE_born = zeros(N,1);
d_CycE_born(6) = 1;

d_CycE_die = zeros(N,1);
d_CycE_die(6) = -1;

% (7: E2F1)
d_RB_E2F1_asso = zeros(N,1);
d_RB_E2F1_asso(7) = -1;

d_RB_E2F1_disso = zeros(N,1);
d_RB_E2F1_disso(7) = 1;


% (8: RBp)
% +
d_RBp_born = zeros(N, 1);
d_RBp_born(8) = 1;

% - 
d_RBp_die = zeros(N,1);
d_RBp_born(8) = -1;    

% (9: p21)

% (10: p21/CycE)

% [11: RB/E2F1]
% [12: RB]
    
% State-change matrix
D = [
     d_p53helper_born, d_p53helper_die, ...     
     d_p53killer_born, d_p53killer_die, d_p53killer_PP, ...
     d_MDM2_born, d_MDM2_p53T, d_MDM2_disso, d_MDM2_die, d_ARF_MDM2_asso, ...
     d_ARF_born, d_ARF_E2F1, d_ARF_disso, d_ARF_die, d_ARF_asso, ...
     d_ARF_MDM2_asso, d_ARF_MDM2_disso, d_ARF_MDM2_dieMDM2, d_ARF_MDM2_dieARF, ...
     d_RB_E2F1_disso, d_RB_E2F1_asso, ...
     d_RBp_born, d_RBp_die, ...
     
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
        a_matrix = [
                    k_sp53, (kp_dp53 + stateVector(3))*stateVector(1), 
                    
                   ];

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