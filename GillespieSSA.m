%% Script for simulating p53 network %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @author: Chao-Ming (Jeremy) Yen
% @model: p53...MDM2...ARF...E2F1...Rb
% @algorithm: Gillespie SSA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
tic
fprintf('----------------------[Start Simulation]----------------------\n');

%% General Setting
M = 1;                                                                      % number of trials
Time_End = 10;                                                              % unit: ?   (end point of time line)  

%% Model parameters
% cell volume
V = 100;

% constant boundary
total_E2F1 = 1 * V;
total_RB = 2 * V;

% Initial Condition
init_p53helper = 0.39 * V;
init_p53killer = 0.39 * V;
init_MDM2 = 2.53 * V;
init_ARF = 2.79 * V;
init_ARF_MDM2 = 1.99 * V;
init_CycE = 3.83 * V; 
init_E2F1 = 0.9 * V;
init_RBp = 1.88 * V;
init_p21 = 1.32 * V;
init_p21_CycE = 1.25 * V;

init_RB_E2F1 = total_E2F1 - init_E2F1;
init_RB = total_RB - init_RBp - init_RB_E2F1;

% assumption
% -----------------------------------
const_DYRK2 = 1 * V;
const_PP = 1 * V;
k_byMDM2 = 1 / V;
% -----------------------------------

% [rate constants]
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
k_sp53 = 0.5 * V;
kp_dp53 = 0.1;
k_pp53 = 1;

% ### increase in virus B-pathway ###
% -----------------------------------
k_dpp53 = 0.5; 
k_dpp53_virus = 0.8;
% -----------------------------------

j_pp53 = 0.1 * V;
j_dpp53 = 0.1 * V;
kp_sMDM2 = 0.02 * V;
kpp_sMDM2 = 0.02;
kp_dMDM2 = 0.1;
k_asAM = 10 / V;
k_dsAM = 2;
kp_sARF = 0.01 * V;
kpp_sARF = 0.3;
k_dARF = 0.1;

kp_sE = 0.01 * V;
kpp_sE = 0.5;
k_dE = 0.12;
k_asp21E = 1 / V;
k_dsp21E = 10;
k_dp21 = 0.2;

k_asRE = 5 / V;

% ### increase in virus A-pathway ###
% -----------------------------------
k_dsRE = 1; 
k_dsRE_virus = 1.2;
% -----------------------------------

k_pRB = 1;
k_dpRB = 0.5;
j_pRB = 0.1 * V;
j_dpRB = 0.1 * V;

kp_sp21 = 0.03 * V;
kpp_sp21 = 0.3;
kppp_sp21 = 0.01;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%%  Preallocation 
[tg,ng] = deal(cell(1,M));                                           

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
sVect = [n_p53helper;  ... 1
         n_p53killer;  ... 2
         n_MDM2;       ... 3
         n_ARF;        ... 4
         n_ARF_MDM2;   ... 5
         n_CycE;       ... 6
         n_E2F1;       ... 7
         n_RBp;        ... 8               
         n_p21;        ... 9
         n_p21_CycE;   ... 10
         n_RB_E2F1;    ... 11 
         n_RB];        ... 12 
%-------------------------------------------
N = length(sVect);

initialState = sVect; % store for future initialization

% Reaction channels(events) = columns of the state-change matrix    
% (1: p53helper) 
d_p53helper_born = zeros(N,1); 
d_p53helper_born(1) = 1;

d_p53helper_die = zeros(N,1);
d_p53helper_die(1) = -1;
    
% (2: p53killer)
d_p53killer_die = zeros(N,1);
d_p53killer_die(2) = -1;

d_p53killer_DYRK2 = zeros(N,1);
d_p53killer_DYRK2(1) = -1;
d_p53killer_DYRK2(2) = 1;

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

d_ARFdie_MDM2 = zeros(N,1);
d_ARFdie_MDM2(3) = 1;
d_ARFdie_MDM2(4) = -1;
d_ARFdie_MDM2(5) = -1;

% (6: CycE)
d_CycE_born = zeros(N,1);
d_CycE_born(6) = 1;

d_CycE_die = zeros(N,1);
d_CycE_die(6) = -1;

% (7: E2F1)
d_RB_E2F1_asso = zeros(N,1);
d_RB_E2F1_asso(7) = -1;
d_RB_E2F1_asso(11) = 1;
d_RB_E2F1_asso(12) = -1;

d_RB_E2F1_disso = zeros(N,1);
d_RB_E2F1_disso(7) = 1;
d_RB_E2F1_asso(11) = -1;
d_RB_E2F1_asso(12) = 1;

% (8: RBp)
d_RB_PP = zeros(N, 1);
d_RB_PP(8) = 1;
d_RB_PP(12) = -1;

d_RBp_dPP = zeros(N,1);
d_RBp_dPP(8) = -1;
d_RBp_dPP(12) = 1;

% (9: p21)
d_p21_born = zeros(N,1);
d_p21_born(9) = 1;

d_p21_die = zeros(N,1);
d_p21_die(9) = -1;

% (10: p21/CycE)
d_p21_CycE_asso = zeros(N,1);
d_p21_CycE_asso(6) = -1;
d_p21_CycE_asso(9) = -1;
d_p21_CycE_asso(10) = 1;

d_p21_CycE_disso = zeros(N,1);
d_p21_CycE_disso(6) = 1;
d_p21_CycE_disso(9) = 1;
d_p21_CycE_disso(10) = -1;

d_p21die_CycE = zeros(N,1);
d_p21die_CycE(6) = 1;
d_p21die_CycE(9) = -1;
d_p21die_CycE(10) = -1;

d_p21_CycEdie = zeros(N,1);
d_p21_CycEdie(6) = -1;
d_p21_CycEdie(9) = 1;
d_p21_CycEdie(10) = -1;

% [11: RB/E2F1]
% N/A

% [12: RB]
% N/A

% State-change matrix
D = [d_p53helper_born, ...1
     d_p53helper_die, ...2                               
     d_p53killer_die, ...3 
     d_p53killer_DYRK2, ...4 
     d_p53killer_PP, ...5              
     d_MDM2_born, ...6 
     d_MDM2_die, ...7                                         
     d_ARF_born, ...8 
     d_ARF_die, ...9                                           
     d_ARF_MDM2_asso, ...10 
     d_ARF_MDM2_disso, ...11 
     d_ARF_MDM2die, ...12 
     d_ARFdie_MDM2, ...13 
     d_CycE_born, ...14 
     d_CycE_die, ...15                                         
     d_RB_E2F1_asso, ...16 
     d_RB_E2F1_disso, ...17                                 
     d_RB_PP, ...18 
     d_RBp_dPP, ...19                                              
     d_p21_born, ...20 
     d_p21_die, ...21                                           
     d_p21_CycE_asso, ...22 
     d_p21_CycE_disso, ...23 
     d_p21die_CycE, ...24 
     d_p21_CycEdie]; % 25 
    

%% Operating Center
for i = 1:M
    
    % [Initialization]
    t = 0;                                                                  % starting time     
    sVect = initialState;    
        
    tg{i} = [];                                                             % for recording time
    ng{i} = [];                                                             % for recording particle number
    
    % [Gillespie SSA]
    % ====================================================================
    while t < Time_End                                                      % [Loop] Gillespie main framework        
        % stochastic sampling
        r1 = rand;                                                          % 1st dice: U(0,1)       
        r2 = rand;                                                          % 2nd dice: U(0,1)

        % update propensity matrix                
        a_matrix = [k_sp53, ...1
                   (kp_dp53 + k_byMDM2 * sVect(3)) * sVect(1), ...2
                   (kp_dp53 + k_byMDM2 * sVect(3)) * sVect(2), ...3
                    k_pp53 * const_DYRK2 * sVect(1) / (j_pp53 + sVect(1)), ...4
                    k_dpp53 * const_PP * sVect(2) / (j_dpp53 + sVect(2)), ...5
                    kp_sMDM2 + kpp_sMDM2 * (sVect(1) + sVect(2)), ...6
                    kp_dMDM2 * sVect(3), ...7
                    kp_sARF + kpp_sARF*sVect(7), ...8
                    k_dARF * sVect(4), ...9
                    k_asAM * sVect(4) * sVect(3), ...10
                    k_dsAM * sVect(5), ...11
                    kp_dMDM2 * sVect(5), ...12
                    k_dARF * sVect(5), ...13
                    kp_sE + kpp_sE * sVect(7), ...14
                    k_dE * sVect(6), ...15
                    k_asRE * sVect(12) * sVect(7), ...16
                    k_dsRE * sVect(11), ...17
                    k_pRB * sVect(6) * sVect(12) / (j_pRB + sVect(12)), ...18
                    k_dpRB * const_PP * sVect(8) / (j_dpRB + sVect(8)), ...19
                    kp_sp21 + kpp_sp21 * sVect(1) + kppp_sp21 * sVect(2), ...20
                    k_dp21 * sVect(9), ...21
                    k_asp21E * sVect(9) * sVect(6), ...22
                    k_dsp21E * sVect(10), ...23
                    k_dE * sVect(10), ...24
                    k_dp21 * sVect(10)]; % 25                    
               

        % [Time & State tracing]
        tg{i} = [tg{i} t];                                                            
        ng{i} = [ng{i} sVect];

        % [Gillespie]
        astr = sum(a_matrix);
        if ~astr 
            disp('bad a_matrix');
            break
        end
        tau = (log(1/r1))/astr;
        u = find(cumsum(a_matrix)>=astr*r2, 1);
        t = t + tau;                                                        % updating time

        % [Update current state]        
        sVect = sVect + D(:,u);

        % [Negative value filters]
        sVect(sVect < 0) = 0;
        
        disp(t);
    end                                                                     % [End] Gillespie main framework
    % ====================================================================
    
    % [update waitbar]
    fprintf('progress: %d %% \r', 100*i/M);

end                                                                         % [End] of one trial 

%% Plot panel
% [Sample trajectory]
% select target to be included in the plot
% ----------------------------------------  
selectVect = [1 2 3 4 12]; 
% ----------------------------------------  
legendDict = {'p53_h_e_l_p_e_r', ...1
              'p53_k_i_l_l_e_r', ...2
              'MDM2',            ...3
              'ARF',             ...4
              'ARF/MDM2',        ...5
              'CycE',            ...6
              'E2F1',            ...7
              'RB_p',            ...8
              'p21',             ...9
              'p21/CycE',        ...10
              'RB/E2F1',         ...11
              'RB',              ...12
              };

legends = cell(length(selectVect), 1);

figure; hold on
for p = 1 : length(selectVect)
    stairs(tg{1},ng{1}(selectVect(p),:));
    legends{p} = legendDict{selectVect(p)};
end

hold off
title('sample trajectort');
xlabel('time');
ylabel('particle numbers');
legend(legends);

%% Miscellaneous
fprintf('-----------------------[End Simulation]-----------------------\n')
toc
