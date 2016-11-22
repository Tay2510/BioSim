clc
clear all;

T_end = 10;
delta_t = 0.01;
N = T_end/delta_t;
timeline = 0:delta_t:(T_end - delta_t);
sigma = 0.001

MDM2_initial = 2.53;
ARF_initial = 2.79;
ARF_MDM2_initial = 1.99;
E2F1_initial = 0.9;
RBp_initial = 1.88;
RB_E2F1_initial = 0;
RB_initial = 0;

p53helper_initial = 0.39;
p53killer_initial = 0.39;
DYRK2_initial = 2;

p53helper = zeros(1, N);
p53killer = zeros(1, N);
MDM2 = zeros(1, N);
ARF = zeros(1, N);
ARF_MDM2 = zeros(1, N);
E2F1 = zeros(1, N);
RBp = zeros(1, N);
RB_E2F1 = zeros(1, N);
RB = zeros(1, N);
p53killer(1) = p53killer_initial;
p53helper(1) = p53helper_initial;
MDM2(1) = MDM2_initial;
ARF(1) = ARF_initial;
ARF_MDM2(1) = ARF_initial;
E2F1(1) = E2F1_initial;
RBp(1) = RBp_initial;
RB_E2F1(1) = RB_E2F1_initial;
RB(1) = RB_initial;
DYRK2(1) = DYRK2_initial;

for i = 2 : N        
    DYRK2(i) = DYRK2(i - 1) + delta_t * (-0.0003 * DYRK2(i - 1));
    p53helper(i) = p53helper(i - 1) + delta_t * (-DYRK2(i - 1) * p53helper(i - 1) / (0.1 + p53helper(i - 1)) + 0.5 - (0.1 + MDM2(i - 1)) * p53helper(i - 1) + 0.5 * p53killer(i - 1) / (0.1 + p53killer(i - 1)));
    p53killer(i) = p53killer(i - 1) + delta_t * (DYRK2(i - 1) * p53helper(i - 1) / (0.1 + p53helper(i - 1)) - 0.5 * p53killer(i - 1) / (p53killer(i - 1)) - (0.1 + MDM2(i - 1) * p53killer(i - 1)));
    MDM2(i) = MDM2(i - 1) + delta_t * (0.02 + 0.3 * (p53killer(i - 1) + p53helper(i - 1)) - 0.1 * MDM2(i - 1) - 10 * ARF(i - 1) * MDM2(i - 1) + 2.1 * ARF_MDM2(i - 1));
    ARF(i) = ARF(i - 1) + delta_t * (0.01 + 0.3 * E2F1(i - 1) - 0.1 * ARF(i - 1) - 10 * ARF(i - 1) * MDM2(i - 1) + 2.1 * ARF_MDM2(i - 1));
    ARF_MDM2(i) = ARF_MDM2(i - 1) + delta_t * (10 * ARF(i - 1) * MDM2(i - 1) - 2 * ARF_MDM2(i - 1) - 0.1 * ARF_MDM2(i - 1) - 0.1 * ARF_MDM2(i - 1));
    E2F1(i) = E2F1(i - 1) + delta_t * (-5 * RB(i - 1) * E2F1(i - 1) + RB_E2F1(i - 1));
    RBp(i) = sigma + RBp(i - 1) + delta_t * (-0.5 * RBp(i - 1) / (0.1 + RBp(i - 1)));
    RB_E2F1(i) = RB_E2F1(i - 1) + delta_t * (1 - E2F1(i - 1));
    RB(i) = RB(i - 1) + delta_t * (2 - RBp(i - 1) - RB_E2F1(i - 1));
end


%% Plot
plot(timeline, p53killer, timeline, p53helper, timeline, MDM2, timeline, ARF, timeline, RB, timeline, RBp);
legend('p53_k_i_l_l_e_r','p53_h_e_l_p_e_r','MDM2','ARF', 'RB', 'RBp');
xlabel('time');
ylabel('concentration');
