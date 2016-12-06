function [ratio2] = helper(kdpp53_change, kdsRE_change)
T_end = 20;
delta_t = 0.001;
N = T_end/delta_t;

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

CycE_initial = 3.83;
p21_initial = 1.32;
p21_CycE_initial = 1.25;

DYRK2 = zeros(1, N);
p53helper = zeros(1, N);
p53killer = zeros(1, N);
MDM2 = zeros(1, N);
ARF = zeros(1, N);
ARF_MDM2 = zeros(1, N);
E2F1 = zeros(1, N);
RBp = zeros(1, N);
RB_E2F1 = zeros(1, N);
RB = zeros(1, N);

CycE = zeros(1, N);
p21_CycE = zeros(1, N);
p21 = zeros(1, N);

CycE(1) = CycE_initial;
p21_CycE(1) = p21_CycE_initial;
p21(1) = p21_initial;

p53killer(1) = p53killer_initial;
p53helper(1) = p53helper_initial;
MDM2(1) = MDM2_initial;
ARF(1) = ARF_initial;
ARF_MDM2(1) = ARF_MDM2_initial;
E2F1(1) = E2F1_initial;
RBp(1) = RBp_initial;
RB_E2F1(1) = RB_E2F1_initial;
RB(1) = RB_initial;
DYRK2(1) = DYRK2_initial;

kdpp53_initial = 0.5;
kdsRE_initial = 1;

kdpp53 = kdpp53_initial;
kdsRE = kdsRE_initial;

for i = 2 : N        
    DYRK2(i) = DYRK2(i - 1) + delta_t * (-0.0003 * DYRK2(i - 1));
    if DYRK2(i) < 0
        DYRK2(i) = 0;
    end
    p53helper(i) = p53helper(i - 1) + delta_t * (-DYRK2(i - 1) * ...
        p53helper(i - 1) / (0.1 + p53helper(i - 1)) + 0.5 - ...
        (0.1 + MDM2(i - 1)) * p53helper(i - 1) + kdpp53 * p53killer(i - 1) ...
        /(0.1 + p53killer(i - 1)));
    if p53helper(i) < 0
       p53helper(i) = 0;
    end
    p53killer(i) = p53killer(i - 1) + delta_t * (DYRK2(i - 1) * ...
        p53helper(i - 1) / (0.1 + p53helper(i - 1)) - kdpp53 * ...
        p53killer(i - 1) / (0.1 + p53killer(i - 1)) - (0.1 + MDM2(i - 1) * ...
        p53killer(i - 1)));
    if p53killer(i) < 0
        p53killer(i) = 0;
    end
    MDM2(i) = MDM2(i - 1) + delta_t * (0.02 + 0.3 * (p53killer(i - 1) ...
        + p53helper(i - 1)) - 0.1 * MDM2(i - 1) - 10 * ARF(i - 1) * ...
        MDM2(i - 1) + 2.1 * ARF_MDM2(i - 1));
    if MDM2(i) < 0
        MDM2(i) = 0;
    end
    ARF(i) = ARF(i - 1) + delta_t * (0.01 + 0.3 * E2F1(i - 1) - 0.1 ...
        * ARF(i - 1) - 10 * ARF(i - 1) * MDM2(i - 1) + 2.1 ...
        * ARF_MDM2(i - 1));
    if ARF(i) < 0
        ARF(i) = 0;
    end
    ARF_MDM2(i) = ARF_MDM2(i - 1) + delta_t * (10 * ARF(i - 1) * ...
        MDM2(i - 1) - 2 * ARF_MDM2(i - 1) - 0.1 * ARF_MDM2(i - 1) - ...
        0.1 * ARF_MDM2(i - 1));
    if ARF_MDM2(i) < 0
        ARF_MDM2(i) = 0;
    end
    E2F1(i) = E2F1(i - 1) + delta_t * (-5 * RB(i - 1) * E2F1(i - 1) + ...
        kdsRE * RB_E2F1(i - 1));
    if E2F1(i) < 0
        E2F1(i) = 0;
    end
    RBp(i) = RBp(i - 1) + delta_t * (CycE(i - 1) * RB(i - 1) /(0.1 + RB(i - 1)) -0.5 * RBp(i - 1) / ...
        (0.1 + RBp(i - 1)));
    if RBp(i) < 0
        RBp(i) = 0;
    end
    CycE(i) = CycE(i - 1) + delta_t * (0.01 + 0.5 * E2F1(i - 1) - 0.12 * CycE(i - 1) - 10 * p21(i - 1) * CycE(i - 1) + 1.2 * p21_CycE(i - 1));
    if CycE(i) < 0
        CycE(i) = 0;
    end
    p21(i) = p21(i - 1) + delta_t * (0.03 + 0.3 * p53helper(i - 1) + 0.01 * p53killer(i - 1) - 0.2 * p21(i - 1) - 10 * p21(i - 1) * CycE(i - 1) + 1.12 * p21_CycE(i - 1));
    if p21(i) < 0
        p21(i) = 0;
    end
    p21_CycE(i) = p21_CycE(i - 1) + delta_t * (10 * p21(i - 1) * CycE(i - 1) - p21_CycE(i - 1) - 0.32 * p21_CycE(i - 1));
    if p21_CycE(i) < 0
        p21_CycE(i) = 0;
    end
    RB_E2F1(i) = 3 - E2F1(i);
    RB(i) = 2 - RBp(i) - RB_E2F1(i);
    if i == N / 2 + 1
       kdpp53 = kdpp53_change;
       kdsRE = kdsRE_change;
    end
end
ratio2 = p53helper(N)/p53killer(N);


