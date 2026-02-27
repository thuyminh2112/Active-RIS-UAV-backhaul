function [AllSolution, Status] = Get_Solution_3_1 (sys, channel,chan, Init_2,r_uav_current, phi_current, p_current, varphi_current, tau_current, slack_3)

rng('default'); 
rng(1);  

AllSolution = cell(7,1);
slack_variable_3 = cell(9,1);

K = sys.K;
N = sys.N; Na = sys.Na; a_max = sys.a_max; I_active =sys.I_active;
Pmax = sys.Pmax;
P_B = sys.P_B;
P_ris = sys.P_ris;
sigma_u = sys.sigma_u;
sigma_r = sys.sigma_r;
r_u = sys.r_u; r_ris = sys.r_ris; r_bs = sys.r_bs;
e0 = sys.e0; e1 = sys.e1; e2= sys.e2;
B = sys.B;

r = sys.R_QoS;

g_B = channel{1}; G = channel{2}; hr = channel{3}; hd = channel{4};
gB_1 = chan{1}; G_0 = chan{2}; hr_1 = chan{3}; hd_1 = chan{4};

%% Slack variable
r0_low_current = slack_3{1}; r1_low_current = slack_3{2}; r0_up_current = slack_3{3}; 
r1_up_current = slack_3{4}; r2_current = slack_3{5}; varrho_low_current = slack_3{6};
varrho_up_current = slack_3{7}; z_current = slack_3{8}; q_current = slack_3{9};

%% Variable definitions
r_uav_opt = sdpvar(3,1,'full', 'real');
r0_low_opt = sdpvar(1,K,'full', 'real');
r1_low_opt = sdpvar(1,1,'full', 'real');
r0_up_opt = sdpvar(1,K,'full', 'real');
r1_up_opt = sdpvar(1,1,'full', 'real');
r2_up_opt = sdpvar(1,1,'full', 'real');
varrho_low_opt = sdpvar(1,K,'full', 'real');
varrho_up_opt = sdpvar(1,K,'full', 'real');
t_opt = sdpvar(1,K,'full','real');

x0 = sdpvar(K,1,'full','real');
x1 = sdpvar(1,1,'full','real');
z_opt = sdpvar(1,K,'full', 'real');
q_opt = sdpvar(1,1,'full', 'real');
p_slack = sdpvar(1,Na,'full', 'real');

F0_low_pow = sdpvar(1,K,'full', 'real');
F_qua_uu = sdpvar(1,K,'full', 'real');
F0_up_pow = sdpvar(1,K,'full', 'real');
F0_low_qua = sdpvar(1,K,'full', 'real');
F1_low_qua = sdpvar(1,K,'full', 'real');
F_low_bil = sdpvar(1,K,'full', 'real');
F_up_bil = sdpvar(1,K,'full', 'real');
%% Constraint
cons = [];

cons = [cons varrho_up_opt >= 0];
cons = [cons varrho_low_opt >= 0];
cons = [cons r0_up_opt >= 0];
cons = [cons r0_low_opt >= 10^(-7)];
% cons = [cons r0_low_opt >= 0];
cons = [cons r1_up_opt >= 0];
% cons = [cons r1_low_opt >= 0];
cons = [cons r1_low_opt >= 10^(-5)];
cons = [cons r2_up_opt >= 0];
cons = [cons q_opt >= 0];
cons = [cons z_opt >= 0];
cons = [cons x0 >= 0];
cons = [cons x1 >= 0];
% cons = [cons t_opt >= 0];

% 
cons = [cons r_uav_opt(1) <= 200];
cons = [cons r_uav_opt(2) <= r_ris(2)-10];
cons = [cons r_uav_opt(1) >= r_bs(1)];
cons = [cons r_uav_opt(2) >= r_bs(2)];
cons = [cons r_uav_opt(3) == 100];

%% (54) - (60)
for ii = 1:K
    %(54)
    F0_low_pow(ii) = get_Fpow(r0_low_opt(ii), -2/e0, r0_low_current(ii));
    cons = [cons norm(r_uav_opt-r_u(:,ii)) + F0_low_pow(ii) <= 0];

    %(56)
    F_qua_uu(ii) = get_Fqua(r_uav_opt, r_u(:,ii), r_uav_current);
    cons = [cons F_qua_uu(ii) <= 0];
    F0_up_pow(ii) = get_Fpow(r0_up_opt(ii), 4/e0, r0_up_current(ii));
    cons = [cons F_qua_uu(ii) <= x0(ii)];
    cons = [cons cone([1, 0.5*(-F0_up_pow(ii)-x0(ii))],0.5*(-F0_up_pow(ii)+x0(ii)))];
    
    %(59)
    c_0(ii) = p_current(ii)*(norm(hd_1(:,ii))^2)*sys.pathloss;
    c_1(ii) = p_current(ii)*(norm(hr(:,ii)'*diag(phi_current)*G_0)^2)*sys.pathloss;
    c_2(ii) = 2*p_current(ii)*real(hd_1(:,ii)'*hr(:,ii)'*diag(phi_current)*G_0)*sys.pathloss;
    
    F0_low_qua(ii) = get_Fqua(r0_low_opt(ii), 0, r0_low_current(ii));
    F1_low_qua(ii) = get_Fqua(r1_low_opt, 0, r1_low_current);       
    if (c_2(ii) > 0)
        F_low_bil(ii) = get_Fbil(r0_low_opt(ii), r1_low_opt, c_2(ii), r0_low_current(ii), r1_low_current);     
    else
        F_low_bil(ii) = get_Fbil(r0_up_opt(ii), r1_up_opt, c_2(ii), r0_up_current(ii), r1_up_current);     
    end
    cons = [cons varrho_low_opt(ii) + c_0(ii)*F0_low_qua(ii) + c_1(ii)*F1_low_qua(ii) + abs(c_2(ii))*F_low_bil(ii) <= 0];

     %(60)
    if(c_2(ii) > 0)
        F_up_bil(ii) = get_Fbil(r0_up_opt(ii), r1_up_opt, c_2(ii), r0_up_current(ii), r1_up_current);   
    else
        F_up_bil(ii) = get_Fbil(r0_low_opt(ii), r1_low_opt, c_2(ii), r0_low_current(ii), r1_low_current);   
    end
    cons = [cons varrho_up_opt(ii) >= c_0(ii)*r0_up_opt(ii)^2 + c_1(ii)*r1_up_opt^2 + abs(c_2(ii))*F_up_bil(ii)];
end

%(55)
F1_low_pow = get_Fpow(r1_low_opt, -2/e1, r1_low_current);
cons = [cons norm(r_uav_opt - r_ris) + F1_low_pow <= 0];

%(57)

F_qua_ur = get_Fqua(r_uav_opt, r_ris, r_uav_current);
cons = [cons F_qua_ur <= 0];
F1_up_pow = get_Fpow(r1_up_opt, 4/e1, r1_up_current);
cons = [cons F_qua_ur <= x1];
cons = [cons cone([1, 0.5*(-F1_up_pow-x1)],0.5*(-F1_up_pow+x1))];

%(58)
F2_pow = get_Fpow(r2_up_opt, -1/2, r2_current);
cons = [cons norm(r_uav_opt) + F2_pow <= 0];

%% (61)
for ii = 1:K
    sigma_k(ii) = (norm(hr(:,ii)'*I_active*diag(phi_current))^2)*sigma_r + sigma_u;
    
    cons = [cons cone([1, 0.5*(varrho_low_opt(ii)/sigma_k(ii)-z_opt(ii))], 0.5*(varrho_low_opt(ii)/sigma_k(ii)+z_opt(ii)))];
end

%% (63)
for ii = 1:K
    Theta_z(ii) = log(1+1/z_current(ii)) + 1/(1+z_current(ii));
    Xi_z(ii) = 1/(z_current(ii)*(z_current(ii)+1));
    
    cons = [cons (1-tau_current)*varphi_current(ii)*B*(Theta_z(ii) - Xi_z(ii)*z_opt(ii)) >= t_opt(ii)];
end

%% (64)
for ii = 1:K
    a(ii) = log(1+varrho_up_current(ii)/sigma_k(ii)) - (varrho_up_current(ii)/sigma_k(ii))/(1+varrho_up_current(ii)/sigma_k(ii));
    b(ii) = 1/(1+varrho_up_current(ii)/sigma_k(ii));
end

% %% (65)
cons = [cons cone([sqrt(sigma_u), 0.5*(q_opt - P_B*gB_1*r2_up_opt)],0.5*(q_opt + P_B*gB_1*r2_up_opt))];

%% (66)
Theta_q = log(1+1/q_current) + 1/(q_current+1);
Xi_q = 1/(q_current*(1+q_current));

%% (67)
sum_RA = 0;
for ii = 1:K
    sum_RA = sum_RA + (1-tau_current)*varphi_current(ii)*(a(ii) + b(ii)*(varrho_up_opt(ii)/sigma_k(ii))) ; 
end
gamma_B = (P_B*g_B)/sigma_u;

%%(68)
for ii = 1:Na
    c_3(ii)= (norm(G_0)^2)*(sys.pathloss)*sum(p_current);

    cons = [cons cone([abs(phi_current(ii))*sqrt(sigma_r), abs(phi_current(ii))*sqrt(c_3(ii))*r1_up_opt, ...
        0.5*(p_slack(ii)-1)], 0.5*(p_slack(ii)+1))];
end

cons = [cons sum(p_slack) <= P_ris];

%% Initialize
init_opt = sdpvar(1,K,'full', 'real');
for ii = 1:K
     if Init_2 == 1   
            cons = [cons tau_current*(Theta_q - Xi_q*q_opt) - sum_RA   >= init_opt(ii)];
            cons = [cons t_opt(ii) - r >= init_opt(ii)];
            cons = [cons init_opt(ii) <= 0];

     else
         cons = [cons sum_RA <= tau_current*(Theta_q - Xi_q*q_opt)];
         cons = [cons t_opt(ii) >= r];
     end
end  

% init_opt = sdpvar(1,1,'full', 'real');
% for ii = 1:K
%      if Init_2 == 1   
%             cons = [cons tau_current*(Theta_q - Xi_q*q_opt) - sum_RA   >= init_opt];
%             cons = [cons t_opt(ii) - 0.5 >= init_opt];
%             cons = [cons init_opt <= 0];
% 
%      else
%          cons = [cons sum_RA <= tau_current*(Theta_q - Xi_q*q_opt)];
%          cons = [cons t_opt(ii) >= 0.5];
%      end
% end 
%% Objective
if Init_2 == 0
   obj = sum(t_opt);
else
    % obj = init_opt;
      obj = sum(init_opt);
end

%% Solver
myops = sdpsettings('solver', 'mosek', 'verbose', 1, 'debug', 1, 'cachesolvers', 1);
diagnotics = solvesdp(cons, -obj, myops)


%% Update
AllSolution{1} = double(obj);
AllSolution{2} = double(r_uav_opt);

slack_variable_3{1} = double(r0_low_opt); 
slack_variable_3{2} = double(r1_low_opt); 
slack_variable_3{3} = double(r0_up_opt); 
slack_variable_3{4} = double(r1_up_opt); 
slack_variable_3{5} = double(r2_up_opt); 
slack_variable_3{6} = double(varrho_low_opt);
slack_variable_3{7} = double(varrho_up_opt); 
slack_variable_3{8} = double(z_opt); 
slack_variable_3{9} = double(q_opt); 

AllSolution{7} = slack_variable_3;

Status       = diagnotics.info;
end