function [AllSolution, Status] = Get_Solution_2 (sys, channel,Init_1,r_uav_current, phi_current, p_current, varphi_current, tau_current, slack_2)

rng('default'); 
rng(1);  

AllSolution = cell(7,1);
slack_variable_2 = cell(6,1);

K = sys.K;
N = sys.N; Na = sys.Na; a_max = sys.a_max; I_active =sys.I_active;
Pmax = sys.Pmax*(10^3);
P_B = sys.P_B*(10^3);
P_ris = sys.P_ris*(10^3);
sigma_u = sys.sigma_u*(10^3);
sigma_r = sys.sigma_r*(10^3);
p_current = p_current*(10^3);

% Pmax = sys.Pmax;
% P_B = sys.P_B;
% P_ris = sys.P_ris;
% sigma_u = sys.sigma_u;
% sigma_r = sys.sigma_r;
% p_current = p_current;

B = sys.B;

r = sys.R_QoS;

g_B = channel{1}; G = channel{2}; hr = channel{3}; hd = channel{4};

%%
mu_current = slack_2{1};
eta_current = slack_2{2};
zeta_current = slack_2{3};
nu_current = slack_2{4};

%% Variable definitions
phi_opt = sdpvar(1,N,'full', 'complex');
t_opt = sdpvar(1,K,'full', 'real');
mu_opt = sdpvar(1,K,'full', 'real');
nu_opt = sdpvar(1,K,'full', 'real');
eta_opt = sdpvar(1,K,'full', 'real');
zeta_opt = sdpvar(1,K,'full', 'real');

%% Constraint
cons = [];

%% (11h)
for ii = 1:Na
    real_part_a(ii) = real(phi_opt(ii));
    imag_part_a(ii) = imag(phi_opt(ii));
    cons = [cons cone([real_part_a(ii), imag_part_a(ii)], a_max)];
    
end

%% (11g)
% for ii = Na+1:N
%    if (sys.check_active == 0)
%         real_part_p(ii) = real(phi_opt(ii));
%         imag_part_p(ii) = imag(phi_opt(ii));
%         cons = [cons cone([real_part_p(ii), imag_part_p(ii)], 1)];
% 
%         real_p_slack(ii) = 2*(real(phi_current(ii))*real(phi_opt(ii))) - norm(real(phi_current(ii)))^2;
%         imag_p_slack(ii) = 2*(imag(phi_current(ii))*imag(phi_opt(ii))) - norm(imag(phi_current(ii)))^2;
%         cons = [cons real_p_slack(ii) + imag_p_slack(ii) >= 0];
%    else
%        cons = [cons phi_opt(ii) == 1];   
%    end
% end
if (sys.check_active == 0)
    for ii = Na+1:N
         real_part_p(ii) = real(phi_opt(ii));
         imag_part_p(ii) = imag(phi_opt(ii));
         cons = [cons cone([real_part_p(ii), imag_part_p(ii)], 1)];
            
         real_p_slack(ii) = 2*(real(phi_current(ii))*real(phi_opt(ii))) - norm(real(phi_current(ii)))^2;
         imag_p_slack(ii) = 2*(imag(phi_current(ii))*imag(phi_opt(ii))) - norm(imag(phi_current(ii)))^2;
         cons = [cons real_p_slack(ii) + imag_p_slack(ii) >= 0];
    end
end


%% (11i)
if Na > 0
    for ii = 1:Na
        xi(ii) = (sigma_r + norm(G(ii,:))^2)*sum(p_current);
    end
    x_opt = sdpvar(1,Na,'full', 'real');
    
    for ii = 1:Na
        cons = [cons cone([sqrt(xi(ii))*phi_opt(ii), 0.5*(x_opt(ii)-1)],0.5*(x_opt(ii)+1))];
    end
    
    cons = [cons sum(x_opt) <= P_ris];
end

for ii = 1:K
    h_slack(ii) = 2*(hd(:,ii)+hr(:,ii)'*diag(phi_current)*G)'*(hd(:,ii)+hr(:,ii)'*diag(phi_opt)*G) - norm(hd(:,ii)+hr(:,ii)'*diag(phi_current)*G)^2;
end

for ii = 1:K
    if Na > 0
        h_1(ii,:) = (hr(:,ii)'*I_active*diag(phi_opt));
        cons = [cons cone([h_1(ii,:)] , nu_opt(ii))];

        cons = [cons cone([sqrt(sigma_r)*nu_opt(ii), sqrt(sigma_u), 0.5*(p_current(ii)*h_slack(ii)-mu_opt(ii))], 0.5*(p_current(ii)*h_slack(ii)+mu_opt(ii)))];
    else
        cons = [cons cone([sqrt(sigma_u), 0.5*(p_current(ii)*h_slack(ii)-mu_opt(ii))], 0.5*(p_current(ii)*h_slack(ii)+mu_opt(ii)))];
    end
end

for ii = 1:K
    B_0(ii) = log(1+1/mu_current(ii))+1/(mu_current(ii)+1);

    B_1(ii) = 1/(mu_current(ii)*(mu_current(ii)+1));

    cons = [cons (1-tau_current)*varphi_current(ii)*B*(B_0(ii) - B_1(ii)*mu_opt(ii)) >= t_opt(ii)];
end

for ii = 1:K
    h_active(ii) = norm(diag(hr(:,ii)')*I_active*phi_current.')^2 + 2*(diag(hr(:,ii)')*I_active*phi_current.').'*(diag(hr(:,ii)')*I_active*phi_opt.' - diag(hr(:,ii)')*I_active*phi_current.');
    cons = [cons h_active(ii)*sigma_r + sigma_u >= eta_opt(ii)];
end

%% (41)
for ii = 1:K
      gamma_A(ii) = (p_current(ii)*norm(hd(:,ii)+hr(:,ii)'*diag(phi_current)*G)^2)/eta_current(ii);
      
      a_gamma(ii) = log(1+gamma_A(ii)) - gamma_A(ii)/(1+gamma_A(ii));
      b_gamma(ii) = 1/(1+gamma_A(ii));
      
      cons = [cons cone([sqrt(b_gamma(ii)*p_current(ii))*(hd(:,ii)+ hr(:,ii)'*diag(phi_opt)*G), ...
        0.5*(eta_opt(ii)-(zeta_opt(ii)-a_gamma(ii)))],0.5*(eta_opt(ii)+(zeta_opt(ii)-a_gamma(ii))))];
end


%% (42c)
gamma_B = (P_B*g_B)/sigma_u;
for ii = 1:K
    sum_zeta(ii) = ((1-tau_current)*varphi_current(ii)*zeta_opt(ii));
end



%%
cons = [cons mu_opt >= 0];
cons = [cons nu_opt >= 0];
cons = [cons eta_opt >= 0];
cons = [cons zeta_opt >= 0];
cons = [cons t_opt >= 0];

%% Initialize
init_opt = sdpvar(1,K,'full', 'real');
 if Init_1 == 1
    for ii = 1:K
        cons = [cons tau_current*log(1+gamma_B) - sum_zeta >= init_opt(ii)];
        cons = [cons t_opt(ii) - r >= init_opt(ii)];
        cons = [cons init_opt(ii) <= 0];
    end 
 else
     cons = [cons sum(sum_zeta) - tau_current*log(1+gamma_B) <= 0];
     for ii = 1:K
        cons = [cons t_opt(ii) >= r];
    end
 end

% init_opt = sdpvar(1,1,'full', 'real');
%  if Init_1 == 1
%     for ii = 1:K
%         cons = [cons tau_current*log(1+gamma_B) - sum_zeta >= init_opt];
%         cons = [cons t_opt(ii) - 0.5 >= init_opt];
%         cons = [cons init_opt <= 0];
%     end 
%  else
%      cons = [cons sum(sum_zeta) - tau_current*log(1+gamma_B) <= 0];
%      for ii = 1:K
%         cons = [cons t_opt(ii) >= 0.5];
%     end
%  end
%% Objective
if Init_1 == 0
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
AllSolution{3} = double(phi_opt);
slack_variable_2{1} = double(mu_opt);
slack_variable_2{2} = double(eta_opt);
slack_variable_2{3} = double(zeta_opt);
slack_variable_2{4} = double(nu_opt);
AllSolution{7} = slack_variable_2;

Status       = diagnotics.info;
end