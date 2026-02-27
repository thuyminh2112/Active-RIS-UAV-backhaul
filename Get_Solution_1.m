function [AllSolution, Status] = Get_Solution_1 (sys, channel, Init, r_uav_current, phi_current, p_current, varphi_current, tau_current, slack_1)

rng('default'); 
rng(1);  

AllSolution = cell(7,1);
slack_variable_1 = cell(6,1);

K = sys.K;
N = sys.N; Na = sys.Na; I_active = sys.I_active;
Pmax = sys.Pmax;
P_ris = sys.P_ris;
P_B = sys.P_B;
sigma_u = sys.sigma_u;
sigma_r = sys.sigma_r;

B = sys.B;

r = sys.R_QoS;

g_B = channel{1}; G = channel{2}; hr = channel{3}; hd = channel{4};

%% Khoi tao cac bien ban dau
v_current = slack_1{1};
varepsilon_current = slack_1{2};
delta_current = slack_1{3};
omega_current = slack_1{4};
u_current = slack_1{5};
beta_current = slack_1{6};

%% Variable definitions
p_opt = sdpvar(1,K,'full', 'real');
varphi_opt = sdpvar(1,K,'full', 'real');
tau_opt = sdpvar(1,1,'full', 'real');
v_opt = sdpvar(1,K,'full', 'real');
varepsilon_opt = sdpvar(1,K,'full', 'real');
delta_opt = sdpvar(1,K,'full', 'real');
omega_opt = sdpvar(1,K,'full', 'real');
u_opt = sdpvar(1,K,'full', 'real');
beta_opt = sdpvar(1,K,'full', 'real');
t_opt = sdpvar(1,K,'full', 'real');

%% Constraint
cons = [];

%% (11b)
for ii = 1:K
   cons = [cons p_opt(ii) >= 0];
end
cons = [cons sum(p_opt) <= Pmax];

%% (11c)
for ii = 1:K
    cons = [cons varphi_opt(ii) <= 1];
    cons = [cons varphi_opt(ii) >= 0];
end

%% (11d)
cons = [cons sum(varphi_opt) <= 1];

%% (11e)
cons = [cons tau_opt >= 0];
cons = [cons tau_opt <= 1];

%% (11i)
if Na > 0
    for ii = 1:Na
        xi(ii) = (sigma_r + norm(G(ii,:))^2)*sum(p_opt);
        sum_abs(ii) = (abs(phi_current(ii)))^2;
        pris_opt(ii) = sum_abs(ii)*xi(ii);
    end
    cons = [cons sum(pris_opt) <= P_ris];
end

%% (20)
A_0 = zeros(1,K);
A_1 = zeros(1,K);
A_2 = zeros(1,K);
for ii = 1:K
    A_0(ii) = (2*log(1+1/omega_current(ii)))/v_current(ii) + 1/(v_current(ii)*(omega_current(ii)+1));
    A_1(ii) = 1/(v_current(ii)*omega_current(ii)*(omega_current(ii)+1));
    A_2(ii) = log(1+1/omega_current(ii))/(v_current(ii)^2);
    
    cons = [cons B*(A_0(ii) - A_1(ii)*omega_opt(ii) - A_2(ii)*v_opt(ii)) >= t_opt(ii)];
     
end

%% (16), (17), (18)
for ii = 1:K
    cons = [cons 2*(delta_current(ii)'*delta_opt(ii)) - (delta_current(ii))^2 >= varepsilon_opt(ii)];
    
    cons = [cons cone([1, 0.5*(v_opt(ii)-varepsilon_opt(ii))], 0.5*(v_opt(ii)+varepsilon_opt(ii)))];           
    
    cons = [cons cone([delta_opt(ii), 0.5*((1-tau_opt)-varphi_opt(ii))], 0.5*((1-tau_opt)+varphi_opt(ii)))];    
end


%% (15)
a = zeros(1,K);
for ii = 1:K
    a(ii) = (norm(hd(:,ii)+hr(:,ii)'*diag(phi_current)*G)^2)/(norm(hr(:,ii)'*I_active*diag(phi_current))^2*sigma_r + sigma_u);
    
    cons = [cons cone([1/sqrt(a(ii)), 0.5*(p_opt(ii)-omega_opt(ii))], 0.5*(p_opt(ii)+omega_opt(ii)))]; 
end


%% (23)
b = zeros(1,K);
c = zeros(1,K);
for ii = 1:K
    b(ii) = sqrt(varphi_current(ii)/(2*(1-tau_current)));
    c(ii) = sqrt((1-tau_current)/(2*varphi_current(ii)));
end

for ii = 1:K
    cons = [cons cone([b(ii)*(1-tau_opt), c(ii)*varphi_opt(ii), 0.5*(u_opt(ii)-1)], 0.5*(u_opt(ii)+1))];
end


%% (25)
a_pk = zeros(1,K);
b_pk = zeros(1,K);
for ii = 1:K
    a_pk(ii) = log(1+a(ii)*p_current(ii)) - (a(ii)*p_current(ii))/(1+a(ii)*p_current(ii));
    b_pk(ii) = a(ii)/(1+a(ii)*p_current(ii));
    t1(ii) = (a_pk(ii) + b_pk(ii)*p_current(ii))/(2*u_current(ii));
    t2(ii) = u_current(ii)/(2*(a_pk(ii)+b_pk(ii)*p_current(ii)));

    cons = [cons cone([sqrt(t1(ii))*u_opt(ii), sqrt(t2(ii))*(a_pk(ii)+b_pk(ii)*p_opt(ii)), 0.5*(beta_opt(ii)-1)],0.5*(beta_opt(ii)+1))];
end


gamma_B = (P_B*g_B)/sigma_u;

%% 
for ii = 1:K
    cons = [cons v_opt(ii) >= 0];
    
    cons = [cons omega_opt(ii) >= 0];
    
    cons = [cons varepsilon_opt(ii) >= 0];
    
    cons = [cons delta_opt(ii) >= 0];

    cons = [cons u_opt(ii) >= 0];
end

% cons = [cons t_opt >= 0];
%% Initialize
init_opt = sdpvar(1,K,'full', 'real');
 if Init == 1
    for ii = 1:K
        cons = [cons tau_opt*log(1+gamma_B) - sum(beta_opt) >= init_opt(ii)];
        cons = [cons t_opt(ii) - r >= init_opt(ii)];
        cons = [cons init_opt(ii) <= 0];
    end  
 else
    cons = [cons sum(beta_opt) <= tau_opt*log(1+gamma_B)];
    for ii = 1:K  
        cons = [cons t_opt(ii) >= r];
    end
 end

% init_opt = sdpvar(1,1,'full', 'real');
%  if Init == 1
%     for ii = 1:K
%         cons = [cons tau_opt*log(1+gamma_B) - sum(beta_opt) >= init_opt];
%         cons = [cons t_opt(ii) - 0.5 >= init_opt];
%         cons = [cons init_opt <= 0];
%     end  
%  else
%     cons = [cons sum(beta_opt) <= tau_opt*log(1+gamma_B)];
%     for ii = 1:K  
%         cons = [cons t_opt(ii) >= 0.5];
%     end
%  end


%% Objective
if Init == 0
    obj = sum(t_opt(:));
else
    % obj = init_opt;
     obj = sum(init_opt);
end

%% Solver
myops = sdpsettings('solver', 'mosek', 'verbose', 1, 'debug', 1, 'cachesolvers', 1);
diagnotics = solvesdp(cons, -obj, myops)

%% Update
AllSolution{1} = double(obj);
AllSolution{2} = double(r_uav_current);
AllSolution{3} = double(phi_current);
AllSolution{4} = double(p_opt);
AllSolution{5} = double(varphi_opt);
AllSolution{6} = double(tau_opt);

slack_variable_1{1} = double(v_opt); 
slack_variable_1{2} = double(varepsilon_opt); 
slack_variable_1{3} = double(delta_opt); 
slack_variable_1{4} = double(omega_opt);
slack_variable_1{5} = double(u_opt); 
slack_variable_1{6} = double(beta_opt);

AllSolution{7} = slack_variable_1;

Status       = diagnotics.info;
end



