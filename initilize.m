function [r_uav_current, phi_current, tau_current, p_current, varphi_current, slack_1, slack_2, slack_3] = initilize(channel, chan, sys)
rng('default'); 
rng(1);  

K = sys.K;
N = sys.N;
sigma_u = sys.sigma_u;
sigma_r = sys.sigma_r;
Pmax = sys.Pmax;
Na = sys.Na; a_max = sys.a_max; I_active = sys.I_active; I_passive = sys.I_passive; NoRIS = sys.NoRIS;

%% Channel
g_B = channel{1}; G = channel{2}; hr = channel{3}; hd = channel{4};

%% Khoi tao bien ban dau
% initialize UAV at the center of the area
r_uav_current = sys.r_uav;

% Phi
phase_active= exp(1i*(rand(1,Na)*2*pi));
abs_active = rand(1,Na)*a_max;
phi_active = zeros(1, Na);
for ii = 1:Na
    phi_active(ii) = abs_active(ii)*phase_active(ii);
end
phi_passive = exp(1i*(rand(1,(N-Na))*2*pi));
phi_current = [phi_active phi_passive];
if NoRIS == 1
    phi_current = 0*rand(1,N);
end

%% Khoi tao bien toi uu
p_current = rand(1,K)*Pmax/K;
varphi_current = rand(1,K)/K;
% varphi_current = ones(1,K)/K;
tau_current = rand(1,1);
% tau_current = (0.1*K)/(sys.B*log(1+(sys.P_B*g_B)/sigma_u));
% tau_current = 0.3;

a = zeros(1,K);
for ii = 1:K
    a(ii) = (norm(hd(:,ii)+ hr(:,ii)'*diag(phi_current)*G)^2)/(norm(hr(:,ii)'*I_active*diag(phi_current))^2*sigma_r + sigma_u);
end

for ii = 1:K
    p_current(ii) = (exp(0.1/((1-tau_current)*varphi_current(ii)*sys.B))-1)/a(ii);
end

%% Slack sub-problem 1
v_current = zeros(1,K);
for ii = 1:K
    v_current(ii) = 1/((1-tau_current)*varphi_current(ii));
end

varepsilon_current = zeros(1,K);
for ii = 1:K
    varepsilon_current(ii) = 1/v_current(ii);
end

delta_current = zeros(1,K);
for ii = 1:K
    delta_current(ii) = sqrt((1-tau_current)*varphi_current(ii));
end

omega_current = zeros(1,K);

for ii = 1:K
    omega_current(ii) = 1/(a(ii)*p_current(ii));
end

u_current = zeros(1,K);
for ii = 1:K
    u_current(ii) = (1-tau_current)*varphi_current(ii);
end


beta_current = zeros(1,K);
for ii = 1:K
   a_p(ii) = log(1+a(ii)*p_current(ii)) - (a(ii)*p_current(ii))/(1+a(ii)*p_current(ii));
   b_p(ii) = a(ii)/(1+a(ii)*p_current(ii));

   beta_current(ii) = u_current(ii)*(a_p(ii)+b_p(ii)*p_current(ii)); 
end

slack_1 = cell(6,1);
slack_1{1} = v_current; slack_1{2} = varepsilon_current; slack_1{3} = delta_current; slack_1{4} = omega_current;
slack_1{5} = u_current; slack_1{6} = beta_current;

%% Slack sub-problem 2
mu_current = zeros(1,K);
eta_current = zeros(1,K);
zeta_current = zeros(1,K);
nu_current = zeros(1,K);

%%
p_current = p_current*(10^3);
sigma_r = sigma_r*(10^3);
sigma_u = sigma_u*(10^3);

for ii = 1:K
    eta_current(ii) = norm(hr(:,ii)'*I_active*diag(phi_current))^2*sigma_r + sigma_u;
    nu_current(ii) = norm(hr(:,ii)'*I_active*diag(phi_current));
    mu_current(ii) = (norm(hr(:,ii)'*I_active*diag(phi_current))^2*sigma_r + sigma_u)/(p_current(ii)*norm(hd(:,ii)+hr(:,ii)'*diag(phi_current)*G)^2);
end
for ii = 1:K  
    gamma_A(ii) = (p_current(ii)*norm(hd(:,ii)+hr(:,ii)'*diag(phi_current)*G)^2)/eta_current(ii);
end

for ii = 1:K
    a_gamma(ii) = log(1+gamma_A(ii)) - gamma_A(ii)/(1+gamma_A(ii));
    b_gamma(ii) = 1/(1+gamma_A(ii));
end

for ii = 1:K
    zeta_current(ii) = a_gamma(ii) + b_gamma(ii)*(p_current(ii)*norm(hd(:,ii)+ hr(:,ii)'*diag(phi_current)*G)^2)/eta_current(ii);
end

slack_2 = cell(6,1);
slack_2{1} = mu_current;
slack_2{2} = eta_current;
slack_2{3} = zeta_current;
slack_2{4} = nu_current;


%% Slack sub-problem 3
gB_1 = chan{1}; G_0 = chan{2}; hr_1 = chan{3}; hd_1 = chan{4};

%%
p_current = p_current*(10^-3);
sigma_r = sys.sigma_r;
sigma_u = sys.sigma_u;

r0_low_current = zeros(1,K);
r1_low_current = zeros(1,1);
r0_up_currrent = zeros(1,K);
r1_up_current = zeros(1,1);
r2_current = zeros(1,1);
varrho_low_current = zeros(1,K);
varrho_up_current = zeros(1,K);
z_current = zeros(1,K);
q_current = zeros(1,1);

c_0 = zeros(1,K); c_1 = zeros(1,K); c_2 = zeros(1,K);
for ii = 1:K
    c_0(ii) = p_current(ii)*(norm(hd_1(:,ii))^2)*sys.pathloss;
    c_1(ii) = p_current(ii)*(norm(hr(:,ii)'*diag(phi_current)*G_0)^2)*sys.pathloss;
    c_2(ii) = 2*p_current(ii)*real(hd_1(:,ii)'*hr(:,ii)'*diag(phi_current)*G_0)*sys.pathloss;
end

for ii = 1:K
    r0_low_current(ii) = norm(sys.r_uav-sys.r_u(:,ii))^(-sys.e0/2);
    r0_up_currrent(ii) = r0_low_current(ii);
end

r1_low_current = norm(sys.r_uav-sys.r_ris)^(-sys.e1/2);
r1_up_current = r1_low_current;
r2_current = norm(sys.r_uav)^(-2);

sigma_k = zeros(1,K);
for ii = 1:K
    varrho_low_current(ii) = c_0(ii)*r0_low_current(ii)^2 + c_1(ii)*r1_low_current^2 + c_2(ii)*r0_low_current(ii)*r1_low_current;
    varrho_up_current(ii) = varrho_low_current(ii);
    sigma_k(ii) = (norm(hr(:,ii)'*I_active*diag(phi_current))^2)*sigma_r + sigma_u;
    z_current(ii) = sigma_k(ii)/varrho_up_current(ii);
end
q_current = sigma_u/(sys.P_B*sys.pathloss*r2_current);

slack_3 = cell(9,1);
slack_3{1} = r0_low_current; slack_3{2} = r1_low_current; slack_3{3} = r0_up_currrent; 
slack_3{4} = r1_up_current; slack_3{5} = r2_current; slack_3{6} = varrho_low_current;
slack_3{7} = varrho_up_current; slack_3{8} = z_current; slack_3{9} = q_current;
end