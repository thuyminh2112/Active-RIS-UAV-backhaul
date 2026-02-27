function sys = config(UE_fix)
% rng('default'); 
% rng(1); 

sys.K_max = 20;
sys.N_max = 100;
sys.Nt_max = 64;
sys.NoMonte= 100;

%% System parameters
sys.K = 4;  % UEs
sys.N = 32; % RIS elements
sys.Nt = 1; % UAV antennas
sys.BS_anten = 1; %BS antennas

%% Large-scale parameters
sys.e0 = 3.2; % UAV-UE pathloss exponent
sys.e1 = 2.0; % UAV-RIS pathloss exponent
sys.e2 = 2.2; % RIS-UE pathloss exponent

sys.kappa = db2pow(10); %Rician factor

%% Pathloss
pathloss_dB = -30;
sys.pathloss = 10^(pathloss_dB/10);

%% Power and noise parameters
% Pmax_dB = 36;
% sys.Pmax = 10^(Pmax_dB/10-3);
sys.Pmax = 1;
P_B_dB = 36;
sys.P_B = 10^(P_B_dB/10-3);

P_ris_dB = 30;
sys.P_ris = 10^(P_ris_dB/10-3);

sigma_u_dBm = -80;
sys.sigma_u = 10^(sigma_u_dBm/10-3);
sys.sigma_r = 1.58*sys.sigma_u;

%% Positions of BS, UAV, RIS and UEs
sys.r_bs = [0; 0; 0]; % position of BS
sys.r_uav = [50; 50; 100];
sys.r_ris = [150; 100; 50]; % position of RIS
filename_channel = strcat('Channel/channel_K-',num2str(sys.K_max), 'N-', num2str(sys.N_max), 'Nt-', num2str(sys.N_max),num2str(sys.NoMonte),'-channels' );
load(filename_channel);
sys.Montein = 1;
sys.r_u = UE_position_cell{sys.Montein}(1:3, 1:sys.K);

%% Bandwidth
sys.B = 10;  %MHz

%% QoS
sys.R_QoS = 0.5;

%% RIS
[N, Na, a_max, I_active, I_passive, check_active] = config_RIS(9, sys.N);
sys.check_active = check_active;
if N > 0
    sys.Na = Na;
    sys.a_max = a_max;
    sys.I_active = I_active;
    sys.I_passive = I_passive;
    sys.NoRIS = 0;
else
    sys.Na = 0;
    sys.a_max = 0;
    sys.I_active = 0;
    sys.I_passive = 0;
    sys.NoRIS = 1;
end
end