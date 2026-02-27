clc;
clear all;
clear classes;
addpath(genpath('./YALMIP-master'));
addpath(genpath('./Mosek'));

% run('./main_resoure_allocation.m');

%% Channel
sys = config(1);
%%--------------------------------------------------------------
K_max = 20;
N_max = 100;
Nt_max = 64;
NoMonte = 100;

filename_channel = strcat('Channel/channel_K-',num2str(K_max), 'N-', num2str(N_max), 'Nt-', num2str(N_max),num2str(NoMonte),'-channels' );

%%--------------------------------------------------------------
load(filename_channel);
Montein = sys.Montein;
gB = gB_cell{Montein};
G = G_cell{Montein}(1:sys.N, 1:sys.Nt);
hr = hr_cell{Montein}(1:sys.N, 1:sys.K);
hd = hd_cell{Montein}(1:sys.Nt, 1:sys.K);

r_uav = sys.r_uav;
channel = cell(4,1);
chan{1} = gB; chan{2} = G; chan{3} = hr; chan{4} = hd;
[gB, G, hr, hd] = update_channel(chan, r_uav);
channel{1} = gB; channel{2} = G; channel{3} = hr; channel{4} = hd;

MMR_opt_chain = [];
r_uav_chain = [];

[MMR_opt_chain, AllSolution, r_uav_chain] = Proposed_Alg_1_3_2(channel,chan, sys);

% if(length(MMR_opt_chain) > 31)
%     MMR_opt_chain = MMR_opt_chain(2:31);
% else
%     MMR_opt_chain = [MMR_opt_chain ones(2,31 - length(MMR_opt_chain)) * MMR_opt_chain(end)];
% end

figure;
plot(1:length(MMR_opt_chain), MMR_opt_chain*(1/log(2)), '-o');
hold on;
 

figure;

plot(r_uav_chain(1,1:length(r_uav_chain)), r_uav_chain(2,1:length(r_uav_chain)), '^', 'Color', 'g');
hold on;

plot(sys.r_ris(1), sys.r_ris(2), '*', 'Color','r');
hold on;
plot(sys.r_bs(1,1), sys.r_bs(2,1), 'rs', 'Color', 'r');
hold on;
plot(sys.r_u(1,:), sys.r_u(2,:), 'ro');

plot(r_uav_chain(1,end), r_uav_chain(2,end), '*', 'Color', 'y');
hold on;

plot(r_uav_chain(1,1), r_uav_chain(2,1), '^', 'Color', 'b');
hold on;

xlabel('X-coordinate (m)');
ylabel('Y-coordinate (m)');

legend('UAV','RIS', 'BS', 'UE');

% RA_NoRIS_conv = MMR_opt_chain;
% RA_NoRIS_ruav = [r_uav_chain(1,end),r_uav_chain(2,end)];
% save('Test/Convergence/P_B_36dBm-Puav_5W-K_4-N_48/P_B_36dBm-Puav_5W-K_4-N_48', 'RA_NoRIS_conv', 'RA_NoRIS_ruav', '-append');