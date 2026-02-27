clc;
clear all;
clear classes;
addpath(genpath('./YALMIP-master'));
addpath(genpath('./Mosek'));

sys = config(1);
%%--------------------------------------------------------------
K_max = 20;
N_max = 100;
Nt_max = 64;
NoMonte = 100;

filename_channel = strcat('Channel/channel_K-',num2str(K_max), 'N-', num2str(N_max), 'Nt-', num2str(N_max),num2str(NoMonte),'-channels' );
load('Channel/channel_K-20N-100Nt-100100-channels.mat', 'gB_cell', 'G_cell', 'hr_cell', 'hd_cell', 'UE_position_cell');

% gB_cell = cell(NoMonte, 1);
% G_cell = cell(NoMonte, 1);
% hr_cell = cell(NoMonte, 1);
% hd_cell = cell(NoMonte, 1);
% UE_position_cell = cell(NoMonte, 1);

for ii = 101:500
    UE_position_cell{ii} = generate_location(K_max);

    [gB, G, hr, hd] = generate_LoS_NLoS(K_max, N_max, Nt_max);

    gB_cell{ii} = gB;
    G_cell{ii} = G;
    hr_cell{ii} = hr;
    hd_cell{ii} = hd;
end

save(filename_channel, 'gB_cell', 'G_cell', 'hr_cell', 'hd_cell', 'UE_position_cell');

% load('Channel/channel_K-20N-100Nt-100100-channels.mat', 'gB_cell', 'G_cell', 'hr_cell', 'hd_cell', 'UE_position_cell');
% 
% for ii = 1:100
%     UE_position_cell{ii} = UE_position_cell{1};
% end
% 
% save('Channel/channel_K-20N-100Nt-100100-channels_1.mat','gB_cell', 'G_cell', 'hr_cell', 'hd_cell', 'UE_position_cell');