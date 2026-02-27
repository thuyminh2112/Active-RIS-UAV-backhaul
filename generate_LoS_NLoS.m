function [gB, G, hr, hd] = generate_LoS_NLoS (K, N, Nt)

sys = config(1);
pathloss = sys.pathloss;
%% BS-UAV (channel power gain)
gB = pathloss;

%% UAV-RIS (Rician fading)
G_LoS = sqrt(sys.kappa/(sys.kappa+1))*LoS_channel(Nt, N, 'UAV-RIS');
G_NLoS = sqrt(1/(sys.kappa+1))*(randn(N, Nt)+ 1i*randn(N, Nt));
G = G_LoS + G_NLoS;

%% RIS-UE (Rician fading)
hr = zeros(N, K);
for k0 = 1:K
    hr_LoS = sqrt(sys.kappa/(sys.kappa+1))*LoS_channel(N, 1, 'RIS-UE');
    hr_NLoS = sqrt(1/(sys.kappa+1))*(randn(1,N)+ 1i*randn(1,N));
    hr(:,k0) = hr_LoS + hr_NLoS; 
end

%% UAV-UE (Rayleigh fading
hd = zeros(1, K);
hd = (randn(Nt, K)+ 1i*randn(Nt, K));

% save channels
% save('LoS_NLoS_channel.mat', 'gB', 'G', 'hr', 'hd');

end