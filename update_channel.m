function [gB, G, hr, hd] = update_channel(chan, r_uav)
    
    sys = config(1);
    K = sys.K; N = sys.N; Nt = sys.Nt;
    pathloss = sys.pathloss;
    e0 = sys.e0; % UAV-UE pathloss exponent
    e1 = sys.e1; % UAV-RIS pathloss exponent
    e2 = sys.e2; % RIS-UE pathloss exponent
    sigma_u = sys.sigma_u;
    sigma_r = sys.sigma_r;
    r_bs = sys.r_bs;
    r_ris = sys.r_ris;
    r_u = sys.r_u;
    gB = zeros(1,1); G = zeros(N,K);hr = zeros(N,Nt);hd = zeros(Nt,K);
    gB_1 = chan{1}; G_1 = chan{2}; hr_1 = chan{3}; hd_1 = chan{4};
    
    
    %% gB -- BS-UAV
    gB = gB_1*norm(r_uav-r_bs)^(-2);
    
    %% G -- UAV-RIS
    PL0 = pathloss*norm(r_uav-r_ris)^(-e1);
    G = sqrt(PL0)*G_1;

    %% hr -- RIS-UE
    for k0 = 1:K
        PL2 = pathloss*norm(r_ris - r_u(:,k0))^(-e2);
        hr(:,k0) = sqrt(PL2)*hr_1(:,k0);
    end

    %% hd -- UAV-UE
    for k0 = 1:K
        PL1 = pathloss*norm(r_uav - r_u(:,k0))^(-e0);
        hd(:,k0) = sqrt(PL1)*hd_1(:,k0);
    end
end