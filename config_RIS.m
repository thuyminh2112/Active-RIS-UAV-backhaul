function [N, Na, a_max, I_active, I_passive, check_active] = config_RIS(RIS_type, N_RIS)
    
    % N, N_RIS: so luong phan tu RIS
    % Na : so luong phan tu active RIS
    % a_max: bien do lon nhat cua active RIS
    check_active = 0;
    %%
    if RIS_type == 0 % No RIS
        N = 0; Na = 0;
    elseif RIS_type == 1 %% Passive RIS
        N = N_RIS; Na = 0;
    elseif RIS_type == 2 %% Hybrid RIS
        N = N_RIS; Na = N_RIS*(1/32);
        % N = N_RIS; Na = 2;
    elseif RIS_type == 3 %% Hybrid RIS
        N = N_RIS; Na = N_RIS*(1/4);
        % N = N_RIS; Na = 4;
    elseif RIS_type == 4 %% Hybrid RIS
        N = N_RIS; Na = N_RIS*(0.375);
        % N = N_RIS; Na = 8;
    elseif RIS_type == 5 %% Hybrid RIS
        N = N_RIS; Na = N_RIS*(1/2);
    elseif RIS_type == 6 %% Hybrid RIS
        N = N_RIS; Na = N_RIS*(0.625);
    elseif RIS_type == 7 %% Hybrid RIS
        N = N_RIS; Na = N_RIS*(3/4);
    elseif RIS_type == 8 %% Hybrid RIS
        N = N_RIS; Na = N_RIS*(0.875);
    else %active RIS
        N = N_RIS; Na = N;
        check_active = 1;
    end

    if Na > 0
        a_max = 20;
    else
        a_max = 1;
    end
        diag_elements = [ones(1, Na), zeros(1, N - Na)];
        I_active = diag(diag_elements);
        I_passive = diag(ones(1, N))-I_active;
end