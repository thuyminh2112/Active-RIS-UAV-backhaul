function [MMR_opt_chain, AllSolution, r_uav_chain] = Proposed_Alg_1_3_2(channel,chan, sys)

rng('default'); 
rng(1); 

[r_uav_current, phi_current, tau_current, p_current, varphi_current, slack_1, slack_2, slack_3] = initilize(channel,chan, sys);
% load('phi_current.mat');
% load('phi_current_1.mat', 'phi_current');

MMR_opt = 0;
MMR_opt_currrent = 0;
MMR_opt_chain = [];
r_uav_chain = [sys.r_uav];
AllSolution = cell(5,1);
%%
Init = 1;
Init_1 = 1;
Init_2 = 1;
iter = 0;
status_1 = 1;
status_2 = 1;
status_3 = 1;
dem1 = 0;
dem2 = 0;
dem3 = 0;
dem_init_1 = 0; dem_init_2 = 0; dem_init_3 = 0;

while(1)
    if (status_1 == 1)
        if Init == 1
           fprintf('Inilizing- Iteration %d -th \n', iter);
        else
           fprintf('Updating_1- Iteration %d -th \n', iter);
        end
        [AllSolution, Status] = Get_Solution_1(sys, channel, Init, r_uav_current, phi_current, p_current, varphi_current, tau_current, slack_1);
        %% Update
        % if (~isempty(findstr(Status,'Successfully')))
            MMR_opt = AllSolution{1};
            p_current = AllSolution{4};
            varphi_current = AllSolution{5};
            tau_current = AllSolution{6};
            slack_1 = AllSolution{7};

            if (MMR_opt >= -1e-5 &&  MMR_opt <= 0.0001)
               Init = 0;
            end
            % if (MMR_opt >= -1e-5 &&  MMR_opt <= 0.0001)
            %    dem_init_1 = dem_init_1 + 1;
            %    if dem_init_1 >= 3
            %     Init = 0;
            %    end
            % end
            dem1 = dem1+1;   
        % end
    end
    if (status_3 == 1 && dem1 >= 1) 
        if Init_2 == 1
           fprintf('Inilizing- Iteration %d -th \n', iter);
        else
           fprintf('Updating_3- Iteration %d -th \n', iter);
        end
        [AllSolution, Status] = Get_Solution_3_1 (sys, channel, chan, Init_2, r_uav_current, phi_current, p_current, varphi_current, tau_current, slack_3);
        %% Update
        if (~isempty(findstr(Status,'Successfully')))
            MMR_opt = AllSolution{1};
            r_uav_current = AllSolution{2};           
            slack_3 = AllSolution{7};
            [gB, G, hr, hd] = update_channel(chan, r_uav_current);
            channel{1} = gB; channel{2} = G; channel{3} = hr; channel{4} = hd;
            if (MMR_opt >= -1e-5 &&  MMR_opt <= 0.0001)
               Init_2 = 0;
            end
            % if (MMR_opt >= -1e-5 &&  MMR_opt <= 0.0001)
            %    dem_init_3 = dem_init_3 + 1;
            %    if dem_init_3 >= 3
            %     Init_2 = 0;
            %    end
            % end            
        end
        dem3 = dem3+1; 
    end    
    if (status_2 == 1 && dem3 >= 1)  
        if Init_1 == 1
           fprintf('Inilizing- Iteration %d -th \n', iter);
        else
           fprintf('Updating_2- Iteration %d -th \n', iter);
        end
        [AllSolution, Status] = Get_Solution_2 (sys, channel, Init_1, r_uav_current, phi_current, p_current, varphi_current, tau_current, slack_2);
        %% Update
        % if (~isempty(findstr(Status,'Successfully')))
            MMR_opt = AllSolution{1};
            phi_current = AllSolution{3};
            slack_2 = AllSolution{7};

            if (MMR_opt >= -1e-5 &&  MMR_opt <= 0.0001)
               Init_1 = 0;
            end
            % if (MMR_opt >= -1e-5 &&  MMR_opt <= 0.0001)
            %    dem_init_2 = dem_init_2 + 1;
            %    if dem_init_2 >= 3
            %     Init_1 = 0;
            %    end
            % end
            dem2 = dem2+1;

            r_uav_chain = [r_uav_chain r_uav_current];
            MMR_opt_chain = [MMR_opt_chain MMR_opt];

            % save('Solution_1.mat', 'p_current', 'varphi_current', 'tau_current', 'phi_current'); 
        % end
    end  

    %% Check convergence
    if((MMR_opt - MMR_opt_currrent <= 1e-3  && iter >= 10) || iter >= 30)
        break;
    end
    if (MMR_opt >= MMR_opt_currrent)
        MMR_opt_currrent = MMR_opt;
    end
    iter = iter + 1;
end

if(length(MMR_opt_chain) > 31)
    MMR_opt_chain = MMR_opt_chain(2:31);
else
    MMR_opt_chain = [MMR_opt_chain(2:end) ones(1,31 - length(MMR_opt_chain)) * MMR_opt_chain(end)];
end

end

