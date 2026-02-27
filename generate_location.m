function [u] = generate_location (K)

% clc;
sys = config(1);

R = 50;
% K = sys.K;
u = zeros(3,K);
%% area of users
x_area = 150; 
y_area = 50;

for k0 = 1:K
    r = rand(1,1)*R;
    theta = rand(1,1)*2*pi;
    u(:,k0) = [x_area+r*cos(theta); y_area+ r*sin(theta); 0];
end

% Save user location to the user_location.mat
% save('user_location.mat', 'u');
% save('user_location_1.mat', 'u');

% load('user_location.mat', 'u');
% Plot user location
% figure;
% plot(sys.r_bs(1,1), sys.r_bs(2,1), 'rs', 'Color', 'r');
% text(sys.r_bs(1,1), sys.r_bs(2,1), ' BS', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
% hold on;
% plot(sys.r_ris(1,1), sys.r_ris(2,1), '^', 'Color','b');
% text(sys.r_ris(1,1), sys.r_ris(2,1), ' RIS', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
% hold on;
% plot(u(1,:), u(2,:), 'ro');
% hold on;
% plot(r_uav_current(1), r_uav_current(2), '*', 'Color', 'm');
% text(r_uav_current(1), r_uav_current(2), ' UAV', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

end