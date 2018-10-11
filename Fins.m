clear
clc

tic

theta = 1 ;                           % Fin offset in degrees - iterate
gamma = 1.4 ;                         % Heat capacity ratio - assumed cst for (dry) air
dcp = 0.1  ;                          % m, distance from Cp of fins to cg. Find better value
A_fin = (0.1+0.278) ;                 % m^2, find better value
I_xx = 4.485;                         % Moment of inertia of filled rocket

dm = 0.1;
M_start = 1.5;
M_fin = 4;
steps = (M_fin-M_start)/dm;          % Needs an integer

M_list = M_start:dm:M_fin;


PX_list = zeros(3,steps);
PY_list = zeros(3,steps);

cnt1 = 0;
cnt2 = 0;


% Situation one, increasing mach for 3 pressures (potential lookup table)
% Function calculates pressures based on oblique shock relations
% Flat plate assumption at given theta is assumed
for p_inf = [1e5, 7e4, 4e4]
    cnt1 = cnt1 + 1;
    cnt2 = 0;
    for M_inf = M_start:dm:M_fin
       
        cnt2 = cnt2 + 1;
        beta_calc = theta_beta_mach_relation(theta,M_inf);                      % Mach-Beta-Theta lookup
        M1_n = M_inf*sind(beta_calc);                                             
        [~, ~, P, ~, M2_n, ~ , ~ ] = flownormalshock(gamma, M1_n);              % Calc for normal component of oblique shock
        PX = p_inf * P;                                                         % Static pressure downstream of the shock 
%       MX = M2_n/sind(beta-theta)                                              % Optionally calculate downstream Mach number                                          
        PX_list(cnt1,cnt2) = PX;
        [~, nu, ~] = flowprandtlmeyer(1.4, M_inf, 'mach');
        nu1 = theta+nu;
        MY = flowprandtlmeyer(1.4,nu1,'nu');
        PY = (((1+((gamma-1)/2)*M_inf)/(1+((gamma-1)/2)*MY^2))^(gamma/(gamma-1)))*p_inf; % Isentropic flow expansion
        PY_list(cnt1,cnt2) = PY;
        
    end
end

DP_list = PX_list-PY_list;
F_fins = 4 * DP_list * A_fin;
T_fins = F_fins * dcp;

figure
plot(M_list,T_fins(1,:)/I_xx)
xlabel('Mach number')
ylabel('Torque [Nm]')
hold on
plot(M_list,T_fins(2,:))
plot(M_list,T_fins(3,:))

legend('p=SL', 'p=3km', 'p=7.2km')

hold off

fprintf('Done \n')

time = toc;
time = vpa(time,2); 

fprintf('%2g seconds to run, pretty damn long still...',time)



