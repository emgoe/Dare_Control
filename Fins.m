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




function beta_calc = theta_beta_mach_relation(T,M)
%% This program accepts 2 of the 3 values: sigma, beta, Mach values, and returns the third
%% NIR KALUSH
%% Implementation of Equation 9.23 from Anderson's Fundamentals of Aerodynamics 
%% ALL RIGHTS RESERVED to Mr. Kalush, University of Maryland, College Park

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NOTE!!!!!! YOU WILL NEED THE SYMBOLIC TOOLBOX TO RUN THIS PROGRAM. Check to see you have it installed  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    LHS = tan(T*pi/180);
    B = 1;
    RHS = inline('2*cot(B*pi/180)*((M*sin(B*pi/180))^2-1)/(M^2*(1.4+cos(2*B*pi/180))+2)','B','M');
    rhs_beta_1_mach_given = RHS(1,M);
    counter = 0;
    while (RHS(B,M) < LHS*0.98 | RHS(B,M) > LHS*1.02)
        if (B > 90 | B < 0 | counter >200)
            break;
        end
        if (RHS(B,M) < LHS*0.70)
            B=B+1;
        elseif(RHS(B,M) < LHS*0.90)
            B=B+0.5;
        elseif(RHS(B,M) < LHS*0.95)
            B=B+0.3;
        elseif(RHS(B,M) < LHS*0.98)
            B=B+0.1;
        end
        
        if(RHS(B,M) > LHS*1.03)
            B=B-0.1;
        elseif(RHS(B,M) > LHS*1.10)
            B = B-0.5;
        end
        counter = counter+1;
        lhs_curr = LHS;
        rhs_curr = RHS(B,M);
    
    end
    
    counter = 0;
    
    while (RHS(B,M) < LHS*0.99 | RHS(B,M) > LHS*1.01)
        if (B > 90 | B < 0 | counter > 100)
            break;
        end
        if (RHS(B,M) < LHS*0.98)
            B=B+0.05;
        elseif(RHS(B,M) < LHS*0.985)
            B=B+0.03;
        elseif(RHS(B,M) < LHS*0.99)
            B=B+0.01;
        end
        
        if(RHS(B,M) > LHS*1.01)
            B=B-0.01;
        elseif(RHS(B,M) > LHS*1.02)
            B = B-0.1;
        end
    
        counter = counter+1;
        lhs_curr = LHS;
        rhs_curr = RHS(B,M);
    end
    
    counter = 0;
    
    while (RHS(B,M) ~= LHS)
        if (B > 90 | B < 0 | counter > 40)
            break;
        end
        
        if (RHS(B,M) < LHS)
            B=B+0.005;
        end
        if(RHS(B,M) > LHS)
            B=B-0.005;
        end
        counter = counter+1;
        lhs_curr = LHS;
        rhs_curr = RHS(B,M);
    end
    
    if (B < 0)
        fprintf('The Beta angle for the given Mach = %2g and Theta = %2g vals is not in bounds (0,90): B = %2g degrees',M,T,B);
    elseif (B > 90)
        fprintf('The Shock is detached from the body \n\n for the given Mach = %2g and Theta = %2g \t B = %2g degrees - a useless quantity',M,T,B);
    else
        beta_calc = B;
    end
end