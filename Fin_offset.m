clear
clc

tic

theta = 5 ;                           % Fin offset in degrees - iterate
gamma = 1.4 ;                         % Heat capacity ratio - assumed cst for (dry) air
dcp = 0.3  ;                          % m, distance from Cp of fins to cg. Find better value
A_fin = (0.1+0.278) ;                 % m^2, find better value
I_xx = 4.485;                         % Moment of inertia of filled rocket
cd = 1;                               % Flat plate Cd of fins

dt = 20/100;


h = 7000;                             % Altitude at M = 1.5 ; Check value
omega = 0;
alpha = 0;

for t = linspace(20,40,101)
    [~, a, p_inf, rho] = atmosisa(h);
    M_inf = 0.15*t - 1.5 ;
    
    beta_calc = theta_beta_mach_relation(theta,M_inf);
    M1_n = M_inf*sind(beta_calc);
    [~, ~, P, ~, M2_n, ~ , ~ ] = flownormalshock(gamma, M1_n);
    PX = p_inf * P;
    
    [~, nu, ~] = flowprandtlmeyer(1.4, M_inf, 'mach');
    nu1 = theta+nu;
    MY = flowprandtlmeyer(1.4,nu1,'nu');
    PY = (((1+((gamma-1)/2)*M_inf)/(1+((gamma-1)/2)*MY^2))^(gamma/(gamma-1)))*p_inf; % Isentropic flow expansion
    
    F_fins = 4 * (PX-PY) * A_fin;
    T_fins = F_fins * dcp;
    alpha = T_fins / I_xx;
    omega = omega + alpha*dt;
    V_rot = omega*dt;
    D_rot = 0.5*rho*V_rot;
    
    T_fins = T_fins - D_rot;
    alpha = T_fins / I_xx;
    omega = omega + alpha*dt;
    
    V_inf = M_inf*a;
    h = h+V_inf*dt
 
end

toc

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