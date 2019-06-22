
clc;
syms t % time [second]

% Parameter Symbols

syms omega % Angular velocity [rad/s]
syms R     % Radius [m]
syms rho   % Density [kg/m3]
syms c     % Chord Length [m]
syms cla   % CL_alpha : Lift curve slope
syms e     % Hinge Offset 
syms m    % mass per length [kg/m]
syms Ib     % Blade Inertia [kg*m2]

% Flapping Symbols

syms beta betaD betaDD  % Flapping derivatives
syms beta0 beta1c beta1s beta1cD beta1cDD beta1sD beta1sDD
syms y                 % Along radius axis

% Angular & Translational Velocities of helicopter

syms u v w p q r

% Non dimensional parameters

syms muX muY muZ 

% Input Symbols

syms theta theta0 theta1c theta1s
syms psi 

psi = omega*t;

 %% Beta Derivatives
beta = beta0 + beta1c * cos(psi) + beta1s * sin(psi);
 betaD = beta1cD * cos(psi) - omega*sin(psi)*beta1c + beta1sD*sin(psi) + omega*cos(psi)*beta1s;
 betaDD = + beta1sDD * sin(psi)  ...
          + beta1cDD * cos(psi)  ...
          + 2.0 * omega * beta1sD * cos(psi)  ...
          - 2.0 * omega * beta1cD * sin(psi)  ...
          - omega^2 * beta1s * sin(psi)  ...
          - omega^2 * beta1c * cos(psi);


% Inflow - Thrust Symbols

syms inflow thrust

% Blade Symbols

syms Ut Up phi

%% Aerodynamic Moment

muX = u / (omega*R);
muY = v / (omega*R);
muZ = w / (omega*R);

theta = theta0 + theta1c * cos(psi) + theta1s * sin(psi);


 % Cross section velocities (CCW rotation)

Ut = omega * y + u * sin(psi) + v * cos(psi);
Up = y * betaD + muX * omega * R * beta * cos(psi) - muZ * omega * R + inflow * omega * R;
 
 % Cross section angles
 
 phi = Up / Ut;
 aoa = theta - phi;
 
 % Lift And Drag 
 
 dL = 0.5 * rho * Ut^2 * c * cla * aoa  ; % * dy
 
 momentAero =int(dL * y , y , e*R , R) ;

 momentAero = expand(momentAero);
 
 %% Trigonometric conversion in order to consider first harmonics only
 
momentAero = subs(momentAero,cos(psi)*sin(psi),sin(2*psi)/2);
momentAero = subs(momentAero,sin(psi)^2,(1-cos(2*psi))/2);
momentAero = subs(momentAero,cos(psi)^2,(1+cos(2*psi))/2);
momentAero = subs(momentAero,sin(psi)^3,(3*sin(psi)-sin(3*psi))/4);
momentAero = subs(momentAero,cos(psi)^3,(3*cos(psi)+cos(3*psi))/4);
momentAero = subs(momentAero,(sin(psi)^2*cos(psi)),((cos(psi)+cos(3*psi))/4));
momentAero = subs(momentAero,(cos(psi)^2*sin(psi)),((sin(psi)+sin(3*psi))/4));

%% Neglecting higher order harmonics

momentAero = subs(momentAero,sin(2*psi),0);
momentAero = subs(momentAero,cos(2*psi),0);
momentAero = subs(momentAero,sin(3*psi),0);
momentAero = subs(momentAero,cos(3*psi),0);
momentAero = simplify(momentAero);

 %% Centrifugal Moment
 dFcf = m*omega^2*y;
 momentCF =int(dFcf*y*beta,y,0,R);
 
 %% Inertial Moment
 
 momentInertial = Ib * betaDD;
 
%% Effective Moments

totalMoment = momentInertial - momentAero + momentCF;
 
totalMoment = expand(totalMoment);
totalMoment = collect(totalMoment,sin(psi));
totalMoment = collect(totalMoment,cos(psi));
totalMomentConstantTerms = subs(subs(totalMoment,cos(psi),0),sin(psi),0);
totalMomentSineTerms = subs((totalMoment - totalMomentConstantTerms),cos(psi),0)/sin(psi);
totalMomentCosineTerms = subs((totalMoment - totalMomentConstantTerms),sin(psi),0)/cos(psi);


 [beta1cDD,beta1sDD] = solve(totalMomentSineTerms,totalMomentCosineTerms,beta1cDD,beta1sDD)
 beta0 = solve(totalMomentConstantTerms,beta0);
 
clc;
disp('%%%');
disp('beta0 = ...');
disp(beta0);
disp('%%%');
disp('beta1cDD = ...');
disp(beta1cDD);
disp('%%%');
disp('beta1sDD = ...');
disp(beta1sDD);
disp('%%%');


% 
% beta0 = simplify(vpa(eval(beta0)));
% beta1cDD = simplify(vpa(eval(beta1cDD)));
% beta1sDD = simplify(vpa(eval(beta1sDD)));
% 
% clc;
% disp('%%%');
% disp('beta0 = ...');
% disp(beta0);
% disp('%%%');
% disp('beta1cDD = ...');
% disp(beta1cDD);
% disp('%%%');
% disp('beta1sDD = ...');
% disp(beta1sDD);
% disp('%%%');
%  
 

 
 
 
 
 



