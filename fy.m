

?çeri?e atla
Gmail ürününü ekran okuyucularla birlikte kullanma
inflow 

28 ileti dizisinden 1.
pehehehe
Gelen Kutusu
x

ali karakaya
Ekler
13 Haz 2019 16:25 (2 gün önce)
Al?c?: ben


Ek alan?

function [Hub_Forces,Hub_Moments,ct_over_solid,inflow_rec,dL_vec,dD_vec,alpha_ang_ave,beta_all,psi,y_out,inf_out] = fcn(U_i,Ang_i,dt,omega_in,theta0,A1,B1,psi_prev,y_data_prev,inflow_data_prev)
%%% Define number of discrete Azimuth data
n = 1;
%%% Define time step of inside integration
h = dt/n;
%%% Define Rotor parameters
omega = omega_in*2*pi/60; %rad/s
ro = 0.002367; %slug/ft^3
a0 = 5.73; %per radian
cd0 = 0.01;
R = 2.80; %ft
chord = 2.6/12; %ft
VT = 624.3; %ft/s
inertia = 2.482159; %slug*ft^2
Lock = ro*a0*chord*R^4/inertia;
C1 = omega^2;
twist = -8*pi/180; %rad/(s/R) 
Force_calc = 0.5*chord*ro;
Psi_calc = omega_in*360/60*h; 
dAzimuth = omega_in*2*pi/60*dt;
hinge = 2/12; % ft 
ib_over_Ib = 1+3/2*hinge/(R-hinge);
delta3 =  0*pi/180; % rad
delta3_cont = -tan(delta3);
hub_spring = 0;
%%% Flapping Solution Vector
y_data = y_data_prev;
%%% Blade section vector
s = (0.70:0.02:2.80);
ds = 0.02; %ft
%%% Construct Flapping Moment Vector
dMy = zeros(1,size(s,2));
dRx_vec = zeros(1,size(s,2));
dRy_vec = zeros(1,size(s,2));
dRz_vec = zeros(1,size(s,2));
inf_ang_ave = zeros(1,size(s,2));
ut_ave = zeros(1,size(s,2));
up_ave = zeros(1,size(s,2));
%%% Define number of blades
blade_num = 16;
actual_blade_num = 4;
Forces = zeros(3,blade_num);
Moments = zeros(3,blade_num);
psi_curr_upd=zeros(blade_num,1);
%%% Auxiliary Variables
alpha_ang_ave = zeros(blade_num,size(s,2)-1);
dL_vec = zeros(blade_num,size(s,2)-1);
dD_vec = zeros(blade_num,size(s,2)-1);
dQ = zeros(blade_num,size(s,2)-1);
%%% Blade Angles initialize
Beta0_sum = 0;
Beta1c_sum = 0;
Beta1s_sum = 0;
%%% Peters-He Vectors
inf_data_next = zeros(15,1);
dlift_vec1 = zeros(1,size(s,2)-1);
dlift_vec2 = zeros(1,size(s,2)-1);
dlift_vec3 = zeros(1,size(s,2)-1);
dmom_vec_cos4 = zeros(1,size(s,2)-1);
dmom_vec_cos5 = zeros(1,size(s,2)-1);
dmom_vec_cos6 = zeros(1,size(s,2)-1);
dmom_vec_cos7 = zeros(1,size(s,2)-1);
dmom_vec_cos8 = zeros(1,size(s,2)-1);
dmom_vec_cos9 = zeros(1,size(s,2)-1);
dmom_vec_sin4 = zeros(1,size(s,2)-1);
dmom_vec_sin5 = zeros(1,size(s,2)-1);
dmom_vec_sin6 = zeros(1,size(s,2)-1);
dmom_vec_sin7 = zeros(1,size(s,2)-1);
dmom_vec_sin8 = zeros(1,size(s,2)-1);
dmom_vec_sin9 = zeros(1,size(s,2)-1);
inflow_rec = zeros(blade_num,(size(s,2)-1));
press11 = zeros(blade_num,1);
press22 = zeros(blade_num,1);
press33 = zeros(blade_num,1); 
press44c = zeros(blade_num,1); 
press55c = zeros(blade_num,1); 
press66c = zeros(blade_num,1);
press77c = zeros(blade_num,1);
press88c = zeros(blade_num,1);
press99c = zeros(blade_num,1);
press44s = zeros(blade_num,1);
press55s = zeros(blade_num,1);
press66s = zeros(blade_num,1);
press77s = zeros(blade_num,1);
press88s = zeros(blade_num,1);
press99s = zeros(blade_num,1);
%%% Define and Get Previous Peters-He Inflow
inf_data = inflow_data_prev; 
 %%%%%%%% Peters-HE constants %%%%%%%%%  
phi1 = sqrt(3); phi2 = sqrt(7); phi3 = sqrt(11); phi4 = sqrt(30)/2; phi5 = sqrt(5);
phi6 = sqrt(210)/4; phi7 = sqrt(1155/8); phi8 = 3/4*sqrt(35); phi9 = sqrt(3465/128);

% M Apperant mass matrices for cosine(2x2) and sine(1x1)
H_pe = [1 0 0 0 0 0 0 0 0; 0 4/9 0 0 0 0 0 0 0; 0 0 (8/15)^2 0 0 0 0 0 0;0 0 0 2/3 0 0 0 0 0;0 0 0 0 16/45 0 0 0 0;...
    0 0 0 0 0 8/15 0 0 0;0 0 0 0 0 0 32/105 0 0;0 0 0 0 0 0 0 48/105 0;0 0 0 0 0 0 0 0 128/315];
mc = 2/pi*H_pe;
mcinv = inv(mc);
ms = 2/pi*[ 2/3 0 0 0 0 0; 0 16/45 0 0 0 0;...
     0 0 8/15 0 0 0; 0 0 0 32/105 0 0; 0 0 0 0 48/105 0; 0 0 0 0 0 128/315];
msinv = inv(ms);
%%%% Inflow Parameters %%%%%%
lamdaind = phi1*inf_data(1);
UU = sqrt(U_i(1)^2 + U_i(2)^2 + U_i(3)^2);
alpha_disc = atan2(-U_i(3),U_i(1));%technical Notes
mu = UU*cos(alpha_disc)/VT;
lamdafree = UU*sin(alpha_disc)/VT;
lamda = lamdafree + lamdaind; 
x = pi/2-atan(abs(lamda/mu)); 
% x = atan(mu/lamda);
VTT = sqrt(mu.^2+lamda.^2); 
vv = (mu.^2 + (lamda+lamdaind)*lamda)/VTT; 
% xx = mu/(VTT + abs(lamda));
xx = tan(abs(x/2)); %check with xx whether they are same or not 
% Flow parameter matricies "V >>>> C"   eqn: 82-83
Ccos = [VTT 0 0 0 0 0 0 0 0;0 vv 0 0 0 0 0 0 0; 0 0 vv 0 0 0 0 0 0; 0 0 0 vv 0 0 0 0 0;...
       0 0 0 0 vv 0 0 0 0;0 0 0 0 0 vv 0 0 0;0 0 0 0 0 0 vv 0 0;0 0 0 0 0 0 0 vv 0;0 0 0 0 0 0 0 0 vv];
Csin = [vv 0 0 0 0 0;0 vv 0 0 0 0; 0 0 vv 0 0 0; 0 0 0 vv 0 0;0 0 0 0 vv 0;0 0 0 0 0 vv];
Taucos1 = [ 0.750000000000000     0.190940653956493    -0.029919597117386    -0.496729413289805      0.000000000000000...
            0.174304172194599    -0.028905076928149     0.000000000000000    -0.025032530918121];
        
Taucos2 = [ 0.190940653956493     0.656250000000000     0.205663227829503    -0.487778600638186     -0.497836949588487...
            0.599071547271275     0.198689554487648    -0.439050920690045     0.172070201652916];
        
Taucos3 = [-0.029919597117386     0.205663227829503     0.644531250000000     0.000000000000000     -0.496420625382602...
            0.187743981885905     0.622676344562828    -0.437801840148508     0.539253532727041];
        
Taucos4 = [ 0.496729413289805     0.487778600638186     0.000000000000000     0.625000000000000     0.191366386154936...
           -0.445278904396397     0.000000000000000     0.168769289021038     0.000000000000000];
        
Taucos5 = [ 0.000000000000000     0.497836949588487     0.496420625382602     0.191366386154936      0.632812500000000...
           -0.454460878748628    -0.479587887133218     0.558088167177687    -0.415335293604671];
       
Taucos6 = [ 0.174304172194599     0.599071547271275     0.187743981885905     0.445278904396397      0.454460878748628...
            0.546875000000000     0.181377918222561    -0.400796821925583     0.157077884866274];
        
Taucos7 = [-0.028905076928149     0.198689554487648     0.622676344562828     0.000000000000000      0.479587887133218...
            0.181377918222561     0.601562500000000    -0.422956760384469     0.520968406964076];
        
Taucos8 = [ 0.000000000000000     0.439050920690045     0.437801840148508     0.168769289021038      0.558088167177687...
            0.400796821925583     0.422956760384469     0.492187500000000    -0.366291299195318];
        
Taucos9 = [-0.025032530918121     0.172070201652916     0.539253532727041     0.000000000000000      0.415335293604671...
            0.157077884866274     0.520968406964076     0.366291299195318     0.451171875000000];
Tau = [Taucos1 ;Taucos2 ;Taucos3 ;Taucos4 ;Taucos5;Taucos6;Taucos7;Taucos8;Taucos9]; 
Lcos1 = [Tau(1,1),Tau(1,2),Tau(1,3),xx*Tau(1,4),xx*Tau(1,5),xx^2*Tau(1,6),xx^2*Tau(1,7),xx^3*Tau(1,8),xx^4*Tau(1,9)];
Lcos2 = [Tau(2,1),Tau(2,2),Tau(2,3),xx*Tau(2,4),xx*Tau(2,5),xx^2*Tau(2,6),xx^2*Tau(2,7),xx^3*Tau(2,8),xx^4*Tau(2,9)];
Lcos3 = [Tau(3,1),Tau(3,2),Tau(3,3),xx*Tau(3,4),xx*Tau(3,5),xx^2*Tau(3,6),xx^2*Tau(3,7),xx^3*Tau(3,8),xx^4*Tau(3,9)];
Lcos4 = [2*xx*Tau(4,1),2*xx*Tau(4,2),2*xx*Tau(4,3),(1-xx^2)*Tau(4,4),(1-xx^2)*Tau(4,5),xx*(1-xx^2)*Tau(4,6),xx*(1-xx^2)*Tau(4,7),xx^2*(1-xx^2)*Tau(4,8),xx^3*(1-xx^2)*Tau(4,9)];
Lcos5 = [2*xx*Tau(5,1),2*xx*Tau(5,2),2*xx*Tau(5,3),(1-xx^2)*Tau(5,4),(1-xx^2)*Tau(5,5),xx*(1-xx^2)*Tau(5,6),xx*(1-xx^2)*Tau(5,7),xx^2*(1-xx^2)*Tau(5,8),xx^3*(1-xx^2)*Tau(5,9)];
Lcos6 = [2*xx^2*Tau(6,1),2*xx^2*Tau(6,2),2*xx^2*Tau(6,3),xx*(1-xx^2)*Tau(6,4),xx*(1-xx^2)*Tau(6,5),(1+xx^4)*Tau(6,6),(1+xx^4)*Tau(6,7),xx*(1+xx^4)*Tau(6,8),xx^2*(1+xx^4)*Tau(6,9)];
Lcos7 = [2*xx^2*Tau(7,1),2*xx^2*Tau(7,2),2*xx^2*Tau(7,3),xx*(1-xx^2)*Tau(7,4),xx*(1-xx^2)*Tau(7,5),(1+xx^4)*Tau(7,6),(1+xx^4)*Tau(7,7),xx*(1+xx^4)*Tau(7,8),xx^2*(1+xx^4)*Tau(7,9)];
Lcos8 = [2*xx^3*Tau(8,1),2*xx^3*Tau(8,2),2*xx^3*Tau(8,3),xx^2*(1-xx^2)*Tau(8,4),xx^2*(1-xx^2)*Tau(8,5),xx*(1+xx^4)*Tau(8,6),xx*(1+xx^4)*Tau(8,7),(1-xx^6)*Tau(8,8),xx*(1-xx^6)*Tau(8,9)];
Lcos9 = [2*xx^4*Tau(9,1),2*xx^4*Tau(9,2),2*xx^4*Tau(9,3),xx^3*(1-xx^2)*Tau(9,4),xx^3*(1-xx^2)*Tau(9,5),xx^2*(1+xx^4)*Tau(9,6),xx^2*(1+xx^4)*Tau(9,7),xx*(1-xx^6)*Tau(9,8),(1+xx^8)*Tau(9,9)];
Lcos = [Lcos1;Lcos2;Lcos3;Lcos4;Lcos5;Lcos6;Lcos7;Lcos8;Lcos9];
Lsin4 = [(1+xx^2)*Tau(4,4),(1+xx^2)*Tau(4,5),xx*(1+xx^2)*Tau(4,6),xx*(1+xx^2)*Tau(4,7),xx^2*(1+xx^2)*Tau(4,8),xx^3*(1+xx^2)*Tau(4,9)];
Lsin5 = [(1+xx^2)*Tau(5,4),(1+xx^2)*Tau(5,5),xx*(1+xx^2)*Tau(5,6),xx*(1+xx^2)*Tau(5,7),xx^2*(1+xx^2)*Tau(5,8),xx^3*(1+xx^2)*Tau(5,9)];
Lsin6 = [xx*(1+xx^2)*Tau(6,4),xx*(1+xx^2)*Tau(6,5),(1-xx^4)*Tau(6,6),(1-xx^4)*Tau(6,7),xx*(1-xx^4)*Tau(6,8),xx^2*(1-xx^4)*Tau(6,9)];
Lsin7 = [xx*(1+xx^2)*Tau(7,4),xx*(1+xx^2)*Tau(7,5),(1-xx^4)*Tau(7,6),(1-xx^4)*Tau(7,7),xx*(1-xx^4)*Tau(7,8),xx^2*(1-xx^4)*Tau(7,9)];
Lsin8 = [xx^2*(1+xx^2)*Tau(8,4),xx^2*(1+xx^2)*Tau(8,5),xx*(1-xx^4)*Tau(8,6),xx*(1-xx^4)*Tau(8,7),(1+xx^6)*Tau(8,8),xx*(1+xx^6)*Tau(8,9)];
Lsin9 = [xx^3*(1+xx^2)*Tau(9,4),xx^3*(1+xx^2)*Tau(9,5),xx^2*(1-xx^4)*Tau(9,6),xx^2*(1-xx^4)*Tau(9,7),xx*(1+xx^6)*Tau(9,8),(1-xx^8)*Tau(9,9)];
Lsin = [Lsin4;Lsin5;Lsin6;Lsin7;Lsin8;Lsin9];
% L-operator ( inflow gain matrix) L^-1 = Ltilde^-1 * C >>>>> D
Dcos = inv(Lcos)*Ccos; 
Dsin = inv(Lsin)*Csin;
%%% Start Outer Loop
%%% Start Outer Loop
for k=1:blade_num
        psi_prev_blade = psi_prev(k);
        %%% Psi position of the blade
        psi_curr = (psi_prev_blade)*pi/180;
        %%% Howlett Method for Flapping
        beta_dot = y_data(3*k)*(sin(dAzimuth)/omega) + y_data(3*k-1)*cos(dAzimuth);
        beta = y_data(3*k-2) + y_data(3*k-1)*(sin(dAzimuth)/omega) + y_data(3*k)*((1-cos(dAzimuth))/omega^2); 
        y_data(3*k-1) = beta_dot;
        y_data(3*k-2) = beta;
        %%% Calculate intermediate varaibles
        COS = cos(psi_curr); COS2 = cos(2*psi_curr);
        COS3 = cos(3*psi_curr); COS4 = cos(4*psi_curr);
        SIN = sin(psi_curr); SIN2 = sin(2*psi_curr);
        SIN3 = sin(3*psi_curr); SIN4 = sin(4*psi_curr);
        ut_cont = U_i(1)*SIN + U_i(2)*COS;
        up_cont = y_data(3*k-2)*(U_i(1)*COS - U_i(2)*SIN) - U_i(3);
        up_ang =(-Ang_i(1)*SIN - Ang_i(2)*COS + y_data(3*k-1));
        pilot_cont = theta0 + A1*COS + B1*SIN + delta3_cont*y_data(3*k-2);
        %%% Calculate radial (Ur), tangential (Ut) and perpendicular (Up)
        ur = U_i(1)*COS - U_i(2)*SIN - y_data(3*k-2)*U_i(3) - 0;
        %%% Start Inner Loop
        for j=1:size(s,2)-1 %%% Trapezoidal Integral Loop to Calculate Lift and Drag
            rbar = s(j)/R;
            rbarnxt = s(j+1)/R;
            rmean = (rbar + rbarnxt)/2;
            inflow_PeHe = phi1*inf_data(1) + phi2*(1-5/2*rbar^2)*inf_data(2)+ phi3*(1-7*rbar^2+63/8*rbar^4)*inf_data(3)...
                          + phi4*rbar*(inf_data(4)*COS+inf_data(10)*SIN)+ phi5*(3*rbar-21/4*rbar^3)*(inf_data(5)*COS + inf_data(11)*SIN)...
                          + phi6*rbar^2*(inf_data(6)*COS2 + inf_data(12)*SIN2) + phi7*(rbar^2-1.5*rbar^4)*(inf_data(7)*COS2 + inf_data(13)*SIN2)...
                          + phi8*rbar^3*(inf_data(8)*COS3 + inf_data(14)*SIN3) + phi9*rbar^4*(inf_data(9)*COS4 + inf_data(15)*SIN4);
                          
                      
            inflow_PeHe_next = phi1*inf_data(1) + phi2*(1-5/2*rbarnxt^2)*inf_data(2)+ phi3*(1-7*rbarnxt^2+63/8*rbarnxt^4)*inf_data(3)...
                              + phi4*rbarnxt*(inf_data(4)*COS+inf_data(10)*SIN)+ phi5*(3*rbarnxt-21/4*rbarnxt^3)*(inf_data(5)*COS + inf_data(11)*SIN)...
                              + phi6*rbarnxt^2*(inf_data(6)*COS2 + inf_data(12)*SIN2) + phi7*(rbarnxt^2-1.5*rbarnxt^4)*(inf_data(7)*COS2 + inf_data(13)*SIN2)...
                              + phi8*rbarnxt^3*(inf_data(8)*COS3 + inf_data(14)*SIN3) + phi9*rbarnxt^4*(inf_data(9)*COS4 + inf_data(15)*SIN4);
                               
            inflow_rec(k,j) = (inflow_PeHe + inflow_PeHe_next)/2*omega*R;
            ut = ut_cont + s(j)*omega -0;
                ut_next = ut_cont + s(j+1)*omega -0;
            up = up_cont + s(j)*up_ang + inflow_PeHe*omega*R;
                up_next = up_cont + s(j+1)*up_ang + inflow_PeHe_next*omega*R;
            inf_ang = atan2(-up,ut);
                inf_ang_next = atan2(-up_next,ut_next);
            alpha = pilot_cont + inf_ang + twist*s(j)/R;
                alpha_next = pilot_cont + inf_ang_next + twist*s(j+1)/R;
            TAS_sqr = ut^2+up^2;%+ur^2;
                TAS_sqr_next = ut_next^2+up_next^2;%+ur^2;
                %% Start of Section Lift and Drag
                    %%% Calculate Section Drag
                       lim1 = 0.212647; 
                       alpha_d = alpha; 
                       alpha_next_d = alpha_next;
                       alpha = abs(alpha); alpha_next = abs(alpha_next);
                    if alpha <= lim1 
                        cd = 0.21359*alpha - 0.002328 + 0.002329;
                    else 
                        cd = 0.2859*alpha^4-1.773*alpha^3 + 2.528*alpha^2 + 0.7592*alpha -0.2162 + 0.002329;
                    end
                    
                    if  alpha_next <= lim1
                        cd_next = 0.21359*alpha_next - 0.002328 + 0.002329;
                    else  
                        cd_next = 0.2859*alpha_next^4-1.773*alpha_next^3 + 2.528*alpha_next^2 + 0.7592*alpha_next -0.2162 + 0.002329;  
                    end
                    %%% End of Section Drag
                    %% Calculate Section Lift
                    aoa = abs(alpha); aoa_next = abs(alpha_next);
                    if aoa >=0 && aoa <= 0.16
                       cl= -1660*(aoa^4)+ 473.2*(aoa^3)- 43.68*(aoa^2)+ 7.738*(aoa)-0.005065;
                    elseif aoa >0.16 && aoa <= 0.24
                       cl=-64080*(aoa^4)+ 50010*(aoa^3)-14640*(aoa^2)+ 1908*(aoa)-92.27; 
                    elseif aoa >0.24 && aoa <= 0.47
                       cl=1074*(aoa^4)-1685*(aoa^3)+988.4 *(aoa^2)-254.1 *(aoa)+24.68;
                    elseif aoa >0.47 && aoa <= 0.61
                       cl=22.53 *(aoa^2)-24.3 *(aoa)+7.456;
                    elseif aoa >0.61 && aoa <= 1.92
                        cl=-0.4341*(aoa^4)+3.026*(aoa^3)-7.863*(aoa^2)+7.341 *(aoa)-1.14;
                    elseif aoa >1.92 && aoa <= 2.71
                       cl=1.964*(aoa^2)-9.426*(aoa)+10.37;
                    elseif aoa >2.71 && aoa <= 2.97
                       cl=-8.478*(aoa^2)+47.75*(aoa)-67.88;
                    elseif aoa >2.97 && aoa <= 3.14
                       cl=33.33*(aoa^2)-198.7*(aoa)+295.3;
                    else
                        cl = 0;
                    end   
                    %%%%%%%%%%%
                    if aoa_next >=0 && aoa_next <= 0.16
                       cl_next= -1660*(aoa_next^4)+ 473.2*(aoa_next^3)- 43.68*(aoa_next^2)+ 7.738*(aoa_next)-0.005065;
                    elseif aoa_next >0.16 && aoa_next <= 0.24
                       cl_next=-64080*(aoa_next^4)+ 50010*(aoa_next^3)-14640*(aoa_next^2)+ 1908*(aoa_next)-92.27; 
                    elseif aoa_next >0.24 && aoa_next <= 0.47
                       cl_next=1074*(aoa_next^4)-1685*(aoa_next^3)+988.4 *(aoa_next^2)-254.1 *(aoa_next)+24.68;
                    elseif aoa_next >0.47 && aoa_next <= 0.61
                       cl_next=22.53 *(aoa_next^2)-24.3 *(aoa_next)+7.456;
                    elseif aoa_next >0.61 && aoa_next <= 1.92
                        cl_next=-0.4341*(aoa_next^4)+3.026*(aoa_next^3)-7.863*(aoa_next^2)+7.341 *(aoa_next)-1.14;
                    elseif aoa_next >1.92 && aoa_next <= 2.71
                       cl_next=1.964*(aoa_next^2)-9.426*(aoa_next)+10.37;
                    elseif aoa_next >2.71 && aoa_next <= 2.97
                       cl_next=-8.478*(aoa_next^2)+47.75*(aoa_next)-67.88;
                    elseif aoa_next >2.97 && aoa_next <= 3.14
                       cl_next=33.33*(aoa_next^2)-198.7*(aoa_next)+295.3;
                    else
                        cl_next = 0;
                    end
                    %%% End of Section Lift
                    if alpha_d < 0
                        cl = -cl;
                    end
                    
                    if alpha_next_d < 0
                        cl_next = -cl_next;
                    end
                %%% End of Section Lift and Drag
                    alpha = alpha_d; 
                    alpha_next = alpha_next_d;
            %% Calculate Differential Lift, Drag and Moment
            dL = Force_calc*TAS_sqr*cl; %%% Use actual data
                dL_next = Force_calc*TAS_sqr_next*cl_next; %%% Use actual data
            dD = Force_calc*TAS_sqr*cd; %%% Use actual data
                dD_next = Force_calc*TAS_sqr_next*cd_next; %%% Use actual data
            
            dL_trape = ds*(dL+dL_next)/2; 
            dD_trape = ds*(dD+dD_next)/2; 
            
            dL_vec(k,j)= dL_trape;
            dD_vec(k,j)= dD_trape;
            
            ut_ave(j) = (ut+ut_next)/2;
            up_ave(j) = (up+up_next)/2;
            inf_ang_ave(j) = (inf_ang+inf_ang_next)/2;
            alpha_ang_ave(k,j) = (alpha+alpha_next)/2*180/pi;

            dRx_vec(j) = y_data(3*k-2)*dL_trape;
            dRy_vec(j) = -dD_trape + inf_ang_ave(j)*dL_trape;
            dRz_vec(j) = -inf_ang_ave(j)*dD_trape -dL_trape;
            dlift_vec1(j) = dL_trape*phi1;
            dlift_vec2(j) = dL_trape*phi2*(1-5/2*rmean^2);
            dlift_vec3(j) = dL_trape*phi3*(1-7*rmean^2+63/8*rmean^4);
            dmom_vec_cos4(j) = dL_trape*phi4*rmean;
            dmom_vec_cos5(j) = dL_trape*phi5*(3*rmean-21/4*rmean^3);
            dmom_vec_cos6(j) = dL_trape*phi6*rmean^2;
            dmom_vec_cos7(j) = dL_trape*phi7*(rmean^2-1.5*rmean^4);
            dmom_vec_cos8(j) = dL_trape*phi8*rmean^3;
            dmom_vec_cos9(j) = dL_trape*phi9*rmean^4;
            dmom_vec_sin4(j) = dL_trape*phi4*rmean;
            dmom_vec_sin5(j) = dL_trape*phi5*(3*rmean-21/4*rmean^3);
            dmom_vec_sin6(j) = dL_trape*phi6*rmean^2;
            dmom_vec_sin7(j) = dL_trape*phi7*(rmean^2-1.5*rmean^4);
            dmom_vec_sin8(j) = dL_trape*phi8*rmean^3;
            dmom_vec_sin9(j) = dL_trape*phi9*rmean^4;
            
            %%% Calculate Flapping Moment
            dMy(j) = (s(j)+ds/2)*(-dL_trape-inf_ang_ave(j)*dD_trape)*1; 
            %%% Calculate Differential Torque
            dQ(k,j) = (s(j)+ds/2)*cos(y_data(3*k-2))*(-dRy_vec(j));
        end
            %%% Output H,Y,T Forces at the end of each simulation time
                Rx = sum(dRx_vec);
                Ry = sum(dRy_vec);
                Rz = sum(dRz_vec);
                H_Force = COS*Rx + SIN*Ry;
                Y_Force = -SIN*Rx + COS*Ry;
                T_Force = Rz;
                Forces(:,k) = [H_Force;Y_Force;T_Force];
            %%% Pressure for Pe-HE %%%
            press11(k) = sum(dlift_vec1); 
            press22(k) = sum(dlift_vec2);
            press33(k) = sum(dlift_vec3); 
            press44c(k) = sum(dmom_vec_cos4)*COS;
            press55c(k) = sum(dmom_vec_cos5)*COS; 
            press66c(k) = sum(dmom_vec_cos6)*COS2; 
            press77c(k) = sum(dmom_vec_cos7)*COS2;
            press88c(k) = sum(dmom_vec_cos8)*COS3; 
            press99c(k) = sum(dmom_vec_cos9)*COS4; 
            press44s(k) = sum(dmom_vec_sin4)*SIN;
            press55s(k) = sum(dmom_vec_sin5)*SIN;
            press66s(k) = sum(dmom_vec_sin6)*SIN2;
            press77s(k) = sum(dmom_vec_sin7)*SIN2;
            press88s(k) = sum(dmom_vec_sin8)*SIN3;
            press99s(k) = sum(dmom_vec_sin9)*SIN4;
            %%% Calculate Additional Moments due to Hinge Offset
                M_Rx = 0;
                M_Ry = hinge*Rz - hub_spring*y_data(3*k-2);
                M_Rz = -hinge*Ry;
                X_Moment = COS*M_Rx + SIN*M_Ry;
                Y_Moment = -SIN*M_Rx + COS*M_Ry;
                Z_Moment = M_Rz;
                Moments(:,k) = [X_Moment;Y_Moment;Z_Moment];    
            %%% Calculate Total Flapping Moment
            M_tot=sum(dMy);
            %%% Howlett Method Continued            
            dot= [0 1;-C1*ib_over_Ib-(hub_spring/inertia) 0]*[y_data(3*k-2);y_data(3*k-1)]+[0;1/inertia]*(-M_tot)+...
                 [0;-2*omega*ib_over_Ib]*(Ang_i(2)*SIN-Ang_i(1)*COS);
            y_data(3*k) = dot(2);
            %%% Actual Psi
            psi_curr_upd(k) = Psi_calc*(1); %degree
            %%% Calculate Tip Path Plane Angles
            Beta0_sum = Beta0_sum + y_data(3*k-2);
            Beta1c_sum = Beta1c_sum + y_data(3*k-2)*cos(psi_curr);
            Beta1s_sum = Beta1s_sum + y_data(3*k-2)*sin(psi_curr);
end
y_out = y_data(:,end); %[beta;beta_dot;beta_ddot];
psi = psi_prev + psi_curr_upd;

Hub_Forces_sum = zeros(3,1);
Hub_Moments_sum = zeros(3,1);
for i=1:size(psi_prev,1)
psi(i,1) = mod(psi(i),360);  
Hub_Forces_sum = Hub_Forces_sum + Forces(:,i);
Hub_Moments_sum = Hub_Moments_sum + Moments(:,i);
end
%%% Calculate Hub Forces
Hub_Forces = (Hub_Forces_sum)/blade_num*actual_blade_num;
%%% Calculate Hub Moments
Hub_Moments = (Hub_Moments_sum)/blade_num*actual_blade_num;        
Hub_Moments(3)= Hub_Moments(3) + sum(sum(dQ))/blade_num*actual_blade_num;
%%% Calculate Total Torque for check
Torque = Hub_Moments(3);
lft=0; yforce = 0; xforce = 0;
for b=1:1:blade_num
    lft = lft + Forces(3,b);
    yforce = yforce + Forces(2,b);
    xforce = xforce + Forces(1,b);
end
        
%%%%Dynamic Inflow %%%%%
%%%% Take Average Pressure %%%%%
Non_dim = pi*ro*VT^2*R^2;
press1 = sum(press11)/blade_num*actual_blade_num/(2*Non_dim);
press2 = sum(press22)/blade_num*actual_blade_num/(2*Non_dim);
press3 = sum(press33)/blade_num*actual_blade_num/(2*Non_dim);
press4c = sum(press44c)/blade_num*actual_blade_num/Non_dim;
press5c = sum(press55c)/blade_num*actual_blade_num/Non_dim;
press6c = sum(press66c)/blade_num*actual_blade_num/Non_dim;
press7c = sum(press77c)/blade_num*actual_blade_num/Non_dim;
press8c = sum(press88c)/blade_num*actual_blade_num/Non_dim;
press9c = sum(press99c)/blade_num*actual_blade_num/Non_dim;
press4s = sum(press44s)/blade_num*actual_blade_num/Non_dim;
press5s = sum(press55s)/blade_num*actual_blade_num/Non_dim;
press6s = sum(press66s)/blade_num*actual_blade_num/Non_dim;
press7s = sum(press77s)/blade_num*actual_blade_num/Non_dim;
press8s = sum(press88s)/blade_num*actual_blade_num/Non_dim;
press9s = sum(press99s)/blade_num*actual_blade_num/Non_dim;
press = [press1;press2;press3;press4c;press5c;press6c;press7c;press8c;press9c];
press_sin = [press4s;press5s;press6s;press7s;press8s;press9s];
%%%% Euler Integration for PETERS-HE %%%%%%
data_cos = [inf_data(1);inf_data(2);inf_data(3);inf_data(4);...
          inf_data(5);inf_data(6);inf_data(7);inf_data(8);inf_data(9)];
data_sin = [inf_data(10);inf_data(11);inf_data(12);inf_data(13);inf_data(14);inf_data(15)];
dot_inf_mc = -mcinv*Dcos*data_cos + 0.5*mcinv*press;
dot_inf_s = -msinv*Dsin*data_sin + 0.5*msinv*press_sin;
%%%% Next Inflow %%%%%
inf_data_next(1) = inf_data(1) + 0.1*dt*omega*dot_inf_mc(1);
inf_data_next(2) = inf_data(2) + 0.1*dt*omega*dot_inf_mc(2);
inf_data_next(3) = inf_data(3) + 0.1*dt*omega*dot_inf_mc(3);
inf_data_next(4) = inf_data(4) + 0.1*dt*omega*dot_inf_mc(4);
inf_data_next(5) = inf_data(5) + 0.1*dt*omega*dot_inf_mc(5);
inf_data_next(6) = inf_data(6) + 0.1*dt*omega*dot_inf_mc(6);
inf_data_next(7) = inf_data(7) + 0.1*dt*omega*dot_inf_mc(7);  
inf_data_next(8) = inf_data(8) + 0.1*dt*omega*dot_inf_mc(8);  
inf_data_next(9) = inf_data(9) + 0.1*dt*omega*dot_inf_mc(9);  
inf_data_next(10) = inf_data(10) + 0.1*dt*omega*dot_inf_s(1);
inf_data_next(11) = inf_data(11) + 0.1*dt*omega*dot_inf_s(2);
inf_data_next(12) = inf_data(12) + 0.1*dt*omega*dot_inf_s(3); 
inf_data_next(13) = inf_data(13) + 0.1*dt*omega*dot_inf_s(4);
inf_data_next(14) = inf_data(14) + 0.1*dt*omega*dot_inf_s(5); 
inf_data_next(15) = inf_data(15) + 0.1*dt*omega*dot_inf_s(6); 
inf_out = inf_data_next(1:15);
%%%Peters_He Integration END %%%%
%%%%% Dynamic Inflow End %%%%%    
ct_over_solid = norm(Hub_Forces)/(ro*pi*R^2*VT^2);  
%%% Output Tip Path Plane Angles
beta0 = Beta0_sum*180/pi/blade_num;
beta1c = Beta1c_sum*180/pi/blade_num;
beta1s = Beta1s_sum*180/pi/blade_num;
beta_all = [beta0 -2*beta1c -2*beta1s];
%%% End of Output Tip Path Plane Angles
        
        
% %%% Analytical Steady-state Solution %%%
% e0 = 0.01; % Equivalent Drag Coeff.
% P = 1; % flap_freq/rotation_freq
% F_omega = 1; % aerodynamic_rotor_speed/inertial_rotor_speed
% BT = 1; % Tiploss Factor
% mu = U_i(1)/(omega*R); % Advance ratio, U_i should be U_a is there exist wind
% F0 = F_omega^2*(BT^4+BT^2*mu^2)/4;
% FT = F_omega^2*(BT^5/5+BT^3*mu^2/6);
% FB1 = F_omega^2*(BT^3*mu/3);
% F_lambda = F_omega^2*(BT^3/3)*(1+e0);
% A0 = F_omega^2*(2*BT^3*mu/3);
% Aa1 = F_omega^2*(BT^4/(4*F_omega)-BT^2*mu^2/8)*(1+e0);
% AT = F_omega^2*(2*BT^4*mu/2);
% AB1 = F_omega^2*(BT^4/(4)+3*BT^2*mu^2/8);
% A_lambda = F_omega^2*(BT^2*mu/2-mu^3/8)*(1+e0);
% Bbeta0 = F_omega^2*(BT^3*mu/3)*(1+e0);
% Bb1 = F_omega^2*(BT^4/(4*F_omega)+BT^2*mu^2/8)*(1+e0);
% BA1 = F_omega^2*(BT^4/(4)+BT^2*mu^2/8);
% 
% analytic_b0 = Lock/2*F0*theta0 +Lock/2*FT*twist +Lock/2*FB1*B1 +Lock/2*F_lambda*(inflow);
% analytic_a1 = A0/Aa1*theta0 + AT/Aa1*twist + AB1/Aa1*B1 + A_lambda/Aa1*(inflow);
% analytic_b1 = Lock/2*Bbeta0/Bb1*F0*theta0 + Lock/2*Bbeta0/Bb1*FT*twist -BA1/Bb1*A1 + Lock/2*Bbeta0/Bb1*FB1*B1 + Lock/2*Bbeta0/Bb1*F_lambda*(inflow);
% beta_analytic =  analytic_b0-cos(psi*pi/180)*analytic_a1-analytic_b1*sin(psi*pi/180);
% beta_dot_analytic = (analytic_a1*sin(psi*pi/180)-analytic_b1*cos(psi*pi/180))*omega;
ali.txt
ali.txt görüntüleniyor.