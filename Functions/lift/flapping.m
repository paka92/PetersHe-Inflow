function [outBeta0, outBeta1c, outBeta1s, outBeta1cD, outBetea1sD] = flapping(GetOrSet, theta0, theta1c, theta1s,RotorRadiusM,RotorSpeedRads, ChordM,ClAlpha,densityKgm3,meanInflowMs,...
    HingeOffset,BladeMassPerLength,u,v,w,dt)
%FLAPPING 1:Get 2:Set
persistent pbeta0;
if isempty(pbeta0)
   pbeta0 = 0.0; 
end
persistent pbeta1cD;
if isempty(pbeta1cD)
   pbeta1cD = 0.0; 
end
persistent pbeta1sD;
if isempty(pbeta1sD)
   pbeta1sD = 0.0; 
end
persistent pbeta1c;
if isempty(pbeta1c)
   pbeta1c = 0.0; 
end
persistent pbeta1s;
if isempty(pbeta1s)
   pbeta1s = 0.0; 
end
R       = RotorRadiusM;
c       = ChordM;
cla     = ClAlpha;
rho     = densityKgm3;
inflow  = meanInflowMs;
Ib      = RotorRadiusM * RotorRadiusM * BladeMassPerLength * 0.5;
e       = HingeOffset;
m       = BladeMassPerLength;
omega   = RotorSpeedRads;

beta1cD = pbeta1cD;
beta1sD = pbeta1sD;
beta1c  = pbeta1c;
beta1s  = pbeta1s;
beta0   = pbeta0;

if(2 == GetOrSet)
beta0 = ...
-((R^3*beta1sD*c*cla*rho*u)/12 + (R^3*beta1cD*c*cla*rho*v)/12 - (R^3*c*cla*omega*rho*w)/6 + (R^4*c*cla*inflow*omega^2*rho)/6 - (R^4*c*cla*omega^2*rho*theta0)/8 - (R^2*c*cla*rho*theta0*u^2)/8 - (R^2*c*cla*rho*theta0*v^2)/8 - (R^4*c*cla*e^3*inflow*omega^2*rho)/6 + (R^4*c*cla*e^4*omega^2*rho*theta0)/8 + (R^2*c*cla*e^2*rho*theta0*u^2)/8 + (R^2*c*cla*e^2*rho*theta0*v^2)/8 + (R^3*beta1s*c*cla*omega*rho*v)/12 - (R^3*c*cla*omega*rho*theta1s*u)/6 - (R^3*c*cla*omega*rho*theta1c*v)/6 - (R^3*beta1sD*c*cla*e^3*rho*u)/12 - (R^3*beta1cD*c*cla*e^3*rho*v)/12 + (R^3*c*cla*e^3*omega*rho*w)/6 - (R^3*beta1s*c*cla*e^3*omega*rho*v)/12 + (R^3*c*cla*e^3*omega*rho*theta1s*u)/6 + (R^3*c*cla*e^3*omega*rho*theta1c*v)/6)/((R^3*m*omega^2)/3 + (R^2*c*cla*rho*u*v)/8 - (R^2*c*cla*e^2*rho*u*v)/8);
%%%
beta1cDD = ...
-(96*Ib*beta1sD*omega - 48*Ib*beta1c*omega^2 + 16*R^3*beta1c*m*omega^2 + 6*R^4*beta1cD*c*cla*omega*rho - 12*R^2*c*cla*rho*v*w + 6*R^4*beta1s*c*cla*omega^2*rho - 6*R^4*c*cla*omega^2*rho*theta1c - 9*R^2*c*cla*rho*theta1c*v^2 - 6*R^4*beta1s*c*cla*e^4*omega^2*rho + 6*R^4*c*cla*e^4*omega^2*rho*theta1c + 9*R^2*c*cla*e^2*rho*theta1c*v^2 + 8*R^3*beta0*c*cla*omega*rho*u + 9*R^2*beta1c*c*cla*rho*u*v + 12*R^3*c*cla*inflow*omega*rho*v - 16*R^3*c*cla*omega*rho*theta0*v - 6*R^4*beta1cD*c*cla*e^4*omega*rho + 12*R^2*c*cla*e^2*rho*v*w - 8*R^3*beta0*c*cla*e^3*omega*rho*u - 9*R^2*beta1c*c*cla*e^2*rho*u*v - 12*R^3*c*cla*e^2*inflow*omega*rho*v + 16*R^3*c*cla*e^3*omega*rho*theta0*v)/(48*Ib);
 
%%%
beta1sDD = ...
(96*Ib*beta1cD*omega + 48*Ib*beta1s*omega^2 - 16*R^3*beta1s*m*omega^2 - 6*R^4*beta1sD*c*cla*omega*rho + 12*R^2*c*cla*rho*u*w + 6*R^4*beta1c*c*cla*omega^2*rho + 6*R^4*c*cla*omega^2*rho*theta1s + 9*R^2*c*cla*rho*theta1s*u^2 - 6*R^4*beta1c*c*cla*e^4*omega^2*rho - 6*R^4*c*cla*e^4*omega^2*rho*theta1s - 9*R^2*c*cla*e^2*rho*theta1s*u^2 - 12*R^3*c*cla*inflow*omega*rho*u + 16*R^3*c*cla*omega*rho*theta0*u + 6*R^4*beta1sD*c*cla*e^4*omega*rho - 12*R^2*c*cla*e^2*rho*u*w + 12*R^3*c*cla*e^2*inflow*omega*rho*u - 16*R^3*c*cla*e^3*omega*rho*theta0*u)/(48*Ib);

% if(beta0 > (20.0 * pi / 180.0))
%     beta0 = (20.0 * pi / 180.0);
% elseif(beta0 < (-20.0 * pi / 180.0))
%     beta0 = (-20.0 * pi / 180.0);
% end
% 
% if(beta1c > (20.0 * pi / 180.0))
%     beta1c = (20.0 * pi / 180.0);
% elseif(beta1c < (-20.0 * pi / 180.0))
%     beta1c = (-20.0 * pi / 180.0);
% end
% 
% if(beta1s > (20.0 * pi / 180.0))
%     beta1s = (20.0 * pi / 180.0);
% elseif(beta1s < (-20.0 * pi / 180.0))
%     beta1s = (-20.0 * pi / 180.0);
% end


%Integrate Betas
beta1cD = beta1cD   + dt * beta1cDD * 0.75;
beta1c  = beta1c    + dt * beta1cD * 0.75;
beta1sD = beta1sD   + dt * beta1sDD * 0.75;
beta1s  = beta1s    + dt * beta1sD * 0.75;
end
    

%
pbeta1cD = beta1cD;
pbeta1sD = beta1sD;
pbeta1c  = beta1c;
pbeta1s  = beta1s;
pbeta0   = beta0;

% Set outputs
outBeta0    = pbeta0;
outBeta1c   = pbeta1c;
outBeta1s   = pbeta1s;
outBeta1cD  = pbeta1cD;
outBetea1sD = pbeta1sD;

end

