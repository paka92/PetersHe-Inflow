function [CosFlowMatrix,SinFlowMatrix] = FlowMatrix(indexCosStates,indexSinStates,uMs,vMs,wMs, meanInflow, tipSpeedMs)
%FLOWMATRIX Summary of this function goes here
%   Detailed explanation goes here
RotorAngleOfAttack = atan2(-wMs,uMs);
AdvanceRatioRotor = sqrt(uMs*uMs + vMs*vMs + wMs*wMs) * cos(RotorAngleOfAttack) / tipSpeedMs; % mu
LambdaFree = sqrt(uMs*uMs + vMs*vMs + wMs*wMs) * sin(RotorAngleOfAttack) / tipSpeedMs; % mu

LambdaInduced = meanInflow;
LambdaTotal = LambdaFree + LambdaInduced;

VTotal = sqrt((AdvanceRatioRotor * AdvanceRatioRotor) + (LambdaTotal * LambdaTotal) );
if(abs(VTotal) < 0.001)
    Vmean = 0.0;
else
    Vmean = ((AdvanceRatioRotor * AdvanceRatioRotor) +  ( LambdaTotal + LambdaInduced) * LambdaTotal) / VTotal;
end
% Place into matricis
cosStateNumber = length(indexCosStates);
sinStateNumber = length(indexSinStates);

%
CosFlowMatrixInternal = eye(cosStateNumber) * Vmean;
CosFlowMatrixInternal(1,1) = VTotal;
CosFlowMatrix = CosFlowMatrixInternal;
%
SinFlowMatrixInternal = eye(sinStateNumber) * Vmean;
SinFlowMatrix = SinFlowMatrixInternal;

end