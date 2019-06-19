function[stateCount,cosStates,sinStates,Mcos,invLc,Msin,invLs] = CreateInflowMatrices(r,k,SkewAngleRad)
% Enter R and J to CreateMatrixNumeric function
% It will create the [superscript,subscript] values for each states in
% cosStates and sinStates cell array.

[cosStates,sinStates] = CreateMatrixNumeric(r,k);
%Get State Numbers
countCosStates = length(cosStates);
countSinStates = length(sinStates);
stateCount = countCosStates + countSinStates;
%Initialize Mass matricie
Mcos    = zeros(countCosStates);
Lcos    = zeros(countCosStates);
GamaCos = zeros(countCosStates);
Msin    = zeros(countSinStates);
Lsin    = zeros(countSinStates);
GamaSin = zeros(countSinStates);

for ii = 1:1:countCosStates
% State Loop    
    for jj = 1:1:countCosStates
    % Create interstate values
    Mcos(ii,jj) = MFunc(cosStates(ii,1),cosStates(ii,2),cosStates(jj,1),cosStates(jj,2));
    Lcos(ii,jj) = GainCSkewFunc(cosStates(ii,1),cosStates(ii,2),cosStates(jj,1),cosStates(jj,2),SkewAngleRad);
    GamaCos(ii,jj) = GamaFunc(cosStates(ii,1),cosStates(ii,2),cosStates(jj,1),cosStates(jj,2));
    end
    
end

for ii = 1:1:countSinStates
% State Loop    
    for jj = 1:1:countSinStates
    % Create interstate values
    Msin(ii,jj) = MFunc(sinStates(ii,1),sinStates(ii,2),sinStates(jj,1),sinStates(jj,2));
    Lsin(ii,jj) = GainSSkewFunc(sinStates(ii,1),sinStates(ii,2),sinStates(jj,1),sinStates(jj,2),SkewAngleRad);
    GamaSin(ii,jj) = GamaFunc(sinStates(ii,1),sinStates(ii,2),sinStates(jj,1),sinStates(jj,2));
    end
    
end
LcValue = Lcos .* GamaCos;
LsValue = Lsin .* GamaSin;
invLc = inv(LcValue);
invLs = inv(LsValue);

% for ii = 1:1:length(c)
%     disp(c{ii});
% end
% 
% for ii = 1:1:length(s)
%     disp(s{ii});
% end
end

