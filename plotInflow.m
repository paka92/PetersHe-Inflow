BladeNumber = 8;
ElementNumber = 10;
HingeOffset = 0.01;
TimeLoc =  200;
[stationMidPoints] = readElementNumber(HingeOffset,ElementNumber);
firstBladeAzimuth = 0.0;
angleBetweenBlades = 2*pi/BladeNumber;
lastBladeAzimuth = 2 * pi - angleBetweenBlades;
Azimuth = [firstBladeAzimuth:angleBetweenBlades:lastBladeAzimuth];

XLoc = zeros(BladeNumber + 1,ElementNumber);
YLoc = zeros(BladeNumber + 1,ElementNumber);
Zval = zeros(BladeNumber + 1,ElementNumber);
        
figure();
for blade=1:1:BladeNumber + 1
    for element = 1:1:10
        bladeValueCall = blade;
        if(blade > BladeNumber)
            bladeValueCall = blade - BladeNumber;
        end
        XLoc(blade,element) = stationMidPoints(element) * cos(Azimuth(bladeValueCall));
        YLoc(blade,element) = stationMidPoints(element) * sin(Azimuth(bladeValueCall));

        Zval(blade,element) = simout.signals.values(bladeValueCall, element, TimeLoc);
    end
end

surf(XLoc,YLoc,Zval);
