function [aRampEnd,aSquareEnd] = simpleCA3Integrator_V2(pStruct,inDur,rampPerc)
% LKW 7/29/21
% Function version of Attractor Model with Adaptation option
% Inputs: 
%   pStruct = struct of parameters from simpleSpaceSearch.m
%   inDur = double; input duration value
%   rampPerc = double; ramp percentage value (0-1)
% Outputs:
%   aRampEnd = final state of ramp afferent pyramidal
%   aSquareEnd = final state of square afferent pyramidal 

% Initialize variables based on pStruct input
T   = 5000;        %Time steps, must be even
W = pStruct.W; H = pStruct.H; Wh = pStruct.Wh; Hh = pStruct.Hh;
tha = pStruct.tha; thh = pStruct.thh; eta = pStruct.eta; Iexcit = pStruct.Iexcit;
% Ionic Currents and related parameters 
thc = pStruct.thc; mu = pStruct.mu; gm = pStruct.gm; om = pStruct.om;
Ek  = -10; Ena = 70; Ecl = 0;	%Reversal potentials of Potassium, Sodium, Chlorine

%Ramp vectors
ARamp = zeros(1,T);   %Initialize input vector
aRamp = zeros(1,T);    %initialize a activity vector
hRamp = zeros(1,T);   %Initialize h activity vector
cRamp = zeros(1,T);   % c activity vector

%Square vectors
ASquare = zeros(1,T);
aSquare = zeros(1,T);
hSquare = zeros(1,T);
cSquare = zeros(1,T);

% %Initialize inputs for afferent input A
ASquare(1:inDur)= Iexcit;    %Square pulse
rampLen         = round(inDur*rampPerc);    %For basic 20% ramp do 0.2

% For holding max current constant
if pStruct.wf_flag == 1
    ARamp(1:inDur)  = Iexcit;
    ARamp(1:rampLen)= linspace(0,Iexcit,rampLen);   %Front Ramp
elseif pStruct.wf_flag == 2
    ARamp(1:inDur)  = Iexcit;
    ARamp(1:rampLen)= linspace(0,Iexcit,rampLen);   %Front Ramp
    ARamp(inDur-rampLen+1:inDur) = linspace(Iexcit,0,rampLen);   %Rear ramp
elseif pStruct.wf_flag == 5
    ARamp(1:inDur)  = Iexcit;
    ARamp(inDur-rampLen+1:inDur) = linspace(Iexcit,0,rampLen);   %Rear ramp
end

% For holding Area Under Curve Constant:
squarea = Iexcit*inDur;     %Area under the square pulse
if pStruct.wf_flag == 3
    IRamp = squarea/(rampLen/2 + inDur - rampLen);  %Calculate Single Ramp Current Max for constant AUC
    ARamp(1:inDur)  = IRamp;
    ARamp(1:rampLen)= linspace(0,IRamp,rampLen);   %Front Ramp
elseif pStruct.wf_flag == 4
    IRamp = squarea/(inDur - rampLen);              %Calculate Double Ramp Current Max for constant AUC
    ARamp(1:inDur)  = IRamp;
    ARamp(1:rampLen)= linspace(0,IRamp,rampLen);   %Front Ramp
    ARamp(inDur-rampLen+1:inDur) = linspace(IRamp,0,rampLen);   %Rear ramp
elseif pStruct.wf_flag == 7
    IRamp = squarea/(rampLen/2 + inDur - rampLen);  %Calculate Single Ramp Current Max for constant AUC
    ARamp(1:inDur)  = IRamp;
    ARamp(inDur-rampLen+1:inDur) = linspace(IRamp,0,rampLen);   %Rear ramp
end

for t=1:T
    if (aRamp(t) - tha)   > 0; aaRamp   = aRamp(t) - tha;   else aaRamp   = 0; end %Threshold activity over 'tha' to 0
    if (aSquare(t) - tha) > 0; aaSquare = aSquare(t) - tha; else aaSquare = 0; end
    if (hRamp(t)-thh)     > 0; hhRamp   = hRamp(t)-thh;     else hhRamp   = 0; end %Threshold activity over 'thh' to 0
    if (hSquare(t)-thh)   > 0; hhSquare = hSquare(t)-thh;   else hhSquare = 0; end
    if (aRamp(t) - thc)   > 0; ccRamp   = aRamp(t)   - thc; else ccRamp   = 0; end %Threshold activity under 'thc' to 0
    if (aSquare(t) - thc) > 0; ccSquare = aSquare(t) - thc; else ccSquare = 0; end

    if pStruct.simTypeFlag == 1
        %Standard Linear variant
        daRamp = ARamp(t) + W*aaRamp - H*hhRamp - eta*aRamp(t); %change in a = Input A(t) + Weight*thesholded activity (a) - Weight*thesholded activity (h) - decay*activity (a)
        dhRamp = Wh*aaRamp - Hh*hhRamp - eta*hRamp(t);          %change in h = Weight*thresholded activity (a) - decay*activity (h)
        daSquare = ASquare(t) + W*aaSquare - H*hhSquare - eta*aSquare(t);
        dhSquare = Wh*aaSquare - Hh*hhSquare - eta*hSquare(t);
    elseif pStruct.simTypeFlag == 2
        %Linear with Adaptation variant
        daRamp = ARamp(t) + W*aaRamp - H*hhRamp - eta*aRamp(t) + mu.*cRamp(t).*(Ek - aRamp(t));
        dcRamp = gm.*ccRamp - om.*cRamp(t);
        dhRamp = Wh*aaRamp - Hh*hhRamp - eta*hRamp(t);
        daSquare = ASquare(t) + W*aaSquare - H*hhSquare - eta*aSquare(t) + mu.*cSquare(t).*(Ek - aSquare(t));
        dcSquare = gm.*ccSquare - om.*cSquare(t);
        dhSquare = Wh*aaSquare - Hh*hhSquare - eta*hSquare(t);
    elseif pStruct.simTypeFlag == 3
        %Ionic potential variant
        daRamp = ARamp(t) + (Ena - aRamp(t)).*W*aaRamp + (Ecl - aRamp(t)).*H*hhRamp - eta*aRamp(t); %change in a = Input A(t) + Weight*thesholded activity (a) - Weight*thesholded activity (h) - decay*activity (a)
        dhRamp = (Ena - hRamp(t)).*Wh*aaRamp + (Ecl - hRamp(t)).*Hh*hhRamp - eta*hRamp(t);          %change in h = Weight*thresholded activity (a) - decay*activity (h)
        daSquare = ASquare(t) + (Ena - aSquare(t)).*W*aaSquare + (Ecl - aSquare(t)).*H*hhSquare - eta*aSquare(t);
        dhSquare = (Ena - hSquare(t)).*Wh*aaSquare + (Ecl - hSquare(t)).*Hh*hhSquare - eta*hSquare(t);
    elseif pStruct.simTypeFlag == 4
        %Ionic potential variant with Adaptation
        daRamp = ARamp(t) + (Ena - aRamp(t)).*W*aaRamp + (Ecl - aRamp(t)).*H*hhRamp - eta*aRamp(t) + mu.*cRamp(t).*(Ek - aRamp(t)); %change in a = Input A(t) + Weight*thesholded activity (a) - Weight*thesholded activity (h) - decay*activity (a)
        dhRamp = (Ena - hRamp(t)).*Wh*aaRamp + (Ecl - hRamp(t)).*Hh*hhRamp - eta*hRamp(t);          %change in h = Weight*thresholded activity (a) - decay*activity (h)
        dcRamp = gm.*ccRamp - om.*cRamp(t);
        daSquare = ASquare(t) + (Ena - aSquare(t)).*W*aaSquare + (Ecl - aSquare(t)).*H*hhSquare - eta*aSquare(t) + mu.*cSquare(t).*(Ek - aSquare(t));
        dhSquare = (Ena - hSquare(t)).*Wh*aaSquare + (Ecl - hSquare(t)).*Hh*hhSquare - eta*hSquare(t);
        dcSquare = gm.*ccSquare - om.*cSquare(t);
    end

    %Neuronal activity update
    aRamp(t+1) = aRamp(t)+daRamp; %Update next time point of a
    hRamp(t+1) = hRamp(t)+dhRamp; %Update next h time point
    aSquare(t+1) = aSquare(t)+daSquare;
    hSquare(t+1) = hSquare(t)+dhSquare;
    
    if pStruct.simTypeFlag == 2 || pStruct.simTypeFlag == 4
        %Adaptation update
        cRamp(t+1) = cRamp(t) + dcRamp;
        cSquare(t+1) = cSquare(t) + dcSquare;
    end

end

aRampEnd = aRamp(end);
aSquareEnd = aSquare(end);

end