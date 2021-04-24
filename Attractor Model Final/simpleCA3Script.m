%%% Sample code for Attractor Model
%LKW 4/24/21 also see Hasselmo et al., J. Neurosci. 1995.
%Plots neural activation "a" over time where a is an autorecurrent CA3 neuron 
%Plots interneuron activation "h" over time.

rampTypeFlag = 1;         %1 = FR IMA; 2 = DR IP; 3 = FR IMA; 4 = DR IP; 5 = BR IMA; 6 = BR IP

W       = 0.016;        %Weight of a to a (recurrency)
H       = 0.06;         %Weight of h to a (feedback inhibition)
Wh      = 0.0042;       %Weight of a to h (feedforward activation of IN)
Hh      = 0.003;        %Inhibition self feedback, try 0.003
tha     = 8;            %Threshold for a
thh     = 8;            %Threshold for h
eta     = 0.01;         %Decay constant
T       = 5000;         %Time steps, must be even
Iexcit  = 0.1;          %Stim strength
inDur   = 500;          %Stim duration
rampPerc= 0.5;          %Ramp percentage 0 (square) to 1

%Ramp vectors
ARamp   = zeros(1,T);   %Initialize input vector
aRamp   = zeros(1,T);   %initialize a activity vector
hRamp   = zeros(1,T);   %Initialize h activity vector
%Square vectors
ASquare = zeros(1,T);
aSquare = zeros(1,T);
hSquare = zeros(1,T);

%Initialize inputs
ASquare(1:inDur)= Iexcit;    %Square pulse
rampLen         = round(inDur*rampPerc);

%For equal max current
if rampTypeFlag == 1
    ARamp(1:inDur)  = Iexcit;
    ARamp(1:rampLen)= linspace(0,Iexcit,rampLen);   %Front Ramp
elseif rampTypeFlag == 2
    ARamp(1:inDur)  = Iexcit;
    ARamp(1:rampLen)= linspace(0,Iexcit,rampLen);   %Front Ramp
    ARamp(inDur-rampLen+1:inDur) = linspace(Iexcit,0,rampLen);   %Rear ramp
elseif rampTypeFlag == 3
    squarea = Iexcit*inDur;     %Area under the square pulse
    IRamp = squarea/(rampLen/2 + inDur - rampLen);  %Calculate Single Ramp Current Max for constant AUC
    ARamp(1:inDur)  = IRamp;
    ARamp(1:rampLen)= linspace(0,IRamp,rampLen);   %Front Ramp
elseif rampTypeFlag == 4
    squarea = Iexcit*inDur;     %Area under the square pulse
    IRamp = squarea/(inDur - rampLen);              %Calculate Double Ramp Current Max for constant AUC
    ARamp(1:inDur)  = IRamp;
    ARamp(1:rampLen)= linspace(0,IRamp,rampLen);   %Front Ramp
    ARamp(inDur-rampLen+1:inDur) = linspace(IRamp,0,rampLen);   %Rear ramp
elseif rampTypeFlag == 5
    ARamp(1:inDur)  = Iexcit;
    ARamp(inDur-rampLen+1:inDur) = linspace(Iexcit,0,rampLen);   %Rear ramp
elseif rampTypeFlag == 6
    squarea = Iexcit*inDur;     %Area under the square pulse
    IRamp = squarea/(rampLen/2 + inDur - rampLen);  %Calculate Single Ramp Current Max for constant AUC
    ARamp(1:inDur)  = IRamp;
    ARamp(inDur-rampLen+1:inDur) = linspace(IRamp,0,rampLen);   %Rear ramp
end

for t=1:T
    if (aRamp(t) - tha)   > 0; aaRamp   = aRamp(t) - tha;   else aaRamp = 0;   end %Threshold activity over 'tha' to 0
    if (aSquare(t) - tha) > 0; aaSquare = aSquare(t) - tha; else aaSquare = 0; end
    if (hRamp(t)-thh)     > 0; hhRamp   = hRamp(t)-thh;     else hhRamp = 0;   end %Threshold activity over 'thh' to 0
    if (hSquare(t)-thh)   > 0; hhSquare = hSquare(t)-thh;   else hhSquare = 0; end
    
    daRamp = ARamp(t) + W*aaRamp - H*hhRamp - eta*aRamp(t); %change in a = Input A(t) + Weight*thesholded activity (a) - Weight*thesholded activity (h) - decay*activity (a)
    dhRamp = Wh*aaRamp - Hh*hhRamp - eta*hRamp(t);          %change in h = Weight*thresholded activity (a) - decay*activity (h)

    daSquare = ASquare(t) + W*aaSquare - H*hhSquare - eta*aSquare(t);
    dhSquare = Wh*aaSquare - Hh*hhSquare - eta*hSquare(t);
    
    aRamp(t+1) = aRamp(t)+daRamp; %Update next time point of a
    hRamp(t+1) = hRamp(t)+dhRamp; %Update next h time point
    aSquare(t+1) = aSquare(t)+daSquare;
    hSquare(t+1) = hSquare(t)+dhSquare;
end

%Observe final activity state of pyr neurons
aRampEnd = aRamp(end);
aSquareEnd = aSquare(end);

set(0,'DefaultLineLineWidth',2)

%% Plotting optional
if max(aRamp) > 0
    maxfig = max(aRamp);
else
    maxfig = 50;
end
maxfig = inf;

figure(); set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.4, 0.28, 0.5]);
axis([1 T 0 maxfig]); axis manual; hold on
figa=gcf; axa=gca; colormap(gray);
set(figa,'Name','a - activation','numbertitle','off');
xlabel('Time (ms)'); ylabel('Activation')
set(gca,'FontSize',20,'fontname','times')

axes(axa); 
plot(aSquare,'r'); 
% plot(ASquare.*50,'-black')
plot(aRamp,'b');    %Plot a in time
% plot(ARamp.*50,'--black'); 
plot(hRamp,'--blue')
plot(hSquare,'--red')
legend('P Square','P Ramp','I Ramp','I Square','Fontsize',16)
% legend('P Square','A Square','P Ramp','A Ramp','I Ramp','I Square')

%Input figure
figure; aae = gca;
set(gcf,'Units','Normalized','OuterPosition', [0.25,0.25,0.28,0.2]);
plot(ASquare,'k-'); hold on;
plot(ARamp,'k--')
legend('A Square', 'A Ramp','FontSize',16)
xlim([0 5000]); ylim([0 max(ARamp)+0.02])
aae.XTick = [];
ylabel('Input')
set(gca,'FontSize',20,'fontname','times')

figure(); set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.25, 0.5, 0.75]);
axis([1 maxfig 0 maxfig]); axis manual; hold on
figc=gcf; axc=gca; colormap(gray);
set(figc,'Name','a vs. h phase plot','numbertitle','off');
xlabel('Activation of Pyramidal'); ylabel('Activation of Interneuron')
set(gca,'FontSize',20)

axes(axc); 
plot(aRamp,hRamp+0.15);   %Plot h by a attractor
plot(aSquare,hSquare);
legend('Ramp','Square')

