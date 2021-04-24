%% Sample Standalone CA3 Pattern Learning Code
% Simulates a learning and test phase for ripple activity in CA3 region
% under varying input waveform conditions
% LKW 4/24/21

%clearvars
%close all

%1 = FR IMA; 2 = DR IMA; 3 = FR IP; 4 = DR IP; 6 = BR IMA; 7 = BR IP
rampTypeFlag        = 1; disp(['Ramp Type ', num2str(rampTypeFlag)]);

cueN                = 1;            %Node index for cue
N                   = 15;           %# Nodes
Ww                  = 0.031;        %Baseline Pyr-Pyr weight strength
Hh                  = 0.05;         %Baseline Weight strength IN-Pyr
Wh                  = 0.05;         %Baseline Weight strength Pyr-IN
HAuto               = 0.003;        %Baseline Weight IN-IN
% tha                 = linspace(2,6,N);
tha                 = 4*ones(1,N);  %Activation Threshold Pyr
% thh                 = linspace(8,4);
thh                 = 4*ones(1,N);  %Activation Threshold IN
eta1                = 0.01;         %Pyr Decay constant
eta2                = 0.01;         %IN Decay constant
T                   = 1500;         %Time steps, must be even
Iexcit1             = 0.5;          %Learn phase afferent input
Iexcit2             = 1;            %Test phase afferent input
achLearn            = 0.10;         %Ach state during learn phase; 0 = high Ach; 1 = low Ach
achTest             = 1;            %Ach state during test phase
e1                  = 0.01;         %Learning Decay Pre synaptic
e2                  = 0.01;         %Learning Decay Post synaptic
e3                  = 0.001;        %Learning rate constant
onsetDelay          = 50;           %Wait time for each new pattern element to arise
learnDur            = 80;           %Duration (ms) of learn pulses
learnOL             = 44;           %Overlap (ms) of learn pattern elements
testDur             = 20;           %Duration (ms) of test cue pulse
rampPerc            = 0.25;         %Percentage of ramp
rampLen1            = round(learnDur*rampPerc);  %Duration ramp length learn phase
rampLen2            = round(testDur*rampPerc);                  %Duration ramp length test phase
AIn1                = zeros(N,T);   %Learn phase afferent input
AIn2                = zeros(N,T);   %Test phase
AIn2(cueN,onsetDelay:onsetDelay+testDur) = Iexcit2;  %Sq Pulse cue test phase
squarea             = Iexcit1*learnDur;     %Area under the square pulse

%Build Stimulation Protocols Based on RampTypeFlag
for i = 1:N
    tmpStart                                    = onsetDelay+1+(i-1)*(learnDur)-(i-1)*(learnOL);
    tmpEnd                                      = tmpStart+learnDur;
    if rampTypeFlag == 1                        %FRC
        AIn1(i,tmpStart:tmpEnd)                 = Iexcit1;
        AIn1(i,tmpStart:tmpStart+rampLen1-1)    = linspace(0,Iexcit1,rampLen1);
    elseif rampTypeFlag == 2                    %DRC
        AIn1(i,tmpStart:tmpEnd)                 = Iexcit1;
        AIn1(i,tmpStart:tmpStart+rampLen1-1)    = linspace(0,Iexcit1,rampLen1);
        AIn1(i,tmpEnd-rampLen1+1:tmpEnd)        = linspace(Iexcit1,0,rampLen1);
    elseif rampTypeFlag == 6
        AIn1(i,tmpStart:tmpEnd)                 = Iexcit1;
        AIn1(i,tmpEnd-rampLen1+1:tmpEnd)        = linspace(Iexcit1,0,rampLen1);
    end
    if rampTypeFlag == 3                        %FRP
        maxI                                    = squarea/(rampLen1/2 + learnDur - rampLen1);  %Calculate Single Ramp Current Max for constant AUC
        AIn1(i,tmpStart:tmpEnd)                 = maxI;
        AIn1(i,tmpStart:tmpStart+rampLen1-1)    = linspace(0,maxI,rampLen1);   %Front Ramp
    elseif rampTypeFlag == 4                    %DRP
        maxI                                    = squarea/(learnDur - rampLen1); %Calculate Double Ramp Current Max for constant AUC
        AIn1(i,tmpStart:tmpEnd)                 = maxI;
        AIn1(i,tmpStart:tmpStart+rampLen1-1)    = linspace(0,maxI,rampLen1);   %Front Ramp
        AIn1(i,tmpEnd-rampLen1+1:tmpEnd)        = linspace(maxI,0,rampLen1);   %Rear ramp
    elseif rampTypeFlag == 7                    %BRP
        maxI                                    = squarea/(rampLen1/2 + learnDur - rampLen1);
        AIn1(i,tmpStart:tmpEnd)                 = maxI;
        AIn1(i,tmpEnd-rampLen1+1:tmpEnd)        = linspace(maxI,0,rampLen1);   %Rear ramp
    end
end

W  = zeros(N,N,T);    %Init wt mat
% W = rand(N,N).*Ww;          %All rand wt mat
AH = zeros(N,N);    %Wt mat of a to h (feedforward activation of IN)
H  = zeros(N,N);    %Weight of h to a (feedback inhibition)

% % Pre-learn wt figure for manuscript
% figure; aac = gca; %set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.2, 0.2, 0.3, 0.5]); hold on;
% imagesc(flipud(W(:,:,1))); axis square; wcb = colorbar; wcb.Label.String = "Weight Strength";
% xlabel("Post-Synaptic CA3"); ylabel("Pre-Synaptic CA3")
% aac.YTick = 1:7:15; yticklabels(linspace(15,1,3));
% aac.XTick = 1:7:15; xticklabels(linspace(1,15,3));
% title("Learn Phase Start Network State")
% set(gca,'fontname','times','FontSize',24)

for i = 1:N         %Build weight mats
    H(i,i)  = Hh;   %Direct feedback IN to Pyr
    AH(i,i) = Wh;   %Direct excitation pyr to IN
end

%Initialize activity vectors for learn/test phases
aLearn = zeros(N,T);
hLearn = zeros(N,T);
aTest  = zeros(N,T);
hTest  = zeros(N,T);

%% Run learn block

for t=1:T
    %Neuron Activity
    for j = 1:N
        if (aLearn(j,t)   - tha(j)) > 0; aaLearn(j)   = aLearn(j,t)   - tha(j); else aaLearn(j)   = 0; end %Threshold activity over 'tha' to 0
        if (hLearn(j,t)   - thh(j)) > 0; hhLearn(j)   = hLearn(j,t)   - thh(j); else hhLearn(j)   = 0; end %Threshold activity over 'thh' to 0
    end
    
    daLearn = AIn1(:,t) + achLearn.*(aaLearn*W(:,:,t))' - (hhLearn*H)' - eta1.*aLearn(:,t); %change in a = Input A(t) + Weight*thesholded activity (a) - Weight*thesholded activity (h) - decay*activity (a)
    dhLearn = (aaLearn*AH)' - HAuto.*hhLearn' - eta2.*hLearn(:,t);          %change in h = Weight*thresholded activity (a) - decay*activity (h)

    aLearn(:,t+1) = aLearn(:,t)+daLearn; %Update next time point of a
    hLearn(:,t+1) = hLearn(:,t)+dhLearn; %Update next h time point
    
    %Weight Update
    dW = e3.*(1-achLearn).*(Ww-W(:,:,t)).*(aaLearn).*(aaLearn');           %Use shunt to impose arbitrary ceiling of Ww

    W(:,:,t+1) = W(:,:,t)+dW;   %Update next Weight time point
end

%% Plot Learn

set(0,'DefaultLineLineWidth',2)
cPyr = autumn(N);
cIN  = jet(N);
cAIn = gray(N+2);

% Plot learn afferent input
figure; aaa = gca;
imagesc(flipud(AIn1));
% title('Afferent Input Pattern');
acb = colorbar; 
aaa.YTick = 1:7:15; yticklabels(linspace(15,1,3));
set(gca,'fontname','times','FontSize',24)

% % Plot learn activity
% figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.5, 0.4]); hold on;
% for i = 1:N
%     plot(aLearn(i,:), 'Color', cPyr(i,:));
%     plot(AIn1(i,:)*5, 'Color', cAIn(i,:));
% end
% xlabel('Time (ms)'); ylabel('Activation'); 
% if rampTypeFlag == 1 && rampPerc == 0
%     legend('Square','Input A Learn','location','southeast')
% elseif rampTypeFlag == 1
%     legend('FR, IMA','location','southeast')
% elseif rampTypeFlag == 2
%     legend('DR, IMA','location','southeast')
% elseif rampTypeFlag == 3
%     legend('FR, IP','location','southeast')
% elseif rampTypeFlag == 4
%     legend('DR, IP','location','southeast')
% elseif rampTypeFlag == 6
%     legend('BR, IMA','location','southeast')
% elseif rampTypeFlag == 7
%     legend('BR, IP','location','southeast')
% end
% % title('Learn Phase')
% set(gca,'fontname','times','FontSize',24);

% Plot Post Learn Weights
figure; aaw = gca; %set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.2, 0.2, 0.3, 0.5]); hold on;
imagesc(flipud(W(:,:,end))); axis square; wcb = colorbar; 
%     xlabel("Post-Synaptic CA3"); ylabel("Pre-Synaptic CA3")
%     wcb.Label.String = "Weight Strength";
aaw.YTick = 1:7:15; yticklabels(linspace(15,1,3));
aaw.XTick = 1:7:15; xticklabels(linspace(1,15,3));
% title("Learn Phase Final Network State")
set(gca,'fontname','times','FontSize',24)

W(:,:,1) = W(:,:,end);       %Set starting time point to learn mat for test phase

%% Run test block

for t=1:T
    %Neuron Activity
    for j = 1:N
        if (aTest(j,t)   - tha(j)) > 0; aaTest(j)   = aTest(j,t)   - tha(j); else aaTest(j)   = 0; end %Threshold activity over 'tha' to 0
        if (hTest(j,t)   - thh(j)) > 0; hhTest(j)   = hTest(j,t)   - thh(j); else hhTest(j)   = 0; end %Threshold activity over 'thh' to 0
    end
    daTest = AIn2(:,t) + achTest.*(aaTest*W(:,:,t))' - (hhTest* H)' - eta1.*aTest(:,t); %change in a = Input A(t) + Weight*thesholded activity (a) - Weight*thesholded activity (h) - decay*activity (a)
    dhTest = (aaTest*AH)' - HAuto.*hhTest' - eta2.*hTest(:,t);          %change in h = Weight*thresholded activity (a) - decay*activity (h)

    aTest(:,t+1) = aTest(:,t)+daTest; %Update next time point of a
    hTest(:,t+1) = hTest(:,t)+dhTest; %Update next h time point
    
    %Weight Update (suppressed by Ach state)
    dW = e3.*(1-achTest).*(Ww-W(:,:,t)).*(aaTest).*(aaTest');           %Use shunt to impose arbitrary ceiling of Ww

    W(:,:,t+1) = W(:,:,t)+dW;   %Update next Weight time point
end

%% Calculate Peaks and Locs of threshold crossing
pksT = zeros(N,1); locsT = zeros(N,1);
actThresh = 10*ones(N,1);
for  i = 1:N
    [pksTmp,locsTmp] = findpeaks(aTest(i,:));
    if ~isempty(locsTmp); pksT(i) = pksTmp(1); locsT(i) = locsTmp(1); end
end
ttt = [];
for j = 1:numel(pksT)     %For each peak in the sim
    indtmp = find(aTest(j,:)>actThresh(j),1);
    if ~isempty(indtmp)
        ttt = [ttt,indtmp];
    end
end
    
%% Plot Test options

% % Test Phase Afferent input figure for manuscript
% figure; aaa = gca;
% imagesc(flipud(AIn2));
% xlabel('Time (ms)'); ylabel('Node #'); title('Afferent Input Pattern');
% aaa.YTick = 1:7:15; yticklabels(linspace(15,1,3));
% acb = colorbar; acb.Label.String = "Afferent Strength";
% set(gca,'fontname','times','FontSize',20)
% breakxaxis([250 1250],0.03);

aat = figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.5, 0.35]);
% subplot(2,1,1)
for i = 1:N
    plot(aTest(i,:),'Color',cPyr(i,:)); hold on;
    plot(AIn2(i,:)*5, 'Color', cAIn(i,:));
end
plot([0,T-1],[10,10],'k--','LineWidth',1)
scatter(ttt,ttt*0+10,'^','k','filled');
% xlabel('Time (ms)'); ylabel('Activation'); 
if rampTypeFlag == 1 && rampPerc == 0
    legend('Square Pyramidal','Test Cue','location','southeast')
elseif rampTypeFlag == 1
    legend('FR IMA Pyramidal','location','southeast')
elseif rampTypeFlag == 2
    legend('DR IMA Pyramidal','location','southeast')
elseif rampTypeFlag == 3
    legend('FR IP Pyramidal','location','southeast')
elseif rampTypeFlag == 4
    legend('DR IP Pyramidal','location','southeast')
elseif rampTypeFlag == 6
    legend('BR IMA Pyramidal','location','southeast')
elseif rampTypeFlag == 7
    legend('BR IP Pyramidal','location','southeast')
end
% title('Test Phase Pyrs')
ylim([-30,25])
set(gca,'fontname','times','FontSize',24);

%% Saves
% disp('saving vars and figs')
% saveDir = 'Your\Save\Dir\Here\';
% OLStr = ['OL',num2str(learnOL), 'ms'];
% if rampTypeFlag == 1 && rampLen1 == 0
%     sBase = ['Square_',OLStr];
%     saveas(aat,[saveDir,sBase,'_TestActivity'],'fig')
%     saveas(aat,[saveDir,sBase,'_TestActivity'],'png')
%     saveas(aat,[saveDir,sBase,'_TestActivity'],'svg')
%     saveas(aaa,[saveDir,sBase,'_Input'],'fig')
%     saveas(aaa,[saveDir,sBase,'_Input'],'png')
%     saveas(aaa,[saveDir,sBase,'_Input'],'svg')
%     saveas(aaw,[saveDir,sBase,'_FinalWeight'],'fig')
%     saveas(aaw,[saveDir,sBase,'_FinalWeight'],'png')
%     saveas(aaw,[saveDir,sBase,'_FinalWeight'],'svg')
% elseif rampTypeFlag == 1
%     sBase = ['FR_IMA_',OLStr];
%     saveas(aat,[saveDir,sBase,'_TestActivity'],'fig')
%     saveas(aat,[saveDir,sBase,'_TestActivity'],'png')
%     saveas(aat,[saveDir,sBase,'_TestActivity'],'svg')
%     saveas(aaa,[saveDir,sBase,'_Input'],'fig')
%     saveas(aaa,[saveDir,sBase,'_Input'],'png')
%     saveas(aaa,[saveDir,sBase,'_Input'],'svg')
%     saveas(aaw,[saveDir,sBase,'_FinalWeight'],'fig')
%     saveas(aaw,[saveDir,sBase,'_FinalWeight'],'png')
%     saveas(aaw,[saveDir,sBase,'_FinalWeight'],'svg')
% elseif rampTypeFlag == 2
%     sBase = ['DR_IMA_',OLStr];
%     saveas(aat,[saveDir,sBase,'_TestActivity'],'fig')
%     saveas(aat,[saveDir,sBase,'_TestActivity'],'png')
%     saveas(aat,[saveDir,sBase,'_TestActivity'],'svg')
%     saveas(aaa,[saveDir,sBase,'_Input'],'fig')
%     saveas(aaa,[saveDir,sBase,'_Input'],'png')
%     saveas(aaa,[saveDir,sBase,'_Input'],'svg')
%     saveas(aaw,[saveDir,sBase,'_FinalWeight'],'fig')
%     saveas(aaw,[saveDir,sBase,'_FinalWeight'],'png')
%     saveas(aaw,[saveDir,sBase,'_FinalWeight'],'svg')
% elseif rampTypeFlag == 6
%     sBase = ['BR_IMA_',OLStr];
%     saveas(aat,[saveDir,sBase,'_TestActivity'],'fig')
%     saveas(aat,[saveDir,sBase,'_TestActivity'],'png')
%     saveas(aat,[saveDir,sBase,'_TestActivity'],'svg')
%     saveas(aaa,[saveDir,sBase,'_Input'],'fig')
%     saveas(aaa,[saveDir,sBase,'_Input'],'png')
%     saveas(aaa,[saveDir,sBase,'_Input'],'svg')
%     saveas(aaw,[saveDir,sBase,'_FinalWeight'],'fig')
%     saveas(aaw,[saveDir,sBase,'_FinalWeight'],'png')
%     saveas(aaw,[saveDir,sBase,'_FinalWeight'],'svg')
% elseif rampTypeFlag == 3
%     sBase = ['FR_IP_',OLStr];
%     saveas(aat,[saveDir,sBase,'_TestActivity'],'fig')
%     saveas(aat,[saveDir,sBase,'_TestActivity'],'png')
%     saveas(aat,[saveDir,sBase,'_TestActivity'],'svg')
%     saveas(aaa,[saveDir,sBase,'_Input'],'fig')
%     saveas(aaa,[saveDir,sBase,'_Input'],'png')
%     saveas(aaa,[saveDir,sBase,'_Input'],'svg')
%     saveas(aaw,[saveDir,sBase,'_FinalWeight'],'fig')
%     saveas(aaw,[saveDir,sBase,'_FinalWeight'],'png')
%     saveas(aaw,[saveDir,sBase,'_FinalWeight'],'svg')
% elseif rampTypeFlag == 4
%     sBase = ['DR_IP_',OLStr];
%     saveas(aat,[saveDir,sBase,'_TestActivity'],'fig')
%     saveas(aat,[saveDir,sBase,'_TestActivity'],'png')
%     saveas(aat,[saveDir,sBase,'_TestActivity'],'svg')
%     saveas(aaa,[saveDir,sBase,'_Input'],'fig')
%     saveas(aaa,[saveDir,sBase,'_Input'],'png')
%     saveas(aaa,[saveDir,sBase,'_Input'],'svg')
%     saveas(aaw,[saveDir,sBase,'_FinalWeight'],'fig')
%     saveas(aaw,[saveDir,sBase,'_FinalWeight'],'png')
%     saveas(aaw,[saveDir,sBase,'_FinalWeight'],'svg')
% elseif rampTypeFlag == 7
%     sBase = ['BR_IP_',OLStr];
%     saveas(aat,[saveDir,sBase,'_TestActivity'],'fig')
%     saveas(aat,[saveDir,sBase,'_TestActivity'],'png')
%     saveas(aat,[saveDir,sBase,'_TestActivity'],'svg')
%     saveas(aaa,[saveDir,sBase,'_Input'],'fig')
%     saveas(aaa,[saveDir,sBase,'_Input'],'png')
%     saveas(aaa,[saveDir,sBase,'_Input'],'svg')
%     saveas(aaw,[saveDir,sBase,'_FinalWeight'],'fig')
%     saveas(aaw,[saveDir,sBase,'_FinalWeight'],'png')
%     saveas(aaw,[saveDir,sBase,'_FinalWeight'],'svg')
% end
