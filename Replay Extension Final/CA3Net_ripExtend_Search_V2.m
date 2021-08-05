%%% Sample code for Replay Extension Model
%LKW 8/2/21
%Relies on ripExtend_fast_V2.m and ripExtend_CohenStats.m
%Searches ripple prolongation parameter space in dimensions of ramp degree
%and input pulse duration
%Plots various metrics and statistics for evaluating network performance
%under varying waveform shape input

clearvars
% close all

pStruct.rampTypeFlag        = 1;        %1 = FR IMA; 2 = DR IMA; 3 = FR IP; 4 = DR IP; 5 = BR IMA; 7 = BR IP;
pStruct.simTypeFlag         = 2;        %1 = Linear; 2 = Linear with Adaptation; 3 = Nonlinear; 4 = Nonlinear with Adaptation

%Setup neuron and weight paramters
pStruct.N                   = 15;               %Nodes per region
N                           = 15;
pStruct.Ww                  = 0.031;            %Weight strength pyr to pyr
pStruct.Hh                  = 0.035;            %Weight strength IN to Pyr
pStruct.Wh                  = 0.05;             %Weight strength pyr to IN
pStruct.HAuto               = 0.003;            %Weight strength IN to IN
pStruct.tha                 = 4*ones(1,N);
pStruct.thh                 = 4*ones(1,N);
pStruct.eta                 = 0.01;             %Decay constant
pStruct.actThresh           = 10*ones(1,N);     %High node firing threshold

%Setup input parameters
pStruct.T                   = 1000;             %Time steps, must be even
pStruct.cueN                = 1;
pStruct.Iexcit1             = 1;                %Cue pulse strength
pStruct.Iexcit2             = 0.09;             %Opto pulse strength
pStruct.inDur1              = 20;               %Duration for cue
pStruct.inDur2              = 100;              %Opto pulse duration
pStruct.rampLen             = pStruct.inDur2*0.0;  %Duration ramp length e.g. 0.0 to 1.0
pStruct.onsetDelay          = 50;               %Wait time to ripple start from sim start
pStruct.stimDelay           = 200;              %Wait time to opto pulse from sim start; try 130 (+50 +20)

% Ionic Currents and related parameters 
pStruct.mu                  = 0.01;      %Ca-dependent K-current
pStruct.gm                  = 0.001;     %Gamma; voltage-dependent Ca-currents
pStruct.om                  = 0.001;     %Omega; constant for diffusion of intracellular Ca
pStruct.thc                 = 4*ones(1,N);

wtBias                      = linspace(0.004,-0.004,N);
% wtBias                      = zeros(1,N);
W                           = zeros(N,N);           %Init wt mat
% W                           = rand(N,N).*Ww;        %All rand wt mat
AH                          = zeros(N,N);           %Wt mat of a to h (feedforward activation of IN)
H                           = zeros(N,N);           %Weight of h to a (feedback inhibition)
for i = 1:N         %Build weight mats
    W(i,i)  = pStruct.Ww+wtBias(i);   %Autorecurrency
    H(i,i)  = pStruct.Hh;   %Direct feedback IN to Pyr
    AH(i,i) = pStruct.Wh;   %Direct excitation pyr to IN
    if i <= N - 1   %Forward 1
        W(i,i+1)  = (pStruct.Ww+wtBias(i))/2;   %Pyr 2 Pyr
%         H(i,i+1)  = Hh/2;   %IN  2 Pyr
%         AH(i,i+1) = Wh/2;   %Pyr 2 IN
    end
    if i <= N - 2  %Forward 2
        W(i,i+2)  = (pStruct.Ww+wtBias(i))/4;   %Pyr 2 Pyr
%         H(i,i+2)  = Hh/4;   %IN  2 Pyr
%         AH(i,i+2) = Wh/4;   %Pyr 2 IN
    end
    if i > 1    %Back 1
%         W(i,i-1)  = Ww/2;   %Pyr 2 Pyr
%         H(i,i-1)  = Hh/4;   %IN  2 Pyr
%         AH(i,i-1) = Wh/2;   %Pyr 2 IN
    end
    if i > 2    %Back 2   
%         W(i,i-2)  = Ww/4;   %Pyr 2 Pyr
%         H(i,i-2)  = Hh/4;   %IN  2 Pyr
%         AH(i,i-2) = Wh/4;   %Pyr 2 IN
    end
end

pStruct.W = W; pStruct.H = H; pStruct.AH = AH;  %Add wtmats to pStruct

%Set up search variables and parameters
nPs                     = 101;
pLim                    = 250;      %250 for inDur2; 500 for stimDelay
pFloor                  = 0;
nRamps                  = 21;
pVect                   = linspace(pFloor,pLim,nPs);    %Vector of parameter values
if mod(pStruct.rampTypeFlag,2) == 1; rampLim = 1; else; rampLim = 0.5; end %Set ramp limit to 1 or 0.5 depending on flag
rampPercs = linspace(0,rampLim,nRamps);

outputs.ds          = zeros(nRamps,nPs);
outputs.cs          = cell(nRamps,nPs);
outputs.means       = zeros(nRamps,nPs);
outputs.SDs         = zeros(nRamps,nPs);
outputs.Ns          = zeros(nRamps,nPs);

%% Gut check above parameters
[actCell,hactCell,inCell] = ripExtend_fast_V2(pStruct);    %Using ramp defined above
aRamp = actCell{1}; aControl = actCell{2};
hRamp = hactCell{1};hControl = hactCell{2};
ARamp = inCell{1};  AControl = inCell{2};
dGut = ripExtend_CohenStats(actCell,hactCell,pStruct);

% Calculate Peaks
pksR = zeros(N,1); locsR = zeros(N,1);
pksC = zeros(N,1); locsC = zeros(N,1);

for  i = 1:N
    [pksTmp,locsTmp] = findpeaks(aRamp(i,:));
    if ~isempty(locsTmp); pksR(i) = pksTmp(1); locsR(i) = locsTmp(1); end
    [pksTmp,locsTmp] = findpeaks(aControl(i,:));
    if ~isempty(locsTmp); pksC(i) = pksTmp(1); locsC(i) = locsTmp(1); end
end
pksC = pksC(pksC>0); locsC = locsC(locsC>0);
locsCell = {locsR,locsC};
pksCell = {pksR,pksC};

%Threshold method
actThresh = 10*ones(1,N); % tha;
tttCell = {};

for i = 1:2
    pksTmp = pksCell{i};        %All the peaks for that sim (e.g. ramp)
    actTmp = actCell{i};        %Pyramidal activation for that sim e.g. NxT
    ttt = [];
    for j = 1:numel(pksTmp)     %For each peak in the sim
        indtmp = find(actTmp(j,:)>actThresh(j),1);
        if ~isempty(indtmp)
            ttt = [ttt,indtmp];
        end
    end
    tttCell(i) = {ttt};
end

cIthI = mean(diff(tttCell{2})); %Get inter-threshold-interval for control

[thaPlot,thax1,thax2] = plotCA3Ripples_fast(actCell,hactCell,inCell,locsCell,pksCell,pStruct.cueN);
set(thaPlot,'defaultLegendAutoUpdate','off');
axCell = {thax1,thax2};
thaLength = 20;
for i = 1:2
    axes(axCell{i});
    tmpTTT = tttCell{i};
    plot([0,999],[10,10],'k--','LineWidth',1)
    for j = 1:length(tmpTTT)
        scatter(tmpTTT(j),actThresh(j),'k^','filled')
    end
end

%% Run Block
for i = 1:nRamps
    pStruct.tmpRamp = rampPercs(i);         %Define ramp percentage
    for j = 1:nPs
        pStruct.inDur2      = pVect(j); 	%2nd search param - inDur
%         pStruct.stimDelay   = pVect(j); %2nd search param - stimDelay
        pStruct.rampLen     = round(pStruct.tmpRamp*pStruct.inDur2);   %Define ramp length
        %Run function
        [actTmp,hactTmp,inTmp] = ripExtend_fast_V2(pStruct);     %Run simulation
        %Calculate stats
        [outputs.ds(i,j),~,~,outputs.Ns(i,j)] = ripExtend_CohenStats(actTmp,hactTmp,pStruct);
    end
end

%% Clean data
outputs.dsNanCor = outputs.ds; tmpMax = max(max(abs(outputs.dsNanCor)));
outputs.dsNanCor(isnan(outputs.dsNanCor)) = tmpMax;  %Convert any nans to artificial P-value 1 for plotting purposes
outputs.meanDiff = abs(outputs.means - cIthI);
% Other metrics
minDLocs = minTimes(outputs.dsNanCor,pVect,nRamps); %Calculate duration of minimum effect size
minDs = min(flipud(abs(outputs.dsNanCor')));    %Calculate minimum effect size

if max(max(outputs.dsNanCor)) > 5
    cbMax = 5;
else
    cbMax = max(max(outputs.dsNanCor));
end
%% Plotting
%Plot map of d-scores
figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.35, 0.65]); aaa = gca;
imagesc(flipud(abs(outputs.dsNanCor')),[0,cbMax]);
aaa.YTick = 1:round((nPs-1)/5):nPs;   %Reversed 0 at top, max at bottom
aaa.XTick = 1:round(nRamps-1)/5:nRamps;
yticklabels(linspace(pLim,pFloor,6));
% xticklabels(0:rampLim*20:rampLim*100);      %For Ramp Percentage Label
xticklabels(0:20:100);
% xlabel('Ramp Percentage'); ylabel('Input Duration (ms)'); 
colormap hot; axis square; hcb = colorbar;
% if mod(pStruct.rampTypeFlag,2) == 1; title("Forward Ramp Opto vs Control");
% else; title("Double Ramp Opto vs Control"); end
% hcb.Label.String = "Cohen's d";
set(aaa,'FontSize',24,'fontname','times')

%Plot map of sequence lengths
figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.35, 0.65]); aaz = gca;
imagesc(flipud(abs(outputs.Ns')),[0,N]);
aaz.YTick = 1:round((nPs-1)/5):nPs;   %Reversed 0 at top, max at bottom
aaz.XTick = 1:round(nRamps-1)/5:nRamps;
yticklabels(linspace(pLim,pFloor,6));
% xticklabels(0:rampLim*20:rampLim*100);      %For Ramp Percentage Label
xticklabels(0:20:100);
% xlabel('Ramp Percentage'); ylabel('Input Duration (ms)'); 
axis square; hcb = colorbar;
% hcb.Label.String = "Sequence Length";
set(aaz,'FontSize',24,'fontname','times')

%% Shuffles
permN = 1000;
z = 1.96;   %confidence level 95% = 1.96, 99% = 2.58
durShufDistro = datasample(flipud(abs(outputs.dsNanCor')),permN,1);
shuf.durShuf = mean(durShufDistro);
seqShufDistro = datasample(flipud(outputs.Ns'),permN,1);
shuf.seqShuf = mean(seqShufDistro);
shuf.sdDurShuffle = std(durShufDistro); 
shuf.sdSeqShuffle = std(seqShufDistro);

% 95% CIs
shuf.upCIdur = shuf.durShuf + z*shuf.sdDurShuffle/sqrt(permN);
shuf.dnCIdur = shuf.durShuf - z*shuf.sdDurShuffle/sqrt(permN);
shuf.upCIseq = shuf.seqShuf + z*shuf.sdSeqShuffle/sqrt(permN);
shuf.dnCIseq = shuf.seqShuf - z*shuf.sdSeqShuffle/sqrt(permN);

% Plot Basic Shuffle Figures
figure;
plot(rampPercs*100,shuf.seqShuf,'k',rampPercs*100,shuf.upCIseq,'k--',rampPercs*100,shuf.dnCIseq,'k--')
ylim([0 15]);
% Minimum effect sizes
figure;
plot(rampPercs*100,minDs)
%% Other variations of describing the ripple extension 

% %Plot map of raw MEANS
% figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.35, 0.65]); aaa = gca;
% imagesc(flipud(abs(outputs.means')),[0,max(max(outputs.means))]);
% aaa.YTick = 1:round((nPs-1)/5):nPs;   %Reversed 0 at top, max at bottom
% aaa.XTick = 1:round(nRamps-1)/5:nRamps;
% yticklabels(linspace(pLim,pFloor,6));
% % xticklabels(0:rampLim*20:rampLim*100);      %For Ramp Percentage Label
% xticklabels(0:20:100);
% xlabel('Ramp Percentage'); ylabel('Input Duration (ms)'); 
% colormap hot; axis square; hcb = colorbar;
% if mod(pStruct.rampTypeFlag,2) == 1; title("Forward Ramp Mean IThI");
% else; title("Double Ramp Opto vs Control"); end
% hcb.Label.String = "Mean IThI (ms)";
% set(aaa,'FontSize',24,'fontname','times')

% %Plot MEAN box plots collapsed for input duration
% figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.5, 0.8]); aad = gca;  %[0.35, 0.39, 0.39, 0.55]
% % boxplot(abs(outputs.ds'),'PlotStyle','compact');
% plot([1,21],[cIthI,cIthI],'--k'); hold on; 
% boxplot(abs(outputs.means'));
% xticklabels(0:rampLim/(nRamps-1)*100:rampLim*100); xtickangle(90)
% ylabel("Raw Mean IThI (ms)"); xlabel("Ramp Percentage");
% if mod(pStruct.rampTypeFlag,2) == 1; title("Forward Ramp Means Collapse Input Duration");
% else; title("Double Ramp Means Collapse Input Duration"); end
% legend('Control Mean IThI')
% set(aad,'FontSize',24,'fontname','times')
% 
% %Plot map of MEAN - CONTROL
% figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.35, 0.65]); aaa = gca;
% imagesc(flipud(abs(outputs.meanDiff')),[0,max(max(outputs.meanDiff))]);
% aaa.YTick = 1:round((nPs-1)/5):nPs;   %Reversed 0 at top, max at bottom
% aaa.XTick = 1:round(nRamps-1)/5:nRamps;
% yticklabels(linspace(pLim,pFloor,6));
% % xticklabels(0:rampLim*20:rampLim*100);      %For Ramp Percentage Label
% xticklabels(0:20:100);
% xlabel('Ramp Percentage'); ylabel('Input Duration (ms)'); 
% colormap hot; axis square; hcb = colorbar;
% if mod(pStruct.rampTypeFlag,2) == 1; title("Forward Ramp Mean - Control IThI");
% else; title("Double Ramp Mean - Control IThI"); end
% hcb.Label.String = "Mean IThI \Delta (ms)";
% set(aaa,'FontSize',24,'fontname','times')

% %Plot MEAN - CONTROL box plots collapsed for input duration
% figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.5, 0.8]); aaf = gca;
% % boxplot(abs(outputs.ds'),'PlotStyle','compact');
% boxplot(abs(outputs.meanDiff'));
% xticklabels(0:rampLim/(nRamps-1)*100:rampLim*100); xtickangle(90)
% ylabel("Mean - Control IthI (ms)"); xlabel("Ramp Percentage");
% if mod(pStruct.rampTypeFlag,2) == 1; title("Forward Ramp Means - Control Collapse Input Duration");
% else; title("Double Ramp Means - Control Collapse Input Duration"); end
% set(aaf,'FontSize',24,'fontname','times')

%% Supplementary Graphs
set(0,'DefaultLineLineWidth',2)
% 
% % Afferent Input 3D plot code
% sqGhost = zeros(1,pStruct.T);
% sqGhost(pStruct.stimDelay+1:pStruct.stimDelay+pStruct.inDur2) = pStruct.Iexcit2;
% figure(); set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.28, 0.35, 0.45]); gaa = gca;
% w3 = waterfall(1:pStruct.T,1:N,ARamp); hold on;
% w3.LineWidth = 3;
% % w3Sq = waterfall(1:pStruct.T,1,sqGhost);
% % w3Sq.LineStyle = '--'; w3Sq.EdgeColor = 'r'; w3Sq.LineWidth = 2;
% % xlabel('Time')
% gaa.YTick = 1:7:15;
% % ylabel('Node'); zlabel('Amplitude')
% set(gca,'FontSize',24,'fontname','times')

% % Combined Ramp Percentage Shuffles - must load in previous data
% patchXs = [1:nRamps,linspace(nRamps,1,nRamps)]; %Vector of x coords;
% figure(); set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.28, 0.365, 0.65]); aac = gca;
% axis square; hold on;
% plot(FR_IMA_shufs.durShuf,'r-');
% plot(DR_IMA_shufs.durShuf,'b-'); 
% plot(BR_IMA_shufs.durShuf,'c-');
% plot(FR_IP_shufs.durShuf,'r--');
% plot(DR_IP_shufs.durShuf,'b--');
% plot(BR_IP_shufs.durShuf,'c--');
% patch(patchXs,[FR_IMA_shufs.dnCIdur,fliplr(FR_IMA_shufs.upCIdur)],'k','EdgeColor','none','FaceAlpha',0.1);
% patch(patchXs,[DR_IMA_shufs.dnCIdur,fliplr(DR_IMA_shufs.upCIdur)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
% patch(patchXs,[BR_IMA_shufs.dnCIdur,fliplr(BR_IMA_shufs.upCIdur)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
% patch(patchXs,[FR_IP_shufs.dnCIdur,fliplr(FR_IP_shufs.upCIdur)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
% patch(patchXs,[DR_IP_shufs.dnCIdur,fliplr(DR_IP_shufs.upCIdur)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
% patch(patchXs,[BR_IP_shufs.dnCIdur,fliplr(BR_IP_shufs.upCIdur)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
% xlabel('Ramp Percentage'); ylabel('Mean Effect Size')
% legend('FR IMA','DR IMA','BR IMA','FR IP','DR IP','BR IP','95% CI','FontSize',16,'location','southwest');
% ylim([0 3]); xlim([1 nRamps]); aac.XTick = linspace(0,nRamps,6); xticklabels(0:20:100)
% set(gca,'FontSize',24,'fontname','times')
% 
% % Combined SeqLen Shuffles - must load in previous data
% patchXs = [1:nRamps,linspace(nRamps,1,nRamps)]; %Vector of x coords;
% figure(); set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.28, 0.365, 0.65]); aac = gca;
% axis square; hold on;
% plot(FR_IMA_shufs.seqShuf,'r-');
% plot(DR_IMA_shufs.seqShuf,'b-'); 
% plot(BR_IMA_shufs.seqShuf,'c-');
% plot(FR_IP_shufs.seqShuf,'r--');
% plot(DR_IP_shufs.seqShuf,'b--');
% plot(BR_IP_shufs.seqShuf,'c--');
% patch(patchXs,[FR_IMA_shufs.dnCIseq,fliplr(FR_IMA_shufs.upCIseq)],'k','EdgeColor','none','FaceAlpha',0.1);
% patch(patchXs,[DR_IMA_shufs.dnCIseq,fliplr(DR_IMA_shufs.upCIseq)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
% patch(patchXs,[BR_IMA_shufs.dnCIseq,fliplr(BR_IMA_shufs.upCIseq)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
% patch(patchXs,[FR_IP_shufs.dnCIseq,fliplr(FR_IP_shufs.upCIseq)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
% patch(patchXs,[DR_IP_shufs.dnCIseq,fliplr(DR_IP_shufs.upCIseq)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
% patch(patchXs,[BR_IP_shufs.dnCIseq,fliplr(BR_IP_shufs.upCIseq)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
% xlabel('Ramp Percentage'); ylabel('Mean Sequence Length')
% legend('FR IMA','DR IMA','BR IMA','FR IP','DR IP','BR IP','95% CI','FontSize',16,'location','southwest');
% ylim([5 15]); xlim([1 nRamps]); aac.XTick = linspace(0,nRamps,6); xticklabels(0:20:100)
% set(gca,'FontSize',24,'fontname','times')
% 
% %Must Load in previously saved minDLoc data
% figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.35, 0.65]); vvv = gca;
% plot(FR_IMA_minLocs,'r-'); hold on;
% plot(DR_IMA_minLocs,'b-');
% plot(BR_IMA_minLocs,'c-');
% plot(FR_IP_minLocs,'r--');
% plot(DR_IP_minLocs-1,'b--');
% plot(BR_IP_minLocs+1,'c--');
% vvv.XTick = 1:round(nRamps-1)/5:nRamps;
% xticklabels(0:20:100);
% xlabel('Ramp Percentage'); ylabel('Input Duration (ms)'); 
% axis square; 
% % title('Input Duration of Minimum Effect Size')
% legend('FR IMA','DR IMA','BR IMA','FR IP','DR IP','BR IP','FontSize',16,'location','northwest');
% ylim([0 101]); xlim([1 nRamps]);
% set(vvv,'FontSize',24,'fontname','times')
% 
% %Must Load in previously saved minD data
% figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.35, 0.65]); vvx = gca;
% plot(FR_IMA_minDs,'r-'); hold on;
% plot(DR_IMA_minDs,'b-');
% plot(BR_IMA_minDs,'c-');
% plot(FR_IP_minDs,'r--');
% plot(DR_IP_minDs,'b--');
% plot(BR_IP_minDs,'c--');
% vvx.XTick = 1:round(nRamps-1)/5:nRamps;
% xticklabels(0:20:100);
% xlabel('Ramp Percentage'); ylabel('Least Temporal Disruption'); 
% axis square; 
% % title('Minimum Effect Size')
% legend('FR IMA','DR IMA','BR IMA','FR IP','DR IP','BR IP','FontSize',16,'location','northwest');
% ylim([0 0.5]); xlim([1 nRamps]);
% set(vvx,'FontSize',24,'fontname','times')

% %Plot Weights
% figure(); set(gcf, 'Colormap', jet) % tightfig();
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.5, 0.75, 0.5]);
% subplot(1,3,1)
% imagesc(W); title('Pyr-Pyr Weight (W)');    axis square
% xlabel('Post-Synaptic Pyr'); ylabel('Pre-Synaptic Pyr')
% set(gca,'FontSize',20)
% subplot(1,3,2)
% imagesc(AH); title("Pyr-IN Weight (W')");   axis square
% xlabel('Post-Synaptic IN'); ylabel('Pre-Synaptic Pyr')
% set(gca,'FontSize',20)
% subplot(1,3,3)
% imagesc(H); title('IN-Pyr Weight (H)');     axis square
% xlabel('Post-Synaptic Pyr'); ylabel('Pre-Synaptic IN')
% set(gca,'FontSize',20)

%% Saves
% disp('saving vars and figs')
% saveDir = 'Your\Dir\Here\';
% if z == 1.96
%     ciStr = '_CI95';
% elseif z == 2.58
%     ciStr = '_CI99';
% end
% 
% if pStruct.rampTypeFlag == 1
%     sBase = ['FR_IMA_', num2str(permN),ciStr];
%     FR_IMA_minLocs  = minDLocs;
%     FR_IMA_minDs    = minDs;
%     FR_IMA_shufs    = shuf;
%     fName = [saveDir, sBase];
%     save(fName,'FR_IMA_minLocs','FR_IMA_minDs','FR_IMA_shufs')
%     saveas(aaa,[saveDir,'FR_IMA_heatmap_ds'],'fig')
%     saveas(aaa,[saveDir,'FR_IMA_heatmap_ds'],'png')
%     saveas(aaa,[saveDir,'FR_IMA_heatmap_ds'],'svg')
%     saveas(aaz,[saveDir,'FR_IMA_heatmap_seqs'],'fig')
%     saveas(aaz,[saveDir,'FR_IMA_heatmap_seqs'],'png')
%     saveas(aaz,[saveDir,'FR_IMA_heatmap_seqs'],'svg')
% elseif pStruct.rampTypeFlag == 2
%     sBase = ['DR_IMA_', num2str(permN),ciStr];
%     DR_IMA_minLocs  = minDLocs;
%     DR_IMA_minDs    = minDs;
%     DR_IMA_shufs    = shuf;
%     fName = [saveDir, sBase];
%     save(fName, 'DR_IMA_minLocs','DR_IMA_minDs','DR_IMA_shufs')
%     saveas(aaa,[saveDir,'DR_IMA_heatmap_ds'],'fig')
%     saveas(aaa,[saveDir,'DR_IMA_heatmap_ds'],'png')
%     saveas(aaa,[saveDir,'DR_IMA_heatmap_ds'],'svg')
%     saveas(aaz,[saveDir,'DR_IMA_heatmap_seqs'],'fig')
%     saveas(aaz,[saveDir,'DR_IMA_heatmap_seqs'],'png')
%     saveas(aaz,[saveDir,'DR_IMA_heatmap_seqs'],'svg')
% elseif pStruct.rampTypeFlag == 5
%     sBase = ['BR_IMA', num2str(permN),ciStr];
%     BR_IMA_minLocs  = minDLocs;
%     BR_IMA_minDs    = minDs;
%     BR_IMA_shufs    = shuf;
%     fName = [saveDir, sBase];
%     save(fName, 'BR_IMA_minLocs','BR_IMA_minDs','BR_IMA_shufs')
%     saveas(aaa,[saveDir,'BR_IMA_heatmap_ds'],'fig')
%     saveas(aaa,[saveDir,'BR_IMA_heatmap_ds'],'png')
%     saveas(aaa,[saveDir,'BR_IMA_heatmap_ds'],'svg')
%     saveas(aaz,[saveDir,'BR_IMA_heatmap_seqs'],'fig')
%     saveas(aaz,[saveDir,'BR_IMA_heatmap_seqs'],'png')
%     saveas(aaz,[saveDir,'BR_IMA_heatmap_seqs'],'svg')
% elseif pStruct.rampTypeFlag == 3
%     sBase = ['FR_IP_', num2str(permN),ciStr];
%     FR_IP_minLocs  = minDLocs;
%     FR_IP_minDs    = minDs;
%     FR_IP_shufs    = shuf;
%     fName = [saveDir, sBase];
%     save(fName, 'FR_IP_minLocs','FR_IP_minDs','FR_IP_shufs')
%     saveas(aaa,[saveDir,'FR_IP_heatmap_ds'],'fig')
%     saveas(aaa,[saveDir,'FR_IP_heatmap_ds'],'png')
%     saveas(aaa,[saveDir,'FR_IP_heatmap_ds'],'svg')
%     saveas(aaz,[saveDir,'FR_IP_heatmap_seqs'],'fig')
%     saveas(aaz,[saveDir,'FR_IP_heatmap_seqs'],'png')
%     saveas(aaz,[saveDir,'FR_IP_heatmap_seqs'],'svg')
% elseif pStruct.rampTypeFlag == 4
%     sBase = ['DR_IP_', num2str(permN),ciStr];
%     DR_IP_minLocs  = minDLocs;
%     DR_IP_minDs    = minDs;
%     DR_IP_shufs    = shuf;
%     fName = [saveDir, sBase];
%     save(fName, 'DR_IP_minLocs','DR_IP_minDs','DR_IP_shufs')
%     saveas(aaa,[saveDir,'DR_IP_heatmap_ds'],'fig')
%     saveas(aaa,[saveDir,'DR_IP_heatmap_ds'],'png')
%     saveas(aaa,[saveDir,'DR_IP_heatmap_ds'],'svg')
%     saveas(aaz,[saveDir,'DR_IP_heatmap_seqs'],'fig')
%     saveas(aaz,[saveDir,'DR_IP_heatmap_seqs'],'png')
%     saveas(aaz,[saveDir,'DR_IP_heatmap_seqs'],'svg')
% elseif pStruct.rampTypeFlag == 7
%     sBase = ['BR_IP_', num2str(permN),ciStr];
%     BR_IP_minLocs  = minDLocs;
%     BR_IP_minDs    = minDs;
%     BR_IP_shufs    = shuf;
%     fName = [saveDir, sBase];
%     save(fName, 'BR_IP_minLocs','BR_IP_minDs','BR_IP_shufs')
%     saveas(aaa,[saveDir,'BR_IP_heatmap_ds'],'fig')
%     saveas(aaa,[saveDir,'BR_IP_heatmap_ds'],'png')
%     saveas(aaa,[saveDir,'BR_IP_heatmap_ds'],'svg')
%     saveas(aaz,[saveDir,'BR_IP_heatmap_seqs'],'fig')
%     saveas(aaz,[saveDir,'BR_IP_heatmap_seqs'],'png')
%     saveas(aaz,[saveDir,'BR_IP_heatmap_seqs'],'svg')
% end
% 
%% Functions

function [minLocs] = minTimes(dsMat,pVector,nRamps)
dsMat = flipud(abs(dsMat'));
[~,dInds] = min(dsMat,[],1);
%     minBin = dsMat == dMins;
indMat = flipud(repmat(pVector,nRamps,1)');
%     minCol = minBin.*indMat;     %Removes all but earliest minimum
minLocs = indMat(dInds);
end



