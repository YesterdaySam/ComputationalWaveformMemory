%%% Sample Code for Learning Model Search Script
%Relies on ca3Seq_learn.m function, learn_CohenStats.m ca3learn_stats.m and
%the plot_CA3LearnAct.m functions
%Searches ripple prolongation parameter space in dimensions of ramp degree
%and input pulse duration
%Plots various metrics and statistics for evaluating network performance
%under varying waveform shape input
%LKW 4/24/21
%%%

% clearvars
load learnControl_wtBiasOff_symmetric_offAxis1_7.mat

pStruct.rampTypeFlag    = 7;    %1 = FRC; 2 = DRC; 3 = FRP; 4 = DRP; 6 = BRC; 7 = BRP
suppGraphFlag = 1;      %Plot extra graphs or not

%Setup neuron and weight paramters
pStruct.N               = 15;       %Nodes per region
N                       = 15;       %For dry run
pStruct.Ww              = 0.031;    %Weight strength pyr to pyr
pStruct.Hh              = 0.05;     %Weight strength IN to Pyr
pStruct.Wh              = 0.05;     %Weight strength pyr to IN
pStruct.HAuto           = 0.003;    %Inhibition self feedback
pStruct.tha             = 4*ones(1,N);
pStruct.thh             = 4*ones(1,N);
pStruct.eta1            = 0.01;     %Decay constant
pStruct.eta2            = 0.01;
pStruct.actThresh       = 10*ones(1,N);
pStruct.e1              = 0.01;     %Presynaptic decay
pStruct.e2              = 0.01;     %Postsynaptic decay
pStruct.e3              = 0.001;    %Learning rate

%Setup input parameters
pStruct.TLearn          = 2000;     %Time steps, must be even
pStruct.TTest           = 1500;
pStruct.cueN            = 1;
pStruct.Iexcit1         = 0.5;      %Strength of learn pulses (Iso current)
pStruct.Iexcit2         = 1;     	%Strength of test cue pulse
pStruct.testDur         = 20;       %Test phase cue duration
pStruct.learnDur        = 80;       %Learn phase pulse duration
pStruct.learnOverlap    = 40;       %Learn phase pulse overlap
pStruct.rampPerc        = 0.50;
pStruct.rampLen1        = round(pStruct.learnDur*pStruct.rampPerc);  %Duration ramp length
pStruct.rampLen2        = round(pStruct.learnDur*pStruct.rampPerc);
pStruct.onsetDelay      = 50;        %Wait time to ripple start from sim start
pStruct.achLearn        = 0.1;       %Acetylcholine state; 0 = high Ach; 1 = low Ach
pStruct.achTest         = 1;

W                           = zeros(N,N,pStruct.TLearn);           %Init wt mat
AH                          = zeros(N,N);           %Wt mat of a to h (feedforward activation of IN)
H                           = zeros(N,N);           %Weight of h to a (feedback inhibition)
for i = 1:N         %Build weight mats
    H(i,i)  = pStruct.Hh;   %Direct feedback IN to Pyr
    AH(i,i) = pStruct.Wh;   %Direct excitation pyr to IN
end

pStruct.W = W; pStruct.H = H; pStruct.AH = AH;  %Add wtmats to pStruct

%Set up search variables and parameters
pLim                    = pStruct.learnDur;      %250 for inDur2
pFloor                  = 0;
nPs                     = pLim - pFloor + 1; %61;
nRamps                  = 21;
pVect                   = linspace(pFloor,pLim,nPs);    %Vector of parameter values
if mod(pStruct.rampTypeFlag,2) == 1; rampLim = 1; else; rampLim = 0.5; end %Set ramp limit to 1 or 0.5 depending on flag
rampPercs = linspace(0,rampLim,nRamps);

outputs.diffs       = cell(nRamps,nPs);
outputs.INdiffs     = cell(nRamps,nPs);
outputs.testPks     = cell(nRamps,nPs);
outputs.testLocs    = cell(nRamps,nPs);
outputs.Ws          = cell(nRamps,nPs);     %Storage of final Wt. Mat snapshot for each simulation
outputs.Ds          = zeros(nRamps,nPs);

%% Gut check single above parameters
gutCheckFlag = 0;
if gutCheckFlag == 1
    [actCell,hactCell,inCell]   = ca3Seq_learn(pStruct);    %Using ramp defined above
    aLearn = actCell{1}; aTest  = actCell{2};
    hLearn = hactCell{1};hTest  = hactCell{2};
    ALearn = inCell{1};  ATest  = inCell{2};
    % dGut = ripExtend_CohenStats(actCell,hactCell,pStruct);
    
    % Calculate Peaks
    pksLearn = zeros(N,1); locsLearn = zeros(N,1);
    pksTest = zeros(N,1); locsTest = zeros(N,1);
    
    for  i = 1:N
        [pksTmp,locsTmp] = findpeaks(aLearn(i,:));
        if ~isempty(locsTmp); pksLearn(i) = pksTmp(1); locsLearn(i) = locsTmp(1); end
        [pksTmp,locsTmp] = findpeaks(aTest(i,:));
        if ~isempty(locsTmp); pksTest(i) = pksTmp(1); locsTest(i) = locsTmp(1); end
    end
    pksTest = pksTest(pksTest>0); locsTest = locsTest(locsTest>0);
    locsCell = {locsLearn,locsTest};
    pksCell = {pksLearn,pksTest};
    
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
    
    [thaPlot,thax1,thax2] = plotCA3Ripples_learn(actCell,hactCell,inCell,locsCell,pksCell,pStruct);
    set(thaPlot,'defaultLegendAutoUpdate','off');
    axCell = {thax1,thax2};
    thaLength = 20;
    for i = 1:2
        axes(axCell{i});
        tmpTTT = tttCell{i};
        plot([0,999],[10,10],'k--','LineWidth',1)
        for j = 1:length(tmpTTT)
            %         plot([tmpTTT(j)-thaLength,tmpTTT(j)+thaLength],[tha(j),tha(j)],'k--','LineWidth',1)
            scatter(tmpTTT(j),actThresh(j),'k^','filled')
        end
    end
    % [rampINplot,rampINax1,rampINax2] = plotCA3INs(actCell,hactCell,inCell,pStruct.cueN,1);
    % [controlINplot,controlINax1,controlINax2] = plotCA3INs(actCell,hactCell,inCell,pStruct.cueN,3);
end

%% Run Block Parameter Search
for i = 1:nRamps
    pStruct.tmpRamp = rampPercs(i);                             %Define ramp percentage
    pStruct.rampLen1= round(pStruct.tmpRamp*pStruct.learnDur);   %Define ramp length
    for j = 1:nPs
        pStruct.learnOverlap = round(pVect(j));                         %Parameter for change
        %Run function and calculate stats
        [actTmp,hactTmp,inTmp,outputs.Ws{i,j}] = ca3Seq_learn(pStruct);
        if max(actTmp{2},[],2) > 100
            outputs.testN(i,j) = 0;
            outputs.Ds(i,j) = NaN;
        else
            [outputs.testPks{i,j},outputs.testLocs{i,j},outputs.diffs{i,j},outputs.INDiffs{i,j}] = ca3learn_stats(pStruct,actTmp,hactTmp);
            outputs.testN(i,j) = numel(outputs.testPks{i,j});
            [outputs.Ds(i,j),~,~] = learn_CohenStats({actTmp{2},aControl},NaN,pStruct);
        end
        %Calculate other outputs numbers
        if ~isempty(outputs.diffs{i,j})
            outputs.meanIThI(i,j) = mean(outputs.diffs{i,j});
%             outputs.meanINIThI(i,j) = mean(outputs.INDiffs{i,j});
%             outputs.sdIThI(i,j) = std(outputs.diffs{i,j});
%             outputs.sdINIThI(i,j) = std(outputs.INDiffs{i,j});
        else
            outputs.meanIThI(i,j) = 0;
%             outputs.meanINIThI(i,j) = 0;
        end
    end
end

%% Map of Effect Size
% Clean ds output
outputs.nanMaxDs = outputs.Ds;
% outputs.cleanDs = outputs.Ds .* nan((outputs.Ds >= 0));
tmpMax = max(max(outputs.nanMaxDs));
% outputs.CleanDs = nan(outputs.CleanDs
outputs.nanMaxDs(isnan(outputs.nanMaxDs)) = tmpMax;  %Convert any nans to max effect size for plotting purposes
outputs.nanZedDs = outputs.Ds;
outputs.nanZedDs(isnan(outputs.nanZedDs)) = 0;      %Convert any nans to 0s for shuffle purposes

figure; % set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.35, 0.65]);
aaz = gca;
imagesc(flipud(abs(outputs.nanMaxDs')),[0,max(max(outputs.nanMaxDs))]);
aaz.YTick = 1:round((nPs-1)/5):nPs;   %Reversed 0 at top, max at bottom
aaz.XTick = 1:round(nRamps-1)/5:nRamps;
yticklabels(linspace(100,0,6))
xticklabels(0:20:100);
% xlabel('Ramp Percentage'); ylabel('OL Duration (ms)'); 
colormap hot; axis square; hcb = colorbar;
% hcb.Label.String = "Cohen's d";
set(aaz,'FontSize',24,'fontname','times')

%% Find lowest OL duration of full replay
minInds = bestSeqN(flipud(outputs.testN'),pVect,nRamps);
minDs   = min(flipud(abs(outputs.Ds')));

%% Shuffles
permN = 500;
z = 1.96;   %confidence level 95% = 1.96, 99% = 2.58
pyrIThIShuf = bootCols(flipud(outputs.meanIThI'),permN);
meanPyrIThIShuf = mean(pyrIThIShuf);
sdPyrIThIShuf = std(pyrIThIShuf);
actNShuf = bootCols(flipud(outputs.testN'),permN);
meanActNShuf = mean(actNShuf);
sdMeanActNShuf = std(actNShuf);
dShuf = bootColDs(flipud(abs(outputs.nanZedDs')),permN);
meanDShuf = mean(dShuf);
sdMeanDShuf = std(dShuf);
% CIs
upCIPyr = meanPyrIThIShuf + z*sdPyrIThIShuf/sqrt(permN);
dnCIPyr = meanPyrIThIShuf - z*sdPyrIThIShuf/sqrt(permN);
upCIN = meanActNShuf + z*sdMeanActNShuf/sqrt(permN);
dnCIN = meanActNShuf - z*sdMeanActNShuf/sqrt(permN);
upCId = meanDShuf + z*sdMeanDShuf/sqrt(permN);
dnCId = meanDShuf - z*sdMeanDShuf/sqrt(permN);

% figure(); % set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.35, 0.65]); 
% % figure(); set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.5, 0.8]);
% subplot(2,2,1); aab = gca;
% tmpMax = max(max(abs(outputs.meanIThI)));
% axis square; hold on;
% plot(rampPercs*100,meanPyrIThIShuf,'k',rampPercs*100,upCIPyr,'k--')      %Plot Duration shuffle and upper CI
% plot(rampPercs*100,dnCIPyr,'k--')
% xlabel('Ramp Percentage'); ylabel("Mean IThI (ms)")
% % legend('OL Shuffle','99% CI');
% ylim([0 tmpMax]); aab.XTick = linspace(0,max(rampPercs*100),6); xticklabels(0:20:100)
% set(gca,'FontSize',24,'fontname','times')
% subplot(2,2,2); aab = gca;
% tmpMax = max(max(abs(outputs.meanIThI)));
% axis square; hold on;
% plot(rampPercs*100,meanDShuf,'k',rampPercs*100,upCId,'k--')      %Plot Duration shuffle and upper CI
% plot(rampPercs*100,dnCId,'k--')
% xlabel('Ramp Percentage'); ylabel("Mean Effect Size (ms)")
% ylim([0 3]);
% aab.XTick = linspace(0,max(rampPercs*100),6); xticklabels(0:20:100);
% set(gca,'FontSize',24,'fontname','times')
% subplot(2,2,3); aac = gca;
% tmpMax = max(max(abs(outputs.testN)));
% axis square; hold on;
% plot(rampPercs*100,meanActNShuf,'k',rampPercs*100,upCIN,'k--')      %Plot Duration shuffle and upper CI
% plot(rampPercs*100,dnCIN,'k--')
% xlabel('Ramp Percentage'); ylabel("Mean Suprathreshold Nodes")
% legend('OL Shuffle','99% CI');
% ylim([0 tmpMax]); aac.XTick = linspace(0,max(rampPercs*100),6); xticklabels(0:20:100)
% set(gca,'FontSize',24,'fontname','times')
% subplot(2,2,4); aad = gca;
% axis square; hold on;
% plot(rampPercs*100, minInds,'k-')
% xlabel('Ramp Percentage'); ylabel('OL dur of Best Seq Length');
% ylim([-10 80]);
% set(gca,'FontSize',24,'fontname','times')

%% Plotting
set(0,'DefaultLineLineWidth',2)

%Plot map of sequence length
seqPlot = figure;  aaa = gca;
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.5, 0.8]);
imagesc(flipud(outputs.testN'),[0,max(max(outputs.testN))]);
aaa.YTick = 1:round((nPs-1)/5):nPs; 
aaa.XTick = 1:round(nRamps-1)/5:nRamps;
% yticklabels(linspace(pLim,pFloor,6));
% xticklabels(0:rampLim*20:rampLim*100);      %For Ramp Percentage Label
yticklabels(linspace(100,0,6))
xticklabels(0:20:100);
% xlabel('Ramp Percentage'); ylabel('Learn Overlap (ms)'); 
colormap jet; axis square; hcb = colorbar;
% if mod(pStruct.rampTypeFlag,2) == 1; title("Sequence Replay Length; Forward Ramp Learn");
% else; title("Sequence Replay Length; Double Ramp Learn"); end
% hcb.Label.String = "# Suprathreshold Test Nodes";
set(aaa,'FontSize',24,'fontname','times')

% %Box plot summaries collapsed row axis
% figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.5, 0.8]); aab = gca;
% % boxplot(abs(outputs.ds'),'PlotStyle','compact');
% boxplot(abs(outputs.ds'));
% xticklabels(0:rampLim/(nRamps-1)*100:rampLim*100); xtickangle(90)
% ylabel("Cohen's d"); xlabel("Ramp Percentage");
% if mod(pStruct.rampTypeFlag,2) == 1; title("Forward Ramp Collapse Input Duration");
% else; title("Double Ramp Collapse Input Duration"); end
% set(aab,'FontSize',24,'fontname','times')
    
% % Plot map of pyramidal mean IthI
% pyrPlot = figure;  aaa = gca;
% % set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.5, 0.8]);
% imagesc(flipud(outputs.meanIThI'),[0,max(max(outputs.meanIThI))]);
% aaa.YTick = 1:round((nPs-1)/5):nPs; 
% aaa.XTick = 1:round(nRamps-1)/5:nRamps;
% % yticklabels(linspace(pLim,pFloor,6));
% % xticklabels(0:rampLim*20:rampLim*100);      %For Ramp Percentage Label
% yticklabels(linspace(100,0,6))
% xticklabels(0:20:100);
% % xlabel('Ramp Percentage'); ylabel('Learn Overlap (ms)'); 
% colormap jet; axis square; hcb = colorbar;
% % if mod(pStruct.rampTypeFlag,2) == 1; title("Replay Mean Pyr IThI; Forward Ramp Learn");
% % else; title("Replay Mean Pyr IThI; Double Ramp Learn"); end
% % hcb.Label.String = "Mean Pyr IThI";
% set(aaa,'FontSize',24,'fontname','times')

% % Plot map of IN mean IthI
% INPlot = figure; aaa = gca;
% % set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.5, 0.8]);
% imagesc(flipud(outputs.meanINIThI'),[0,max(max(outputs.meanINIThI))]);
% aaa.YTick = 1:round((nPs-1)/5):nPs; 
% aaa.XTick = 1:round(nRamps-1)/5:nRamps;
% yticklabels(linspace(pLim,pFloor,6));
% xticklabels(0:rampLim*20:rampLim*100);      %For Ramp Percentage Label
% % xlabel('Ramp Percentage'); 
% % ylabel('Learn Overlap (ms)'); 
% colormap jet; axis square; hcb = colorbar;
% % if mod(pStruct.rampTypeFlag,2) == 1; title("Replay Mean IN IThI; Forward Ramp Learn");
% % else; title("Replay Mean IN IThI; Double Ramp Learn"); end
% % hcb.Label.String = "Mean IN IThI";
% set(aaa,'FontSize',24,'fontname','times')

% if pStruct.rampTypeFlag == 1; sname = '\forwardRamp_';
% elseif pStruct.rampTypeFlag == 2; sname = '\doubleRamp_';
% elseif pStruct.rampTypeFlag == 3; sname = '\isoForwardRamp_';
% else; sname = '\isoDoubleRamp_'; 
% end
% dirname = 'C:\Users\Kvothe\Documents\Research\Code\CA3 Region Code\CA3Net Versions\Learning Code\Figures\Weight Shunt Search';

% saveas(seqPlot,strcat(dirname,sname,'searchLearnOverlap_learnDur80_learnAch_0_10_testN'))
% saveas(seqPlot,strcat(dirname,sname,'searchOL_learnDur80_learnAch_0_10_testN.png'))
% saveas(pyrPlot,strcat(dirname,sname,'searchLearnOverlap_learnDur80_learnAch_0_10_testIThI'))
% saveas(pyrPlot,strcat(dirname,sname,'searchOL_learnDur80_learnAch_0_10_testIThI.png'))
% saveas(INPlot,strcat(dirname,sname,'searchLearnOverlap_learnDur80_learnAch_0_10_testINIThI'))
% saveas(INPlot,strcat(dirname,sname,'searchOL_learnDur80_learnAch_0_10_testINIThI.png'))

%% Supplementary Graphs

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

% % Afferent 3D plot code
% sqGhost = zeros(1,pStruct.TLearn);
% sqGhost(pStruct.onsetDelay+1:pStruct.onsetDelay+pStruct.learnDur) = pStruct.Iexcit1;
% figure(); set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.28, 0.35, 0.45]); gaa = gca;
% w3 = waterfall(1:pStruct.TLearn,1:N,ALearn); hold on;
% w3.LineWidth = 3;
% w3Sq = waterfall(1:pStruct.TLearn,1,sqGhost);
% w3Sq.LineStyle = '--'; w3Sq.EdgeColor = 'r'; w3Sq.LineWidth = 2;
% gaa.YTick = 1:7:15;
% zlim([0 1]); xlim([0 1200])
% xlabel('Time'); ylabel('Node','position',[-148.7513525348877,6.348030863114843,-0.099984467105344]); zlabel('Amplitude')
% view([-13.2608328513592 45.9028169014085]);
% set(gca,'FontSize',24,'fontname','times')

if suppGraphFlag == 1
    load('C:\Users\Kvothe\Documents\Research\Code\CA3 Region Code\Pattern Learning Model\learnVars\DRC_sumStats_wtBiasOff_offAxis1_7.mat')
    load('C:\Users\Kvothe\Documents\Research\Code\CA3 Region Code\Pattern Learning Model\learnVars\DRP_sumStats_wtBiasOff_offAxis1_7.mat')
    load('C:\Users\Kvothe\Documents\Research\Code\CA3 Region Code\Pattern Learning Model\learnVars\FRC_sumStats_wtBiasOff_offAxis1_7.mat')
    load('C:\Users\Kvothe\Documents\Research\Code\CA3 Region Code\Pattern Learning Model\learnVars\FRP_sumStats_wtBiasOff_offAxis1_7.mat')
    load('C:\Users\Kvothe\Documents\Research\Code\CA3 Region Code\Pattern Learning Model\learnVars\BRC_sumStats_wtBiasOff_offAxis1_7.mat')
    load('C:\Users\Kvothe\Documents\Research\Code\CA3 Region Code\Pattern Learning Model\learnVars\BRP_sumStats_wtBiasOff_offAxis1_7.mat')
    
    % Combined mean IThI x Ramp shuffle - must load in previous data
    patchXs = [1:nRamps,linspace(nRamps,1,nRamps)]; %Vector of x coords;
    figure(); set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.28, 0.365, 0.65]); aac = gca;
    axis square; hold on;
    plot(FRC_meanIThI,'r-');
    plot(DRC_meanIThI,'b-');
    plot(BRC_meanIThI,'c-');
    plot(FRP_meanIThI,'r--');
    plot(DRP_meanIThI,'b--');
    plot(BRP_meanIThI,'c--');
    patch(patchXs,[FRC_meanIThIdn,fliplr(FRC_meanIThIup)],'k','EdgeColor','none','FaceAlpha',0.1);
    patch(patchXs,[DRC_meanIThIdn,fliplr(DRC_meanIThIup)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
    patch(patchXs,[FRP_meanIThIdn,fliplr(FRP_meanIThIup)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
    patch(patchXs,[DRP_meanIThIdn,fliplr(DRP_meanIThIup)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
    patch(patchXs,[BRC_meanIThIdn,fliplr(BRC_meanIThIup)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
    patch(patchXs,[BRP_meanIThIdn,fliplr(BRP_meanIThIup)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
    xlabel('Ramp Percentage'); ylabel('Mean IThI')
    % legend('FR Equal Max','DR Equal Max','FR Equal Power','DR Equal Power','95% CI','FontSize',16,'location','northwest');
    legend('FR IMA','DR IMA','BR IMA','FR IP','DR IP','BR IP','95% CI','FontSize',16,'location','northwest');
    ylim([25 35]);xlim([1 nRamps]); aac.XTick = linspace(0,nRamps,6); xticklabels(0:20:100)
    set(gca,'FontSize',24,'fontname','times')
    
    % Combined min D x Ramp shuffle - must load in previous data
    figure(); set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.28, 0.365, 0.65]); aac = gca;
    axis square; hold on;
    plot(FRC_minDs,'r-');
    plot(DRC_minDs,'b-');
    plot(BRC_minDs,'c-');
    plot(FRP_minDs,'r--');
    plot(DRP_minDs,'b--');
    plot(BRP_minDs,'c--');
    xlabel('Ramp Percentage'); ylabel('Min Effect Size')
    % legend('FR Equal Max','DR Equal Max','FR Equal Power','DR Equal Power','FontSize',16,'location','northwest');
    legend('FR IMA','DR IMA','BR IMA','FR IP','DR IP','BR IP','FontSize',16,'location','northwest');
    ylim([0 1]); xlim([1 nRamps]); aac.XTick = linspace(0,nRamps,6); xticklabels(0:20:100)
    set(gca,'FontSize',24,'fontname','times')
    
    % Combined mean D x Ramp shuffle - must load in previous data
    patchXs = [1:nRamps,linspace(nRamps,1,nRamps)]; %Vector of x coords;
    figure(); set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.28, 0.365, 0.65]); aac = gca;
    axis square; hold on;
    plot(FRC_meanD,'r-');
    plot(DRC_meanD,'b-');
    plot(BRC_meanD,'c-');
    plot(FRP_meanD,'r--');
    plot(DRP_meanD,'b--');
    plot(BRP_meanD,'c--');
    patch(patchXs,[FRC_meanDdn,fliplr(FRC_meanDup)],'k','EdgeColor','none','FaceAlpha',0.1);
    patch(patchXs,[DRC_meanDdn,fliplr(DRC_meanDup)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
    patch(patchXs,[FRP_meanDdn,fliplr(FRP_meanDup)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
    patch(patchXs,[DRP_meanDdn,fliplr(DRP_meanDup)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
    patch(patchXs,[BRC_meanDdn,fliplr(BRC_meanDup)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
    patch(patchXs,[BRP_meanDdn,fliplr(BRP_meanDup)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
    xlabel('Ramp Percentage'); ylabel('Mean Effect Size')
    % legend('FR Equal Max','DR Equal Max','FR Equal Power','DR Equal Power','95% CI','FontSize',16,'location','northwest');
    legend('FR IMA','DR IMA','BR IMA','FR IP','DR IP','BR IP','95% CI','FontSize',16,'location','northwest');
    ylim([0 3]); xlim([1 nRamps]); aac.XTick = linspace(0,nRamps,6); xticklabels(0:20:100)
    set(gca,'FontSize',24,'fontname','times')
    
    % Combined mean Seq Length x Ramp shuffle - must load in previous data
    patchXs = [1:nRamps,linspace(nRamps,1,nRamps)]; %Vector of x coords;
    figure(); set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.28, 0.365, 0.65]); aac = gca;
    axis square; hold on;
    plot(FRC_meanN,'r-');
    plot(DRC_meanN,'b-');
    plot(BRC_meanN,'c-');
    plot(FRP_meanN,'r--');
    plot(DRP_meanN,'b--');
    plot(BRP_meanN,'c--');
    patch(patchXs,[FRC_meanNdn,fliplr(FRC_meanNup)],'k','EdgeColor','none','FaceAlpha',0.1);
    patch(patchXs,[DRC_meanNdn,fliplr(DRC_meanNup)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
    patch(patchXs,[FRP_meanNdn,fliplr(FRP_meanNup)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
    patch(patchXs,[DRP_meanNdn,fliplr(DRP_meanNup)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
    patch(patchXs,[BRC_meanNdn,fliplr(BRC_meanNup)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
    patch(patchXs,[BRP_meanNdn,fliplr(BRP_meanNup)],'k','EdgeColor','none','FaceAlpha',0.1,'HandleVisibility','off');
    xlabel('Ramp Percentage'); ylabel('Mean Sequence Length')
    % legend('FR Equal Max','DR Equal Max','FR Equal Power','DR Equal Power','95% CI','FontSize',16,'location','northwest');
    legend('FR IMA','DR IMA','BR IMA','FR IP','DR IP','BR IP','95% CI','FontSize',16,'location','northwest');
    ylim([0 15]); xlim([1 nRamps]); aac.XTick = linspace(0,nRamps,6); xticklabels(0:20:100)
    set(gca,'FontSize',24,'fontname','times')
end
%% Functions

function [shufDistro] = bootCols(sampMat,iters)
%Equal weight for any nonzero, non 1 value, 0 weight for 0s,1s
for i = 1:size(sampMat,2)   %For number of columns
    tmpWts = sampMat(:,i) > 1;
    tmpWts = ones(size(sampMat,1),1).*tmpWts;   %Set wt 0 for any 1 or 0
    shufDistro(:,i) = datasample(sampMat(:,i),iters,1,'Weights',tmpWts);
end
end

function [shufDistro] = bootColDs(sampMat,iters)
%Equal weight for any nonzero, non 1 value, 0 weight for 0s,1s
for i = 1:size(sampMat,2)   %For number of columns
    tmpWts = sampMat(:,i) > 0;
    tmpWts = ones(size(sampMat,1),1).*tmpWts;   %Set wt 0 for any 1 or 0
    shufDistro(:,i) = datasample(sampMat(:,i),iters,1,'Weights',tmpWts);
end
end

function [bestNLocs] = bestSeqN(nMat,pVector,nRamps)
% Looks for first index of best seq length == N
indMat = flipud(repmat(pVector,nRamps,1)');

bestNInds = zeros(1,nRamps);
for i = 1:nRamps
    bestNInds(i) = find(nMat(:,i) == 15,1,'last');
    if bestNInds(i) == 0
        bestNInds(i) = find(nMat(:,i) == max(nMat(:,i)),1,'first');
    end
end
bestNLocs = indMat(bestNInds);
end


