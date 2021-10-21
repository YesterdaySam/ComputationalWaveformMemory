%%% Sample code for Attractor Model
% LKW 7/29/21 also see Hasselmo et al., J. Neurosci. 1995.
% Searches attractor outcome across dimensions of ramp degree and pulse
% length.
% Incorporates simple adaptation current using simulated intracellular Calcium dynamics
% Relies on simpleCA3Integrator.m

% clear
pStruct.wf_flag = 1;         %1 = FR IMA; 2 = DR IMA; 3 = FR IP; 4 = DR IP; 5 = BR IMA; 7 = BR IP
pStruct.simTypeFlag = 1;     %1 = Linear; 2 = Linear with Adaptation; 3 = Nonlinear; 4 = Nonlinear with Adaptation

pStruct.W       = 0.016;     %Weight of a to a (recurrency)
pStruct.H       = 0.06;      %Weight of h to a (feedback inhibition)
pStruct.Wh      = 0.0042;    %Weight of a to h (feedforward activation of IN)
pStruct.Hh      = 0.003;     %Inhibition self feedback
pStruct.tha     = 8;         %Threshold for a
pStruct.thh     = 8;         %Threshold for h
pStruct.thc     = 8;         %Threshold for voltage-dep ca-currents
pStruct.mu      = 0.01;      %Ca-dependent K-current
pStruct.gm      = 0.001;     %Gamma; voltage-dependent Ca-currents
pStruct.om      = 0.001;     %Omega; constant for diffusion of intracellular Ca
pStruct.eta     = 0.01;      %Decay constant
pStruct.Iexcit  = 0.1;       %Stim strength
pStruct.inDur   = 500;       %Stim duration

nPs     = 401;     %Number of steps of parameter
nXs     = 21;      %Number of ramp degrees
pLim    = 2000;    %Upper limit of dur search parameter
pFloor  = 1;       %Floor of dur search parameter space
pVect   = round(linspace(pFloor,pLim,nPs));    %Vector of dur parameter values
if mod(pStruct.wf_flag,2) == 1; rampCap = 1; else; rampCap = 0.5; end
rampPercs = linspace(0,rampCap,nXs);    %Vector of percentages
rampOuts = zeros(nXs,nPs);  %Mat of attractor outputs
squareOuts = zeros(1,nPs);  %for vector of square outputs

%% Run Space Search Simulation 

for j = 1:nXs
    for i = 1:nPs   %nPs or 399
%         [rampOuts(j,i),squareOuts(i)] = simpleCA3Integrator(pStruct,pVect(i),rampPercs(j));
        [rampOuts(j,i),squareOuts(i)] = simpleCA3Integrator_V2(pStruct,pVect(i),rampPercs(j));
%         [rampOuts(i),squareOuts(i)] = simpleCA3Integrator(pStruct,inDur,0.2);
    end
end

rampHeur = sum(round(rampOuts) > 0, 2)';    %Counts number of attractor states reached at each level of ramp
rampHeurNorm = rampHeur/nPs*100;

flipRampOuts = flipud(rampOuts');
set(0,'DefaultLineLineWidth',2)

%% Plot heatmap

%Attractor state image x-axis Ramp percentage, y-axis input duration
% figure(); set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.4, 0.28, 0.5]); aaa = gca;  %Optimized Figure size
figure(); set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.27, 0.5]); aaa = gca;    %Normalized Figure size
imagesc(flipRampOuts,[-10,33]);
aaa.YTick = 1:round(nPs-1)/5:nPs; aaa.XTick = 1:round(nXs-1)/5:nXs;
yticklabels(fliplr(0:pLim/5:pLim));
% xticklabels(0:rampCap*20:rampCap*100)     %For variable 50/100% ramp across stim classes
xticklabels(0:20:100)   %For straight 100% ramp across all sim classes
xlabel('Ramp Percentage');ylabel('Input Duration (ms)'); colormap gray; axis square;
set(gca,'FontSize',24,'fontname','times')

%% Perform Shuffles
permN = 1000;
nullShuf = ca3BootstrapAll(flipRampOuts,permN);
[durShuf,sdDurShuffle] = ca3BootstrapCol(flipRampOuts,permN);
% CIs
z = 1.96;   %confidence level 95% = 1.96, 99% = 2.58
upCIdur = durShuf + z*sdDurShuffle/sqrt(permN);
dnCIdur = durShuf - z*sdDurShuffle/sqrt(permN);

minLocs = minTimes(flipRampOuts);

figure(); set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.28, 0.42, 0.65]); aab = gca;
% figure(); set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.5, 0.8]);
% plot(rampPercs*100,rampHeurNorm,'k');                             %Plot actual duration sum
axis square; hold on;
plot(rampPercs*100,durShuf*100,'k',rampPercs*100,upCIdur*100,'k--')      %Plot Duration shuffle and upper CI
plot(rampPercs*100,dnCIdur*100,'k--')      %Plot lower CIs for legend purposes
xlabel('Ramp Percentage'); ylabel({'Likelihood of any input duration';'To reach Attractor'})
legend('Duration Shuffle','CI');
ylim([0 100]); aab.XTick = linspace(0,max(rampPercs*100),6); xticklabels(0:20:100)
set(gca,'FontSize',24,'fontname','times')

%% For generating ramp examples figures

rampsMat        = zeros(nXs,1000);
for i = 1:nXs
    inDur = pStruct.inDur;
    Iexcit = pStruct.Iexcit;
    squarea = Iexcit*inDur;     %Area under the square pulse
    rampLen         = round(pStruct.inDur*rampPercs(i));    %For basic 20% ramp do 0.2
    offset = 150+(i-1)*10;
    if pStruct.wf_flag == 1
        rampsMat(i,offset+1:offset+inDur)   = Iexcit;
        rampsMat(i,offset+1:offset+rampLen) = linspace(0,Iexcit,rampLen);   %Front Ramp
    elseif pStruct.wf_flag == 2
        rampsMat(i,offset+1:offset+inDur)   = Iexcit;
        rampsMat(i,offset+1:offset+rampLen) = linspace(0,Iexcit,rampLen);   %Front Ramp
        rampsMat(i,offset+1+inDur-rampLen:offset+inDur) = linspace(Iexcit,0,rampLen);   %Rear ramp
    elseif pStruct.wf_flag == 3
        IRamp = squarea/(rampLen/2 + inDur - rampLen);  %Calculate Single Ramp Current Max for constant AUC
        rampsMat(i,offset+1:offset+inDur)  = IRamp;
        rampsMat(i,offset+1:offset+rampLen)= linspace(0,IRamp,rampLen);   %Front Ramp
    elseif pStruct.wf_flag == 4
        IRamp = squarea/(inDur - rampLen);              %Calculate Double Ramp Current Max for constant AUC
        rampsMat(i,offset+1:offset+inDur)  = IRamp;
        rampsMat(i,offset+1:offset+rampLen)= linspace(0,IRamp,rampLen);   %Front Ramp
        rampsMat(i,offset+1+inDur-rampLen:offset+inDur) = linspace(IRamp,0,rampLen);   %Rear ramp
    elseif pStruct.wf_flag == 5
        rampsMat(i,offset+1:offset+inDur)   = Iexcit;
        rampsMat(i,offset+1+inDur-rampLen:offset+inDur) = linspace(Iexcit,0,rampLen);   %Rear ramp
    elseif pStruct.wf_flag == 7
        IRamp = squarea/(rampLen/2 + inDur - rampLen);  %Calculate Single Ramp Current Max for constant AUC   
        rampsMat(i,offset+1:offset+inDur)  = IRamp;
        rampsMat(i,offset+1+inDur-rampLen:offset+inDur) = linspace(IRamp,0,rampLen);   %Rear ramp
        
    end
end
rampsMat(nXs+1,offset+1:offset+inDur) = Iexcit;

c = gray(nXs+5);

figure(); set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.27, 0.36]); hold on;
for i = 1:nXs
    if i == 1
        plot(rampsMat(i,:)+(nXs-i+1)*0.005,'r')
    elseif i == 11
        plot(rampsMat(i,:)+(nXs-i+1)*0.005,'b')
    else
        plot(rampsMat(i,:)+(nXs-i+1)*0.005,'Color',c(nXs-i+1,:))
    end
end
plot(rampsMat(nXs+1,:)+(nXs-i+1)*0.005,'r--')
xlabel('Time (ms)'); ylabel('Input Strength')
ylim([0,max(max(rampsMat))+0.025])
set(gca,'FontSize',24,'fontname','times')

%% Combined Shuffle Demo Plot
% % Must load in previously saved vectors from shuffle
% patchXs = [1:21,linspace(21,1,21)]; %Vector of x coords;
% figure(); set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.28, 0.365, 0.65]); aac = gca;
% axis square; hold on;
% plot(durShufSFC,'r-');
% plot(durShufDRC,'b-');
% plot(durShufBRC,'c-');
% plot(durShufSFP,'r--');
% plot(durShufDRP,'b--');
% plot(durShufBRP,'c--');
% patch(patchXs,[durShufSFCDn,fliplr(durShufSFCUp)],'k','EdgeColor','none','FaceAlpha',0.1);
% patch(patchXs,[durShufDRCDn,fliplr(durShufDRCUp)],'k','EdgeColor','none','FaceAlpha',0.1);
% patch(patchXs,[durShufBRCDn,fliplr(durShufBRCUp)],'k','EdgeColor','none','FaceAlpha',0.1);
% patch(patchXs,[durShufSFPDn,fliplr(durShufSFPUp)],'k','EdgeColor','none','FaceAlpha',0.1);
% patch(patchXs,[durShufDRPDn,fliplr(durShufDRPUp)],'k','EdgeColor','none','FaceAlpha',0.1);
% patch(patchXs,[durShufBRPDn,fliplr(durShufBRPUp)],'k','EdgeColor','none','FaceAlpha',0.1);
% xlabel('Ramp Percentage');  %ylabel({'Likelihood of any input';'duration to reach Attractor'})
% ylabel('Attractor Likelihood');
% legend('FR IMA','DR IMA','BR IMA','FR IP','DR IP','BR IP','95% CI','FontSize',16,'Location','southwest');
% ylim([0 1]); xlim([1 21]); aac.XTick = linspace(0,21,6); xticklabels(0:20:100)
% set(gca,'FontSize',24,'fontname','times')
% 
% % Plot first attractor states must load in previously saved vector variable
% figure(); set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.28, 0.365, 0.65]); baa = gca;
% axis square; hold on;
% plot(minLocsSFC,'r-');
% plot(minLocsDRC,'b-');
% plot(minLocsBRC,'c-');
% plot(minLocsSFP,'r--');
% plot(minLocsDRP,'b--');
% plot(minLocsBRP,'c--')
% baa.YTick = 1:round(nPs-1)/5:nPs; baa.XTick = 1:round(nXs-1)/5:nXs;
% yticklabels(0:pLim/5:pLim);
% ylim([1 nPs]); xlim([1 nXs]);
% xticklabels(0:20:100)   %For straight 100% ramp across all sim classes
% xlabel('Ramp Percentage');ylabel('Input Duration (ms)'); colormap gray; axis square;
% legend('FR IMA','DR IMA','BR IMA','FR IP','DR IP','BR IP','FontSize',16,'Location','northwest');
% set(gca,'FontSize',24,'fontname','times')

%% Functions

function [shufSums] = ca3BootstrapAll(sampMat,iters)
nYs = size(sampMat,1);
nXs = size(sampMat,2);
attracVect = reshape(sampMat,nYs*nXs,1);
shufSums = zeros(1,nXs);
for i = 1:nXs
    shufDistro = datasample(attracVect,iters);
    tmpSum = sum(round(shufDistro) > 0);
    shufSums(i) = tmpSum/iters*100;
end
end

function [shufAve,shufStd] = ca3BootstrapCol(sampMat,iters)
shufDistro = datasample(sampMat,iters,1);   %Sample the columns iters number of times
shufStd = std(round(shufDistro) > 0);    %Count any nonzero attractor states
shufAve = mean(round(shufDistro) > 0);    %Normalize against iters and make percentage
end

function [minLocs] = minTimes(attMat)
[~,dInds] = max(flipud(round(attMat)),[],1);
tmpInds = dInds == 1;
minLocs = dInds;
minLocs(tmpInds) = [];
end
