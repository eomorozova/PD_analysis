
% A function to calculate and plot a phase response curve 

function PRC(V,I,Fs)

% inputs:
% V membrane potential of the neuron
% I injected current (containing positive pulses)
% sampling frequency

if nargin<3
    Fs=10000;
end

%% find a start of the positive pulses

Ithresh=0.5;
idx_pulse=find(diff(I>Ithresh)>0);

%% calculate phase response curves

thresh=-43; % threshold to find the beginning of each burst
Vcross=[]; delta_phase = []; stim_phase=[]; stim_time=[];

for i = 1:length(idx_pulse)
    
    % pick window around the pulse [t1,t2] seconds
    t1 = 2; t2 = 8;
    
    V1 = V(idx_pulse(i)-t1*Fs:idx_pulse(i)+t2*Fs);
    I1 = I(idx_pulse(i)-t1*Fs:idx_pulse(i)+t2*Fs);
    
    % find the beginning of each burst
    Vcross{i} = find(diff(smooth(V1,100)'>thresh)>0)./Fs; %(in seconds)
    
    ISI = diff(Vcross{i}); %interburst intervals
    
    % subtract t1 to recenter the traces relative to the pulse
    st = Vcross{i}-t1;
    
    idx=find(st>0); IdxAfter(i)=idx(1); % find the time of the fist burst after the pulse
    idx=find(st<0); IdxBefore(i)=idx(end); % last burst before the pulse
    
    P0(i)=median(ISI(end-5:end)); % period of the last 5 bursts after stimulus (unperturbed period)
    
    P1(i)=ISI(IdxAfter(i)-1); % new period after the pulse
    
    delta_phase(i)=(P1(i)-P0(i))/P0(i); % shift in phase
    
    stim_phase(i) = abs(2-Vcross{i}(IdxBefore(i)))/P0(i); % phase of the pulse
    stim_time(i) = abs(2-Vcross{i}(IdxBefore(i))); % time of the pulse in seconds
    
end

%% summary plot

clf

% raster plot aligned to the beginning of the burst before the pulse
%(each dot corresponds to the beginning of a burst)
subplot(2,3,1)

[~,idx]=sort(stim_time); % arrange the trials based on the time of the pulse

k=1;
for i = idx
    
    plot(stim_time(i),k,'.','markersize',10,'color','r'), hold on
    
    for j = 1:8
        plot(Vcross{i}(IdxBefore(i)+j-1)-Vcross{i}(IdxBefore(i)),k,'.','markersize',10,'color','k'), hold on
    end
    k=k+1;
    
end

% plot unperturbed burst period with vertical lines
dalim=get(gca,'YLim');
for ii=1:8
    h=line([median(P0)*ii median(P0)*ii],dalim,'Color','r','linewidth',1.5);
end
ylim([0 length(idx)]); xlim([-0.5 5])
ylabel('Trial #'); xlabel('Time, sec')
set(gca,'Fontsize',14);


% distribution of baseline periods
subplot(2,3,4)
histogram(P0,'facecolor','k')
xlabel('Period, sec')
rand_pase_shift=std(P0)/ median(P0);
text(0.8,60,['mean(P0)=',num2str(round(mean(P0),2)),' sec'],'fontsize',12)
text(0.8,50,['std(P0)=',num2str(round(std(P0),2)),' sec'],'fontsize',12)
text(0.8,40,['chance \delta \theta=',num2str(round(std(P0)/mean(P0),2))],'fontsize',12)
%xlim([0.5 1])
title('Distribution of PD neuron burst periods')
dalim=get(gca,'YLim');
h=line([median(P0) median(P0)],dalim,'Color',[1 0 0],'LineStyle','--','linewidth',1);
set(gca,'Fontsize',14)

% phase responce curve
subplot(2,3,2)
plot(stim_phase,delta_phase,'.','markersize',10,'color','k')
xlabel('Stimulus Phase'); ylabel('Phase shift')
title('PRC')
h = line([0.5 0.5], [-1 1]); set(h,'color','k','linestyle','--','linewidth',1), hold on
h = line([0 1], [0 0]); set(h,'color','k','linestyle','--','linewidth',1)
ylim([-1 1]); xlim([0 1])
set(gca,'Fontsize',14)


% plot a waveform normilized to the period
subplot(2,3,5)
i=3;
V_1=V(Vcross{i}(IdxBefore(i)+2)*Fs:Vcross{i}(IdxBefore(i)+3)*Fs);
x=1/length(V_1):1/length(V_1):1;

plot(x,smooth(V_1,10),'k','linewidth',1.5)

ylabel('V_M, mV'); xlabel('Stimulus Phase');
set(gca,'Fontsize',16)


% plot example traces with different pulse phases
for ii=1:3
    subplot(3,3,ii*3)
    i=idx(ii*50);
    V1=V(idx_pulse(i)-t1*Fs:idx_pulse(i)+5*Fs);
    I1=I(idx_pulse(i)-t1*Fs:idx_pulse(i)+5*Fs);
    
    x=1/Fs:1/Fs:1*length(V1)/Fs;
    plot(x,smooth(V1,100),'k','linewidth',1)
    hold on, plot(x,smooth(I1,100)-55,'b','linewidth',1)
    hold on, plot(Vcross{i},-45,'.','markersize',20)
    
    dalim=get(gca,'YLim');
    % beginning of the burst before the pulse
    line([Vcross{i}(IdxBefore(i)) Vcross{i}(IdxBefore(i))],dalim,'Color','c','LineStyle','-','Linewidth',1.5);
    % location of whether the burst would be without the stimulus
    line([Vcross{i}(IdxBefore(i))+P0(i) Vcross{i}(IdxBefore(i))+P0(i)],dalim,'Color','r','LineStyle','-','Linewidth',1.5);
    
    text(0.1,-60,['P1=',num2str(round(P1(i),2)),' s'],'fontsize',12)
    text(1,-60,['P0=',num2str(round(P0(i),2)),' s'],'fontsize',12)
    text(2.6,-60,['\Delta\theta=',num2str(round(delta_phase(i),2))],'fontsize',12)
    text(4,-60,['stim \theta=',num2str(round(stim_phase(i),2))],'fontsize',12)
    h=mmyplothorzline(-45); set(h,'color','k','linewidth',1)
    xlim([1 5]); ylim([-65 -24])
    ylabel('V_M, mV'); xlabel('Time, sec')
    set(gca,'Fontsize',14)
end

