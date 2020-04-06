function out =  plotOutputsStrokeMotionNoMotion()
%for motion=0:2
%imagesToAnalyze={'ROI817.P3.fig','ROI910.P8.fig','ROI817.P8.fig','ROI910.P9.fig'};
figureName='ROI910.p4.fig';
%experimentNamesToAnalyze={'8.17.P3','9.10.P8','8.17.P8','9.10.P9'};
experimentName='9.10.P4';

filename='output';

%analyze motion, no motion or all
%motion=0; % largest no motion section
%motion=1; % largest motion section
%motion=2; % whole movie
motion=1;


correlationType='xcorr'; %xcorr,corr,FCA
folderName='strokeData';

linkageParam=1.2;
percentageSyncRequired=0.5;
windowSize=1;


%get all xyposition from fig
openfig(fullfile(folderName,figureName));
allCellLables=findobj(gca,'Type','text');
allCellLables=allCellLables(2:end);
xyPositions=zeros(length(allCellLables),2);
for i=1:size(xyPositions,1)
    xyPositions(i,:)=allCellLables(i).Position(1:2);
end
xyPositions=xyPositions(end:-1:1,:);
xyPositions=xyPositions-20;
%remove excluded cells
excludedCellLabel=fullfile(folderName,strcat(experimentName,'.Excluded.mat'));
excluded_cells=load(excludedCellLabel);
excluded_cells=excluded_cells.excluded_cells;

try
    spikeCourseLabel=fullfile(folderName,strcat(experimentName,'.Deconvdata.mat'));
    spikeCourse=load(spikeCourseLabel);

catch
    spikeCourseLabel=fullfile(folderName,strcat(experimentName,'.Deconv.mat'));
    spikeCourse=load(spikeCourseLabel);

end
spikeCourse=spikeCourse.DeconvMat;

try
    motionLabel=fullfile(folderName,strcat(experimentName,'.Motion.mat'));
    motionData=load(motionLabel);
catch
    motionLabel=fullfile(folderName,strcat(experimentName,'.Motion1.mat'));
    motionData=load(motionLabel);
end


try
    motionData=motionData.test;
catch
    motionData=motionData.motion(:,1);
end

motionRegions=bwconncomp(motionData);
numPixels=cellfun(@numel,motionRegions.PixelIdxList);
[biggestMotion, idxOfMotion] = max(numPixels);
motionFrameList=motionRegions.PixelIdxList{idxOfMotion};

noMotionData=~motionData;
noMotionRegions=bwconncomp(noMotionData);
numPixels=cellfun(@numel,noMotionRegions.PixelIdxList);
[biggestNoMotion, idxOfNoMotion] = max(numPixels);
noMotionFrameList=noMotionRegions.PixelIdxList{idxOfNoMotion};

totalTime=min(biggestMotion,biggestNoMotion);
motionFrameList=motionFrameList(1:totalTime);
noMotionFrameList=noMotionFrameList(1:totalTime);
oldexperimentName=experimentName;
if(motion==1)
    spikeCourse=spikeCourse(motionFrameList,:);
    experimentName=strcat(experimentName,'-motion-');
elseif(motion==0)
    spikeCourse=spikeCourse(noMotionFrameList,:);
    experimentName=strcat(experimentName,'-noMotion-');
end


%added to remove single spiking neurons from xcorr
if(strcmp(correlationType,'xcorr'))
    for i=1:size(spikeCourse,2)
        thisSpikeCourse=(spikeCourse(:,i));
        peaksAbove=bwconncomp(thisSpikeCourse);
        spikeCount=peaksAbove.NumObjects;
        if(spikeCount<2)
            excluded_cells(end+1)=i;
        end
    end
end

neuronIDs=[1:size(xyPositions,1)];
exciteNeurons=zeros(size(neuronIDs));
try
    exciteIDLabels=fullfile(folderName,strcat(oldexperimentName,'.Pyramidals.mat'));
    exciteIDs=load(exciteIDLabels);
catch
    exciteIDLabels=fullfile(folderName,strcat(oldexperimentName,'.Excitatory.mat'));
    exciteIDs=load(exciteIDLabels);
end
exciteIDs=exciteIDs.ind_pyramidal;
exciteNeurons(exciteIDs)=1;

neuronIDs(excluded_cells)=[];
spikeCourse(:,excluded_cells)=[];
xyPositions(excluded_cells,:)=[];
exciteNeurons(excluded_cells)=[];
validNeuronIds=neuronIDs;

%new figure with neuron positions
figure;

hold on;
for i=1:length(neuronIDs)
    thisNeuronID=neuronIDs(i);
    if(any(exciteIDs==thisNeuronID))
        color=[0,1,0];
    else
        color=[1,0,0];
    end
    plot(xyPositions(i,1),xyPositions(i,2),'o','color',color,'MarkerSize',15,'linewidth',5);
    text(xyPositions(i,1)+10,xyPositions(i,2),num2str(thisNeuronID),'color',color);
end
set(gca,'color','k');

numberOfNeurons=length(validNeuronIds);
timeCourseLength=size(spikeCourse,1);
allRawTraces=zeros(numberOfNeurons,timeCourseLength);
allRawTracesNormalized=zeros(numberOfNeurons,timeCourseLength);

meanInterPeakTimePerNeuron=zeros(numberOfNeurons,1);
allRawInterpeakTimes={};
allInterpeakTimesWithoutClusters={};
arrayOfBustines=zeros(numberOfNeurons,1);
arrayOfMemoryWithCluster=zeros(numberOfNeurons,1);
arrayOfMemoryWithoutCluster=zeros(numberOfNeurons,1);
for i=1:length(validNeuronIds)
    thisID=validNeuronIds(i);
    allRawTraces(i,:)=spikeCourse(:,i)';
    allRawTracesNormalized(i,:)=allRawTraces(i,:)./max(allRawTraces(i,:));
    thisSpikeCourse=spikeCourse(:,i)';
    thisNeuronPeakTimes=[];
    lastPeak=0;
    for k=1:length(thisSpikeCourse)
        if(thisSpikeCourse(k)>0)
            if(lastPeak~=0)
                thisNeuronPeakTimes=[thisNeuronPeakTimes,k-lastPeak];
            end
            lastPeak=k;
        end
    end
    if(length(thisNeuronPeakTimes)>2)
        allRawInterpeakTimes{i}=thisNeuronPeakTimes;
        allInterpeakTimesWithoutClusters{i}=thisNeuronPeakTimes(thisNeuronPeakTimes>1);
        meanInterPeakTimeThisNeuron=mean(thisNeuronPeakTimes);
        meanInterPeakTimePerNeuron(i)=meanInterPeakTimeThisNeuron;
        stdInterPeakTimeThisNeuron=std(thisNeuronPeakTimes);
        arrayOfBustines(i)=(stdInterPeakTimeThisNeuron-meanInterPeakTimeThisNeuron)/...
            (stdInterPeakTimeThisNeuron+meanInterPeakTimeThisNeuron);
        [arrayOfMemoryWithCluster(i),~]= corr(thisNeuronPeakTimes(1:end-1)',thisNeuronPeakTimes(2:end)','type', 'Spearman');
        if(length(allInterpeakTimesWithoutClusters{i})>2)
            [arrayOfMemoryWithoutCluster(i),~]= corr(allInterpeakTimesWithoutClusters{i}(1:end-1)',allInterpeakTimesWithoutClusters{i}(2:end)','type', 'Spearman');
        end

    else
        arrayOfBustines(i)=NaN;
        arrayOfMemoryWithCluster(i)=NaN;
        arrayOfMemoryWithoutCluster(i)=NaN;
    end

end

%%%uncomment to make zstack image

%%%%%% calculate the overall interpeak times for all neurons 

%loop through the time course 
populationSpikes=zeros(size(spikeCourse,1),1);
populationInterpeakTimes=[];
lastSpikeTime=0;
for i=1:size(spikeCourse,1)
    spikeAtTimepoint=0;
    for j=1:size(spikeCourse,2)
        if(spikeCourse(i,j)>0)
            spikeAtTimepoint=1;
        end
    end
    populationSpikes(i)=spikeAtTimepoint;
    if(spikeAtTimepoint==1)
        timeSinceLastSpike=i-lastSpikeTime;
        populationInterpeakTimes(end+1)=timeSinceLastSpike;
        lastSpikeTime=i;
    end
end

meanInterPeakTimePopulation=mean(populationInterpeakTimes);
stdInterPeakTimePopulation=std(populationInterpeakTimes);
burstinessPopulation=(stdInterPeakTimePopulation-meanInterPeakTimePopulation)/...
(stdInterPeakTimePopulation+meanInterPeakTimePopulation);
memoryPopulation= corr(populationInterpeakTimes(1:end-1)',populationInterpeakTimes(2:end)','type', 'Spearman');



figure();
histogram(arrayOfBustines);
xlabel('Burstiness Per Neuron');
ylabel('count')

figure();
histogram(arrayOfMemoryWithCluster);
xlabel('Memory Per Neuron');
ylabel('count')


%fig1=figure;
figure();
stackedplot(allRawTracesNormalized',3,1);
grid off;
colormap('parula');
ylabel('time');
view(-90,85);
saveas(gcf,strcat(experimentName,'-traces.png'));
figure();
stackedplot(allRawTracesNormalized',5,1);
xlabel('time');
saveas(gcf,strcat(experimentName,'-topviewtraces.png'));
figure();
histogram(meanInterPeakTimePerNeuron((meanInterPeakTimePerNeuron>0)));
xlabel('mean interpeak times per neuron (s)');
ylabel('count')
saveas(gcf,strcat(experimentName,'-interpeak.png'))

R=zeros(size(allRawTracesNormalized,1));
if(strcmp(correlationType,'corr'))
    R=abs(corrcoef(allRawTracesNormalized'));
elseif(strcmp(correlationType,'xcorr'))
    for i=1:size(R,1)
        for j=1:size(R,1)
            R(i,j)=max(xcorr(allRawTracesNormalized(i,:),allRawTracesNormalized(j,:),'coeff'));
            if(isnan(R(i,j)))
                pause;
            end
        end
    end
elseif(strcmp(correlationType,'fca'))
    fileID=fopen('FCAInputFile.txt','w');
    for i=1:size(allRawTracesNormalized,1)
        thisTimeCourse=spikeCourse(:,i)';
        spikeTimes=find(thisTimeCourse>0);
        for j=1:length(spikeTimes)
            fprintf(fileID,'%d %d\n',i,spikeTimes(j));
        end

    end
    fclose(fileID);
    %data_file_in,start,finish,data_file_out,measure_type,surrogate_type,num_surr,binsize,jit,width
    [R,distMatrix]=FCA('FCAInputFile.txt',1,size(allRawTracesNormalized,2),'FCAOutput','mean','jitter_uniform',10000,0,10,0);    
end

if(strcmp(correlationType,'fca'))
    clustering=dlmread('FCAOutputclust.txt',',');
    distances=dlmread('FCAOutputdist.txt',',');
    tree=[clustering,distances];
    tree(:,3)=cumsum(tree(:,3));
    figure;
    colormap(gcf,jet) 
    [H,T,Outperm]=dendrogram(tree,0,'ColorThreshold',13);
    
else
    value=sum(R(:)>0.3);
    value = value - size(R,1);
    totalSquares=(length(R(:))-size(R,1))/2;
    value=value/2;
    percentPositive=(value/totalSquares)*100
    figure();
    mySpikes= spikeCourse';
    mySpikes(mySpikes>0)=1;
    h=surf(mySpikes);
    set(h,'edgecolor','none');
    view(2);
    set(gcf,'color','w');
    saveas(gcf,strcat(experimentName,'-spikes.png'));

    %Percent of >30% correlated spikes
    R2=abs(corrcoef(mySpikes'));
    value2=sum(R2(:)>0.3);
    value2 = value2 - size(R2,1);
    totalSquares=(length(R2(:))-size(R2,1))/2;
    value2=value2/2;
    percentPositive3=(value2/totalSquares)*100
    
    cgo= clustergram(R,'colormap','jet','ShowDendrogram',true,'symmetric',false);
    set(cgo,'Linkage','complete','Dendrogram',linkageParam)
    set(gcf,'color','w');
    caxis([0,1]);

    figure;
    plot(cgo);
    saveas(gcf,strcat(experimentName,'-clustergram.png'));

    figure;
    tree=linkage(R,'complete');
    figure;
    colormap(gcf,jet) 
    [H,T,Outperm]=dendrogram(tree,0,'Reorder',optimalleaforder(tree,pdist(R)),'ColorThreshold',linkageParam);
end

set(H,'LineWidth',3)
colorArray=zeros(3,length(H));
set(gcf,'color','w');
daspect([4 1 1]);
box on;
set(gca,'XTick',[]);
set(gca,'YTick',[]);
saveas(gcf,strcat(experimentName,'-dendrogram.png'));

for i=1:length(H)
    thisLine=H(i);
    thisColor=get(H(i),'Color');
    thisYData=get(H(i),'YData');
    thisXData=get(H(i),'XData');
    for j=1:length(thisYData)
        if thisYData(j)==0
            thisXIndex=thisXData(j);
            leafIndex=Outperm(thisXIndex);
            colorArray(:,leafIndex)=thisColor;
        end
    end
end
[uniqueColors,~,colorIndex]=unique(colorArray','rows');
clusters=cell(1,length(uniqueColors));
for i=1:length(colorIndex)
    thisCluster=colorIndex(i);
    if(uniqueColors(colorIndex(i),:)==[0,0,0])
        continue;
    end
    if isempty(clusters{thisCluster})
        clusters{thisCluster}=i;
    else
        clusters{thisCluster}(end+1)=i;
    end
end


nonZeroElements=[];
for i=1:length(clusters)
    if(~isempty(clusters{i})&&length(clusters{i})>2)
        nonZeroElements=[nonZeroElements,i];
    end
end

clusters=clusters(nonZeroElements);
uniqueColors(max(uniqueColors,[],2)==0,:)=[];
%let's plot the clusters spatially back on the tiff image
figure;


set(gca,'color','k');
spatialExtentOfClusters=zeros(1,length(clusters));
totalSpatialExtentOfClusters=zeros(1,length(clusters));
neuronsInCluster=zeros(1,length(clusters));
pyramidalInCluster=zeros(1,length(clusters));
interneuronsInCluster=zeros(1,length(clusters));
clusterColors=zeros(3,length(clusters));
averagePairwiseDistances=zeros(1,length(clusters));
maximumPairwiseDistances=zeros(1,length(clusters));
averageActivationsPerFrame=zeros(1,length(clusters));
timeCourseLength=size(spikeCourse,1);
totalActivationsInCluster=zeros(1,length(clusters));
averageInterpeakTimes=zeros(1,length(clusters));
distanceBetweenPyramidal=zeros(1,length(clusters));
distanceBetweenInterneurons=zeros(1,length(clusters));

figure;
hold on;
colorArray=zeros(length(clusters),3);
for i=1:length(clusters)
    neuronsInThisCluster=clusters{i};
    clusterCenters=zeros(2,length(neuronsInThisCluster));
    activationsInCluster=0;
    thisClusterInterpeakTimes=[];
    for j=1:length(neuronsInThisCluster)
        thisNeuron=neuronsInThisCluster(j);
        
        plot(xyPositions(thisNeuron,1),xyPositions(thisNeuron,2),'o','color',uniqueColors(i,:),'lineWidth',5,'MarkerSize',10);
        text(xyPositions(thisNeuron,1)+13,xyPositions(thisNeuron,2),num2str(thisNeuron),'color',uniqueColors(i,:));
        clusterColors(:,i)=uniqueColors(i,:)';
        clusterCenters(:,j)=xyPositions(thisNeuron,:)';
        colorArray(thisNeuron,:)=clusterColors(:,i)';
        thisNeuronTimeCourse=spikeCourse(:,thisNeuron)';
        spikesInThisTimecourse=thisNeuronTimeCourse>0;
        totalSpikesInThisTimescourse=sum(spikesInThisTimecourse);
        activationsInCluster=activationsInCluster+totalSpikesInThisTimescourse;
        lastPeakTime=0;
        thisNeuronInterpeakTimes=[];
        for k=1:length(spikesInThisTimecourse)
            if spikesInThisTimecourse(k)==1
                if(lastPeakTime==0)
                    lastPeakTime=k;
                else
                   thisNeuronInterpeakTimes(end+1)=k-lastPeakTime;
                   lastPeakTime=k;
                end
            end
        end
        thisClusterInterpeakTimes=[thisClusterInterpeakTimes,thisNeuronInterpeakTimes];
    end
    xPoints=clusterCenters(1,:);
    yPoints=clusterCenters(2,:);
    [boundaryPoints,boundaryArea]=boundary(xPoints',yPoints');
    totalSpatialExtentOfClusters(i)=boundaryArea;
    neuronsInCluster(i)=length(neuronsInThisCluster);
    
    pyramidalInCluster(i)=sum(exciteNeurons(neuronsInThisCluster));
    interneuronsInCluster(i)=neuronsInCluster(i)-pyramidalInCluster(i);
    
    if pyramidalInCluster(i)>=2
        indexesOfExciteNeurons=find(exciteNeurons);
        indexesOfExciteNeuronsInThisCluster=intersect(neuronsInThisCluster,indexesOfExciteNeurons);        
        xyPositionsOfExciteNeurons=xyPositions(indexesOfExciteNeuronsInThisCluster,:);
        distanceBetweenPyramidal(i)=mean(pdist(xyPositionsOfExciteNeurons));
    end
    if interneuronsInCluster(i)>=2
        indexesOfInhNeurons=find(~exciteNeurons);
        indexesOfInhNeuronsInThisCluster=intersect(neuronsInThisCluster,indexesOfInhNeurons);        
        xyPositionsOfInhNeurons=xyPositions(indexesOfInhNeuronsInThisCluster,:);
        distanceBetweenInterneurons(i)=mean(pdist(xyPositionsOfInhNeurons));
    end
    
    spatialExtentOfCluster=boundaryArea/length(neuronsInThisCluster);
    spatialExtentOfClusters(i)=spatialExtentOfCluster;
    patch('XData',xPoints(boundaryPoints),'YData',yPoints(boundaryPoints),'FaceColor',uniqueColors(i,:),'FaceAlpha',.15,'EdgeColor',uniqueColors(i,:),'LineWidth',2);
    averagePairwiseDistance=mean(pdist([xPoints',yPoints']));
    averagePairwiseDistances(i)=averagePairwiseDistance;
    maximumPairwiseDistance=max((pdist([xPoints',yPoints'])));
    maximumPairwiseDistances(i)=maximumPairwiseDistance;
    totalActivationsInCluster(i)=activationsInCluster;
    averageActivationsPerFrame(i)=activationsInCluster/(timeCourseLength*length(neuronsInThisCluster));
    averageInterpeakTimesInCluster=mean(thisClusterInterpeakTimes);
    averageInterpeakTimes(i)=averageInterpeakTimesInCluster;
end
set(gcf,'color','k');
hold on;
for i=1:length(neuronIDs)
    thisNeuronID=neuronIDs(i);
    if(exciteNeurons(i)==1)
        color=[0,1,0];
    else
        color=[1,0,0];
    end
    plot(xyPositions(i,1),xyPositions(i,2),'o','color',color,'MarkerSize',4,'linewidth',5);    
end
fig=gcf;
fig.InvertHardcopy = 'off';
saveas(gcf,strcat(experimentName,'-clusterLocations.png'));


%peaks visualized
%saveas(gcf,strcat(num2str(fileName),'-heatmap.png'))
%pairwise cross correlation between all neurons
%correlation matrix
%plot as heatmap
%Below code counts the number of neurons "S" that are spiking at any given
%time.  PercentSpikeCorr is the percentage of the total that are spiking.
%NumPos is the number of neurons correlated at 1, 10, 20, 30% etc.  This
%also calculates the averge spike rate (Avg).
M=mean(mySpikes,2);
Smean=sum(M);
Si=size(mySpikes,1)
Avg=((Smean/Si)/0.3)
S=sum(mySpikes);
PercentSpikeCorr=S/Si;
NumPos1=sum(PercentSpikeCorr(:)>0.01);
NumPos10=sum(PercentSpikeCorr(:)>0.1);
%This is to calculate burtsiness by the Schleiss and Smith methodology, M
%is defined above
StDev=std(mySpikes,0,2);
A1=StDev-M;
B1=StDev+M;
Burstiness=A1./B1;
for k=1:size(mySpikes,1)
    idx{k}=find(mySpikes(k,:)>0);
end
B2=cellfun(@diff,idx,'UniformOutput',false);
B=cell2mat(B2);
C=cellfun(@(v)v(1),idx);
D=horzcat(B,C);


%figure out how many synchronization events there are.
synchronizationThresh=percentageSyncRequired;

figure;
subplot(6,1,[1,2]);
h=surf(mySpikes);
set(h,'linestyle','none');
view(2);
set(gcf,'color','w');
set(gca,'yTick',[]);
set(gca,'xTick',[]);


subplot(6,1,3);
synchronization=zeros(1,size(mySpikes,2));
for i=1:length(synchronization)
    spikesInWindow=sum(sum(mySpikes(:,max(1,i-windowSize):min(i+windowSize,size(mySpikes,2)))));
    spikesForEachFrame=sum(mySpikes(:,max(1,i-windowSize):min(i+windowSize,size(mySpikes,2))),1);
    maxInWindowSyncInWindow=max(spikesForEachFrame)/size(mySpikes,1);
    RelevantFrames=mySpikes(:,max(1,i-windowSize):min(i+windowSize,size(mySpikes,2)));
    totalFramesInWindow=length(RelevantFrames);
    synchronization(i)=maxInWindowSyncInWindow;
end
plot(synchronization,'k','lineWidth',2);
hold on;
plot([0,length(synchronization)],[synchronizationThresh,synchronizationThresh],'r','lineWidth',1);
set(gca,'yTick',[0,1]);
ylim([0,1]);
set(gca,'xTick',[]);
framesSynced=sum(synchronization>synchronizationThresh);
peaksAbove=bwconncomp(synchronization>synchronizationThresh);
SyncEvents=peaksAbove.NumObjects;
xlabel(sprintf('frames synced:%d, num of sync events: %d',framesSynced,SyncEvents));

subplot(6,1,4);
h=surf([synchronization;synchronization]);
set(h,'linestyle','none');
view(2);
set(gca,'yTick',[]);
set(gca,'xTick',[]);
caxis([0,1]);

subplot(6,1,[5,6]);
stackedplot(allRawTracesNormalized',3,1);
grid off;
colormap('parula');
ylabel('time');
view(-90,85);
set(gca,'yTick',[1:1000:8000]);
set(gca,'yTickLabel',[0:1000:8000]);
set(gca,'xTick',[]);
set(gca,'zTick',[]);
saveas(gcf,strcat(experimentName,'-snycEvents.png'));

figure;
meanAmplitudeOfNeuron=mean(allRawTraces,2);
try
    scatter(meanAmplitudeOfNeuron,meanInterPeakTimePerNeuron,35,colorArray,'filled');
catch
    scatter(meanAmplitudeOfNeuron,meanInterPeakTimePerNeuron,35,'filled');
end
set(gcf,'color','w');
box on;
xlabel('mean amplitude');
ylabel('mean interpeak time');
r = corrcoef(meanAmplitudeOfNeuron,meanInterPeakTimePerNeuron);
title(strcat('R=',num2str(r(1,2))));
saveas(gcf,strcat(experimentName,'-ampFreqCorr.png'));


csvwrite(strcat(experimentName,'-myfileBursty.csv'),Burstiness);
csvwrite(strcat(experimentName,'-myfileR.csv'),R);
csvwrite(strcat(experimentName,'-myfileR2.csv'),R2);
csvwrite(strcat(experimentName,'-myfileNormtrace.csv'),allRawTracesNormalized);
csvwrite(strcat(experimentName,'-myfileSpikes.csv'),mySpikes);

csvwrite(strcat(experimentName,'-meatInterpeakTimePerNeuron.csv'),meanInterPeakTimePerNeuron);
csvwrite(strcat(experimentName,'-burstiness.csv'),arrayOfBustines);
csvwrite(strcat(experimentName,'-memory.csv'),arrayOfMemoryWithCluster);
csvwrite(strcat(experimentName,'-burstinessAndMemoryPopulation.csv'),[burstinessPopulation,memoryPopulation]);
csvwrite(strcat(experimentName,'-averagePairwiseDistances.csv'),averagePairwiseDistances);
csvwrite(strcat(experimentName,'-spatialExtentOfClusters.csv'),spatialExtentOfClusters);
csvwrite(strcat(experimentName,'-maximumPairwiseDistances.csv'),maximumPairwiseDistances);
csvwrite(strcat(experimentName,'-averageActivationsPerFrame.csv'),averageActivationsPerFrame);
csvwrite(strcat(experimentName,'-clusterColors.csv'),clusterColors);
csvwrite(strcat(experimentName,'-totalSpatialExtentOfClusters.csv'),totalSpatialExtentOfClusters);
csvwrite(strcat(experimentName,'-neuronsInCluster.csv'),neuronsInCluster);
csvwrite(strcat(experimentName,'-totalActivationsInCluster.csv'),totalActivationsInCluster);
csvwrite(strcat(experimentName,'-averageInterpeakTimes.csv'),averageInterpeakTimes);
    
csvwrite(strcat(experimentName,'-pyramidal.csv'),pyramidalInCluster);
csvwrite(strcat(experimentName,'-interneurons.csv'),interneuronsInCluster);
csvwrite(strcat(experimentName,'-distanceBetweenPyramidal.csv'),distanceBetweenPyramidal);
csvwrite(strcat(experimentName,'-distanceBetweenInterneurons.csv'),distanceBetweenInterneurons);

end
