clear all;
clc;

%%remove silence
allFiles = 'allList.txt';
fid = fopen(allFiles);
myData = textscan(fid, '%s');
fclose(fid);

audioPath = '/Users/yucongwang/Documents/EE214/project/new/';
myFiles = myData{1};
for iter = 1 : length(myFiles)
    audioName = strcat(audioPath, myFiles{iter});
    [audioFileVoiced, Fs] = detectVoiced(myFiles{iter});
    audiowrite(audioName, audioFileVoiced{1}, Fs);
    if(mod(iter,10)==0)
        disp(['Rewrite ',num2str(iter),' of ',num2str(length(myFiles)),' files.']);
    end
end

%%
clear all;
clc;
% Define lists
allFiles = 'allList.txt';
% trainList = 'trainCleanList.txt';
% testList = 'testCleanList.txt';
tic

% Extract features
featureDict = containers.Map;
fid = fopen(allFiles);
myData = textscan(fid,'%s');
fclose(fid);
myFiles = myData{1};
for(i = 1:length(myFiles))
    [snd,fs] = audioread(myFiles{i});
%      mfccs = msf_mfcc(snd,fs,'nfilt',40,'ncep',12);
%      lpcs = msf_lpc(snd,fs,'order',10);
%      feature = [mfccs,lpcs];
     feature=PNCC(snd,fs);
     
    featureDict(myFiles{i}) = feature';
    featureSet{i} = feature';
    if(mod(i,10)==0)
        disp(['Completed ',num2str(i),' of ',num2str(length(myFiles)),' files.']);
    end
end
%%
% Train the classifier
trainList = 'trainCleanList.txt';
clc;
fid = fopen(trainList);
myData = textscan(fid,'%s %s %f');
fclose(fid);
fileList1 = myData{1};
fileList2 = myData{2};
labels = myData{3};

for(i = 1:length(labels))
     dist{i} = dtw(featureDict(fileList1{i}),featureDict(fileList2{i}));
end
%%
% Train threshold
DMat = cell2mat(dist);
[~,threshold] = compute_eer(DMat,labels);
%%
% Test the classifier
testList = 'trialsEvaluation.txt';

fid = fopen(testList);
myData = textscan(fid,'%s %s');
fclose(fid);
fileList1 = myData{1};
fileList2 = myData{2};
%testlabels = myData{3};
% scores = zeros(length(labels),1);
% bCell = cell(length(testlabels),1);

for(i = 1:length(fileList1))
     scores(i) = dtw(featureDict(fileList1{i}),featureDict(fileList2{i}));
     prediction(i) = (scores(i)<= threshold);
     
end

predictions = double(prediction)';

% 
% FPR = sum(~testlabels & predictions)/sum(~testlabels);
% FNR = sum(testlabels & ~predictions)/sum(testlabels);
% disp(['The false positive rate is ',num2str(FPR*100),'%.'])
% disp(['The false negative rate is ',num2str(FNR*100),'%.'])

toc

%%
fid = fopen('result.txt','w');
fprintf(fid,'%f \n',predictions);
fclose(fid);
%%


result=[predictions,testlabels];
roc(result);