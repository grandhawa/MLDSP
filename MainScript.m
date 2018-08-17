% ML-DSP
% the main program that user needs to run for the results.
% default dataSet is Primates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Gurjit Singh Randhawa  %
% Department of Computer Science,%
% Western University, Canada     %
% email: grandha8@uwo.ca         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear workspace
close all;
clear all;
clc ;
set(0, 'DefaultFigureRenderer', 'painters');
% can be uncommented for MATLAB2017 or earlier versions
%digits(4);

%read fasta files
dataSet = 'Primates';
selectedFolder = dataSet;
fprintf('Reading sequences .... \n');
[AcNmb, Seq, numberOfClusters, clusterNames, pointsPerCluster] = readFasta(dataSet);
totalSeq = length(Seq);

%calculate length stats
[maxLen, minLen, meanLen, medLen] = lengthCalc(Seq);

%numerical sequences, length normalized by median length
%code can be modified to use other length stats for length normalization
mLen = medLen;
nmValSH=cell(1,totalSeq);
f=cell(1,totalSeq);
lg=cell(1,totalSeq);

fprintf('Generating numerical sequences, applying DFT, computing magnitude spectra .... \n');
parfor a = 1:totalSeq
    %using Purin/Pyramidine representation by default
    %change method call for other representations
    ns = numMappingPP(Seq{a});
    %change "medLen" to other length stat for length normalization
    I = mLen-length(ns);
    if(I>0)
        nsTemp = wextend('1','asym',ns,I);
        nsNew = nsTemp((I+1):length(nsTemp));
    elseif(I<0)
        nsNew=ns(1:mLen);
    else
        nsNew = ns;
    end
    nmValSH{a} = nsNew;
    %fourier transform
    f{a} = fft(nsNew);
    %magnitude spectra
    lg{a} = abs(f{a});   
end

%distance calculation by Pearson correlation coefficient
% change 'cor' to 'euc' for Euclidean
fprintf('Computing Distance matrix .... \n');
fm=cell2mat(lg(:));
disMat = f_dis(fm,'cor',0,1);

%Multi-dimensional Scaling
fprintf('Performing Multi-dimensional scaling .... \n');
[Y,eigvals] = cmdscale(disMat);

%3D  plot
fprintf('Generating 3D plot .... \n');
index=1;
counter=1;
Cluster = zeros(1,totalSeq);
for i=1:totalSeq   
    Cluster(i)=index;
    if(counter==pointsPerCluster{index})
        index=index+1;
        counter=0;
    end
    counter= counter+1;
end
uniqueClusters  = unique(Cluster);
cmap = distinguishable_colors(numberOfClusters);
hf = figure;
hold on;
for h=1:numberOfClusters
    cIndex = Cluster == uniqueClusters(h);
    plot3(Y(cIndex,1),Y(cIndex,2),Y(cIndex,3),'.','markersize', 15, 'Color',cmap(h,:),'DisplayName',clusterNames{h});
end
view(3), axis vis3d, box on, datacursormode on
xlabel('x'), ylabel('y'), zlabel('z')
tname = strcat(selectedFolder,' (',int2str(totalSeq),' Sequences',')');
title(tname)
hdt = datacursormode(hf);
set(hdt,'UpdateFcn',{@myupdatefcn,Y,AcNmb})
legend('show');


%Phylogenetic Tree
% plotted for datasets with atmost 50 sequences
% parameter can be changed in following statement 
if(totalSeq<=50)
    fprintf('Creating Phylogenetic Tree .... \n');
    UPGMAtree = seqlinkage(disMat,'UPGMA',AcNmb);
    plot(UPGMAtree);
end

%Classification Code
clear a;
a=[];
for i=1:numberOfClusters
    for j=1:pointsPerCluster{i}
        a=[a; i];
    end
end
ATestlg = [disMat a];
rng(15,'twister');

alabels = a;
fprintf('Performing classification .... \n');
folds=10;
if (totalSeq<folds)
    folds = totalSeq;
end
[accuracy, avg_accuracy, clNames] = classificationCode(disMat,alabels, folds, totalSeq);
acc = [accuracy avg_accuracy];
s.ClassifierModel=cellstr(clNames.');
s.Accuracy=cell2mat(acc).';
ClassificationAccuracyScores = struct2table(s)

%new sequence classification example
%uncomment the following line and replace the accesion number 
%newSeqClassify('NC_028718.1', mLen, disMat, lg, alabels, clusterNames)

fprintf('**** Processing completed ****\n');

