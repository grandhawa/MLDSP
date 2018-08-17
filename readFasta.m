function [AcNmb, Seq, numberOfClusters, clusterNames, pointsPerCluster] = readFasta(dataSet)
% read fasta files in folder DataBase/'dataSet'
% make subFolders for each cluster and place respective fasta sequences insides
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Gurjit Singh Randhawa  %
% Department of Computer Science,%
% Western University, Canada     %
% email: grandha8@uwo.ca         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    path = pwd;
    dbPath = strcat(path,'\','DataBase');
    cd(dbPath);
    tStr = strcat(dataSet,'.mat');
    if ~isempty(dir(tStr))
        load(tStr);
    else
        dbPath = strcat(dbPath,'\',dataSet);
        folderInfo = dir (dbPath);
        Seq = {};
        AcNmb = {};    

        if(isempty(folderInfo))
             error('Error : DataSet does not exist.')
        else
%             numberOfClusters = length(folderInfo)-2;
%             clusterNames = cell(1,numberOfClusters);
%             pointsPerCluster = cell(1,numberOfClusters);
            cd(dbPath); 
            for i=3:length(folderInfo)
               fnm=folderInfo(i).name; 
               subFolderPath = strcat(dbPath,'\',fnm);
               subFolderInfo = dir(subFolderPath);
               result = subFolderInfo.isdir;
               if ~(result)
                    flag = downloadFasta( fnm );
               end
            end
            
            numberOfClusters = 0;
            index = 1;
            folderInfo = dir (dbPath);
            for i=3:length(folderInfo)
                cd(dbPath);                
                subFolderPath = strcat(dbPath,'\',folderInfo(i).name);
                subFolderInfo = dir(subFolderPath);
                result = subFolderInfo.isdir;
                if ~(result)
                    continue;
                end
                cd(subFolderPath);
                clusterNames{index} =  folderInfo(i).name;
                pts = length(subFolderInfo)-2;
                pointsPerCluster{index} = pts;
                seqTemp = cell(1,pts);
                acTemp = cell(1,pts);
                for j=3:length(subFolderInfo)
                    [Header, Sequence] = fastaread(subFolderInfo(j).name);    
                    Sequence = regexprep(Sequence,'[^A,^C, ^G, ^T]','','ignorecase');
                    seqTemp{j-2} = Sequence;
                    acTemp{j-2} = Header; 
                end
                Seq = [Seq seqTemp];
                AcNmb = [AcNmb acTemp]; 
                numberOfClusters = numberOfClusters+1;
                index = index+1;
            end 
            if(numberOfClusters<2)
                error('Error : DataSet should have atleast two clusters.')            
            end
        end
    end
    cd(path);
end

