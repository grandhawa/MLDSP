function [ sFlag ] = downloadFasta( fname )
%function to download .fasta files 
%input is .txt file with comma seperated accession numbers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Gurjit Singh Randhawa  %
% Department of Computer Science,%
% Western University, Canada     %
% email: grandha8@uwo.ca         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    warning('off','all');
    sFlag = 0;
    path = pwd;
    info  = dir(fname);
    if info.bytes == 0
         error('Error : Empty txt file , add accession numbers or delete the empty file.')
    end    
    txt = fileread(fname);
    aNmList = strsplit(txt,',');
    len = length(aNmList);
    fnm = strsplit(fname,'.');
    fInfo= dir(fnm{1});
    if ~isempty(fInfo)
        return
    end
    mkdir(fnm{1});
    cd(fnm{1});
    for a = 1:len
        seq = getgenbank(aNmList{a},'FileFormat', 'fasta');
        sname = strcat(aNmList{a},'.fasta');
        fastawrite(sname,seq);    
    end  
    sFlag=1;
    cd(path);
end

