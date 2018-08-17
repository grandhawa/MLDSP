function [ numSeq ] = numMappingCodons( sq)
%function for Codon representation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Gurjit Singh Randhawa  %
% Department of Computer Science,%
% Western University, Canada     %
% email: grandha8@uwo.ca         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    len = length(sq);  
    numSeq = zeros(1,len,'double');
    codons = {'TTT','TTC','TTA','TTG','CTT','CTC','CTA','CTG','TCT','TCC','TCA','TCG','AGT','AGC','TAT','TAC',
              'TAA','TAG','TGA','TGT','TGC','TGG','CCT','CCC','CCA','CCG','CAT','CAC','CAA','CAG','CGT','CGC',
              'CGA','CGG','AGA','AGG','ATT','ATC','ATA','ATG','ACT','ACC','ACA','ACG','AAT','AAC','AAA','AAG',
              'GTT','GTC','GTA','GTG','GCT','GCC','GCA','GCG','GAT','GAC','GAA','GAG','GGT','GGC','GGA','GGG'};
    
    for K = 1:len
            if(K<len-1)
                t = sq(K:K+2);
            elseif(K==len-1)
                t = strcat(sq(K:K+1),sq(1:1));
            else
                t = strcat(sq(K:K),sq(1:2));
            end
            tp = find(strcmpi(codons, t))-1;
            numSeq(K) = tp;    
    end       
end

