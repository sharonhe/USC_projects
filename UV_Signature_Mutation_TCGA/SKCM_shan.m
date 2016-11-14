clear all
close all
clc

tdfread('SKCM.txt');
[Patients,IA,IC] = unique(Tumor_Sample_Barcode,'rows');

fileID = fopen('mutCount.txt','w');
fprintf(fileID,'Barcode\tACO\tAGO\tATO\tCAO\tCGO\tCTO\tGAO\tGCO\tGTO\tTAO\tTCO\tTGO\tACX\tAGX\tATX\tCAX\tCGX\tCTX\tGAX\tGCX\tGTX\tTAX\tTCX\tTGX\n');
fclose(fileID);

for i = 1:size(Patients,1)
    
    clear PatientMutation;
    clear PatientAjacentBases;
    
    PatientMutation(:,1) = Reference_Allele((IC == i) & (Variant_Type(:,1) == 'S'));
    
    PatientMutation(:,2) = Tumor_Seq_Allele2((IC == i) & (Variant_Type(:,1) == 'S'));
    
    PatientAjacentBases = Ajacent_Bases((IC == i) & (Variant_Type(:,1) == 'S'),:);
    
    for j = 1:size(PatientMutation,1)
        
        if (PatientMutation(j,1) == 'C')
            if (PatientAjacentBases(j,:) == 'ACG'|...
               PatientAjacentBases(j,:) == 'ACA'|...
               PatientAjacentBases(j,:) == 'GCG'|...
               PatientAjacentBases(j,:) == 'GCA')
               PatientMutation(j,3) = 'X';
            else
                PatientMutation(j,3) = 'O';
            end
          
            
        elseif (PatientMutation(j,1) == 'G')
            if (PatientAjacentBases(j,:) == 'TGC'|...
               PatientAjacentBases(j,:) == 'TGT'|...
               PatientAjacentBases(j,:) == 'CGC'|...
               PatientAjacentBases(j,:) == 'CGT')
               PatientMutation(j,3) = 'X';
            else
                PatientMutation(j,3) = 'O';
            end
            
        
        elseif (PatientMutation(j,1) == 'A')
            if (PatientAjacentBases(j,:) == 'TAC'|...
               PatientAjacentBases(j,:) == 'TAT'|...
               PatientAjacentBases(j,:) == 'CAC'|...
               PatientAjacentBases(j,:) == 'CAT')
               PatientMutation(j,3) = 'X';
            else
                PatientMutation(j,3) = 'O';
            end
            
        elseif (PatientMutation(j,1) == 'T')
            if (PatientAjacentBases(j,:) == 'ATG'|...
               PatientAjacentBases(j,:) == 'ATA'|...
               PatientAjacentBases(j,:) == 'GTG'|...
               PatientAjacentBases(j,:) == 'GTA')
               PatientMutation(j,3) = 'X';
            else
                PatientMutation(j,3) = 'O';
            end
        end
        
    end
    
    [Mutations,mA,mC] = unique(PatientMutation,'rows');
    
    for j = 1:size(Mutations,1)
        MutationCount(j,1) = sum(mC == j);
    end
    
    % ACO
    ACO = MutationCount(Mutations(:,1) == 'A' & Mutations(:,2) == 'C' & Mutations(:,3) == 'O');
    % AGO
    AGO = MutationCount(Mutations(:,1) == 'A' & Mutations(:,2) == 'G' & Mutations(:,3) == 'O');
    % ATO
    ATO = MutationCount(Mutations(:,1) == 'A' & Mutations(:,2) == 'T' & Mutations(:,3) == 'O');
    % ACX
    ACX = MutationCount(Mutations(:,1) == 'A' & Mutations(:,2) == 'C' & Mutations(:,3) == 'X');
    % AGX
    AGX = MutationCount(Mutations(:,1) == 'A' & Mutations(:,2) == 'G' & Mutations(:,3) == 'X');
    % ATX
    ATX = MutationCount(Mutations(:,1) == 'A' & Mutations(:,2) == 'T' & Mutations(:,3) == 'X');
    % CAO
    CAO = MutationCount(Mutations(:,1) == 'C' & Mutations(:,2) == 'A' & Mutations(:,3) == 'O');
    % CGO
    CGO = MutationCount(Mutations(:,1) == 'C' & Mutations(:,2) == 'G' & Mutations(:,3) == 'O');
    % CTO
    CTO = MutationCount(Mutations(:,1) == 'C' & Mutations(:,2) == 'T' & Mutations(:,3) == 'O');
    % CAX
    CAX = MutationCount(Mutations(:,1) == 'C' & Mutations(:,2) == 'A' & Mutations(:,3) == 'X');
    % CGX
    CGX = MutationCount(Mutations(:,1) == 'C' & Mutations(:,2) == 'G' & Mutations(:,3) == 'X');
    % CTX
    CTX = MutationCount(Mutations(:,1) == 'C' & Mutations(:,2) == 'T' & Mutations(:,3) == 'X');
    % GAO
    GAO = MutationCount(Mutations(:,1) == 'G' & Mutations(:,2) == 'A' & Mutations(:,3) == 'O');
    % GCO
    GCO = MutationCount(Mutations(:,1) == 'G' & Mutations(:,2) == 'C' & Mutations(:,3) == 'O');
    % GTO
    GTO = MutationCount(Mutations(:,1) == 'G' & Mutations(:,2) == 'T' & Mutations(:,3) == 'O');
    % GAX
    GAX = MutationCount(Mutations(:,1) == 'G' & Mutations(:,2) == 'A' & Mutations(:,3) == 'X');
    % GCX
    GCX = MutationCount(Mutations(:,1) == 'G' & Mutations(:,2) == 'C' & Mutations(:,3) == 'X');
    % GTX
    GTX = MutationCount(Mutations(:,1) == 'G' & Mutations(:,2) == 'T' & Mutations(:,3) == 'X');
    % TAO
    TAO = MutationCount(Mutations(:,1) == 'T' & Mutations(:,2) == 'A' & Mutations(:,3) == 'O');
    % TCO
    TCO = MutationCount(Mutations(:,1) == 'T' & Mutations(:,2) == 'C' & Mutations(:,3) == 'O');
    % TGO
    TGO = MutationCount(Mutations(:,1) == 'T' & Mutations(:,2) == 'G' & Mutations(:,3) == 'O');
    % TAX
    TAX = MutationCount(Mutations(:,1) == 'T' & Mutations(:,2) == 'A' & Mutations(:,3) == 'X');
    % TCX
    TCX = MutationCount(Mutations(:,1) == 'T' & Mutations(:,2) == 'C' & Mutations(:,3) == 'X');
    % TGX
    TGX = MutationCount(Mutations(:,1) == 'T' & Mutations(:,2) == 'G' & Mutations(:,3) == 'X');
    
    % Write this patient's stats to file
    fileID = fopen('mutCount.txt','a');
    fprintf(fileID,'%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n',Tumor_Sample_Barcode(IA(i),1:12),ACO,AGO,ATO,CAO,CGO,CTO,GAO,GCO,GTO,TAO,TCO,TGO,ACX,AGX,ATX,CAX,CGX,CTX,GAX,GCX,GTX,TAX,TCX,TGX);
    fclose(fileID);
    
end