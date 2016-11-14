clear all

mutPercentage = tdfread('mut_percent.txt');
clinical = tdfread('Clinical.txt');

fileID = fopen('MatchClinical.txt','w');
fprintf(fileID,'Barcode\tT\tN\tM\tStage\tSurvival\tStatus\n');
fclose(fileID);

for i = 1:size(mutPercentage.Barcode,1)
    
    for j = 1:size(clinical.PatientID,1)
        if(clinical.PatientID(j,:) == mutPercentage.Barcode(i,:))
          fileID = fopen('MatchClinical.txt','a');
          fprintf(fileID,'%s\t%s\t%s\t%s\t%s\t%s\t%s\n',clinical.PatientID(j,:),...
              clinical.T(j,:),clinical.N(j,:),clinical.M(j,:),...
              clinical.Stage(j,:),clinical.Survival(j,:),clinical.Status(j,:));
          fclose(fileID);
            
       
        end
        
    end
    
    
end