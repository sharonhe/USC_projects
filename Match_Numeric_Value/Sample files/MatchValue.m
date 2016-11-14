clear all

mutCount = tdfread('mut_percent.txt');
mRNA = tdfread('Additional_mRNA.txt');

MatchedValue = zeros(size(mutCount.Barcode,1),15);

for i = 1:size(mutCount.Barcode,1)
    
    for j = 1:size(mRNA.COMMON,1)
        if(mRNA.COMMON(j,1:12) == mutCount.Barcode(i,1:12))
          MatchedValue(i,1) = mRNA.ATG9A(j,1);
          MatchedValue(i,2) = mRNA.BECN1(j,1);
          MatchedValue(i,3) = mRNA.CUL4A(j,1);
          MatchedValue(i,4) = mRNA.DDB1(j,1);
          MatchedValue(i,5) = mRNA.DDB2(j,1);
          MatchedValue(i,6) = mRNA.ERCC2(j,1);
          MatchedValue(i,7) = mRNA.ERCC3(j,1);
          MatchedValue(i,8) = mRNA.ERCC4(j,1);
          MatchedValue(i,9) = mRNA.ERCC5(j,1);
          MatchedValue(i,10) = mRNA.PIK3C3(j,1);
          MatchedValue(i,11) = mRNA.PRKDC(j,1);
          MatchedValue(i,12) = mRNA.XPA(j,1);
          MatchedValue(i,13) = mRNA.XPC(j,1);
          MatchedValue(i,14) = mRNA.XRCC5(j,1);
          MatchedValue(i,15) = mRNA.XRCC6(j,1);
        end
        
    end
    
    
end
