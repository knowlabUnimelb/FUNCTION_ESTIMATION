clear all
clc

load parameterClusters

for i = 1:3
    [r, ~] = corrcoef(covParms(:,:,i));
    temp = triu(circshift(triu(r), [0 -1])); 
    temp(:,end) = []; 
    temp = temp(:);
    temp(temp == 0) = []; 
    
    meanr(1,i) = mean(temp); 
    stdr(1,i)  = std(temp);
    [h(1,i),p(1,i),ci(:,i,1),stats(1,i)]=ttest(temp, 0);
    
    
%     correctedParms(:,:,i) = covParms(:,:,i) - repmat(mean(covParms(:,:,i)), size(covParms, 1), 1) + mean(reshape(covParms(:,:,i), numel(covParms(:,:,i)), 1));
%     [r2, ~] = corrcoef(correctedParms(:,:,i));
%     temp = triu(circshift(triu(r2), [0 -1])); 
%     temp(:,end) = []; 
%     temp = temp(:);
%     temp(temp == 0) = []; 
%     
%     meanr(2,i) = mean(temp); 
%     stdr(2,i)  = std(temp);
%     [h(2,i),p(2,i),ci(:,i,2),stats(2,i)]=ttest(temp, 0);
end