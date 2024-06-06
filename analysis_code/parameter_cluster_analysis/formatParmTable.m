function [Nc, v] = formatParmTable(clusters, clusterData, idx)


cnt = aggregate([clusters{idx},  clusterData{idx}], 1, [2 3 4], @count,1)';
Nc = numel(cnt(cnt>10));
N = sum(cnt(cnt > 10));
v  = N * cnt./N .* (1 - cnt./N);
m   = exp(aggregate([clusters{idx},  clusterData{idx}], 1, [2 3 4], @mean,1)');
s   = exp(aggregate([clusters{idx},  clusterData{idx}], 1, [2 3 4], @std,1)');

for i = 1:numel(cnt)
    if cnt(i) > 2
        [r, p] = corrcoef(clusterData{idx}(clusters{idx} == i, :));
        
        corrs(i,:) = [r(1,2), r(1,3), r(2,3)];
        ps(i,:)    = [p(1,2), p(1,3), p(2,3)];
    else
        corrs(i,:) = nan(1,3);
        ps(i,:) = nan(1,3);
    end
end

str1 = repmat('%10d\t', 1, numel(cnt)); str1(end) = 'n';
str2 = repmat('%10.2f (%3.2f)\t', 1, numel(cnt)); str2(end) = 'n';
str3 = repmat('%2.2f, %2.2f, %2.2f\t', 1, numel(cnt)); str3(end) = 'n';

fprintf(str1, cnt)
for i = 1:3
    fprintf(str2, weave(m(i,:),s(i,:)))
end
corrs = corrs';
fprintf(str3, corrs(:)')

disp(' ')
disp('---------------------------------------------------')
disp(ps)

