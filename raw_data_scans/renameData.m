clear all
clc

dataFolder = fullfile(pwd, 'Other Scanned Data');
oldSubjectId = {'A1' , 'A2' , 'A3' , 'A4' , 'A5' , 'A6' , 'A7' , 'A8' , 'A9' , 'A10', 'A11'};
newSubjectId = {'171', '172', '173', '174', '175', '176', '177', '178', '179', '180', '181'};   

for j = 2:numel(oldSubjectId)
    currentFolder = pwd;
    cd(fullfile(dataFolder, oldSubjectId{j}))
    folderContents = dir(fullfile(pwd, sprintf('%s*', oldSubjectId{j})));
    
    for i = 1:size(folderContents, 1);
        oldname = folderContents(i).name;
        [pre, suf] = strtok(folderContents(i).name, '_');
        newname = sprintf('%s%s', newSubjectId{j}, suf);
        
        movefile(fullfile(pwd, oldname), fullfile(pwd, newname))
        
    end
end

cd(currentFolder)
