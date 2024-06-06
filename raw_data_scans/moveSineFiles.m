clear all
clc

folders = {'Scanned Data', 'Other Scanned Data'};

for i = 1:numel(folders)
    fd = dir(fullfile(pwd, folders{i})); fd(1:2) = [];
    for j = 1:numel(fd)
        if fd(j).isdir
            [success, message, messid] =...
              copyfile(fullfile(pwd, folders{i}, fd(j).name, '*_9.JPG'),...
                     'C:\Documents and Settings\littled\My Documents\My Dropbox\Dan''s Experiments\2011 Sem2 Social Functions');
                 
            [success, message, messid] =...
              copyfile(fullfile(pwd, folders{i}, fd(j).name, '*_21.JPG'),...
                     'C:\Documents and Settings\littled\My Documents\My Dropbox\Dan''s Experiments\2011 Sem2 Social Functions');
        end
    end
end