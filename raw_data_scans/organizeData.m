clear all
clc

datafolder = 'Other Scanned Data'; % 'Data'; 
dn = dir(fullfile(pwd, datafolder)); % read directory contents
% each folder is a subject

data = [];
for i = 3:numel(dn)
    if dn(i).isdir % if folder
        fn = dir(fullfile(pwd, datafolder, dn(i).name, '*.txt'));
        subj = str2double(dn(i).name);
        disp(sprintf('==========================\nNow processing subject %d\n==========================\n', subj))
        for j = 1:numel(fn)
            fid = fopen(fullfile(pwd, datafolder, dn(i).name, fn(j).name));
            xy = textscan(fid, '%f%f', 'HeaderLines', 1, 'Delimiter', ',');
            fclose(fid);
            
            [pref1, suff1] = strtok(fn(j).name, '_'); 
            [pref2, suff2] = strtok(suff1, '.');
            [pref3, suff3] = strtok(pref2, '_');
            pref2 = deblank(pref2);
            pref3 = fliplr(pref3); %pref3(1) = [];
            if isletter(pref3(1))
                pref3(1) = []; 
            end
            qn = str2num(fliplr(pref3));
            
            if j == 1; 
                oldqn = qn; 
                disp(sprintf('Now processing question %d', qn))
            else
                if oldqn ~= qn; 
                    disp(sprintf('Now processing question %d', qn))
                    oldqn = qn;
                end
            end
            
            data = [data;... 
                repmat(subj, size(xy{1}, 1), 1),... % subject number
                repmat(qn, size(xy{1}, 1), 1),... % question number
                xy{1},... % x
                xy{2}];   % y
        end
    end
end