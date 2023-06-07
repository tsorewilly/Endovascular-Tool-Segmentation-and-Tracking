clear; clc;
warning off;
FolderDir = 'G:\Documents\My Postdoc Work\Publications\IJPRAI\Heneghan\';

%Load Image Frames    
Lists=dir(fullfile(FolderDir, '*.png'));
imageLists = {Lists.name};

for id = 1:length(imageLists)
    % Get the file name 
    [~, oldName,ext] = fileparts(Lists(id).name);
    newName = strcat(extractBefore(oldName, '1_bw'), '8',ext);
    movefile([FolderDir, Lists(id).name], [FolderDir, newName]); 
end