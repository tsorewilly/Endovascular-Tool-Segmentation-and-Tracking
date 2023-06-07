clear; clc; 
prntFlda = 'J:\20200715-02\New folder\';
Lists=dir(fullfile(prntFlda));
FldList = {Lists.name};


for fld=3:length(FldList)
    FileList=dir(fullfile([prntFlda, '\', FldList{fld}], '*dcm'));
    Files = {FileList.name};
    for dsFile=1:length(Files)
        fNAme = ['F:\20200715\DSA\',Files{dsFile}];
        if ~isfile(fNAme)
            copyfile([prntFlda, '\', FldList{fld}, '\', Files{dsFile}], fNAme)
        else
            disp([fNAme, ' already exists']);
        end
    end
end

