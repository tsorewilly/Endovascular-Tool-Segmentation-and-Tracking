%function ProcessData(FolderDir)
%clear; 
%[images, headers] = dicomfolder(FolderDir);
warning off;

%[fileEMG,  = uigetfile('*.acq', 'Select Signal To Process', FolderDir);
%FolderDir = uigetdir('E:\');
FolderDir ='E:\Rabbit\2021.05.11\';

%Load Image Frames
Lists=dir(fullfile(FolderDir, '*.ima'));
imageLists = {Lists.name};
img_NUM = length(imageLists);
LagLead=1; incr = 5; Start=1;%While incr of 11 is 1 Sec Image Data Frames, LagLead can be optimized using rawViewer.exe
 %%
tracker=1; oldXcoord=1; oldYcoord=1;
BgFrame = Start;
se = strel('diamond',7); %size(se.Neighborhood)
se_wire=strel('line',17,100); 
n = 1; %ceil(6.*rand(1,1)) % n = neighborhood

Options.FrangiScaleRange=[7 7];

%Define the background   
bgFrame = dicomread(fullfile(FolderDir, imageLists{BgFrame}));
figure(101); clf; imshow(bgFrame,[]); title('Guidewire Tracking Frames');
Temp=zeros(size(bgFrame));
% hWaitBar = waitbar(0, 'Reading the CT Images');

for f = LagLead : 1 : img_NUM-LagLead
    startTime = tic;
    dsFrame = dicomread(fullfile(FolderDir, imageLists{f}));
    figure(1002); clf; imshow(dsFrame,[]); title(['Guidewire Navigation in Frame ', num2str(f)]);
    pause(0.2)
end 