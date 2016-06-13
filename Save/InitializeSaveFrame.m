%Initialize SaveFrame

FileLocation='C:\Users\thatmathguy\Documents\MCTP\MATLAB';%where to create new folder
NewFolder='saveframetest';%name of new folder to be created
Title='plot';%name prefix for saved files
if exist(strcat(FileLocation,'\',NewFolder),'file')
    rmdir(strcat(FileLocation,'\',NewFolder),'s');
end
mkdir(FileLocation,NewFolder)%creates new folder