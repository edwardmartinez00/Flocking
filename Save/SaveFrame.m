%SaveFrame
%Saves frames for compilation into video
%kstep, step number
function SaveFrame(FileLocation,NewFolder,Title,kstep)
    FileName=strcat(FileLocation,'\',NewFolder,'\',Title,num2str(kstep));
    print(FileName,'-djpeg')
end