function [latVecTif,lonVecTif,tifData] = readTifAndTfw(tifFileName)
tifData = 0;
lonVecTif = 0;
latVecTif = 0;
try
    tifData = imread(tifFileName);
catch
    msgbox('No TIF File Found','TIF Plot Error','error');
end
tfwFileName = tifFileName;
tfwFileName(end-2:end) = 'tfw';
try
    tfwData = load(tfwFileName);
    lonVecTif = tfwData(5):tfwData(1):tfwData(5)+(size(tifData,2)-1)*tfwData(1).';
    latVecTif = (tfwData(6):tfwData(4):tfwData(6)+(size(tifData,1)-1)*tfwData(4));
catch
    msgbox('No TFW File Found','TIF Plot Error','error');
end

