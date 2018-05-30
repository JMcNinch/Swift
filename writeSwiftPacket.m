function writeSwiftPacket(fid,UTCHH,UTCMM,UTCSS,UTChh,lat,lon,track,...
                            bearingZeroError, bitPerSize, compassInvalid,...
                            spokeSamples, rangeCellsDiv2, rangeCellSize_mm ,...
                        spokeSequenceNumber, spokeAzimuth, spokeCompass,...
                        spokeLengthBytes, trueNorth, intensity, isgood)
%% File Format
%4 UTCHH,UTCMM,UTXSS,UTChh (int8)
%3 Lat,Lon,Track, (float)
%3 BearingError,BitsPerSamp,CompassInv, (int8)
%8 Samples,CellsDiv2,  CellSize,Sequence,Azimuth,Compass, SpokeLength,TrueNorth (int 16)
%1024 Data (int8)
%1 isgood (int8)

%numElements = 1043;

for n = 1%:size(data,1)
%     fwrite(fid,data(n,1:4),'int8');
%     fwrite(fid,data(n,5:7),'float');
%     fwrite(fid,data(n,8:10),'int8');
%     fwrite(fid,data(n,11:18),'int16');
%     fwrite(fid,data(n,19:(18+1024)),'int8');
%     fwrite(fid,data(n,1043),'int8');
     fwrite(fid,[UTCHH UTCMM UTCSS UTChh],'int8');
    fwrite(fid,[lat lon track],'float');
    fwrite(fid,[bearingZeroError bitPerSize compassInvalid],'int8');
    fwrite(fid,[spokeSamples rangeCellsDiv2 rangeCellSize_mm spokeSequenceNumber spokeAzimuth spokeCompass ...
                spokeLengthBytes trueNorth],'int16');
    fwrite(fid,intensity,'int8');
    fwrite(fid,isgood,'int8');
end

