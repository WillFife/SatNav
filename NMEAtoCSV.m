%  Original by: Robert B Harris
%  Modified for Matlab by: Hector Escobar
%  Last Modified: September 1, 2020 by Vishnuu Mallik
%  Modified: September 16, 2010. 
%  Modifications: Addapted for Matlab. Fixed bug where sod didn't show up.
%  Number of satellites tracked added. 
%
%  Modified: September 1, 2020. 
%  Modifications: According to 
%  https://www.gpsworld.com/what-exactly-is-gps-nmea-data/
%  latitude is expressed in DDMM.MMMMM format and longitude is expressed in
%  DDDMM.MMMMM format, aka degrees and decimal minutes. The last
%  modification of this script assumed that longitude was expressed in
%  DDMM.MMMMM format, causing an error in the output.
%
%
%  [sod, lats, longs, hts, sats] = readNMEApos(fname)
%
%     Reads NMEA messages recorded to a text file, and transforms the GPGGA
%     (position fix) messages into four vectors:
%
%        sod   -- Seconds of day.
%        lats  -- Latitude, with decimal fractions, in whole degrees.
%                 North is positive.
%        longs -- Longitude, same units as latitude. East is positive.
%        hts   -- Height above ellipsoid, in meters. 
%        sats  -- Number of satellites in view

function [sod, lats, longs, hts, sats] = readNMEAposV3(fname)
fid = fopen(fname,'rt');

numRead = 0;
numGood = 0;
numGoodPos = 0;

% First read -- count number of spaces

tline=fgetl(fid);
while ischar(tline)
   
   numRead = numRead+1;
   
   % Extract the recorded checksum
   csumidx = findstr(tline,'*');
   ocsum = hex2dec(tline(end-1:end));
   
   % Comput expected checksum
   ccsum = 0;
   for k=2:(csumidx-1)
     ccsum = bitxor(ccsum,abs(tline(k)));
   end;

   if (ocsum==ccsum)
      numGood = numGood+1;
      if (tline(2:6)=='GPGGA')
         numGoodPos = numGoodPos+1;
      end
   end

   tline = fgetl(fid);   
end

fclose(fid);

fprintf('Read %d lines of data\n', numRead);
fprintf('%d passed checksum.\n', numGood);
fprintf('%d good position fixes.\n', numGoodPos);


% Second read -- store GPGGA information
sod = zeros(numGoodPos,1);
lats = sod;
longs = sod;
hts = sod;
sats=sod;
fid = fopen(fname);
tline=fgetl(fid);

idx = 0;

while ischar(tline)
   
   % Extract the recorded checksum
   csumidx = findstr(tline,'*');
   ocsum = hex2dec(tline(end-1:end));
   
   % Comput expected checksum
   ccsum = 0;
   for k=2:(csumidx-1)
     ccsum = bitxor(ccsum,abs(tline(k)));
   end;

   if (ocsum==ccsum)
      if (tline(2:6)=='GPGGA')
         idx = idx+1;
         cidx = findstr(tline,',');
         sodStr = tline((cidx(1)+1):(cidx(2)-1));
         latStr = tline((cidx(2)+1):(cidx(3)-1));
         latDirStr = tline(cidx(3)+1);
         lonStr = tline((cidx(4)+1):(cidx(5)-1));
         lonDirStr = tline(cidx(5)+1);
         htStr = tline((cidx(9)+1):(cidx(10)-1));
         sats(idx)=str2num(tline((cidx(7)+1):(cidx(8)-1)));
         sod(idx) = str2num(sodStr(1:2))*3600 + str2num(sodStr(3:4))*60+str2num(sodStr(5:(length(sodStr)-4)));

         latSgn=1;
         if (latDirStr=='S')
            latSgn=-1;
         end
  
         lonSgn=1;
         if (lonDirStr=='W')
            lonSgn=-1;
         end
         
           
         lats(idx) = latSgn*(str2num(latStr(1:2))+str2num(latStr(3:length(latStr)))/60); % DDMM.MMMMM format
         %longs(idx) =
         %%lonSgn*(str2num(lonStr(1:2))+str2num(lonStr(3:length(lonStr)))/60);
         %%DDMM.MMMMM format 
         longs(idx) = lonSgn*(str2num(lonStr(1:3))+str2num(lonStr(4:length(lonStr)))/60); %DDDMM.MMMMM format
         hts(idx) = str2num(htStr);
      end
   end

   tline = fgetl(fid);   
end

fclose(fid);

arr_out      = zeros(length(lats), 5);
arr_out(:,1) = sod;
arr_out(:,2) = lats;
arr_out(:,3) = longs;
arr_out(:,4) = hts;
arr_out(:,5) = sats;
T = array2table(arr_out);
T.Properties.VariableNames(1:5) = {'Sec_of_day {s}','latitude {deg}','longitude {deg}','h_ellipsiod {m}','num_sats'};
filesplit = split(fname, ".");
filen = filesplit(1);
ext      = '.csv';
dir      = '/Users/williamfife/Projects/dev/SatNav/data/';
filename = strcat(dir, filen, ext);
writetable(T, filename{1,1});