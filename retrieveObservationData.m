function [channel] = ...
    retrieveObservationData(gpsWeek,gpsSec, utctime, markerId,obsDirectory)
% retrieveObservationData : Retrieve from the NOAA (NGS) CORS data repository
%                           the GNSS observation data (i.e., receiver
%                           observables) for the specified CORS reference
%                           site.  A day's worth of data are retrieved which
%                           span the time interval covering the specified GPS
%                           time.  Only retrieval and parsing of uncompressed
%                           RINEX 2.11 files is currently supported (those
%                           having filenames of the form {ssss}{ddd}0.{yy}o).
%                           NOAA (NGS) deletes these uncompressed files after
%                           one year, leaving only the Hatanaka-compressed
%                           RINEX 2.11 files (those having filenames of the
%                           form {ssss}{ddd}0.{yy}d), which can be
%                           uncompressed using the teqc utility.  For more
%                           information on the NOAA (NGS) data and sites, see
%                           http://www.ngs.noaa.gov/CORS/
%
%
% INPUTS
%
% gpsWeek ------------- GPS week number.
%
% gpsSec -------------- GPS second of week.
%
% markerId ------------ 4-alphanumeric-character marker identifier (e.g.,
%                       txau).
%
% obsDirectory -------- Path to directory where the RINEX observation file
%                       retrieved from the CORS directory and the converted
%                       channel.mat file will be stored.  Leave blank to store
%                       in navsol/obsFiles (if on path) or in working
%                       directory (otherwise).
%
%
% OUTPUTS
%
% channel ------------- Matrix formatted according to the format described
%                       in channeldef.txt.  Some entries of the matrix are
%                       left blank because there is no corresponding
%                       quantity in the retrieved RINEX observation file.
%
%+------------------------------------------------------------------------------+
% References:
%
%
% Author: Todd Humphreys
%+==============================================================================+

% Choose default obsDirectory if none specified
if(nargin == 3)
    % if navsol is in path, default obsDirectory is 'navsol/obsFiles'
    % otherwise, use current directory (don't use './obsFiles' because
    % it may not exist)
    obsDirectory = '.';
    p = strsplit(path,pathsep);
    for ii = 1:length(p)
        if(~isempty(strfind(p{ii},'navsol')))
            obsDirectory = fullfile(p{ii},'obsFiles');
            break;
        end
    end
end

% Get UTC time corresponding to gpsWeek and gpsSec.
tFps = datenum(utctime, 'yyyy-mm-dd HH:MM:SS');
dv = datevec(tFps);
tFpsYear = dv(1); tFpsHour = dv(4); tFpsMin = dv(5);
t0Year = datenum(tFpsYear,1,1,0,0,0);
doy = floor(tFps - t0Year) + 1;
markerId = lower(markerId);
fn = sprintf('%4s%03d0.%02do',markerId,doy,mod(tFpsYear,100));
s = sprintf('ftp://www.ngs.noaa.gov/cors/rinex/%4d/%03d/%4s/',...
    tFpsYear,doy,markerId);
s = [s fn];
% Download the file if not already present
fid = fopen([obsDirectory '/' fn]);
if(fid == -1)
    try
        disp(['Retrieving RINEX observation file ' fn ' from NGS ftp site ...']);
        fnn = [obsDirectory '/' fn '.gz'];
        syscmd = sprintf('cd %s && wget %s.gz', obsDirectory, s);
        system(syscmd);
        gunzip(fnn);
        delete(fnn);
        disp('Observation file successfully retrieved.');
    catch mException
        try
            urlwrite([s '.gz'],fnn);
            gunzip(fnn);
            delete(fnn);
            disp('Observation file successfully retrieved.');
        catch mException
            disp('');
            disp('Failed to retrieve observation file.');
            disp(['Error description: ' mException.identifier]);
            disp('You may wish to download the file manually; it can be found at:');
            disp([s '.gz']);
            disp(['After downloading, unzip the file using the Matlab gunzip command, ' ...
                'then place in the directory ' obsDirectory ' and run your code again.']);
            error(mException.message);
        end
    end
else
    fclose(fid);
end
disp('Converting RINEX 2.11 data to channel.mat format ...');
%fn = 'txau280a.19o';
fname = [obsDirectory '/' fn];
[channel] = rinObs2channel(fname);
channel = channel';
writematrix(channel', [obsDirectory 'channel.csv']);
save([obsDirectory '/channel.mat'], 'channel');
disp(['Data successfully converted and saved to ' obsDirectory '/channel.mat']);
channel = channel';

% Create and store a dummy navsol.mat file with zeros for the receiver clock
% offset (column 7) for use in CDGNSS processing.  Setting c*deltR to zero
% assumes that the RINEX CORS data receiver timestamps are very close (within
% a few ns) of true time, which seems to be the case in practice.
columnDefinitions;
wkInvalid = 9999;
iidum = find(channel(:,ortWkCol) ~= wkInvalid);
tRw = channel(iidum,ortWkCol);
tRws = channel(iidum,ortWholeSecCol);
tRfs = channel(iidum,ortFracSecCol);
tRs = tRws + tRfs;
[tRsU,II,JJ] = unique(tRs);
navsol = zeros(length(tRsU),12);
navsol(:,ns_ortWkCol) = tRw(II);
navsol(:,ns_ortWholeSecCol) = tRws(II);
navsol(:,ns_ortFracSecCol) = tRfs(II);
navsol = navsol';
writematrix(navsol', [obsDirectory 'navsol.csv']);
save([obsDirectory '/navsol.mat'], 'navsol');
disp(['Associated navsol.mat with c*deltR = 0 created and saved to ' ...
    obsDirectory '/navsol.mat']);