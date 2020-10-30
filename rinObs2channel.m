function [channel] = rinObs2channel(filename)
% rinObs2channel : Parse RINEX 2.11 observation file and convert it to a
%                  partially-filled channel.mat file of the format described
%                  in channeldef.txt.
%               
% INPUTS
% filename --- Name of RINEX observation file including path
% 
% OUTPUTS
% channel --- Array in channel.mat format containing various data for each
%             observed satellite at each epoch:
%
%   Col 1 ---- RRT week number (LEFT AS ZERO)
%   Col 2 ---- RRT seconds of week (LEFT AS ZERO)
%   Col 3 ---- ORT week number
%   Col 4 ---- ORT whole seconds of week
%   Col 5 ---- ORT fractional seconds
%   Col 6 ---- Apparent Doppler frequency in Hz (LEFT AS ZERO)
%   Col 7 ---- Beat carrier phase in cycles
%   Col 8 ---- Pseudorange in meters
%   Col 9 ---- Carrier-to-noise ratio (C/N0) in dB-Hz
%   Col 10 --- Flag indicating whether (1) or not (0) the pseudorange and
%              carrier phase measurements are valid 
%   Col 11 --- A flag indicating whether (1) or not (0) an anomaly has been
%              detected in PLL's phase tracking 
%   Col 12 --- Channel status indicator 
%   Col 13 --- Signal type 
%   Col 14 --- Transmitter identification number (TXID)
%
%+------------------------------------------------------------------------+
% References:
%
% Author: Deep Mukherji, Matthew Murrian
%+========================================================================+

warning('off','MATLAB:nonIntegerTruncatedInConversionToChar');
ortWkCol = 3; ortWholeSecCol = 4; ortFracSecCol = 5;
thetahatCol = 7; prCol = 8; C_N0Col = 9; prValidFlagCol = 10;
statusCol = 12; signalTypeCol = 13; txIdCol = 14;
GPS_L1_CA = 0; GPS_L2_CL = 2; STATUS_DATA_LOCK = 6;

rinex_stream = fopen(filename);

if rinex_stream == -1
  error(['rinObs2channel: Unable to open ' filename]);
else
  disp(['Parsing observation data in ' filename '...']);
end

str_in = fgetllen(rinex_stream);

% Pre-allocate and initialize various items
gps_zero_time = datenum(1980,1,6,00,00,00);
epoch_satellites_type = char(zeros(100,1));
epoch_satellites_num = zeros(100,1);
observation_types = char(zeros(100,2));
interval = 1; % Default assumption for preallocation below.
kk = 1;

% Read information from RINEX header
while(cast(str_in(1),'uint8') ~= 0)
  % Check for version tag
  in_header = strfind(str_in,'RINEX VERSION / TYPE');
  
  % Read version tag
  if in_header == 61
    format_version = str2double(str_in(1:9));
  end
  
  % Check for observation types tag
  in_header = strfind(str_in,'# / TYPES OF OBSERV');
  
  % Read observation types tag
  if in_header == 61
    % Read number of observation types to expect
    number_observation_types = str2double(str_in(1:6));
    % Trim out only the observation types into a sub-string
    str_in = str_in(7:60);
    
    % Iterate to load all observation types
    for ii=1:number_observation_types
      % Load a continuation line if on observation type 10, 19, etc
      if (mod(ii,9) == 1) && (ii > 1)
        str_in = fgetllen(rinex_stream);
        
        in_header = strfind(str_in,'# / TYPES OF OBSERV');
        
        % Exit this inner loop if line wasn't a continuation
        if in_header ~= 61
          continue;
        end
        
        % Trim waste characters from beginning
        str_in = str_in(7:60);
      end
      
      % Extract observation type
      observation_types(ii,:) = str_in(5:6);
      % Trim observation type for next interation
      str_in = str_in(7:end);
    end
  end
  
  % Check for interval tag
  in_header = strfind(str_in,'INTERVAL');
  
  if in_header == 61
    interval = sscanf(str_in,'%10f');
  end
  
  % Check for first observation time tag
  in_header = strfind(str_in,'TIME OF FIRST OBS');
  
  % Load first observation time
  if in_header == 61
    first_obs.year = str2double(str_in(1:6));
    first_obs.month = str2double(str_in(7:12));
    first_obs.day = str2double(str_in(13:18));
    first_obs.hour = str2double(str_in(19:24));
    first_obs.min = str2double(str_in(25:30));
    first_obs.sec = str2double(str_in(31:43));
  end
  
  % Check for end of header
  in_header = strfind(str_in,'END OF HEADER');
  
  % If end of header, read the next line and exit outside while(...) loop
  if in_header == 61
    str_in = fgetllen(rinex_stream);
    break;
  end
  
  % If we made it here, read a new line and try again
  str_in = fgetllen(rinex_stream);
end

% Pre-initialize observation indexes to one past the actual list. This
% observation value will always refer to a 'Null' observation for
% short-circuiting.
idxC1 = number_observation_types + 1; idxL1 = number_observation_types + 1;
idxL2 = number_observation_types + 1; 
idxP2 = number_observation_types + 1; idxS1 = number_observation_types + 1;
idxS2 = number_observation_types + 1;

% Pre-initialize observation storage vectors to NaN. This should not need
% to be 'cleared' every iteration since blank observations convert to NaN.
current_observation = NaN(number_observation_types+1,1);

% Assign index values to observation types for ease-of-use
for ll=1:number_observation_types
  switch observation_types(ll,:)
    case 'C1'
      idxC1 = ll;
    case 'L1'
      idxL1 = ll;
    case 'L2'
      idxL2 = ll;
    case 'P2'
      idxP2 = ll;
    case 'S1'
      idxS1 = ll;
    case 'S2'
      idxS2 = ll;
  end
end

% Ensure we're reading the correct version of RINEX file
if format_version ~= 2.11
  error(['rinObs2channel: RINEX version 2.11 required. File reports version ' ...
         format_version]);
end

epoch_hour_last = -1;

% Assign observation values which are not expected to change
channel_epoch_template = zeros(1,14);
channel_epoch_template(1,prValidFlagCol) = 1;
channel_epoch_template(1,statusCol) = STATUS_DATA_LOCK;
channel_epoch_template(1,signalTypeCol) = GPS_L1_CA;

% Pre-allocate space for 24 hours of samples, for 20 satellites, 2
% frequencies per, every (interval) seconds.
channel = zeros(ceil(24*60*60*20*2/interval),14);

% Begin to ingest observations
while(cast(str_in(1),'uint8') ~= 0)
  % First line of each Epoch should contain the following fields.
  epoch_tag = sscanf(str_in,' %2d %2d %2d %2d %2d%11f  %1d%3d');
  
  % If any field is lacking, it's because we need to skip past header
  % information and comment lines.
  if length(epoch_tag) ~= 8
    str_in = fgetllen(rinex_stream);
    continue;
  elseif (epoch_tag(7) > 1) % Indicates header information likely follows
    str_in = fgetllen(rinex_stream);
    continue;
  end

  epoch_year = epoch_tag(1);
  epoch_month = epoch_tag(2);
  epoch_day = epoch_tag(3);
  epoch_hour = epoch_tag(4);
  epoch_min = epoch_tag(5);
  epoch_sec = epoch_tag(6);
  epoch_number_satellites = epoch_tag(8);
  
  epoch_sec_whole = floor(epoch_sec);
  epoch_sec_frac = epoch_sec - epoch_sec_whole;

  % Provide some indication of action so that a lock-up is not suspected
  if epoch_hour_last ~= epoch_hour
    epoch_hour_last = epoch_hour;
    disp(['Processing hour ' num2str(epoch_hour)]);
  end
  
  % Generate a long version of the year that allows proper accounting if
  % the RINEX file passes over a year boundary.
  epoch_year_long = first_obs.year + (epoch_year - mod(first_obs.year,100));

  % Generate a scalar value referring to the number of whole seconds at the
  % current epoch since the GPS reference epoch.
  delta_n = datenum(epoch_year_long,epoch_month,epoch_day,epoch_hour,...
                    epoch_min,epoch_sec_whole) - gps_zero_time;
  
  % Produce familiar values of gpsWeek and gpsSec for the current epoch.
  gpsWeek = floor(delta_n/7);
  gpsSecWhole = round((delta_n - gpsWeek*7)*86400);
  gpsSecFrac  = epoch_sec_frac;
  
  % Advance to satellite list for epoch
  str_loc = 33;
  
  % Update the observation values which change only every epoch
  channel_epoch_template(1,ortWkCol) = gpsWeek;
  channel_epoch_template(1,ortWholeSecCol) = gpsSecWhole;
  channel_epoch_template(1,ortFracSecCol) = gpsSecFrac;
  
  % Iterate to load all observation types
  for ii=1:epoch_number_satellites
    % Load a continuation line if over a multiple of 12 satellites
    if (mod(ii,12) == 1) && (ii > 1)
      str_in = fgetllen(rinex_stream);
      str_loc = 33;
    end
    
    % Extract satellite records specifier for current epoch
    if str_in(str_loc) ~= ' '
      epoch_satellites_type(ii,1) = str_in(str_loc);
    else
      % Express the default for a blank explicitly as G
      epoch_satellites_type(ii,1) = 'G';
    end
    
    epoch_satellites_num(ii,1) = sscanf(str_in(str_loc + 1:str_loc + 2),'%2d');
    
    % Trim observation type for next interation
    str_loc = str_loc + 3;
  end
  
  % Iterate to load all observations for this epoch
  for ii=1:epoch_number_satellites
    str_in = fgetllen(rinex_stream);
    str_loc = 1;
    
    for jj=1:number_observation_types
      if (mod(jj,5) == 1) && (jj > 1)
        str_in = fgetllen(rinex_stream);
        str_loc = 1;
      end
      
      parse_in = sscanf(str_in(str_loc:str_loc+13),'%14f%');
      
      if ~isempty(parse_in)
        current_observation(jj) = parse_in(1);
      else
        current_observation(jj) = NaN;
      end
      
      str_loc = str_loc + 16;
    end
    
    % Only process the data if system is GPS
    if (epoch_satellites_type(ii) == 'G' && ...
        ~(isnan(current_observation(idxL1)) || isnan(current_observation(idxC1))))
      % Set txID value for observation template
      channel_epoch_template(1,txIdCol) = epoch_satellites_num(ii,1);
      
      % Record observables for L1 observation.
      channel_epoch_template(1,thetahatCol) = current_observation(idxL1);
      channel_epoch_template(1,prCol) = current_observation(idxC1);
      if(isnan(current_observation(idxS1)))
        channel_epoch_template(1,C_N0Col) = 0;
      else
        channel_epoch_template(1,C_N0Col) = current_observation(idxS1);
      end
      
      % Copy observation template to current observation. Observation
      % values were copied to the template initially in order to
      % capitalize on data locality. The template and current
      % observations are likely on stack where channel is most 
      % certaintly on the heap.
      channel(kk,:) = channel_epoch_template;

      % Record an L2 observation if (inclusive) either L2 or P2
      % observation exists for the present epoch
      if ~isnan(current_observation(idxL2)) || ~isnan(current_observation(idxP2))
        kk = kk + 1;
        
        % Copy L2 observations to the template. Similar
        % implications for data locality as discussed above.
        channel_epoch_template(1,thetahatCol) = current_observation(idxL2);
        channel_epoch_template(1,prCol) = current_observation(idxP2);
        channel_epoch_template(1,C_N0Col) = current_observation(idxS2);
        
        % Copy observation template to current observation.
        channel(kk,:) = channel_epoch_template;
        
        % This is assigned directly to channel so I don't have to
        % toggle it back in the template.
        channel(kk,signalTypeCol) = GPS_L2_CL;
      end
      
      kk = kk + 1;
    end        
  end
  
  str_in = fgetllen(rinex_stream);
end

% Trim the channel matrix prior to returning.
channel = channel(1:kk-1,:);

% Pad the retrieved line to (length_in) number of characters (i.e., 80 for the
% RINEX format). This eliminates a number of conditionals above in the case
% that the last observation of a line is blank (thereby shortening the line
% and risking a SIGSEGV during parsing).
%
% In the case that an EOF is reached, this function will return an 80-length
% character array where the first indexed position is equal to a 0 (binary 0,
% not char '0'). Catch this condition by: cast(ret_str(1),'uint8') == 0. A
% warning will also be issued during this truncation of a -1 (return from
% fgets at EOF) to 0. Despite this unclean method for the single EOF event, it
% is a much faster way to handle the 99.999% of non-EOF calls to this
% function.
function [ret_str] = fgetllen(file_stream)
line_in = fgets(file_stream);
ret_str = repmat(' ',1,80);
ret_str(1:length(line_in)) = line_in(1:length(line_in));
return;

warning('on','MATLAB:nonIntegerTruncatedInConversionToChar');          

