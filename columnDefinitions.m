% columnDefinitions.m
%
% See definition files in trunk/doc/


% channel data
rrtWkCol = 1;
rrtSecCol = 2;
ortWkCol = 3;
ortWholeSecCol = 4;
ortFracSecCol = 5;
fpllCol = 6;
thetahatCol = 7;
prCol = 8;
C_N0Col = 9;
prValidFlagCol = 10;
phaseErrorFlagCol = 11;
statusCol = 12;
signalTypeCol = 13;
txIdCol = 14;

% Signal Types
GPS_L1_CA = 0;
GPS_L2_CM = 1;
GPS_L2_CL = 2;
GPS_L2_CLM = 3;
GPS_L1_ALT = 4;
CDMA_UHF_PILOT_I = 5;
CDMA_UHF_PILOT_Q = 6;
CDMA_UHF_PILOT_IQ = 9;
WAAS_L1_I = 13;
WAAS_L1_I_ALT1 = 14;
% Channel status 
STATUS_NULL = 0;
STATUS_ALLOCATED = 1;
STATUS_ACQUIRED = 2;
STATUS_SYMBOL_LOCK = 3;
STATUS_FREQ_LOCK = 4;
STATUS_PHASE_LOCK = 5; 
STATUS_DATA_LOCK = 6;

% Aliases 
svIdCol = 14;

% navsol data
ns_ortWkCol = 1; 
ns_ortWholeSecCol = 2;
ns_ortFracSecCol = 3;
ecefXCol = 4;
ecefYCol = 5;
ecefZCol = 6;
deltRxCol = 7;
ecefXdotCol = 8;
ecefYdotCol = 9;
ecefZdotCol = 10;
deltRxdotCol = 11;
solFlagCol = 12;

% iq data 
iq_rrtWkCol = 1;
iq_rrtSecCol = 2;
iq_ortWkCol = 3;
iq_ortWholeSecCol = 4;
iq_ortFracSecCol = 5;
iq_beatCarrPhaseCol = 6;
iq_iCol = 7;
iq_qCol = 8;
iq_dbCol = 9;
iq_signalTypeCol = 10;
iq_txIdCol = 11;


% scint data
sc_ortWkCol = 1;
sc_ortWholeSecCol = 2;
sc_ortFracSecCol = 3;
sc_measIntCol = 4;
sc_s4WholeCol = 5;
sc_sigPhiWholeCol = 8;
sc_tau0Col = 11;
sc_sprCol = 12;
sc_refIndicatorCol = 13;
sc_signalTypeCol = 14;
sc_txIdCol = 15;


% tx data
tx_ortWkCol = 1;
tx_ortWholeSecCol = 2;
tx_ortFracSecCol = 3;
tx_azCol = 4;
tx_elCol = 5;
tx_healthCol = 6;
tx_systemCol = 7;
tx_txIdCol = 8;

% uBlox columns
uBloxSvIdCol = 6;
uBloxRecTimeSecondsCol = 1;
uBloxPrCol = 3;
uBloxThetahatCol = 2;
uBloxSolMillisecondsCol = 1;
uBloxSolNanosecondsCol = 2;
uBloxWeekNumberCol = 8;
uBloxCN0Col = 4;



