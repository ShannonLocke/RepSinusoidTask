function [] = run_dataAnalysisPipeline()
% This script will run the entire data analysis pipeline for the 
% Repeated-Sinusoids task in the Metacognition for Eye Tracking project. 
% (Locke, Goettker, Gegenfurtner, & Mamassian, 2021)
% 
% Data: https://osf.io/j9asn/
% Code: https://github.com/ShannonLocke/RepSinusoidTask
%
% Created by SML Dec 2021
% License: CC-By Attribution 4.0 International

%% PREAMBLE:

% Switchboard:
testDataYN = false; % analyse test data or pilot data

% Ensure correct folder structure:
% <== DO THIS!!!

% Find the participant IDs of available data:
if testDataYN
    dataFilePattern = '../data/s*';
    resDir = 'results/';
else
    dataFilePattern = '../data_pilot/s*';
    resDir = 'results_pilot/';
end
sOnFile = dir(dataFilePattern);
nSs = length(sOnFile);
for nn = 1:nSs
    getFolderName = sOnFile(nn).name;
    all_sID(nn) = str2num(getFolderName(2:end));
end

% If option selected, limit to only unprocessed Ss data:
onlyNewDataYN = false;
if onlyNewDataYN
    dataFilePattern = 'output_data/processed/processedEyeData_s*.mat';
    sOnFile = dir(dataFilePattern);
    nSs = length(sOnFile);
    for nn = 1:nSs
        sel_sID = sOnFile(nn).name;
        sel_sID = str2num(sel_sID(13:15));
        all_sID(all_sID == sel_sID) = [];
    end
end
% <== ISSUE: WON'T PROCESSING ONLY NEW AFFECT THE DATA SUMMARY?

%% Extract the raw data:
disp('1/7. EXTRACTING THE RAW DATA...')
prep_extractRawData(all_sID, testDataYN, resDir)

%% Process the raw data:
disp('2/7. PROCESSING THE RAW DATA...')
prep_processData(all_sID, resDir)

%% Compute data summary:
disp('3/7. COMPUTING DATA SUMMARY...')
prep_summaryData(all_sID, testDataYN, resDir)

%% H1 analysis:
disp('4/7. TESTING HYPOTHESIS 1 NOW...')
H1_metacogSensitivity(resDir)

%% H2 analysis:
disp('5/7. TESTING HYPOTHESIS 2 NOW...')
H2_repeatedTrajectories(resDir)

%% H3 analysis:
disp('6/7. TESTING HYPOTHESIS 3 NOW...')
H3_temporalMetacogSensitivity(resDir)

%% H4 analysis:
disp('7/7. TESTING HYPOTHESIS 4 NOW...')
H4_repeatedSessions(resDir)

end