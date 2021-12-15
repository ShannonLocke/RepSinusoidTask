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

% Find the participant IDs of available data:
all_sID = [];

% Extract the raw data:
disp('...EXTRACTING THE RAW DATA...')
prep_extractRawData(all_sID)

% Process the raw data:
disp('...PROCESSING THE RAW DATA...')
prep_processData(all_sID)

% Compute data summary:
disp('...COMPUTING DATA SUMMARY...')
prep_summaryData(all_sID)

% H1 analysis:
disp('...TESTING HYPOTHESIS 1 NOW...')
H1_metacogSensitivity

% H2 analysis:
disp('...TESTING HYPOTHESIS 2 NOW...')
H2_repeatedTrajectories

% H3 analysis:
disp('...TESTING HYPOTHESIS 3 NOW...')
H3_temporalMetacognitiveSensitivity

% H4 analysis:
disp('...TESTING HYPOTHESIS 4 NOW...')
H4_repeatedSessions

end