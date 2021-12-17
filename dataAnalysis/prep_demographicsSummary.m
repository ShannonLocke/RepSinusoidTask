function [] = prep_demographicsSummary()
% This script will report the demographic information for the participants 
% of the Repeated-Sinusoids task (Locke, Goettker, Gegenfurtner, & Mamassian, 2021).
%
% Data: https://osf.io/j9asn/
% Code: https://github.com/ShannonLocke/RepSinusoidTask
%
% Created by SML Dec 2021
% License: CC-By Attribution 4.0 International

% Source file information:
dataPath = '../data/participantInfo/';
dataFilePattern = 'participantInfoCard_SPC_*.mat';
dataOnFile = dir([dataPath, dataFilePattern]);
nSs = length(dataOnFile);

% Get participant age:
age = NaN([1,nSs]);
lefthanded = NaN([1,nSs]);
% female = NaN([1,nSs]);
for nn = 1:nSs
    fname = [dataPath, dataOnFile(nn).name];
    load(fname,'participantInfo')
    age(nn) = str2double(participantInfo.Age);
    lefthanded(nn) = strcmp(participantInfo.Handedness,'L');
    % female(nn) = strcmp(participantInfo.Gender,'F');
end

% Compute mean and range:
disp(['There are ' num2str(length(age)) ' participants in this study.'])
disp(['The mean age is ' num2str(mean(age)) ' years old.'])
disp(['The age range is ' num2str(min(age)) '--' num2str(max(age)) ' years old.'])
% disp(['The number of females participants is: ' num2str(sum(female))])
disp(['The number of left-handed participants is: ' num2str(sum(lefthanded))])

end