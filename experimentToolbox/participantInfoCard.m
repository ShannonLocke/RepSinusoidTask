function [participantInfo,saveFile] = participantInfoCard(expName,fieldNames,prompt,cardToScrubID)
%
% Created by SML March 2021

% Defaults:
if nargin < 4
    cardToScrubID = [];
    if nargin < 3
        prompt = {'First Name (no spaces)','Last Name (no spaces)','Age',...
            'Vision (1=normal, 2=corrected)','Handedness ("L", "R")','ID number (3 digit)'};
        if nargin < 2
            fieldNames = {'FirstName','LastName','Age','Vision','Handedness','sID'};
            if nargin < 1
                expName = 'myExp';
            end
        end
    end
end
if isempty(expName); expName = 'myExp'; end
if isempty(fieldNames); fieldNames = {'FirstName','LastName','Age','Vision','Handedness','sID'}; end
if isempty(prompt); prompt = {'First Name (no spaces)','Last Name (no spaces)','Age', 'Vision (1=normal, 2=corrected)','Handedness ("L", "R")','ID number (3 digit)'}; end

% Scrub participant ID from card when experiment complete:
if ~isempty(cardToScrubID)
    participantInfo = cardToScrubID;
    participantInfo.sID = [];
    participantInfo.finishExp = getTimeStamp;
    subjectSpecificpattern = [upper(participantInfo.LastName) '_' lower(participantInfo.FirstName) '_' participantInfo.startExp];
    fname = ['data/participantInfoCard_' expName '_' subjectSpecificpattern '.mat'];
    save(fname,'participantInfo');
    saveFile = [];
    return
end

% Make data directory if none exists:
% if ~isfolder('data'); mkdir('data'); end

% Check for exisiting participant info cards:
dataFilePattern = ['data/participantInfoCard_' expName '_*.mat'];
cardsOnFile = dir(dataFilePattern);
if ~isempty(cardsOnFile) % participant info cars found
    makeNewEntryYN = false;
    if length(cardsOnFile) == 1 % only 1 card, load it
        fname = cardsOnFile(1).name;
        load(['data/' fname],'participantInfo')
    else % more than 1, select between cards
        msg_multCards = 'Multiple participant info cards found! Select which you would like to load or click cancel to make new.';
        ff = msgbox(msg_multCards,'','help');
        uiwait(ff);
        fname = uigetfile('data/participantInfoCard_*.mat');
        if fname
            load(['data/' fname],'participantInfo')
        else
            makeNewEntryYN = true;
        end
    end
else % no cards, must make new card
    makeNewEntryYN = true;
end

% Get participant info from input dialogue:
dim = [1,25];
if makeNewEntryYN
    dlgtitle = 'Enter new:';
    participantInfo = inputdlg(prompt,dlgtitle,dim);
    participantInfo = cell2struct(participantInfo,fieldNames);
else
    dlgtitle = 'Confirm info:';
    definput = struct2cell(participantInfo); definput(end) = [];
    participantInfo_new = inputdlg(prompt,dlgtitle,dim,definput);
    if ~all(cellfun(@isequal, definput, participantInfo_new)) % participant did not confirmed info (i.e., made any changes)
        participantInfo = cell2struct(participantInfo_new,fieldNames);
        makeNewEntryYN = true;
    end
end

% Save new entry if one was made:
if makeNewEntryYN
    participantInfo.startExp = getTimeStamp;
    subjectSpecificpattern = [upper(participantInfo.LastName) '_' lower(participantInfo.FirstName) '_' participantInfo.startExp];
    fname = ['data/participantInfoCard_' expName '_' subjectSpecificpattern '.mat'];
    save(fname,'participantInfo');
end

% Determine the filename for storing participant data:
saveFile = ['data/results_' expName '_s' participantInfo.sID '_' participantInfo.startExp '.mat'];

end