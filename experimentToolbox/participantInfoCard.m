function [participantInfo,saveFile] = participantInfoCard(expName,fieldNames,prompt,cardToScrubID)
%
% Created by SML March 2021, Updated by SML November 2021

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

% Define data path (and create one if it doesn't exist):
dataPath = 'data/participantInfo/';
if ~isfolder(dataPath); mkdir(dataPath); end

% Scrub participant ID from card when experiment complete:
if ~isempty(cardToScrubID)
    participantInfo = cardToScrubID;
    participantInfo.sID = [];
    getTime = getTimeStamp;
    participantInfo.finishExp = getTime(1:10); % only store date (exact time too identifying)
    subjectSpecificpattern = [upper(participantInfo.LastName) '_' lower(participantInfo.FirstName)];
    fname = ['participantInfoCard_' expName '_' subjectSpecificpattern '.mat'];
    save([dataPath, fname],'participantInfo');
    saveFile = [];
    return
end

% Check for exisiting participant info cards:
dataFilePattern = [dataPath 'participantInfoCard_' expName '_*.mat'];
cardsOnFile = dir(dataFilePattern);
if ~isempty(cardsOnFile) % participant info cards found
    makeNewEntryYN = false;
    if ~isempty(cardsOnFile) % at least 1 card, select between cards
        msg_multCards = 'Participant info cards found! Select which you would like to load or click cancel to make a new card.';
        ff = msgbox(msg_multCards,'','help');
        uiwait(ff);
        fname = uigetfile([dataPath 'participantInfoCard_*.mat']);
        if fname
            load([dataPath fname],'participantInfo')
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
    if isempty(participantInfo); error('No participant was specified.'); end
    participantInfo = cell2struct(participantInfo,fieldNames);
else
    dlgtitle = 'Confirm info:';
    definput = struct2cell(participantInfo); definput(end) = [];
    participantInfo_new = inputdlg(prompt,dlgtitle,dim,definput);
    if isempty(participantInfo_new) % experimenter selected cancel, will make new participant
        makeNewEntryYN = true;
        msg_multCards = 'You selected cancel. Please provide the details for the new participant or press cancel again to quit the setup.';
        ff = msgbox(msg_multCards,'','help');
        uiwait(ff);
        dlgtitle = 'Enter new:';
        participantInfo = inputdlg(prompt,dlgtitle,dim);
        if isempty(participantInfo); error('No participant was specified.'); end
        participantInfo = cell2struct(participantInfo,fieldNames);
    elseif ~all(cellfun(@isequal, definput, participantInfo_new)) % experimenter made any changes to the info
        participantInfo = cell2struct(participantInfo_new,fieldNames);
        msg_multCards = 'You have made changes to this entry. Creating a new participant.';
        ff = msgbox(msg_multCards,'','help');
        uiwait(ff);
        makeNewEntryYN = true;
    end
end

% Save new entry if one was made:
if makeNewEntryYN
    getTime = getTimeStamp;
    participantInfo.startExp = getTime(1:10); % only store date (exact time too identifying)
    subjectSpecificpattern = [upper(participantInfo.LastName) '_' lower(participantInfo.FirstName)];
    fname = ['participantInfoCard_' expName '_' subjectSpecificpattern '.mat'];
    save([dataPath, fname],'participantInfo');
end

% Determine the filename for storing participant data:
dataPath = ['data/s' num2str(participantInfo.sID) '/'];
if ~isfolder(dataPath); mkdir(dataPath); end
saveFile = [dataPath 'results_' expName '_s' participantInfo.sID '_' participantInfo.startExp '.mat'];

end