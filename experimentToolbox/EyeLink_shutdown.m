% Send end-of-block event trigger:
Eyelink('message', 'Block_End');

% Close the file:
Eyelink('CloseFile');
disp(['Receiving data file: ' edfFileName])
status = Eyelink('ReceiveFile', [], edfFilePath, 1);
if status > 0
    disp(['ReceiveFile status: ' num2str(status)])
end

% Shutdown the EyeLink::
Eyelink('Shutdown');

