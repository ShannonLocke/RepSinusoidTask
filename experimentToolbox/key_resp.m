function [key] = key_resp(KbIdx)

while KbCheck(KbIdx); end % Wait until all keys are released.

keyIsDown = false;

while ~keyIsDown
    % Check the state of the keyboard.
    [ keyIsDown, ~, keyCode ] = KbCheck(KbIdx);
    
    % If the user is pressing a key, then display its code number and name.
    if keyIsDown
        
        % Note that we use find(keyCode) because keyCode is an array.
        % See 'help KbCheck'
        %fprintf('You pressed key %i which is %s\n', find(keyCode), KbName(keyCode));
        key = find(keyCode);
        % If the user holds down a key, KbCheck will report multiple events.
        % To condense multiple 'keyDown' events into a single event, we wait until all
        % keys have been released.
        while KbCheck(-1); end
    end
    
end

key = key(1); % select only the first key that was pressed instances where more than one key was pressed.

end