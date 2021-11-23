% Re-calibrate if required:
if mod(trial,calibEvery) == 0
    EyelinkDoTrackerSetup(el);
end

% Drift correction:
Eyelink('command','set_idle_mode');
WaitSecs(0.05);
drift_success = EyelinkDoDriftCorrection(el, pos00(1), pos00(2), 1, 1);
if drift_success ~= 1; disp(drift_success); end
Eyelink('command','set_idle_mode');
WaitSecs(0.05);
Eyelink('StartRecording');

% Start recording:
Eyelink('command','mark_playback_start');

% Record trial number:
Eyelink('message', ['TrialID' num2str(trial)]);
WaitSecs(ISI/2);

% Record stimulus start:
Eyelink('message','TRIAL_STIMSTART');