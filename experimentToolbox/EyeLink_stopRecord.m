% Stop recording:
Eyelink('Message', 'TRIAL_END');
WaitSecs(0.1);
Eyelink('StopRecording');
error = Eyelink('CheckRecording');
fprintf('Stop Recording: %d; ',error);
Eyelink('SetOfflineMode');
WaitSecs(0.1);