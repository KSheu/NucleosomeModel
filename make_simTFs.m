%% ------ simulated TF vectors - osc vs non-osc
START_TIME =0;
END_TIME = 480;
% END_TIME = 900;
length_pulse = [3.75 .75];
num_pulses = [1 5];
for p = 1:length(length_pulse)
    
    pulse_on = length_pulse(p); % (hrs)
    pulse_off = .75; % (hrs)
    pulse_num = num_pulses(p); % (hrs)
    tf = zeros(1, END_TIME-START_TIME);
    dt = 1;

    for i = 1:pulse_num
        start_idx = (i-1)*(60*(pulse_on+pulse_off))/dt+60*pulse_off/dt;
        tf(round(start_idx):round(start_idx+60*pulse_on/dt)) = 1;
    end
    plot(tf);
    if p == 1
        save('nfkb_nonoscillatory', 'tf');
    else
        save('nfkb_oscillatory', 'tf');
    end
end


%% ------ simulated TF vectors - high amp. vs low amp.
START_TIME =0;
END_TIME = 480;

length_pulse = [3.75 .75];
num_pulses = [1 5];
for p = 1:length(length_pulse)
    p=1;
    pulse_on = length_pulse(p); % (hrs)
    pulse_off = .75; % (hrs)
    pulse_num = num_pulses(p); % (hrs)
    tf = zeros(1, END_TIME-START_TIME);
    dt = 1;

    for i = 1:pulse_num
        start_idx = (i-1)*(60*(pulse_on+pulse_off))/dt+60*pulse_off/dt;
        tf(round(start_idx):round(start_idx+60*pulse_on/dt)) = 10;
    end
    plot(tf);
    ylim([0 11]);
    if p == 1
        save('nfkb_nonoscillatory_hiamp', 'tf');
    else
        save('nfkb_oscillatory_hiamp', 'tf');
    end
end

%% ---- simulated tf vectors - total activity
START_TIME =0;
% END_TIME = 480;
END_TIME = 900;

length_pulse = .75;
num_pulses = 10;
for p = 1:length(length_pulse)
    
    pulse_on = length_pulse(p); % (hrs)
    pulse_off = .75; % (hrs)
    pulse_num = num_pulses(p); % (hrs)
    tf = zeros(1, END_TIME-START_TIME);
    dt = 1;

    for i = 1:pulse_num
        start_idx = (i-1)*(60*(pulse_on+pulse_off))/dt+60*pulse_off/dt;
        tf(round(start_idx):round(start_idx+60*pulse_on/dt)) = 1;
    end
    plot(tf);
    ylim([0 1.2]);
    save('nfkb_oscillatory_2xtotalactivity', 'tf');
    
end
%% ---- simulated tf vectors - total activity long persistent
START_TIME =0;
% END_TIME = 480;
END_TIME = 900;

length_pulse = 15;
num_pulses = 1;
for p = 1:length(length_pulse)
    
    pulse_on = length_pulse(p); % (hrs)
    pulse_off = .75; % (hrs)
    pulse_num = num_pulses(p); % (hrs)
    tf = zeros(1, END_TIME-START_TIME);
    dt = 1;

    for i = 1:pulse_num
        start_idx = (i-1)*(60*(pulse_on+pulse_off))/dt+60*pulse_off/dt;
        tf(round(start_idx):round(start_idx+60*pulse_on/dt)) = 0.5;
    end
    plot(tf);
    ylim([0 1.2]);
    save('nfkb_persistent_2xtotalactivity', 'tf');
    
end

