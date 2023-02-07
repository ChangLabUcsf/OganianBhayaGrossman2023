function f0=addF0(sent, corpus, datapath)
    sentName=sent.name;  
    
    soundir='stim_info/';
    
    fid = fopen([datapath '/' soundir corpus 'Sounds/' ...
        sentName '_pitch_int.txt']);    
    pitchtext = textscan(fid, '%s%s%s%s', 'Delimiter', '	');
    fclose(fid);
    pitchtext = cell2mat(cellfun(@(x) str2double(x), ...
        pitchtext{3}, 'UniformOutput', false))';

    % clip last row
    pitchtext = pitchtext (1:end-1);

    % make all NaNs into 0s
    pitchtext(isnan(pitchtext ))=0;

    % pad with zeros
    befpad = floor((length(sent.stress)-length(pitchtext))/2);
    aftpad = ceil((length(sent.stress)-length(pitchtext))/2);
    f0=[zeros(1, befpad) pitchtext  zeros(1, aftpad)];
end