function config=read_parse_expt_config(configFN)
if ~isfile(configFN)
    fprintf('ERROR: configuration file %s does not exist.\n',configFN);
    config=[];
    return;
end

fid1=fopen(configFN,'r');

a=textscan(fid1,'%s','delimiter','\n');
a=a{1};
fclose(fid1);

config=struct;
% items={'SUBJECT_ID', 'SUBJECT_GENDER', 'SUBJECT_DOB', 'SHIFT_DIRECTION','SHIFT_RATIO',...
%        'MOUTH_MIC_DIST','SPL_TARGET','SPL_RANGE','VOWEL_LEN_TARG',...
%        'VOWEL_LEN_RANGE','PRE_WORDS','TRAIN_WORDS','TEST_WORDS',...
%     'PRE_REPS','PRACT1_REPS','PRACT2_REPS','START_REPS','RAMP_REPS',...
%     'STAY_REPS','STAY2_REPS','END_REPS','TEST1_REPS','TEST2_REPS','TEST3_REPS',...
%     'SAMPLING_RATE','FRAME_SIZE','DOWNSAMP_FACT','TRIAL_LEN',...
%     'NATURAL_REPS','NATURAL_WORDS'};
items={ 'SUBJECT_ID', 'SUBJECT_GENDER', 'SUBJECT_DOB', 'SUBJECT_GROUP', ...
        'DATA_DIR', ...
        'BACKGROUND_COLOR', 'TEXT_COLOR', ...
        'PITCH_SHIFT_DOWN', 'PITCH_SHIFT_UP', ...
        'MOUTH_MIC_DIST', 'SPL_TARGET', 'SPL_RANGE', ...
        'VOWEL_LEN_TARG', 'VOWEL_LEN_RANGE', ...
        'USE_SPL_TARGET', 'USE_VOWEL_LEN_TARGET', ...
        'NUM_RUNS', 'TRIALS_PER_RUN', ...
        'MRI_TR', 'MRI_TA', 'STIM_DELAY', 'REC_DELAY', 'REC_LEN', ...
        'SAMPLING_RATE', 'FRAME_SIZE', 'DOWNSAMP_FACT', ...
        'MASK_NOISE_LV'};

for i1=1:numel(items)
    item=items{i1};
    
%     if isequal(item,'NATURAL_WORDS') || isequal(item,'NATURAL_REPS')
%         pause(0);
%     end
       
    bFound=0;
    for j1=1:numel(a)
        if ~isempty(strfind(a{j1},item))
            bFound=1;
            break;
        end
    end

    dba=deblank(a{j1});
    
    if bFound==1 && ~isequal(dba(1),'%')
        str=a{j1};
        if iscell(str)
            str=str{1};
        end
        str=strrep(str,item,'');
        if ~isempty(strfind(str,'%'))
            idxp=strfind(str,'%');
            str=str(1:idxp(1)-1);
        end
    else
        fprintf('WARNING: item %s not specified.\n',item);
        config.(item)=[];
        continue;
    end
    
    if isequal(item,'SUBJECT_ID') || isequal(item, 'SUBJECT_GROUP') 
        str=strtrim(str);
        config.(item)=str;
    elseif isequal(item,'SUBJECT_GENDER')
        str=strtrim(str);
        config.(item)=lower(str);
    elseif isequal(item,'SUBJECT_DOB')
        str=strtrim(str);
        config.(item)=str;
        config.EXPT_DATE=date;
        config.SUBJECT_AGE=(datenum(config.EXPT_DATE)-datenum(config.SUBJECT_DOB))/365.2442;
    elseif isequal(item, 'DATA_DIR') || isequal(item, 'BACKGROUND_COLOR') || isequal(item, 'TEXT_COLOR') 
        str=strtrim(str);
        config.(item) = str;
    elseif isequal(item,'SHIFT_DIRECTION') || isequal(item, 'SHIFT_DIRECTION_SUST') 
        str=strtrim(str);
        if isequal(lower(str),'down') || isequal(lower(str),'f1down') || isequal(lower(str),'f1_down')
            str='F1Down';
        elseif isequal(lower(str),'up') || isequal(lower(str),'f1up') || isequal(lower(str),'f1_up')
            str='F1Up';
        end
        config.(item)=str;
    elseif isequal(item,'SHIFT_RATIO') || isequal(item,'MOUTH_MIC_DIST') || ...
            isequal(item,'SPL_TARGET') || isequal(item,'SPL_RANGE') || ...
            isequal(item,'VOWEL_LEN_TARG') || isequal(item,'VOWEL_LEN_RANGE') || ...
            isequal(item,'USE_SPL_TARGET') || isequal(item,'USE_VOWEL_LEN_TARGET') || ...
            isequal(item,'PRE_REPS') || isequal(item,'PRACT1_REPS') || isequal(item,'PRACT2_REPS') || ...
            isequal(item,'NUM_RUNS') || isequal(item, 'TRIALS_PER_RUN') || ...
            isequal(item,'MRI_TR') || isequal(item, 'MRI_TA') || ... 
            isequal(item,'STIM_DELAY') || isequal(item, 'REC_DELAY') || isequal(item, 'REC_LEN') || ... 
            isequal(item,'START_REPS') || isequal(item,'RAMP_REPS') || ...
            isequal(item,'STAY_REPS')  || isequal(item,'STAY2_REPS')  || isequal(item,'END_REPS') || ...
            isequal(item,'TEST1_REPS')  || isequal(item,'TEST2_REPS') || isequal(item,'TEST3_REPS') || ...
            isequal(item,'SAMPLING_RATE') || isequal(item,'FRAME_SIZE') || ...
            isequal(item,'DOWNSAMP_FACT') || isequal(item,'TRIAL_LEN') || isequal(item,'TRIAL_LEN_MAX') || ...
            isequal(item,'MASK_NOISE_LV') || ...
            isequal(item,'NATURAL_REPS') || ...
            isequal(item, 'PITCH_SHIFT_DOWN') || isequal(item, 'PITCH_SHIFT_UP') || ...
            isequal(item, 'RAND_BLOCKS') || isequal(item, 'RAND_TRIALS_PER_BLOCK') || ...
            isequal(item, 'RAND_LOWER_TRIALS_PER_BLOCK') || isequal(item, 'RAND_HIGHER_TRIALS_PER_BLOCK') || ...
            isequal(item, 'SUST_START_REPS') || isequal(item, 'SUST_RAMP_REPS') || ...
            isequal(item, 'SUST_STAY_REPS') || isequal(item, 'SUST_END_REPS')
        str=strtrim(str);
        config.(item)=str2num(str);
%     elseif isequal(item,'PRE_WORDS') || isequal(item,'TEST_WORDS') || isequal(item,'TRAIN_WORDS') || isequal(item,'NATURAL_WORDS') || ...
%            isequal(item, 'RAND_WORDS') || isequal(item, 'SUST_WORDS')
%         str=strtrim(str);
%         idxs=strfind(str,' ');
%         words=cell(1,0);
%         for k1=1:numel(idxs)+1
%             if k1==1
%                 words{end+1}=str(1:idxs(k1)-1);
%             elseif k1==numel(idxs)+1
%                 words{end+1}=str(idxs(k1-1)+1:end);
%             else
%                 words{end+1}=str(idxs(k1-1)+1:idxs(k1)-1);
%             end
%         end
%         config.(item)=words;
    end
end
return