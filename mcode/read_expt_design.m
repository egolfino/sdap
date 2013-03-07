function expt_design = read_expt_design(xlsfn, nruns, tpr)
%%
% xlsfn ............ Input experiment design xls file
% nruns ............ Number of runs
% tpr   ............ Number of trials in each run
%%
if ~isfile(xlsfn)
    error('Cannot experiment design xls file: %s', xlsfn);
end
[N, T] = xlsread(xlsfn);

expt_design = struct;

for i1 = 1 : nruns
    runScript = struct;
    runScript.nTrials = tpr;
    runScript.trialTypes = nan(1, tpr);
    runScript.stimUtters = cell(1, tpr);
    
    col0 = mod(i1 - 1, 4) * 2 + 1;
    row0 = floor((i1 - 1) / 4) * 42 + 2;
%     fprintf('i1 = %d, row0 = %d, col0 = %d\n', i1, row0, col0);

    for i2 = 1 : tpr
        txt1 = lower(deblank(T{row0 + i2 - 1, col0}));
        txt2 = T{row0 + i2 - 1, col0 + 1};
        
        if isequal(txt1, 'normal')
            runScript.trialTypes(i2) = 1; 
        elseif isequal(txt1, 'pitch perturbed down')
            runScript.trialTypes(i2) = 2;
        elseif isequal(txt1, 'pitch perturbed up')
            runScript.trialTypes(i2) = 3;
        elseif isequal(txt1, 'noise masked')
            runScript.trialTypes(i2) = 4;
        elseif isequal(txt1, 'baseline')
            runScript.trialTypes(i2) = 5;
        else
            error('Unrecognized trial type in %s: %s', xlsfn, txt1);
        end
        
%         if runScript.trialTypes(i2) >= 1 && runScript.trialTypes(i2) <= 5
        if runScript.trialTypes(i2) == 5
            for j1 = 1 : numel(txt2)
                if txt2(j1) >= 256
                    txt2(j1) = '#';
                end
            end
        end
        runScript.stimUtters{i2} = txt2;
%         end
    end
    
    expt_design.(['run', num2str(i1)]) = runScript;
end

return