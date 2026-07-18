% main_export_submission_RS_CMSA_ESII_v9_A6.m
% Export the final CEC 2026 submission folder.
%
% The TR requires 960 CSV files, one for each PID/PIN/D:
%   "pid05 pin04 dim20.csv"
% Each file has D solution columns and one fitness column.
%
% IMPORTANT:
% Select ONE global reporting rule based on your validation experiment.
% Do not choose different rules per PID/PIN/run using post-hoc RPR/F1.

clear; clc;
addpath(pwd);

% ---------------- Final submission settings ----------------
pidList = 1:16;
pinList = 1:15;
dimList = [2 5 10 20];

% Default keeps current accepted branch. After A6 validation, you may set:
selectedRuleName = 'A6_medium'; %or another fixed global rule if it wins.
% selectedRuleName = 'A3_relaxed';

seedBase = 500000;   % deterministic final-run seed base
verbose = false;

submissionFolder = fullfile(pwd, ['Submission_S_CARD_CMSA_' selectedRuleName]);
if ~exist(submissionFolder, 'dir')
    mkdir(submissionFolder);
end

logRows = {};
rowNo = 0;
totalJobs = numel(pidList)*numel(pinList)*numel(dimList);
jobNo = 0;

for ip = 1:numel(pidList)
    pid = pidList(ip);
    for ii = 1:numel(pinList)
        pin = pinList(ii);
        for id = 1:numel(dimList)
            D = dimList(id);
            jobNo = jobNo + 1;

            seedNo = seedBase + 10000*pid + 100*pin + D;
            fprintf('\nSubmission job %d/%d: PID=%d, PIN=%d, D=%d, seed=%d\n', ...
                jobNo, totalJobs, pid, pin, D, seedNo);

            status = "OK";
            errMsg = "";

            try
                problem = ProblemMM(pid, pin, D);
                problem.form();

                result = S_CARD_CMSA(problem, seedNo, verbose);
                [X, f, ruleFound] = local_select_report(result, selectedRuleName);

                if ~ruleFound
                    warning('Selected rule %s not found. Falling back to A3_relaxed/default final_pop.', selectedRuleName);
                    X = result.final_pop;
                    f = result.final_f;
                end

                % Ensure correct shape and finite output.
                if isempty(X)
                    X = zeros(0,D);
                    f = zeros(0,1);
                end
                f = f(:);
                out = [X, f];

                fileName = sprintf('pid%02d pin%02d dim%d.csv', pid, pin, D);
                filePath = fullfile(submissionFolder, fileName);
                writematrix(out, filePath);

                rowNo = rowNo + 1;
                logRows(rowNo,:) = {pid, pin, D, seedNo, selectedRuleName, size(X,1), ...
                    result.problemInfo.used_eval, result.problemInfo.max_eval, ...
                    result.rpr_selected, result.score_selected, status, errMsg}; %#ok<SAGROW>

            catch ME
                status = "ERROR";
                errMsg = string(ME.message);
                fprintf(2, '\nERROR PID=%d PIN=%d D=%d: %s\n', pid, pin, D, getReport(ME, 'extended', 'hyperlinks', 'off'));

                rowNo = rowNo + 1;
                logRows(rowNo,:) = {pid, pin, D, seedNo, selectedRuleName, NaN, NaN, NaN, NaN, NaN, status, errMsg}; %#ok<SAGROW>
            end

            logTbl = cell2table(logRows, 'VariableNames', { ...
                'PID','PIN','Dim','Seed','SelectedRule','Nsol','UsedEval','MaxEval', ...
                'InternalRPR_DefaultA3','InternalScore_DefaultA3','Status','ErrorMessage'});
            save(fullfile(submissionFolder, 'submission_log.mat'), 'logTbl');
            writetable(logTbl, fullfile(submissionFolder, 'submission_log.csv'));
        end
    end
end

fprintf('\nSubmission export completed. Folder:\n%s\n', submissionFolder);

function [X, f, found] = local_select_report(result, selectedRuleName)
    found = false;
    X = [];
    f = [];

    % Default accepted branch.
    if strcmp(selectedRuleName, 'A3_relaxed')
        X = result.final_pop_A3;
        f = result.final_f_A3;
        found = true;
        return;
    elseif strcmp(selectedRuleName, 'A1_strict')
        X = result.final_pop_A1;
        f = result.final_f_A1;
        found = true;
        return;
    end

    % A6 extra reporting rules.
    if isfield(result, 'A6') && isfield(result.A6, 'names')
        idx = find(strcmp(result.A6.names, selectedRuleName), 1, 'first');
        if ~isempty(idx)
            fr = result.A6.reports{idx};
            X = fr.solution;
            f = fr.value;
            found = true;
        end
    end
end
