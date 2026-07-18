function result = RS_CMSA_ESII_v9_A6(problem, seedNo, verbose)
% RS_CMSA_ESII_v9_A6: A6 = A3b + official-score-aware final reporting diagnostics
% -------------------------------------------------------------------------
% MATLAB RS-CMSA-ESII A3b version for the IEEE CEC 2026
% Multimodal Optimization benchmark (ProblemMM.m).
%
% A6 = A3b plus multiple fixed final-reporting rules for RPR/F1/official-score evaluation.
% The secondary archive stores the best solution from every restart and is
% used only in the final reporting stage. It does not change sampling,
% taboo regions, restart behavior, covariance adaptation, or archive updates.
%
% Usage:
%   problem = ProblemMM(pid, instance_no, dim);
%   problem.form();
%   result = RS_CMSA_ESII_v9_A6(problem, seedNo, true);
%
% Input:
%   problem : formed ProblemMM object. The objective is assumed minimization.
%   seedNo  : MATLAB RNG seed.
%   verbose : true/false for restart progress display.
%
% Output fields:
%   result.archive      : internal RS-CMSA-ESII archive
%   result.process      : process/state information
%   result.final_pop    : final reported solutions = archive.solution
%   result.final_f      : function values for final_pop = archive.value
%   result.rpr          : robust peak ratio, if final_pop is non-empty
%   result.i_pr         : per-minimum robust peak score, if available
%   result.usedTimeSec  : wall-clock time
%   result.problemInfo  : pid, instance, dim, max_eval, used_eval
%
% Notes:
%   1) This version keeps the A0 optimizer unchanged and computes both
%      A1 cleanup and A3 refined cleanup on exactly the same archive/candidate set.
%   2) This code does NOT use problem.minima.X during optimization. Minima
%      data are used only after optimization to calculate RPR.
%   3) Keep all vectors as row vectors and all populations as N-by-D matrices.
% -------------------------------------------------------------------------

    if nargin < 2 || isempty(seedNo)
        seedNo = 1;
    end
    if nargin < 3 || isempty(verbose)
        verbose = true;
    end

    rng(seedNo, 'twister');
    startTime = tic;

    if isnan(problem.max_eval)
        error('ProblemMM object is not formed. Please call problem.form() before running the optimizer.');
    end
    if problem.used_eval > 0
        warning(['problem.used_eval is already > 0. For a clean run, create a fresh ', ...
                 'ProblemMM object and call problem.form() before optimization.']);
    end

    opt     = rs_options(problem);
    process = rs_process_init(opt, problem);
    archive = rs_archive_init(problem);
    candidateArchive = rs_candidate_archive_init(problem);
    diagnostics = rs_diagnostics_init(problem, seedNo);

    while process.usedEvalTillRestart < problem.max_eval
        restart = rs_restart_init(process, opt, problem);
        subpop  = rs_initialize_subpop(archive, process, opt, problem, restart);

        [restart, subpop, archive] = rs_run_one_restart(subpop, archive, process, opt, problem, restart);

        % A1 addition: passively store the best solution of this restart.
        % This archive is NOT used by the optimizer during the run.
        candidateArchive = rs_candidate_archive_append(candidateArchive, subpop, restart, process, archive, problem);

        if restart.terminationFlag == -1
            diagnostics = rs_diagnostics_update(diagnostics, restart, subpop, archive, process, NaN);
            break;
        end

        [archive, problem] = rs_archive_update(subpop, restart, process, opt, problem, archive);
        repeatedHitsThisRestart = sum(archive.hitTimesThisRestart);
        diagnostics = rs_diagnostics_update(diagnostics, restart, subpop, archive, process, repeatedHitsThisRestart);
        process = rs_process_update(restart, archive, process, opt, problem);

        if verbose
            fprintf('restartNo = %d, usedEval = %d (%.2f%%), archiveSize = %d, termFlag = %d, condC = %.3e\n', ...
                process.restartNo-1, process.usedEvalTillRestart, ...
                100*process.usedEvalTillRestart/problem.max_eval, archive.size, restart.subpopTermFlag, diagnostics.restartLog.condC(end));
        end

        if problem.used_eval >= problem.max_eval
            break;
        end
    end

    % A3b final reporting: compute BOTH A1 cleanup and A3 refined cleanup
    % on the exact same search trajectory, primary archive, and secondary
    % candidate archive. These steps do not consume extra FEs and do not use
    % the known true minima.
    [final_pop_A1, final_f_A1, finalReportA1] = rs_final_report_secondary_archive_A1(archive, candidateArchive, opt, problem);
    [final_pop_A3, final_f_A3, finalReportA3] = rs_final_report_refined_cleanup_A3(archive, candidateArchive, opt, problem);

    result.archive          = archive;              % original A0 primary archive
    result.candidateArchive = candidateArchive;     % passive secondary candidate archive
    result.finalReportA1    = finalReportA1;        % v2/A1 cleanup report
    result.finalReportA3    = finalReportA3;        % v4/A3 refined cleanup report
    result.finalReport      = finalReportA3;        % default report = A3 refined cleanup
    result.process          = process;
    result.primary_final_pop = archive.solution;
    result.primary_final_f   = archive.value;
    result.final_pop_A1     = final_pop_A1;
    result.final_f_A1       = final_f_A1;
    result.final_pop_A3     = final_pop_A3;
    result.final_f_A3       = final_f_A3;
    result.final_pop        = final_pop_A3;         % default final output = A3
    result.final_f          = final_f_A3;
    result.usedTimeSec      = toc(startTime);
    result.problemInfo = struct('pid', problem.pid, 'instance_no', problem.instance_no, ...
                                'dim', problem.dim, 'max_eval', problem.max_eval, ...
                                'used_eval', problem.used_eval, 'seedNo', seedNo);
    result.diagnostics = rs_diagnostics_finalize(diagnostics, archive, candidateArchive, ...
        finalReportA3, result.problemInfo, result.usedTimeSec);

    % Extra same-run cleanup diagnostics.
    result.diagnostics.summary.nReportedOptima_A1 = size(finalReportA1.solution,1);
    result.diagnostics.summary.nReportedOptima_A3 = size(finalReportA3.solution,1);
    result.diagnostics.summary.nCandidateUsed_A1 = finalReportA1.nCandidateUsed;
    result.diagnostics.summary.nCandidateUsed_A3 = finalReportA3.nCandidateUsed;
    result.diagnostics.summary.nCandidateAdded_A1 = finalReportA1.nCandidateAddedToReport;
    result.diagnostics.summary.nCandidateAdded_A3 = finalReportA3.nCandidateAddedToReport;
    result.diagnostics.summary.nCandidateReplaced_A1 = finalReportA1.nCandidateReplacedExisting;
    result.diagnostics.summary.nCandidateReplaced_A3 = finalReportA3.nCandidateReplacedExisting;

    % Post-processing only: robust peak ratio using official UtilityMethod.
    % result.rpr_A1_cleanup and result.rpr_A3_cleanup are directly comparable
    % because they come from the same run and the same candidate pool.
    ftol0 = [1e-5, 1];
    result.rpr_A1_cleanup = NaN;
    result.i_pr_A1_cleanup = [];
    result.rpr_A3_cleanup = NaN;
    result.i_pr_A3_cleanup = [];
    result.rpr  = NaN;     % default = A3 cleanup
    result.i_pr = [];
    result.primary_rpr  = NaN;
    result.primary_i_pr = [];

    if ~isempty(result.final_pop_A1)
        try
            [result.rpr_A1_cleanup, result.i_pr_A1_cleanup] = UtilityMethod.calc_robust_peak_ratio( ...
                result.final_pop_A1, result.final_f_A1, problem.minima.X, problem.minima.f, ftol0);
        catch ME
            warning('A1 cleanup RPR calculation failed: %s', ME.message);
        end
    end
    if ~isempty(result.final_pop_A3)
        try
            [result.rpr_A3_cleanup, result.i_pr_A3_cleanup] = UtilityMethod.calc_robust_peak_ratio( ...
                result.final_pop_A3, result.final_f_A3, problem.minima.X, problem.minima.f, ftol0);
        catch ME
            warning('A3 cleanup RPR calculation failed: %s', ME.message);
        end
    end
    if ~isempty(result.primary_final_pop)
        try
            [result.primary_rpr, result.primary_i_pr] = UtilityMethod.calc_robust_peak_ratio( ...
                result.primary_final_pop, result.primary_final_f, problem.minima.X, problem.minima.f, ftol0);
        catch ME
            warning('A0 primary RPR calculation failed: %s', ME.message);
        end
    end

    result.rpr = result.rpr_A3_cleanup;
    result.i_pr = result.i_pr_A3_cleanup;
    result.deltaRPR_A3_vs_A1 = result.rpr_A3_cleanup - result.rpr_A1_cleanup;
    result.deltaRPR_A1_vs_Primary = result.rpr_A1_cleanup - result.primary_rpr;
    result.deltaRPR_A3_vs_Primary = result.rpr_A3_cleanup - result.primary_rpr;
    result.diagnostics.summary.rpr_A1_cleanup = result.rpr_A1_cleanup;
    result.diagnostics.summary.rpr_A3_cleanup = result.rpr_A3_cleanup;
    result.diagnostics.summary.rpr_A0_primary = result.primary_rpr;
    result.diagnostics.summary.deltaRPR_A3_vs_A1 = result.deltaRPR_A3_vs_A1;
    result.diagnostics.summary.deltaRPR_A1_vs_Primary = result.deltaRPR_A1_vs_Primary;
    result.diagnostics.summary.deltaRPR_A3_vs_Primary = result.deltaRPR_A3_vs_Primary;

    % A6 addition: evaluate several fixed final-reporting rules on exactly
    % the same already-evaluated candidate pool. This section is for
    % development/validation of the official CEC 2026 score:
    %   score = 0.5*(RPR + F1),
    % where precision = RPR * NGM/Nsol and recall = RPR.
    %
    % These calculations are POST-PROCESSING ONLY. They do not affect the
    % search trajectory and do not consume function evaluations.
    result.A6 = rs_evaluate_A6_reporting_rules(archive, candidateArchive, opt, problem, ftol0);

    % Safe default: keep the accepted A3 report as the optimizer's default
    % final output. After the A6 validation experiment, a fixed rule can be
    % selected globally for the final submission runner.
    result.finalReportSelected = result.finalReportA3;
    result.final_pop_selected  = result.final_pop_A3;
    result.final_f_selected    = result.final_f_A3;
    result.rpr_selected        = result.rpr_A3_cleanup;
    result.f1_selected         = rs_f1_from_rpr_nsol(result.rpr_A3_cleanup, size(result.final_pop_A3,1), problem.n_minima);
    result.score_selected      = 0.5 * (result.rpr_selected + result.f1_selected);

    result.diagnostics.summary.f1_A1_cleanup = rs_f1_from_rpr_nsol(result.rpr_A1_cleanup, size(result.final_pop_A1,1), problem.n_minima);
    result.diagnostics.summary.f1_A3_cleanup = result.f1_selected;
    result.diagnostics.summary.score_A1_cleanup = 0.5 * (result.rpr_A1_cleanup + result.diagnostics.summary.f1_A1_cleanup);
    result.diagnostics.summary.score_A3_cleanup = result.score_selected;
end


% =========================================================================
% DIAGNOSTICS
% =========================================================================
function diagnostics = rs_diagnostics_init(problem, seedNo)
    diagnostics.meta.pid = problem.pid;
    diagnostics.meta.instance_no = problem.instance_no;
    diagnostics.meta.dim = problem.dim;
    diagnostics.meta.seedNo = seedNo;

    diagnostics.restartLog.restartNo = zeros(0,1);
    diagnostics.restartLog.termFlag = zeros(0,1);
    diagnostics.restartLog.subpopTermFlag = zeros(0,1);
    diagnostics.restartLog.usedEvalEvolve = zeros(0,1);
    diagnostics.restartLog.usedEvalMerge = zeros(0,1);
    diagnostics.restartLog.usedEvalArchive = zeros(0,1);
    diagnostics.restartLog.usedEvalTotal = zeros(0,1);
    diagnostics.restartLog.archiveSize = zeros(0,1);
    diagnostics.restartLog.repeatedHits = zeros(0,1);
    diagnostics.restartLog.bestVal = zeros(0,1);
    diagnostics.restartLog.condC = zeros(0,1);
end

function diagnostics = rs_diagnostics_update(diagnostics, restart, subpop, archive, process, repeatedHitsThisRestart)
    if nargin < 6 || isempty(repeatedHitsThisRestart) || isnan(repeatedHitsThisRestart)
        repeatedHitsThisRestart = 0;
    end

    if isfield(subpop, 'mutProfile') && isfield(subpop.mutProfile, 'stretch') && ~isempty(subpop.mutProfile.stretch)
        condC = (max(subpop.mutProfile.stretch) / max(min(subpop.mutProfile.stretch), realmin))^2;
    else
        condC = NaN;
    end

    usedEvalArchive = archive.usedEval;
    usedEvalTotal = restart.usedEvalEvolve + restart.usedEvalMerge + usedEvalArchive;

    diagnostics.restartLog.restartNo(end+1,1) = process.restartNo;
    diagnostics.restartLog.termFlag(end+1,1) = restart.terminationFlag;
    diagnostics.restartLog.subpopTermFlag(end+1,1) = restart.subpopTermFlag;
    diagnostics.restartLog.usedEvalEvolve(end+1,1) = restart.usedEvalEvolve;
    diagnostics.restartLog.usedEvalMerge(end+1,1) = restart.usedEvalMerge;
    diagnostics.restartLog.usedEvalArchive(end+1,1) = usedEvalArchive;
    diagnostics.restartLog.usedEvalTotal(end+1,1) = usedEvalTotal;
    diagnostics.restartLog.archiveSize(end+1,1) = archive.size;
    diagnostics.restartLog.repeatedHits(end+1,1) = repeatedHitsThisRestart;
    diagnostics.restartLog.bestVal(end+1,1) = subpop.bestVal;
    diagnostics.restartLog.condC(end+1,1) = condC;
end

function diagnostics = rs_diagnostics_finalize(diagnostics, archive, candidateArchive, finalReport, problemInfo, usedTimeSec)
    log = diagnostics.restartLog;
    nRestarts = numel(log.restartNo);

    condVals = log.condC(isfinite(log.condC));
    if isempty(condVals)
        condMin = NaN; condMean = NaN; condMedian = NaN; condMax = NaN; condStd = NaN;
    else
        condMin = min(condVals);
        condMean = mean(condVals);
        condMedian = median(condVals);
        condMax = max(condVals);
        condStd = std(condVals);
    end

    if isempty(log.usedEvalTotal)
        avgEvalPerRestart = NaN;
    else
        avgEvalPerRestart = mean(log.usedEvalTotal);
    end

    diagnostics.summary = struct();
    diagnostics.summary.pid = problemInfo.pid;
    diagnostics.summary.instance_no = problemInfo.instance_no;
    diagnostics.summary.dim = problemInfo.dim;
    diagnostics.summary.seedNo = problemInfo.seedNo;
    diagnostics.summary.nReportedOptima = size(finalReport.solution,1);
    diagnostics.summary.nArchiveMembers = archive.size;
    diagnostics.summary.nPrimaryReportedOptima = size(archive.solution,1);
    diagnostics.summary.nCandidateArchiveRaw = candidateArchive.size;
    diagnostics.summary.nCandidateArchiveFinite = finalReport.nCandidateFinite;
    diagnostics.summary.nCandidateArchiveFiltered = finalReport.nCandidateFiltered;
    diagnostics.summary.nCandidateArchiveUsed = finalReport.nCandidateUsed;
    diagnostics.summary.nCandidateAddedToReport = finalReport.nCandidateAddedToReport;
    diagnostics.summary.nCandidateReplacedExisting = finalReport.nCandidateReplacedExisting;
    diagnostics.summary.nFinalReportedOptima = size(finalReport.solution,1);
    diagnostics.summary.nFinalCleanupExactDuplicates = getfield_safe(finalReport, 'nExactDuplicateMerged', NaN);
    diagnostics.summary.nFinalCleanupNearDuplicates = getfield_safe(finalReport, 'nNearDuplicateMerged', NaN);
    diagnostics.summary.nFinalCleanupTrimmed = getfield_safe(finalReport, 'nTrimmedByMaxReport', NaN);
    diagnostics.summary.finalCleanupStrategyCode = getfield_safe(finalReport, 'strategyCode', NaN);
    diagnostics.summary.nRestarts = nRestarts;
    diagnostics.summary.nMergeTerminations = sum(log.subpopTermFlag == 3);
    diagnostics.summary.nLocalFailureTerminations = sum(log.subpopTermFlag == 4);
    diagnostics.summary.nRepeatedArchivedHits = sum(log.repeatedHits);
    diagnostics.summary.finalRepeatedArchiveHits = sum(archive.hitTimesSoFar);
    diagnostics.summary.avgEvalPerRestart = avgEvalPerRestart;
    diagnostics.summary.condC_min = condMin;
    diagnostics.summary.condC_mean = condMean;
    diagnostics.summary.condC_median = condMedian;
    diagnostics.summary.condC_max = condMax;
    diagnostics.summary.condC_std = condStd;
    diagnostics.summary.usedEval = problemInfo.used_eval;
    diagnostics.summary.maxEval = problemInfo.max_eval;
    diagnostics.summary.usedTimeSec = usedTimeSec;
end

function val = getfield_safe(s, fieldName, defaultVal)
    if isstruct(s) && isfield(s, fieldName)
        val = s.(fieldName);
    else
        val = defaultVal;
    end
end

% =========================================================================
% OPTIONS
% =========================================================================
function opt = rs_options(problem)
    D = problem.dim;

    opt.archiving.hillVallBudget   = 10;
    opt.archiving.iniNormTabDis    = 1;
    opt.archiving.newNormTabDisPrc = 10;
    opt.archiving.targetNewNicheFr = 0.5;
    opt.archiving.targetGlobFr     = 0.5;
    opt.archiving.tauNormTabDis    = sqrt(1/D);
    opt.archiving.tolFunArch       = 1e-5;
    opt.archiving.neighborSize     = 5;

    opt.coreSearch.algorithm          = 'CMSA';
    opt.coreSearch.targetNumSubpop    = 1;
    opt.coreSearch.iniSubpopSizeCoeff = 6;
    opt.coreSearch.finSubpopSizeCoeff = 6;
    opt.coreSearch.iniSigCoeff        = 2;
    opt.coreSearch.muToPopSizeRatio   = 0.2;
    opt.coreSearch.maxIniSigma        = 0.3;
    opt.coreSearch.objValUpLimit      = 1e100;
    opt.coreSearch.tauSigmaCoeff      = 0.5;
    opt.coreSearch.eltRatio           = 0.1;
    opt.coreSearch.tauCovCoeff        = 1;
    opt.coreSearch.sigmaUpdateBiasImp = 1;

    opt.niching.criticTabooThresh = 0.01;
    opt.niching.redCoeff          = 0.99^(1/D);
    opt.niching.iniR0IncFac       = 1.04;
    opt.niching.maxRejectIni      = 100;

    opt.stopCr.tolHistFun = 1e-6;
    opt.stopCr.maxCondC   = 1e14;
    opt.stopCr.maxIterPar = [100, 50];
    opt.stopCr.tolX       = 1e-12;
    opt.stopCr.stagPar    = [120, 0.2, 30];

    opt.stopCr.merge.threshold       = 0.5;
    opt.stopCr.merge.windowSizeCoeff = 0.1;
    opt.stopCr.merge.chkIntervalCoeff= 0.1;
    opt.stopCr.merge.maxEval         = 10;
    opt.stopCr.merge.addedConst      = 1;

    opt.stopCr.localConverge.tolCoeff            = 0.04;
    opt.stopCr.localConverge.windowSizeCoeff     = 0.5;
    opt.stopCr.localConverge.tabooCriticUpLimit  = 0.01;

    % A1 secondary candidate archive options retained.
    % A3 changes only the final cleanup/reporting policy. It does not alter
    % sampling, taboo regions, restarts, covariance adaptation, or FEs.
    opt.secondaryArchive.enable          = true;
    opt.secondaryArchive.tolFunWindow    = 1.0;
    opt.secondaryArchive.dupNormDistCoeff = 0.005; % retained for compatibility; A3 uses finalCleanup thresholds
    opt.secondaryArchive.maxCandidateMultiplier = 10; % retained for compatibility
    opt.secondaryArchive.keepOnlyFeasibleFinite = true;

    % A3 refined final cleanup/reporting options.
    % Rationale: the CEC 2026 robust peak ratio does not penalize extra
    % reported points; it takes the best reported value assigned to each
    % true minimum during post-processing. Therefore, A3 uses a less
    % aggressive cleanup than v2/A1, preserving more spatially distinct
    % evaluated candidates while removing only exact/very-near duplicates.
    opt.finalCleanup.enable = true;
    opt.finalCleanup.strategyCode = 3;
    opt.finalCleanup.valueWindow = 1.0;          % aligned with upper ftol in RPR calculation
    opt.finalCleanup.exactDupNormDistCoeff = 1e-10;
    opt.finalCleanup.nearDupNormDistCoeff  = 5e-4;
    opt.finalCleanup.maxReportMultiplier = 100; % generous cap: 100*n_minima
    opt.finalCleanup.minMaxReport = 500;
    opt.finalCleanup.keepPrimaryAlways = true;
    opt.finalCleanup.keepOnlyFeasibleFinite = true;
end

% =========================================================================
% PROCESS
% =========================================================================
function process = rs_process_init(opt, problem)
    D = problem.dim;
    process.restartNo  = 0;
    process.subpopSize = max(1, floor(opt.coreSearch.finSubpopSizeCoeff * sqrt(D)));
    process.mu         = max(1, floor(0.5 + process.subpopSize * opt.coreSearch.muToPopSizeRatio));
    W = log(1 + process.mu) - log(1:process.mu);
    process.recWeights = W / sum(W);
    process.muEff      = sum(process.recWeights)^2 / sum(process.recWeights.^2);
    process.coreSpecStrPar = rs_spec_str_par_CMSA(process, opt, problem);
    process.defNormTabDis  = opt.archiving.iniNormTabDis;
    process.usedEvalTillRestart = 0;
    process.bestValTillRestart  = inf;
    process.iniR0 = sqrt(D)/2;
end

function core = rs_spec_str_par_CMSA(process, opt, problem)
    D = problem.dim;
    core.numElt   = ceil(opt.coreSearch.eltRatio * process.subpopSize);
    core.tauCov   = 1 + D*(1+D)/(2*process.muEff*opt.coreSearch.tauCovCoeff);
    core.tauSigma = sqrt(0.5*opt.coreSearch.tauSigmaCoeff/D);
end

function process = rs_process_update(restart, archive, process, opt, problem)
    process.restartNo = process.restartNo + 1;

    usedEvalThisRestart = archive.usedEval + restart.usedEvalEvolve + restart.usedEvalMerge;
    usedEvalSoFar       = process.usedEvalTillRestart + usedEvalThisRestart;

    coeff = opt.coreSearch.iniSubpopSizeCoeff * ...
        (opt.coreSearch.finSubpopSizeCoeff/opt.coreSearch.iniSubpopSizeCoeff)^(usedEvalSoFar/problem.max_eval);
    process.subpopSize = max(1, floor(coeff * sqrt(problem.dim)));
    process.mu         = max(1, floor(0.5 + process.subpopSize * opt.coreSearch.muToPopSizeRatio));

    W = log(1 + process.mu) - log(1:process.mu);
    process.recWeights = W / sum(W);
    process.muEff      = sum(process.recWeights)^2 / sum(process.recWeights.^2);
    process.coreSpecStrPar = rs_spec_str_par_CMSA(process, opt, problem);

    if archive.size > 0
        process.defNormTabDis = prctile(archive.normTabDis, opt.archiving.newNormTabDisPrc);
    else
        process.defNormTabDis = opt.archiving.iniNormTabDis;
    end

    process.usedEvalTillRestart = process.usedEvalTillRestart + usedEvalThisRestart;
    process.bestValTillRestart  = min([restart.subpopBestVal, process.bestValTillRestart]);
    process.iniR0 = min(restart.recIniR0 * opt.niching.iniR0IncFac, 0.5*sqrt(problem.dim));
end

% =========================================================================
% ARCHIVE
% =========================================================================
function archive = rs_archive_init(problem)
    D = problem.dim;
    archive.solution = zeros(0,D);
    archive.value    = zeros(0,1);
    archive.normTabDis = zeros(0,1);
    archive.foundEval  = zeros(0,1);
    archive.usedEval   = 0;
    archive.hitTimesSoFar = zeros(0,1);
    archive.hitTimesThisRestart = zeros(0,1);
    archive.size = 0;
    archive.foundTime = zeros(0,1);
    archive.numCallF  = zeros(0,1);

    archive.dummyArchive.solution   = zeros(0,D);
    archive.dummyArchive.value      = zeros(0,1);
    archive.dummyArchive.foundEval  = zeros(0,1);
    archive.dummyArchive.numCallF   = zeros(0,1);
    archive.dummyArchive.foundTime  = zeros(0,1);
    archive.dummyArchive.actionCode = zeros(0,1);
end

function [archive, problem] = rs_archive_update(subpop, restart, process, opt, problem, archive)
    archive.usedEval = 0;
    archive.hitTimesThisRestart = zeros(archive.size,1);
    bestValSoFar = min([restart.subpopBestVal, process.bestValTillRestart]);

    % Remove archive entries that are no longer desirable.
    if archive.size > 0
        keepIt = (archive.value - opt.archiving.tolFunArch) < bestValSoFar;
        discardInd = find(~keepIt);
        archive = rs_dummy_archive_append(archive, -1, discardInd, problem);

        archive.solution = archive.solution(keepIt,:);
        archive.value    = archive.value(keepIt);
        archive.normTabDis = archive.normTabDis(keepIt);
        archive.foundEval  = archive.foundEval(keepIt);
        archive.hitTimesThisRestart = archive.hitTimesThisRestart(keepIt);
        archive.hitTimesSoFar = archive.hitTimesSoFar(keepIt);
        archive.numCallF  = archive.numCallF(keepIt);
        archive.foundTime = archive.foundTime(keepIt);
        archive.size = numel(archive.value);
    end

    Ndesirable = 0;
    isNew = false;
    chkEvolved  = restart.iterNo > 1;
    chkIsGlobal = (subpop.bestVal - opt.archiving.tolFunArch) <= bestValSoFar;

    if chkEvolved && chkIsGlobal
        Ndesirable = Ndesirable + 1;
        if archive.size == 0
            isNew = true;
            archive = rs_archive_append(archive, subpop.bestSol, subpop.bestVal, restart, process, opt, problem);
            archive = rs_dummy_archive_append(archive, 1, archive.size, problem);
        else
            [isNew, matchArchNo, usedEval] = rs_archive_is_new_basin(archive, subpop.bestSol, subpop.bestVal, opt, problem);
            archive.usedEval = archive.usedEval + usedEval;

            if ~isNew
                archive.hitTimesSoFar(matchArchNo) = archive.hitTimesSoFar(matchArchNo) + 1;
                archive.hitTimesThisRestart(matchArchNo) = archive.hitTimesThisRestart(matchArchNo) + 1;

                if subpop.bestVal < (archive.value(matchArchNo) - opt.stopCr.tolHistFun)
                    archive = rs_dummy_archive_append(archive, -1, matchArchNo, problem);
                    archive.value(matchArchNo) = subpop.bestVal;
                    archive.solution(matchArchNo,:) = subpop.bestSol;
                    % Static MMO version: usedEvalChangeDetect = 0.
                    archive.foundEval(matchArchNo) = archive.usedEval + restart.usedEvalMerge + ...
                        restart.usedEvalEvolve + process.usedEvalTillRestart;
                    archive.foundTime(matchArchNo) = 0;
                    archive.numCallF(matchArchNo)  = problem.used_eval;
                    archive = rs_dummy_archive_append(archive, 1, matchArchNo, problem);
                end
            else
                archive = rs_archive_append(archive, subpop.bestSol, subpop.bestVal, restart, process, opt, problem);
                archive = rs_dummy_archive_append(archive, 1, archive.size, problem);
            end
        end
    end

    % Adapt normalized taboo distances.
    if archive.size > 0
        if Ndesirable == 0
            archive.normTabDis = archive.normTabDis .* ...
                exp(-opt.archiving.tauNormTabDis * opt.archiving.targetGlobFr / archive.size);
        elseif isNew
            % Do nothing.
        else
            repDiff = archive.hitTimesThisRestart;
            if any(archive.hitTimesThisRestart == 0)
                if archive.size > 1
                    repDiff(archive.hitTimesThisRestart == 0) = ...
                        -(1 - opt.archiving.targetNewNicheFr)/(archive.size - 1);
                else
                    repDiff(archive.hitTimesThisRestart == 0) = 0;
                end
            end
            archive.normTabDis = archive.normTabDis .* exp(opt.archiving.tauNormTabDis * repDiff);
        end
    end
end

function archive = rs_archive_append(archive, sol, val, restart, process, opt, problem)
    archive.solution = [archive.solution; sol(:)']; %#ok<AGROW>
    archive.value    = [archive.value; val]; %#ok<AGROW>
    archive.normTabDis = [archive.normTabDis; process.defNormTabDis]; %#ok<AGROW>
    archive.foundEval  = [archive.foundEval; archive.usedEval + restart.usedEvalMerge + ...
        restart.usedEvalEvolve + process.usedEvalTillRestart]; %#ok<AGROW>
    archive.size = archive.size + 1;
    archive.hitTimesThisRestart = [archive.hitTimesThisRestart; 0]; %#ok<AGROW>
    archive.hitTimesSoFar = [archive.hitTimesSoFar; 0]; %#ok<AGROW>
    archive.foundTime = [archive.foundTime; 0]; %#ok<AGROW>
    archive.numCallF  = [archive.numCallF; problem.used_eval]; %#ok<AGROW>
end

function archive = rs_dummy_archive_append(archive, action, index, problem)
    if isempty(index)
        return;
    end
    index = index(:);
    index(index < 1 | index > archive.size) = [];
    if isempty(index)
        return;
    end
    archive.dummyArchive.solution = [archive.dummyArchive.solution; archive.solution(index,:)]; %#ok<AGROW>
    archive.dummyArchive.value    = [archive.dummyArchive.value; archive.value(index)]; %#ok<AGROW>
    archive.dummyArchive.foundEval= [archive.dummyArchive.foundEval; archive.foundEval(index)]; %#ok<AGROW>
    archive.dummyArchive.numCallF = [archive.dummyArchive.numCallF; archive.numCallF(index)]; %#ok<AGROW>
    archive.dummyArchive.foundTime= [archive.dummyArchive.foundTime; archive.foundTime(index)]; %#ok<AGROW>
    archive.dummyArchive.actionCode = [archive.dummyArchive.actionCode; action*ones(numel(index),1)]; %#ok<AGROW>
end

function [isNew, sameArchInd, totalUsedEval] = rs_archive_is_new_basin(archive, x, f, opt, problem)
    totalUsedEval = 0;
    sameArchInd = NaN;
    isNew = true;

    if archive.size == 0
        return;
    end

    dis = sqrt(sum((archive.solution - x(:)').^2, 2));
    [~, candidInd] = sort(dis, 'ascend');
    candidInd = candidInd(1:min(opt.archiving.neighborSize, numel(candidInd)));

    for ii = 1:numel(candidInd)
        archNo = candidInd(ii);
        usedEval = 0;
        basinIsNewWrtThisArchive = false;
        while usedEval < opt.archiving.hillVallBudget && problem.used_eval < problem.max_eval
            r = 0.8*rand() + 0.1;
            testX = archive.solution(archNo,:) + r*(x(:)' - archive.solution(archNo,:));
            testF = problem.func_eval(testX);
            totalUsedEval = totalUsedEval + 1;
            usedEval = usedEval + 1;
            if testF > (max(f, archive.value(archNo)) + opt.stopCr.tolHistFun)
                basinIsNewWrtThisArchive = true;
                break;
            end
        end
        if ~basinIsNewWrtThisArchive
            isNew = false;
            sameArchInd = archNo;
            break;
        end
    end
end


% =========================================================================
% SECONDARY CANDIDATE ARCHIVE AND FINAL REPORTING (A1 ADDITION)
% =========================================================================
function candidateArchive = rs_candidate_archive_init(problem)
    D = problem.dim;
    candidateArchive.solution = zeros(0,D);
    candidateArchive.value = zeros(0,1);
    candidateArchive.restartNo = zeros(0,1);
    candidateArchive.termFlag = zeros(0,1);
    candidateArchive.subpopTermFlag = zeros(0,1);
    candidateArchive.foundEval = zeros(0,1);
    candidateArchive.usedEvalInternal = zeros(0,1);
    candidateArchive.condC = zeros(0,1);
    candidateArchive.archiveSizeAtCapture = zeros(0,1);
    candidateArchive.sourceCode = zeros(0,1);
    candidateArchive.size = 0;
end

function candidateArchive = rs_candidate_archive_append(candidateArchive, subpop, restart, process, archive, problem)
    % Passive capture: store the best solution produced by the restart.
    % It is intentionally not used for future sampling/taboo/restart decisions.
    if ~isfield(subpop, 'bestSol') || isempty(subpop.bestSol) || ~isfinite(subpop.bestVal)
        return;
    end
    if any(~isfinite(subpop.bestSol))
        return;
    end

    condC = NaN;
    if isfield(subpop, 'mutProfile') && isfield(subpop.mutProfile, 'stretch') && ~isempty(subpop.mutProfile.stretch)
        condC = (max(subpop.mutProfile.stretch) / max(min(subpop.mutProfile.stretch), realmin))^2;
    end

    % Source code follows the subpopulation termination flag:
    % 1/2: convergence, 3: merge, 4: local-failure prediction,
    % -1: stagnation, -2: covariance condition failure, 0: budget interruption.
    sourceCode = subpop.terminationFlag;

    candidateArchive.solution = [candidateArchive.solution; subpop.bestSol(:)']; %#ok<AGROW>
    candidateArchive.value = [candidateArchive.value; subpop.bestVal]; %#ok<AGROW>
    candidateArchive.restartNo = [candidateArchive.restartNo; process.restartNo]; %#ok<AGROW>
    candidateArchive.termFlag = [candidateArchive.termFlag; restart.terminationFlag]; %#ok<AGROW>
    candidateArchive.subpopTermFlag = [candidateArchive.subpopTermFlag; subpop.terminationFlag]; %#ok<AGROW>
    candidateArchive.foundEval = [candidateArchive.foundEval; problem.used_eval]; %#ok<AGROW>
    candidateArchive.usedEvalInternal = [candidateArchive.usedEvalInternal; process.usedEvalTillRestart + restart.usedEvalEvolve + restart.usedEvalMerge]; %#ok<AGROW>
    candidateArchive.condC = [candidateArchive.condC; condC]; %#ok<AGROW>
    candidateArchive.archiveSizeAtCapture = [candidateArchive.archiveSizeAtCapture; archive.size]; %#ok<AGROW>
    candidateArchive.sourceCode = [candidateArchive.sourceCode; sourceCode]; %#ok<AGROW>
    candidateArchive.size = candidateArchive.size + 1;
end

function [final_pop, final_f, finalReport] = rs_final_report_secondary_archive_A1(archive, candidateArchive, opt, problem)
    % Build final reported set as a cleaned union of the A0 archive and the
    % passive candidate archive. This function does NOT evaluate new points.
    D = problem.dim;
    primaryX = archive.solution;
    primaryF = archive.value(:);
    nPrimary = size(primaryX,1);

    candX = candidateArchive.solution;
    candF = candidateArchive.value(:);
    finiteCand = false(size(candF));
    if ~isempty(candF)
        finiteCand = isfinite(candF) & all(isfinite(candX),2);
        if opt.secondaryArchive.keepOnlyFeasibleFinite
            finiteCand = finiteCand & (candF < opt.coreSearch.objValUpLimit);
        end
    end
    nCandidateFinite = sum(finiteCand);

    % Value filter: only candidates that can potentially contribute to the
    % robust peak ratio with ftol upper bound = 1.0 are retained.
    allKnownF = [primaryF; candF(finiteCand)];
    if isempty(allKnownF)
        final_pop = zeros(0,D);
        final_f = zeros(0,1);
        finalReport = rs_empty_final_report(D);
        return;
    end
    bestKnownF = min(allKnownF);
    candKeep = finiteCand & (candF <= bestKnownF + opt.secondaryArchive.tolFunWindow);
    candInd = find(candKeep);
    nCandidateFiltered = numel(candInd);

    % Keep only the best bounded number of secondary candidates to avoid an
    % excessively large final report.
    maxCand = opt.secondaryArchive.maxCandidateMultiplier * max(1, problem.n_minima);
    if ~isempty(candInd)
        [~, ord] = sort(candF(candInd), 'ascend');
        candInd = candInd(ord);
        candInd = candInd(1:min(maxCand, numel(candInd)));
    end
    nCandidateUsed = numel(candInd);

    % Order: primary archive first, then secondary candidates by value.
    combX = [primaryX; candX(candInd,:)];
    combF = [primaryF; candF(candInd)];
    combSource = [ones(nPrimary,1); 2*ones(nCandidateUsed,1)];
    combOrigIndex = [(1:nPrimary)'; candInd(:)];

    if isempty(combF)
        final_pop = zeros(0,D);
        final_f = zeros(0,1);
        finalReport = rs_empty_final_report(D);
        finalReport.nCandidateFinite = nCandidateFinite;
        finalReport.nCandidateFiltered = nCandidateFiltered;
        finalReport.nCandidateUsed = nCandidateUsed;
        return;
    end

    % Primary entries are sorted by value, candidates by value after them.
    primaryOrder = find(combSource == 1);
    candOrder = find(combSource == 2);
    [~, po] = sort(combF(primaryOrder), 'ascend');
    [~, co] = sort(combF(candOrder), 'ascend');
    order = [primaryOrder(po); candOrder(co)];

    dupThresh = opt.secondaryArchive.dupNormDistCoeff * sqrt(D);
    keptX = zeros(0,D);
    keptF = zeros(0,1);
    keptSource = zeros(0,1);
    keptOrigIndex = zeros(0,1);
    nCandidateAdded = 0;
    nCandidateReplaced = 0;

    for kk = 1:numel(order)
        ii = order(kk);
        x = combX(ii,:);
        f = combF(ii);
        src = combSource(ii);
        origIdx = combOrigIndex(ii);

        if isempty(keptF)
            keptX = x;
            keptF = f;
            keptSource = src;
            keptOrigIndex = origIdx;
            if src == 2
                nCandidateAdded = nCandidateAdded + 1;
            end
            continue;
        end

        xNorm = (x - problem.low_bound) ./ (problem.up_bound - problem.low_bound);
        keptNorm = (keptX - problem.low_bound) ./ (problem.up_bound - problem.low_bound);
        dis = sqrt(sum((keptNorm - xNorm).^2, 2));
        [minDis, minInd] = min(dis);

        if minDis <= dupThresh
            % Same reported basin/point according to distance cleanup.
            % Keep the better approximation. This can replace a primary
            % archive member only for final reporting; the primary archive
            % stored in result.archive remains unchanged.
            if f < keptF(minInd)
                if src == 2
                    nCandidateReplaced = nCandidateReplaced + 1;
                end
                keptX(minInd,:) = x;
                keptF(minInd) = f;
                keptSource(minInd) = src;
                keptOrigIndex(minInd) = origIdx;
            end
        else
            keptX = [keptX; x]; %#ok<AGROW>
            keptF = [keptF; f]; %#ok<AGROW>
            keptSource = [keptSource; src]; %#ok<AGROW>
            keptOrigIndex = [keptOrigIndex; origIdx]; %#ok<AGROW>
            if src == 2
                nCandidateAdded = nCandidateAdded + 1;
            end
        end
    end

    final_pop = keptX;
    final_f = keptF;

    finalReport.solution = final_pop;
    finalReport.value = final_f;
    finalReport.source = keptSource; % 1 = primary archive, 2 = secondary candidate archive
    finalReport.originalIndex = keptOrigIndex;
    finalReport.bestKnownF = bestKnownF;
    finalReport.valueWindow = opt.secondaryArchive.tolFunWindow;
    finalReport.dupNormDist = dupThresh;
    finalReport.nPrimary = nPrimary;
    finalReport.nCandidateRaw = candidateArchive.size;
    finalReport.nCandidateFinite = nCandidateFinite;
    finalReport.nCandidateFiltered = nCandidateFiltered;
    finalReport.nCandidateUsed = nCandidateUsed;
    finalReport.nCandidateAddedToReport = nCandidateAdded;
    finalReport.nCandidateReplacedExisting = nCandidateReplaced;
end


function [final_pop, final_f, finalReport] = rs_final_report_refined_cleanup_A3(archive, candidateArchive, opt, problem)
    % A3 final reporting refinement.
    % ---------------------------------------------------------------
    % This function uses only already evaluated points from:
    %   1) the primary RS-CMSA-ESII archive, and
    %   2) the passive secondary candidate archive.
    % It does not call problem.func_eval(), does not consume extra FEs,
    % and does not use problem.minima.X or problem.minima.f.
    %
    % Difference from v2/A1:
    %   - v2 used a relatively strong distance duplicate filter and kept at
    %     most 10*n_minima secondary candidates.
    %   - A3 keeps substantially more candidate points and removes only
    %     exact/very-near duplicates. This is safer for CEC 2026 RPR because
    %     extra reported points are not penalized; for each true optimum,
    %     post-processing uses the best reported value assigned to it.
    % ---------------------------------------------------------------

    D = problem.dim;
    primaryX = archive.solution;
    primaryF = archive.value(:);
    nPrimary = size(primaryX,1);

    candX = candidateArchive.solution;
    candF = candidateArchive.value(:);

    finalReport = rs_empty_final_report(D);
    finalReport.strategyCode = opt.finalCleanup.strategyCode;

    finiteCand = false(size(candF));
    if ~isempty(candF)
        finiteCand = isfinite(candF) & all(isfinite(candX),2);
        if opt.finalCleanup.keepOnlyFeasibleFinite
            finiteCand = finiteCand & (candF < opt.coreSearch.objValUpLimit);
        end
    end
    nCandidateFinite = sum(finiteCand);

    allKnownF = [primaryF; candF(finiteCand)];
    if isempty(allKnownF)
        final_pop = zeros(0,D);
        final_f = zeros(0,1);
        finalReport.nCandidateRaw = candidateArchive.size;
        finalReport.nCandidateFinite = nCandidateFinite;
        return;
    end

    bestKnownF = min(allKnownF);

    % Value-window filter. A candidate outside this window cannot receive
    % credit under the final RPR tolerance if bestKnownF is close to f*. If
    % bestKnownF is slightly above f*, the window remains conservative.
    candKeep = finiteCand & (candF <= bestKnownF + opt.finalCleanup.valueWindow);
    candInd = find(candKeep);
    nCandidateFiltered = numel(candInd);

    if ~isempty(candInd)
        [~, ord] = sort(candF(candInd), 'ascend');
        candInd = candInd(ord);
    end

    maxReport = max(opt.finalCleanup.minMaxReport, opt.finalCleanup.maxReportMultiplier * max(1, problem.n_minima));
    maxCand = max(0, maxReport - nPrimary);
    nTrimmedByMaxReport = 0;
    if numel(candInd) > maxCand
        nTrimmedByMaxReport = numel(candInd) - maxCand;
        candInd = candInd(1:maxCand);
    end
    nCandidateUsed = numel(candInd);

    combX = [primaryX; candX(candInd,:)];
    combF = [primaryF; candF(candInd)];
    combSource = [ones(nPrimary,1); 2*ones(nCandidateUsed,1)];
    combOrigIndex = [(1:nPrimary)'; candInd(:)];

    if isempty(combF)
        final_pop = zeros(0,D);
        final_f = zeros(0,1);
        finalReport.nCandidateRaw = candidateArchive.size;
        finalReport.nCandidateFinite = nCandidateFinite;
        finalReport.nCandidateFiltered = nCandidateFiltered;
        finalReport.nCandidateUsed = nCandidateUsed;
        finalReport.nTrimmedByMaxReport = nTrimmedByMaxReport;
        return;
    end

    [~, order] = sort(combF, 'ascend');

    exactDupThresh = max(realmin, opt.finalCleanup.exactDupNormDistCoeff * sqrt(D));
    nearDupThresh  = max(exactDupThresh, opt.finalCleanup.nearDupNormDistCoeff * sqrt(D));

    keptX = zeros(0,D);
    keptF = zeros(0,1);
    keptSource = zeros(0,1);
    keptOrigIndex = zeros(0,1);

    nCandidateAdded = 0;
    nCandidateReplaced = 0;
    nExactDuplicateMerged = 0;
    nNearDuplicateMerged = 0;

    lb = problem.low_bound;
    ub = problem.up_bound;
    range = ub - lb;
    range(range == 0) = 1;

    for kk = 1:numel(order)
        ii = order(kk);
        x = combX(ii,:);
        f = combF(ii);
        src = combSource(ii);
        origIdx = combOrigIndex(ii);

        if isempty(keptF)
            keptX = x;
            keptF = f;
            keptSource = src;
            keptOrigIndex = origIdx;
            if src == 2
                nCandidateAdded = nCandidateAdded + 1;
            end
            continue;
        end

        xNorm = (x - lb) ./ range;
        keptNorm = (keptX - lb) ./ range;
        dis = sqrt(sum((keptNorm - xNorm).^2, 2));
        [minDis, minInd] = min(dis);

        if minDis <= exactDupThresh
            nExactDuplicateMerged = nExactDuplicateMerged + 1;
            if f < keptF(minInd)
                if src == 2
                    nCandidateReplaced = nCandidateReplaced + 1;
                end
                keptX(minInd,:) = x;
                keptF(minInd) = f;
                keptSource(minInd) = src;
                keptOrigIndex(minInd) = origIdx;
            end
        elseif minDis <= nearDupThresh
            nNearDuplicateMerged = nNearDuplicateMerged + 1;
            if f < keptF(minInd)
                if src == 2
                    nCandidateReplaced = nCandidateReplaced + 1;
                end
                keptX(minInd,:) = x;
                keptF(minInd) = f;
                keptSource(minInd) = src;
                keptOrigIndex(minInd) = origIdx;
            end
        else
            keptX = [keptX; x]; %#ok<AGROW>
            keptF = [keptF; f]; %#ok<AGROW>
            keptSource = [keptSource; src]; %#ok<AGROW>
            keptOrigIndex = [keptOrigIndex; origIdx]; %#ok<AGROW>
            if src == 2
                nCandidateAdded = nCandidateAdded + 1;
            end
        end
    end

    final_pop = keptX;
    final_f = keptF;

    finalReport.solution = final_pop;
    finalReport.value = final_f;
    finalReport.source = keptSource;
    finalReport.originalIndex = keptOrigIndex;
    finalReport.bestKnownF = bestKnownF;
    finalReport.valueWindow = opt.finalCleanup.valueWindow;
    finalReport.dupNormDist = nearDupThresh;
    finalReport.exactDupNormDist = exactDupThresh;
    finalReport.nearDupNormDist = nearDupThresh;
    finalReport.nPrimary = nPrimary;
    finalReport.nCandidateRaw = candidateArchive.size;
    finalReport.nCandidateFinite = nCandidateFinite;
    finalReport.nCandidateFiltered = nCandidateFiltered;
    finalReport.nCandidateUsed = nCandidateUsed;
    finalReport.nCandidateAddedToReport = nCandidateAdded;
    finalReport.nCandidateReplacedExisting = nCandidateReplaced;
    finalReport.nExactDuplicateMerged = nExactDuplicateMerged;
    finalReport.nNearDuplicateMerged = nNearDuplicateMerged;
    finalReport.nTrimmedByMaxReport = nTrimmedByMaxReport;
    finalReport.strategyCode = opt.finalCleanup.strategyCode;
end


function A6 = rs_evaluate_A6_reporting_rules(archive, candidateArchive, opt, problem, ftol0)
    % Evaluate multiple fixed cleanup/reporting rules from the same candidate pool.
    % This is a diagnostic/development tool for choosing one GLOBAL final
    % reporting policy. It must not be used to choose a rule separately for
    % each competition instance, because that would use post-hoc performance
    % information.

    cfgs = rs_A6_reporting_configs(problem);
    n = numel(cfgs);

    A6 = struct();
    A6.names = cell(n,1);
    A6.metrics = repmat(struct( ...
        'name','', 'strategyCode',NaN, 'nsol',NaN, ...
        'rpr',NaN, 'precision',NaN, 'f1',NaN, 'officialScore',NaN, ...
        'nCandidateFinite',NaN, 'nCandidateFiltered',NaN, 'nCandidateUsed',NaN, ...
        'nCandidateAdded',NaN, 'nCandidateReplaced',NaN, ...
        'nExactDuplicateMerged',NaN, 'nNearDuplicateMerged',NaN, 'nTrimmedByMaxReport',NaN), n, 1);
    A6.reports = cell(n,1);

    for k = 1:n
        cfg = cfgs(k);

        % Use the exact accepted implementations for A1 and A3 so the
        % reported A6 diagnostics are directly comparable with earlier runs.
        if strcmp(cfg.name, 'A1_strict')
            [fp, ff, fr] = rs_final_report_secondary_archive_A1(archive, candidateArchive, opt, problem);
            fr.strategyName = cfg.name;
            fr.strategyCode = cfg.strategyCode;
        elseif strcmp(cfg.name, 'A3_relaxed')
            [fp, ff, fr] = rs_final_report_refined_cleanup_A3(archive, candidateArchive, opt, problem);
            fr.strategyName = cfg.name;
            fr.strategyCode = cfg.strategyCode;
        else
            [fp, ff, fr] = rs_final_report_custom_cleanup_A6(archive, candidateArchive, opt, problem, cfg);
        end

        rpr = NaN;
        ipr = [];
        if ~isempty(fp)
            try
                [rpr, ipr] = UtilityMethod.calc_robust_peak_ratio(fp, ff, problem.minima.X, problem.minima.f, ftol0); %#ok<NASGU>
            catch ME
                warning('A6 reporting rule %s RPR calculation failed: %s', cfg.name, ME.message);
            end
        end

        nsol = size(fp,1);
        [precision, f1, score] = rs_submission_metrics_from_rpr(rpr, nsol, problem.n_minima);

        A6.names{k} = cfg.name;
        A6.reports{k} = fr;
        A6.metrics(k).name = cfg.name;
        A6.metrics(k).strategyCode = cfg.strategyCode;
        A6.metrics(k).nsol = nsol;
        A6.metrics(k).rpr = rpr;
        A6.metrics(k).precision = precision;
        A6.metrics(k).f1 = f1;
        A6.metrics(k).officialScore = score;
        A6.metrics(k).nCandidateFinite = getfield_safe(fr, 'nCandidateFinite', NaN);
        A6.metrics(k).nCandidateFiltered = getfield_safe(fr, 'nCandidateFiltered', NaN);
        A6.metrics(k).nCandidateUsed = getfield_safe(fr, 'nCandidateUsed', NaN);
        A6.metrics(k).nCandidateAdded = getfield_safe(fr, 'nCandidateAddedToReport', NaN);
        A6.metrics(k).nCandidateReplaced = getfield_safe(fr, 'nCandidateReplacedExisting', NaN);
        A6.metrics(k).nExactDuplicateMerged = getfield_safe(fr, 'nExactDuplicateMerged', NaN);
        A6.metrics(k).nNearDuplicateMerged = getfield_safe(fr, 'nNearDuplicateMerged', NaN);
        A6.metrics(k).nTrimmedByMaxReport = getfield_safe(fr, 'nTrimmedByMaxReport', NaN);
    end

    % Development-only best rule for this run. Do not use this per-run best
    % selection for final competition submission. Use the batch results to
    % choose one global rule.
    scores = arrayfun(@(s) s.officialScore, A6.metrics);
    [A6.bestScoreDevelopmentOnly, A6.bestIndexDevelopmentOnly] = max(scores);
    if isempty(A6.bestIndexDevelopmentOnly) || isnan(A6.bestScoreDevelopmentOnly)
        A6.bestIndexDevelopmentOnly = NaN;
        A6.bestNameDevelopmentOnly = '';
    else
        A6.bestNameDevelopmentOnly = A6.metrics(A6.bestIndexDevelopmentOnly).name;
    end
end

function cfgs = rs_A6_reporting_configs(problem)
    % Fixed, global reporting-rule candidates. These are not PID-specific.
    % Some use problem.n_minima only as a generous cap; rules marked
    % primaryPlus do not depend on NGM for capping.
    unusedD = problem.dim; %#ok<NASGU>

    cfgs = repmat(struct( ...
        'name','', 'strategyCode',NaN, ...
        'valueWindow',1.0, ...
        'exactDupNormDistCoeff',1e-10, ...
        'nearDupNormDistCoeff',5e-4, ...
        'maxReportMultiplier',100, ...
        'minMaxReport',500, ...
        'maxReportMode','nMinima', ...
        'primaryCapMultiplier',1.25, ...
        'primaryCapAdd',5, ...
        'maxAbsReport',500, ...
        'keepOnlyFeasibleFinite',true), 6, 1);

    % Rule 1: accepted A1/v2 strict cleanup.
    cfgs(1).name = 'A1_strict';
    cfgs(1).strategyCode = 1;
    cfgs(1).valueWindow = 1.0;
    cfgs(1).nearDupNormDistCoeff = 0.005;
    cfgs(1).exactDupNormDistCoeff = 1e-10;
    cfgs(1).maxReportMultiplier = 10;
    cfgs(1).minMaxReport = 0;
    cfgs(1).maxReportMode = 'nMinima';

    % Rule 2: accepted A3b relaxed cleanup.
    cfgs(2).name = 'A3_relaxed';
    cfgs(2).strategyCode = 3;
    cfgs(2).valueWindow = 1.0;
    cfgs(2).nearDupNormDistCoeff = 5e-4;
    cfgs(2).exactDupNormDistCoeff = 1e-10;
    cfgs(2).maxReportMultiplier = 100;
    cfgs(2).minMaxReport = 500;
    cfgs(2).maxReportMode = 'nMinima';

    % Rule 3: medium cleanup, intended to balance RPR and precision/F1.
    cfgs(3).name = 'A6_medium';
    cfgs(3).strategyCode = 61;
    cfgs(3).valueWindow = 0.75;
    cfgs(3).nearDupNormDistCoeff = 1e-3;
    cfgs(3).exactDupNormDistCoeff = 1e-10;
    cfgs(3).maxReportMultiplier = 50;
    cfgs(3).minMaxReport = 300;
    cfgs(3).maxReportMode = 'nMinima';

    % Rule 4: value-filtered relaxed cleanup.
    cfgs(4).name = 'A6_value05';
    cfgs(4).strategyCode = 62;
    cfgs(4).valueWindow = 0.5;
    cfgs(4).nearDupNormDistCoeff = 5e-4;
    cfgs(4).exactDupNormDistCoeff = 1e-10;
    cfgs(4).maxReportMultiplier = 100;
    cfgs(4).minMaxReport = 500;
    cfgs(4).maxReportMode = 'nMinima';

    % Rule 5: precision-aware cap relative to primary archive size; does
    % not use the true number of global optima as a cap.
    cfgs(5).name = 'A6_primaryCap125';
    cfgs(5).strategyCode = 63;
    cfgs(5).valueWindow = 1.0;
    cfgs(5).nearDupNormDistCoeff = 5e-4;
    cfgs(5).exactDupNormDistCoeff = 1e-10;
    cfgs(5).primaryCapMultiplier = 1.25;
    cfgs(5).primaryCapAdd = 5;
    cfgs(5).maxAbsReport = 500;
    cfgs(5).maxReportMode = 'primaryPlus';

    % Rule 6: slightly stricter distance with still-generous reporting.
    cfgs(6).name = 'A6_dist002';
    cfgs(6).strategyCode = 64;
    cfgs(6).valueWindow = 1.0;
    cfgs(6).nearDupNormDistCoeff = 0.002;
    cfgs(6).exactDupNormDistCoeff = 1e-10;
    cfgs(6).maxReportMultiplier = 30;
    cfgs(6).minMaxReport = 300;
    cfgs(6).maxReportMode = 'nMinima';
end

function [precision, f1, score] = rs_submission_metrics_from_rpr(rpr, nsol, nGlobalMinima)
    if isnan(rpr) || isempty(rpr) || nsol <= 0 || nGlobalMinima <= 0
        precision = NaN;
        f1 = NaN;
        score = NaN;
        return;
    end
    precision = rpr * nGlobalMinima / nsol;
    precision = max(0, min(1, precision));
    recall = rpr;
    if precision + recall <= 0
        f1 = 0;
    else
        f1 = 2 * precision * recall / (precision + recall);
    end
    score = 0.5 * (rpr + f1);
end

function f1 = rs_f1_from_rpr_nsol(rpr, nsol, nGlobalMinima)
    [~, f1, ~] = rs_submission_metrics_from_rpr(rpr, nsol, nGlobalMinima);
end

function [final_pop, final_f, finalReport] = rs_final_report_custom_cleanup_A6(archive, candidateArchive, opt, problem, cfg)
    % Generic A6 final cleanup/reporting rule.
    % It uses only already evaluated primary archive and passive candidate
    % archive points. It does not call problem.func_eval() and does not use
    % known minima locations.

    D = problem.dim;
    primaryX = archive.solution;
    primaryF = archive.value(:);
    nPrimary = size(primaryX,1);

    candX = candidateArchive.solution;
    candF = candidateArchive.value(:);

    finalReport = rs_empty_final_report(D);
    finalReport.strategyCode = cfg.strategyCode;
    finalReport.strategyName = cfg.name;

    finiteCand = false(size(candF));
    if ~isempty(candF)
        finiteCand = isfinite(candF) & all(isfinite(candX),2);
        if cfg.keepOnlyFeasibleFinite
            finiteCand = finiteCand & (candF < opt.coreSearch.objValUpLimit);
        end
    end
    nCandidateFinite = sum(finiteCand);

    allKnownF = [primaryF; candF(finiteCand)];
    if isempty(allKnownF)
        final_pop = zeros(0,D);
        final_f = zeros(0,1);
        finalReport.nCandidateRaw = candidateArchive.size;
        finalReport.nCandidateFinite = nCandidateFinite;
        return;
    end
    bestKnownF = min(allKnownF);

    candKeep = finiteCand & (candF <= bestKnownF + cfg.valueWindow);
    candInd = find(candKeep);
    nCandidateFiltered = numel(candInd);

    if ~isempty(candInd)
        [~, ord] = sort(candF(candInd), 'ascend');
        candInd = candInd(ord);
    end

    switch lower(cfg.maxReportMode)
        case 'primaryplus'
            maxReport = max(nPrimary, ceil(cfg.primaryCapMultiplier * max(1,nPrimary)) + cfg.primaryCapAdd);
            maxReport = min(maxReport, cfg.maxAbsReport);
        otherwise
            maxReport = max(cfg.minMaxReport, cfg.maxReportMultiplier * max(1, problem.n_minima));
    end

    maxCand = max(0, maxReport - nPrimary);
    nTrimmedByMaxReport = 0;
    if numel(candInd) > maxCand
        nTrimmedByMaxReport = numel(candInd) - maxCand;
        candInd = candInd(1:maxCand);
    end
    nCandidateUsed = numel(candInd);

    combX = [primaryX; candX(candInd,:)];
    combF = [primaryF; candF(candInd)];
    combSource = [ones(nPrimary,1); 2*ones(nCandidateUsed,1)];
    combOrigIndex = [(1:nPrimary)'; candInd(:)];

    if isempty(combF)
        final_pop = zeros(0,D);
        final_f = zeros(0,1);
        finalReport.nCandidateRaw = candidateArchive.size;
        finalReport.nCandidateFinite = nCandidateFinite;
        finalReport.nCandidateFiltered = nCandidateFiltered;
        finalReport.nCandidateUsed = nCandidateUsed;
        finalReport.nTrimmedByMaxReport = nTrimmedByMaxReport;
        return;
    end

    [~, order] = sort(combF, 'ascend');

    exactDupThresh = max(realmin, cfg.exactDupNormDistCoeff * sqrt(D));
    nearDupThresh  = max(exactDupThresh, cfg.nearDupNormDistCoeff * sqrt(D));

    keptX = zeros(0,D);
    keptF = zeros(0,1);
    keptSource = zeros(0,1);
    keptOrigIndex = zeros(0,1);

    nCandidateAdded = 0;
    nCandidateReplaced = 0;
    nExactDuplicateMerged = 0;
    nNearDuplicateMerged = 0;

    lb = problem.low_bound;
    ub = problem.up_bound;
    range = ub - lb;
    range(range == 0) = 1;

    for kk = 1:numel(order)
        ii = order(kk);
        x = combX(ii,:);
        f = combF(ii);
        src = combSource(ii);
        origIdx = combOrigIndex(ii);

        if isempty(keptF)
            keptX = x;
            keptF = f;
            keptSource = src;
            keptOrigIndex = origIdx;
            if src == 2
                nCandidateAdded = nCandidateAdded + 1;
            end
            continue;
        end

        xNorm = (x - lb) ./ range;
        keptNorm = (keptX - lb) ./ range;
        dis = sqrt(sum((keptNorm - xNorm).^2, 2));
        [minDis, minInd] = min(dis);

        if minDis <= exactDupThresh
            nExactDuplicateMerged = nExactDuplicateMerged + 1;
            if f < keptF(minInd)
                if src == 2
                    nCandidateReplaced = nCandidateReplaced + 1;
                end
                keptX(minInd,:) = x;
                keptF(minInd) = f;
                keptSource(minInd) = src;
                keptOrigIndex(minInd) = origIdx;
            end
        elseif minDis <= nearDupThresh
            nNearDuplicateMerged = nNearDuplicateMerged + 1;
            if f < keptF(minInd)
                if src == 2
                    nCandidateReplaced = nCandidateReplaced + 1;
                end
                keptX(minInd,:) = x;
                keptF(minInd) = f;
                keptSource(minInd) = src;
                keptOrigIndex(minInd) = origIdx;
            end
        else
            keptX = [keptX; x]; %#ok<AGROW>
            keptF = [keptF; f]; %#ok<AGROW>
            keptSource = [keptSource; src]; %#ok<AGROW>
            keptOrigIndex = [keptOrigIndex; origIdx]; %#ok<AGROW>
            if src == 2
                nCandidateAdded = nCandidateAdded + 1;
            end
        end
    end

    final_pop = keptX;
    final_f = keptF;

    finalReport.solution = final_pop;
    finalReport.value = final_f;
    finalReport.source = keptSource;
    finalReport.originalIndex = keptOrigIndex;
    finalReport.bestKnownF = bestKnownF;
    finalReport.valueWindow = cfg.valueWindow;
    finalReport.dupNormDist = nearDupThresh;
    finalReport.exactDupNormDist = exactDupThresh;
    finalReport.nearDupNormDist = nearDupThresh;
    finalReport.nPrimary = nPrimary;
    finalReport.nCandidateRaw = candidateArchive.size;
    finalReport.nCandidateFinite = nCandidateFinite;
    finalReport.nCandidateFiltered = nCandidateFiltered;
    finalReport.nCandidateUsed = nCandidateUsed;
    finalReport.nCandidateAddedToReport = nCandidateAdded;
    finalReport.nCandidateReplacedExisting = nCandidateReplaced;
    finalReport.nExactDuplicateMerged = nExactDuplicateMerged;
    finalReport.nNearDuplicateMerged = nNearDuplicateMerged;
    finalReport.nTrimmedByMaxReport = nTrimmedByMaxReport;
    finalReport.strategyCode = cfg.strategyCode;
    finalReport.strategyName = cfg.name;
end


function finalReport = rs_empty_final_report(D)
    finalReport.solution = zeros(0,D);
    finalReport.value = zeros(0,1);
    finalReport.source = zeros(0,1);
    finalReport.originalIndex = zeros(0,1);
    finalReport.bestKnownF = NaN;
    finalReport.valueWindow = NaN;
    finalReport.dupNormDist = NaN;
    finalReport.exactDupNormDist = NaN;
    finalReport.nearDupNormDist = NaN;
    finalReport.nPrimary = 0;
    finalReport.nCandidateRaw = 0;
    finalReport.nCandidateFinite = 0;
    finalReport.nCandidateFiltered = 0;
    finalReport.nCandidateUsed = 0;
    finalReport.nCandidateAddedToReport = 0;
    finalReport.nCandidateReplacedExisting = 0;
    finalReport.nExactDuplicateMerged = 0;
    finalReport.nNearDuplicateMerged = 0;
    finalReport.nTrimmedByMaxReport = 0;
    finalReport.strategyCode = NaN;
end

% =========================================================================
% RESTART
% =========================================================================
function restart = rs_restart_init(process, opt, problem)
    restart.stagSize = floor(opt.stopCr.stagPar(1) + opt.stopCr.stagPar(3)*problem.dim/process.subpopSize);
    restart.tolHistSize = floor(10 + 30.0*problem.dim/process.subpopSize);
    restart.usedEvalEvolve = 0;
    restart.usedEvalMerge  = 0;
    restart.iterNo = 0;
    restart.terminationFlag = 0;
    restart.subpopTermFlag  = 0;
    restart.subpopBestVal   = inf;
    restart.recIniR0 = process.iniR0;
end

function subpop = rs_initialize_subpop(archive, process, opt, problem, restart) %#ok<INUSD>
    D = problem.dim;
    archSolRescaled = zeros(archive.size, D);
    for k = 1:archive.size
        archSolRescaled(k,:) = (archive.solution(k,:) - problem.low_bound) ./ (problem.up_bound - problem.low_bound);
    end

    numReject = 0;
    R0 = process.iniR0;
    wasSuccess = false;

    while ~wasSuccess
        X = rand(1,D);
        chkDis = true;
        if archive.size > 0
            dis2dis = sqrt(sum((archSolRescaled - X).^2, 2));
            chkDis = all(dis2dis > archive.normTabDis * R0);
        end

        if chkDis
            wasSuccess = true;
            recIniR0 = R0;
        else
            numReject = numReject + 1;
        end

        if numReject > opt.niching.maxRejectIni
            numReject = 0;
            R0 = R0 * opt.niching.redCoeff;
        end
    end

    center = X .* (problem.up_bound - problem.low_bound) + problem.low_bound;
    smean  = min(opt.coreSearch.maxIniSigma, R0 * opt.coreSearch.iniSigCoeff);
    stretch = problem.up_bound - problem.low_bound;

    subpop = rs_subpop_init(center, smean, stretch, process.subpopSize);
    subpop.recIniR0_from_init = recIniR0;
end

function [restart, subpop, archive] = rs_run_one_restart(subpop, archive, process, opt, problem, restart)
    restart.recIniR0 = subpop.recIniR0_from_init;
    while restart.terminationFlag == 0
        restart.iterNo = restart.iterNo + 1;
        restart.stagSize = floor(opt.stopCr.stagPar(1) + ...
            floor(opt.stopCr.stagPar(2)*restart.iterNo + opt.stopCr.stagPar(3)*problem.dim/process.subpopSize));

        subpop = rs_update_taboo_region(subpop, archive, opt, problem);
        subpop = rs_update_merge_check(subpop, archive, opt, problem);
        [subpop, restart] = rs_evolve(subpop, restart, archive, process, opt, problem);
        [subpop, restart] = rs_update_term_flag(subpop, restart, archive, process, opt, problem);

        restart.subpopBestVal = subpop.bestVal;
        restart.subpopTermFlag = subpop.terminationFlag;

        remainEvalAfter = problem.max_eval - (process.usedEvalTillRestart + ...
            restart.usedEvalEvolve + restart.usedEvalMerge + process.subpopSize);
        reqEvalForDetectMult = min(archive.size, opt.archiving.neighborSize) * opt.archiving.hillVallBudget;

        if remainEvalAfter < reqEvalForDetectMult
            restart.terminationFlag = -1;
            break;
        end

        if restart.terminationFlag == 0 && restart.subpopTermFlag ~= 0
            restart.terminationFlag = 1;
        end

        if problem.used_eval >= problem.max_eval
            restart.terminationFlag = -1;
            break;
        end
    end
end

% =========================================================================
% SUBPOPULATION
% =========================================================================
function subpop = rs_subpop_init(center, smean, stretch, popSize)
    D = numel(center);
    subpop.center = center(:)';
    subpop.mutProfile.smean   = smean;
    subpop.mutProfile.stretch = stretch(:)';
    subpop.mutProfile.C       = diag(subpop.mutProfile.stretch.^2);
    subpop.mutProfile.Cinv    = diag(subpop.mutProfile.stretch.^(-2));
    subpop.mutProfile.rotMat  = eye(D);

    subpop.samples = rs_sampling_init(popSize, D);
    subpop.bestSol = zeros(0,D);
    subpop.bestVal = inf;

    subpop.tabooRegion.center = zeros(0,D);
    subpop.tabooRegion.normTabDis = zeros(0,1);
    subpop.tabooRegion.criticality = zeros(0,1);
    subpop.tabooRegion.criticInd = zeros(0,1);

    subpop.mergeCheck.bestCandidArchIndHist = zeros(0,1);
    subpop.mergeCheck.bestCandidArchMergeabilityHist = zeros(0,1);
    subpop.mergeCheck.candidArchCountHist = zeros(0,1);
    subpop.mergeCheck.matchArchInd = NaN;
    subpop.mergeCheck.checkAfterIterNo = 0;
    subpop.mergeCheck.mergeAtUsedEval = inf;

    subpop.localConvergeCheck.stopAtUsedEval = inf;

    subpop.bestValNonEliteHist = zeros(0,1);
    subpop.medValNonEliteHist  = zeros(0,1);
    subpop.maxCriticalityHist  = zeros(0,1);
    subpop.terminationFlag = 0;
    subpop.usedEval = 0;
    subpop.iterNo = 0;

    subpop.elite.sol = zeros(0,D);
    subpop.elite.val = zeros(0,1);
    subpop.elite.s   = zeros(0,1);
    subpop.elite.Z   = zeros(0,D);
    subpop.elite.wasRepaired = false(0,1);
end

function samples = rs_sampling_init(popSize, D)
    samples.s = zeros(popSize,1);
    samples.X = zeros(popSize,D);
    samples.Z = zeros(popSize,D);
    samples.f = inf(popSize,1);
    samples.isFeas = false(popSize,1);
    samples.wasRepaired = false(popSize,1);
    samples.argsortNoElite = [];
    samples.argsortWithElite = [];
end

function normDis = rs_calc_norm_dis(subpop, x1, x2, disMetric)
    if strcmpi(disMetric, 'Mahalanobis')
        dx = x1(:)' - x2(:)';
        normDis = sqrt(max(dx * subpop.mutProfile.Cinv * dx', 0)) / subpop.mutProfile.smean;
    elseif strcmpi(disMetric, 'Euclidean')
        dx = x1(:)' - x2(:)';
        mean_str = exp(mean(log(subpop.mutProfile.stretch)));
        normDis = sqrt(sum(dx.^2)) / (subpop.mutProfile.smean * mean_str);
    else
        error('Unknown distance metric.');
    end
end

function tabAccept = rs_is_taboo_acceptable(subpop, sample, tempRedRatio)
    tabAccept = true;
    for ii = 1:numel(subpop.tabooRegion.criticInd)
        tabInd = subpop.tabooRegion.criticInd(ii);
        normDis = rs_calc_norm_dis(subpop, sample, subpop.tabooRegion.center(tabInd,:), 'Mahalanobis');
        tabAccept = normDis >= (subpop.tabooRegion.normTabDis(tabInd) * tempRedRatio);
        if ~tabAccept
            break;
        end
    end
end

function subpop = rs_update_taboo_region(subpop, archive, opt, problem)
    maxCount = archive.size;
    centers = zeros(maxCount, problem.dim);
    normTab = zeros(maxCount,1);
    count = 0;

    for k = 1:maxCount
        if archive.value(k) < subpop.bestVal
            count = count + 1;
            centers(count,:) = archive.solution(k,:);
            normTab(count) = archive.normTabDis(k);
        end
    end

    subpop.tabooRegion.center = centers(1:count,:);
    subpop.tabooRegion.normTabDis = normTab(1:count);
    subpop.tabooRegion.criticality = zeros(count,1);
    subpop.tabooRegion.criticInd = zeros(0,1);

    if count > 0
        L = zeros(count,1);
        for k = 1:count
            L(k) = rs_calc_norm_dis(subpop, subpop.tabooRegion.center(k,:), subpop.center, 'Mahalanobis');
        end
        intU = L + subpop.tabooRegion.normTabDis;
        intL = L - subpop.tabooRegion.normTabDis;
        subpop.tabooRegion.criticality = rs_normcdf(intU) - rs_normcdf(intL);
        [~, crInd] = sort(subpop.tabooRegion.criticality, 'descend');
        Ncr = sum(subpop.tabooRegion.criticality > opt.niching.criticTabooThresh);
        subpop.tabooRegion.criticInd = crInd(1:Ncr);
    end

    subpop.maxCriticalityHist = [subpop.maxCriticalityHist; max([subpop.tabooRegion.criticality; 0])]; %#ok<AGROW>
end

function subpop = rs_update_merge_check(subpop, archive, opt, problem)
    if archive.size > 0
        L = zeros(archive.size,1);
        for k = 1:archive.size
            L(k) = rs_calc_norm_dis(subpop, archive.solution(k,:), subpop.center, 'Mahalanobis');
        end
        Mergeability = (archive.normTabDis + opt.stopCr.merge.addedConst) ./ max(L, realmin);
        [maxMergeability, indMax] = max(Mergeability);
        candidArchCount = sum(Mergeability > opt.stopCr.merge.threshold);
    else
        indMax = -1;
        maxMergeability = 0;
        candidArchCount = 0;
    end

    subpop.mergeCheck.bestCandidArchIndHist = [subpop.mergeCheck.bestCandidArchIndHist; indMax]; %#ok<AGROW>
    subpop.mergeCheck.bestCandidArchMergeabilityHist = [subpop.mergeCheck.bestCandidArchMergeabilityHist; maxMergeability]; %#ok<AGROW>
    subpop.mergeCheck.candidArchCountHist = [subpop.mergeCheck.candidArchCountHist; candidArchCount]; %#ok<AGROW>
end

function [subpop, restart] = rs_evolve(subpop, restart, archive, process, opt, problem)
    [subpop, restart] = rs_sample_and_eval(subpop, restart, archive, process, opt, problem);
    subpop = rs_select(subpop, process, opt);
    subpop = rs_recombine_CMSA(subpop, process, opt, problem);
end

function [subpop, restart] = rs_sample_and_eval(subpop, restart, archive, process, opt, problem) %#ok<INUSD>
    tempRedRatio = 1;
    solNo = 1;
    D = problem.dim;
    popSize = process.subpopSize;

    subpop.samples = rs_sampling_init(popSize, D);

    while solNo <= popSize
        subpop.samples.s(solNo) = subpop.mutProfile.smean * exp(randn() * process.coreSpecStrPar.tauSigma);
        zcol = subpop.mutProfile.rotMat * (subpop.mutProfile.stretch(:) .* randn(D,1));
        subpop.samples.Z(solNo,:) = zcol(:)';
        subpop.samples.X(solNo,:) = subpop.center + subpop.samples.s(solNo) * subpop.samples.Z(solNo,:);
        [subpop, solNo] = rs_repair_infeas(subpop, solNo, problem);

        acceptIt = rs_is_taboo_acceptable(subpop, subpop.samples.X(solNo,:), tempRedRatio);
        if acceptIt
            solNo = solNo + 1;
        else
            tempRedRatio = tempRedRatio * opt.niching.redCoeff;
        end
    end

    for k = 1:popSize
        penU = subpop.samples.X(k,:) - problem.up_bound;
        penU = penU .* (penU > 0);
        penL = problem.low_bound - subpop.samples.X(k,:);
        penL = penL .* (penL > 0);
        penUL = sum(penU + penL);
        subpop.samples.isFeas(k) = ~(penUL > 0);
        if ~subpop.samples.isFeas(k)
            subpop.samples.f(k) = opt.coreSearch.objValUpLimit * (1 + penUL);
        else
            if problem.used_eval < problem.max_eval
                subpop.samples.f(k) = problem.func_eval(subpop.samples.X(k,:));
                restart.usedEvalEvolve = restart.usedEvalEvolve + 1;
                subpop.usedEval = subpop.usedEval + 1;
            else
                subpop.samples.f(k) = opt.coreSearch.objValUpLimit;
            end
        end
    end

    subpop.iterNo = subpop.iterNo + 1;
    subpop.bestValNonEliteHist = [subpop.bestValNonEliteHist; min(subpop.samples.f)]; %#ok<AGROW>
    subpop.medValNonEliteHist  = [subpop.medValNonEliteHist; median(subpop.samples.f)]; %#ok<AGROW>

    if numel(subpop.bestValNonEliteHist) > restart.stagSize
        subpop.bestValNonEliteHist = subpop.bestValNonEliteHist(end-restart.stagSize+1:end);
        subpop.medValNonEliteHist  = subpop.medValNonEliteHist(end-restart.stagSize+1:end);
    end
end

function [subpop, solNo] = rs_repair_infeas(subpop, solNo, problem)
    x = subpop.samples.X(solNo,:);
    if any(problem.up_bound < x) || any(problem.low_bound > x)
        subpop.samples.wasRepaired(solNo) = true;

        relocItU = find(x > problem.up_bound);
        if ~isempty(relocItU)
            tmpUp  = problem.up_bound(relocItU);
            tmpLow = 2*subpop.center(relocItU) - problem.up_bound(relocItU);
            tmpLow = max(tmpLow, problem.low_bound(relocItU));
            x(relocItU) = tmpLow + rand(1,numel(relocItU)) .* (tmpUp - tmpLow);
        end

        relocItL = find(x < problem.low_bound);
        if ~isempty(relocItL)
            tmpLow = problem.low_bound(relocItL);
            tmpUp  = 2*subpop.center(relocItL) - problem.low_bound(relocItL);
            tmpUp  = min(tmpUp, problem.up_bound(relocItL));
            x(relocItL) = tmpLow + rand(1,numel(relocItL)) .* (tmpUp - tmpLow);
        end

        subpop.samples.X(solNo,:) = x;
        subpop.samples.Z(solNo,:) = (subpop.samples.X(solNo,:) - subpop.center) / subpop.samples.s(solNo);
    end
end

function subpop = rs_select(subpop, process, opt)
    [~, subpop.samples.argsortNoElite] = sort(subpop.samples.f, 'ascend');

    if strcmp(opt.coreSearch.algorithm, 'CMSA') && process.coreSpecStrPar.numElt > 0 && ~isempty(subpop.elite.val)
        appendIt = true(numel(subpop.elite.val),1);
        for eltNo = 1:numel(subpop.elite.val)
            appendIt(eltNo) = rs_is_taboo_acceptable(subpop, subpop.elite.sol(eltNo,:), 1);
        end
        subpop.samples.X = [subpop.samples.X; subpop.elite.sol(appendIt,:)];
        subpop.samples.Z = [subpop.samples.Z; subpop.elite.Z(appendIt,:)];
        subpop.samples.s = [subpop.samples.s; subpop.elite.s(appendIt)];
        subpop.samples.f = [subpop.samples.f; subpop.elite.val(appendIt)];
        subpop.samples.wasRepaired = [subpop.samples.wasRepaired; subpop.elite.wasRepaired(appendIt)];
    end

    [~, subpop.samples.argsortWithElite] = sort(subpop.samples.f, 'ascend');
end

function subpop = rs_recombine_CMSA(subpop, process, opt, problem)
    oldCenter = subpop.center;
    ind = subpop.samples.argsortWithElite;

    muInd = ind(1:process.mu);
    subpop.center = process.recWeights * subpop.samples.X(muInd,:);
    subpop.center = min(max(subpop.center, problem.low_bound), problem.up_bound);

    subpop.bestSol = subpop.samples.X(ind(1),:);
    subpop.bestVal = subpop.samples.f(ind(1));

    oldSmean = subpop.mutProfile.smean;
    correctionTerm = (rs_geomean(subpop.samples.s) / oldSmean)^opt.coreSearch.sigmaUpdateBiasImp;
    subpop.mutProfile.smean = exp(sum(process.recWeights(:) .* log(max(subpop.samples.s(muInd), realmin)))) / correctionTerm;

    suggC = zeros(problem.dim, problem.dim);
    for parNo = 1:process.mu
        z = (subpop.samples.X(muInd(parNo),:) - oldCenter) / subpop.samples.s(muInd(parNo));
        suggC = suggC + process.recWeights(parNo) * (z(:) * z(:)');
    end
    cc = 1 / process.coreSpecStrPar.tauCov;
    newC = (1-cc) * subpop.mutProfile.C + cc * suggC;
    subpop.mutProfile.C = 0.5 * (newC + newC');

    [V, E] = eig(subpop.mutProfile.C, 'vector');
    E = real(E);
    E(E < realmin) = realmin;
    subpop.mutProfile.rotMat = real(V);
    subpop.mutProfile.stretch = sqrt(E(:))';
    subpop.mutProfile.Cinv = subpop.mutProfile.rotMat * diag(subpop.mutProfile.stretch.^(-2)) * subpop.mutProfile.rotMat';

    if process.coreSpecStrPar.numElt > 0
        [~, surviveInd] = sort(subpop.samples.f, 'ascend');
        limit1 = sum(subpop.samples.isFeas) + numel(subpop.elite.val);
        limit2 = process.coreSpecStrPar.numElt;
        actEltNum = floor(0.5 + min(limit1, limit2));
        actEltNum = min(actEltNum, numel(surviveInd));
        if actEltNum > 0
            keep = surviveInd(1:actEltNum);
            subpop.elite.sol = subpop.samples.X(keep,:);
            subpop.elite.val = subpop.samples.f(keep);
            subpop.elite.Z   = subpop.samples.Z(keep,:);
            subpop.elite.s   = subpop.samples.s(keep);
            subpop.elite.wasRepaired = subpop.samples.wasRepaired(keep);
        else
            D = problem.dim;
            subpop.elite.sol = zeros(0,D);
            subpop.elite.val = zeros(0,1);
            subpop.elite.Z   = zeros(0,D);
            subpop.elite.s   = zeros(0,1);
            subpop.elite.wasRepaired = false(0,1);
        end
    end
end

function [subpop, restart] = rs_update_term_flag(subpop, restart, archive, process, opt, problem)
    condC = (max(subpop.mutProfile.stretch) / max(min(subpop.mutProfile.stretch), realmin))^2;
    if condC > opt.stopCr.maxCondC
        subpop.terminationFlag = -2;
    end

    if subpop.iterNo >= restart.stagSize && subpop.terminationFlag == 0
        histN = numel(subpop.bestValNonEliteHist);
        if histN >= restart.stagSize && histN >= 20
            firstInd = 1:20;
            lastInd  = histN-19:histN;
            minImpBest = median(subpop.bestValNonEliteHist(lastInd)) - median(subpop.bestValNonEliteHist(firstInd));
            minImpMed  = median(subpop.medValNonEliteHist(lastInd))  - median(subpop.medValNonEliteHist(firstInd));
            if min(minImpBest, minImpMed) > 0
                subpop.terminationFlag = -1;
            end
        end
    end

    if subpop.iterNo >= restart.tolHistSize && subpop.terminationFlag == 0
        recent = subpop.bestValNonEliteHist(end-restart.tolHistSize+1:end);
        maxDiff = max(recent) - min(recent);
        if maxDiff < opt.stopCr.tolHistFun
            subpop.terminationFlag = 1;
        end
    end

    if (max(subpop.mutProfile.stretch) * subpop.mutProfile.smean) < opt.stopCr.tolX && subpop.terminationFlag == 0
        subpop.terminationFlag = 2;
    end

    % Merge predictor.
    if subpop.terminationFlag == 0
        windowCount = floor(1 + opt.stopCr.merge.windowSizeCoeff * restart.tolHistSize);
        st = max(1, numel(subpop.mergeCheck.candidArchCountHist) - windowCount + 1);
        chk1 = subpop.iterNo >= (opt.stopCr.merge.windowSizeCoeff * restart.tolHistSize);
        chk2 = ~isempty(subpop.mergeCheck.candidArchCountHist) && ...
               all(subpop.mergeCheck.candidArchCountHist(st:end) == 1);
        chk3 = subpop.bestVal < opt.coreSearch.objValUpLimit;
        chk4 = subpop.iterNo > subpop.mergeCheck.checkAfterIterNo;
        chk5 = isinf(subpop.mergeCheck.mergeAtUsedEval);

        if chk1 && chk2 && chk3 && chk4 && chk5
            candidateIndex = subpop.mergeCheck.bestCandidArchIndHist(end);
            if candidateIndex >= 1 && candidateIndex <= archive.size
                usedEval = 0;
                isNew = false;
                rgrid = linspace(0.5/opt.stopCr.merge.maxEval, ...
                    1 - 0.5/opt.stopCr.merge.maxEval, opt.stopCr.merge.maxEval);

                endX1 = subpop.bestSol;
                endF1 = subpop.bestVal;
                endX2 = archive.solution(candidateIndex,:);
                endF2 = archive.value(candidateIndex);

                while usedEval < opt.stopCr.merge.maxEval && problem.used_eval < problem.max_eval
                    usedEval = usedEval + 1;
                    testX = rgrid(usedEval) * endX1 + (1-rgrid(usedEval)) * endX2;
                    testF = problem.func_eval(testX);
                    subpop.usedEval = subpop.usedEval + 1;
                    restart.usedEvalMerge = restart.usedEvalMerge + 1;
                    if testF > (max(endF1, endF2) + opt.stopCr.tolHistFun)
                        isNew = true;
                        break;
                    end
                end

                if ~isNew
                    if isnan(subpop.mergeCheck.matchArchInd)
                        subpop.mergeCheck.matchArchInd = candidateIndex;
                        subpop.mergeCheck.mergeAtUsedEval = subpop.usedEval;
                        subpop.terminationFlag = 3;
                        subpop.bestSol = archive.solution(candidateIndex,:);
                        subpop.bestVal = archive.value(candidateIndex);
                    end
                else
                    chkInterval = opt.stopCr.merge.chkIntervalCoeff * restart.tolHistSize;
                    subpop.mergeCheck.checkAfterIterNo = chkInterval + subpop.iterNo;
                end
            end
        end
    end

    % Local convergence predictor.
    if subpop.terminationFlag == 0
        ws = floor(2 + opt.stopCr.localConverge.windowSizeCoeff * restart.tolHistSize);
        histN = numel(subpop.bestValNonEliteHist);
        st = max(1, histN - ws + 1);
        stCrit = max(1, numel(subpop.maxCriticalityHist) - ws + 1);

        chk1 = subpop.iterNo > (2 + max(opt.stopCr.localConverge.windowSizeCoeff * restart.tolHistSize, ws));
        chk2 = archive.size > 0;
        chk3 = ~isempty(subpop.maxCriticalityHist) && ...
               max(subpop.maxCriticalityHist(stCrit:end)) < opt.stopCr.localConverge.tabooCriticUpLimit;
        chk4 = subpop.bestVal < opt.coreSearch.objValUpLimit;
        chk5 = isinf(subpop.localConvergeCheck.stopAtUsedEval);

        if chk1 && chk2 && chk3 && chk4 && chk5 && histN >= 2
            recent = subpop.bestValNonEliteHist(st:end);
            meanDiff = mean(abs(diff(recent)));
            willBeLocal = meanDiff < (opt.stopCr.localConverge.tolCoeff * ...
                (subpop.bestVal - max(archive.value) - opt.archiving.tolFunArch));
            if willBeLocal
                subpop.localConvergeCheck.stopAtUsedEval = subpop.usedEval;
                subpop.terminationFlag = 4;
            end
        end
    end
end

% =========================================================================
% NUMERICAL HELPERS
% =========================================================================
function y = rs_normcdf(x)
    y = 0.5 * (1 + erf(x ./ sqrt(2)));
end

function g = rs_geomean(x)
    x = x(:);
    g = exp(mean(log(max(x, realmin))));
end
