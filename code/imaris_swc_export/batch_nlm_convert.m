function summary = batch_nlm_convert(inputDir, outputDir, converterPath, options)
%BATCH_NLM_CONVERT Convert all SWC files in a folder with NLMorphologyConverter.
%
% SUMMARY = BATCH_NLM_CONVERT(INPUTDIR, OUTPUTDIR, CONVERTERPATH) converts
% each .swc file in INPUTDIR and writes the result to OUTPUTDIR.
%
% OPTIONAL EXTERNAL DEPENDENCY:
%   NLMorphologyConverter. The local installer name suggests version 0.9.0,
%   but the installed binary and its provenance must be verified before use.
%   Pass its executable path explicitly.

    arguments
        inputDir (1,1) string
        outputDir (1,1) string
        converterPath (1,1) string
        options.OutputFormat (1,1) string = "SWC"
        options.OverwriteExisting (1,1) logical = false
    end

    if ~isfolder(inputDir)
        error("batch_nlm_convert:MissingInputDirectory", ...
            "Input directory does not exist: %s", inputDir);
    end
    if ~isfile(converterPath)
        error("batch_nlm_convert:MissingConverter", ...
            "Converter executable does not exist: %s", converterPath);
    end
    if isempty(regexp(options.OutputFormat, "^[A-Za-z0-9_-]+$", "once"))
        error("batch_nlm_convert:InvalidOutputFormat", ...
            "OutputFormat may contain only letters, numbers, underscores and hyphens.");
    end
    if ~isfolder(outputDir)
        mkdir(outputDir);
    end

    swcFiles = dir(fullfile(inputDir, "*.swc"));
    if isempty(swcFiles)
        warning("batch_nlm_convert:NoFiles", ...
            "No SWC files were found in: %s", inputDir);
        summary = struct("converted", 0, "failed", 0, "skipped", 0);
        return;
    end

    fprintf("NLM conversion started\n");
    fprintf("SWC files: %d\n", numel(swcFiles));
    fprintf("Output directory: %s\n\n", outputDir);

    converted = 0;
    failed = 0;
    skipped = 0;

    for k = 1:numel(swcFiles)
        inPath = fullfile(swcFiles(k).folder, swcFiles(k).name);
        outPath = fullfile(outputDir, swcFiles(k).name);

        if ~options.OverwriteExisting && isfile(outPath)
            skipped = skipped + 1;
            fprintf("[%d/%d] SKIP %s\n", k, numel(swcFiles), swcFiles(k).name);
            continue;
        end

        fprintf("[%d/%d] Converting %s\n", k, numel(swcFiles), swcFiles(k).name);
        command = sprintf('"%s" "%s" "%s" %s', ...
            converterPath, inPath, outPath, options.OutputFormat);
        [status, commandOutput] = system(command);

        if status == 0
            converted = converted + 1;
        else
            failed = failed + 1;
            warning("batch_nlm_convert:ConversionFailed", ...
                "Conversion failed for %s (status %d): %s", ...
                swcFiles(k).name, status, strtrim(commandOutput));
        end
    end

    summary = struct("converted", converted, "failed", failed, "skipped", skipped);
    fprintf("\nNLM conversion complete: converted=%d, failed=%d, skipped=%d\n", ...
        converted, failed, skipped);
end
