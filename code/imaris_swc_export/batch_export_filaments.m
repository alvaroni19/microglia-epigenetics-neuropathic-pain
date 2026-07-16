function summary = batch_export_filaments(inputDir, swcOutputDir, options)
%BATCH_EXPORT_FILAMENTS Export Imaris filament objects from .ims files to SWC.
%
% SUMMARY = BATCH_EXPORT_FILAMENTS(INPUTDIR, SWCOUTPUTDIR) processes every
% .ims file in INPUTDIR and writes its filament objects as SWC files.
%
% REQUIRED EXTERNAL DEPENDENCY:
%   ImarisReader by Peter Beemiller, commit
%   9b71e5ca3b570c52698e9494ba339eb13c971de2 (MIT license), must be
%   obtained separately and added to the MATLAB path:
%   https://github.com/PeterBeemiller/ImarisReader
%
% This file contains all custom helper functions needed to transform the
% Filaments objects returned by ImarisReader into SWC files.

    arguments
        inputDir (1,1) string
        swcOutputDir (1,1) string
        options.MinNodes (1,1) double {mustBeInteger,mustBeNonnegative} = 1
        options.Recursive (1,1) logical = false
        options.RunNLM (1,1) logical = false
        options.NLMOutputDir (1,1) string = ""
        options.ConverterPath (1,1) string = ""
        options.OutputFormat (1,1) string = "SWC"
        options.OverwriteExistingNLM (1,1) logical = false
    end

    if ~isfolder(inputDir)
        error("batch_export_filaments:MissingInputDirectory", ...
            "Input directory does not exist: %s", inputDir);
    end
    if ~isfolder(swcOutputDir)
        mkdir(swcOutputDir);
    end

    if options.Recursive
        imsFiles = dir(fullfile(inputDir, "**", "*.ims"));
    else
        imsFiles = dir(fullfile(inputDir, "*.ims"));
    end

    if isempty(imsFiles)
        warning("batch_export_filaments:NoFiles", ...
            "No .ims files were found in: %s", inputDir);
        summary = struct("imsFiles", 0, "exported", 0, "omitted", 0, ...
            "failedFiles", 0, "nlm", []);
        return;
    end

    fprintf("Imaris-to-SWC export started\n");
    fprintf("IMS files: %d\n", numel(imsFiles));
    fprintf("Output directory: %s\n", swcOutputDir);
    fprintf("Minimum nodes per filament: %d\n\n", options.MinNodes);

    totalExported = 0;
    totalOmitted = 0;
    failedFiles = 0;

    for k = 1:numel(imsFiles)
        imsPath = string(fullfile(imsFiles(k).folder, imsFiles(k).name));
        fprintf("[%d/%d] Processing %s\n", k, numel(imsFiles), imsFiles(k).name);

        try
            fileObj = ImarisReader(imsPath);
            cleanupObj = onCleanup(@() delete(fileObj));
            [exported, omitted] = export_all_filaments_to_swc( ...
                fileObj, swcOutputDir, options.MinNodes, imsPath);
            totalExported = totalExported + exported;
            totalOmitted = totalOmitted + omitted;
            clear cleanupObj fileObj;
        catch exception
            failedFiles = failedFiles + 1;
            warning("batch_export_filaments:FileFailed", ...
                "Failed to process %s: %s", imsFiles(k).name, exception.message);
        end
    end

    nlmSummary = [];
    if options.RunNLM
        if strlength(options.ConverterPath) == 0
            error("batch_export_filaments:MissingConverter", ...
                "ConverterPath is required when RunNLM is true.");
        end
        nlmOutputDir = options.NLMOutputDir;
        if strlength(nlmOutputDir) == 0
            nlmOutputDir = swcOutputDir + "_nlm";
        end
        nlmSummary = batch_nlm_convert(swcOutputDir, nlmOutputDir, ...
            options.ConverterPath, OutputFormat=options.OutputFormat, ...
            OverwriteExisting=options.OverwriteExistingNLM);
    end

    summary = struct("imsFiles", numel(imsFiles), "exported", totalExported, ...
        "omitted", totalOmitted, "failedFiles", failedFiles, "nlm", nlmSummary);
    fprintf("\nExport complete: exported=%d, omitted=%d, failed IMS files=%d\n", ...
        totalExported, totalOmitted, failedFiles);
end

function [totalExported, totalOmitted] = export_all_filaments_to_swc( ...
        fileObj, outDir, minNodes, imsFileName)
%EXPORT_ALL_FILAMENTS_TO_SWC Export every Imaris Filaments object to SWC.

    filamentsList = fileObj.Filaments;
    if isempty(filamentsList)
        error("batch_export_filaments:NoFilaments", ...
            "The selected Imaris file does not contain Filaments objects.");
    end

    totalExported = 0;
    totalOmitted = 0;
    for filamentsIndex = 1:numel(filamentsList)
        try
            [exported, omitted] = export_filaments_object_to_swc( ...
                fileObj, filamentsIndex, outDir, minNodes, imsFileName);
            totalExported = totalExported + exported;
            totalOmitted = totalOmitted + omitted;
        catch exception
            warning("batch_export_filaments:FilamentsObjectFailed", ...
                "Failed to process Filaments object %d: %s", ...
                filamentsIndex, exception.message);
        end
    end
end

function [exportedCount, omittedCount] = export_filaments_object_to_swc( ...
        fileObj, filamentsIndex, outDir, minNodes, imsFileName)
%EXPORT_FILAMENTS_OBJECT_TO_SWC Export one Imaris Filaments object.

    filamentsList = fileObj.Filaments;
    if filamentsIndex > numel(filamentsList)
        error("batch_export_filaments:IndexOutOfRange", ...
            "Filaments index %d exceeds the number of available objects (%d).", ...
            filamentsIndex, numel(filamentsList));
    end

    filamentObject = filamentsList(filamentsIndex);
    imsBase = sanitize_name(resolve_ims_base_name(imsFileName));
    objectName = sanitize_name(string(filamentObject.Name));
    if strlength(objectName) == 0
        objectName = "Filaments_" + string(filamentsIndex);
    end

    numberOfFilaments = filamentObject.NumberOfFilaments;
    filamentIDs = [];
    if numberOfFilaments > 1
        try
            filamentIDs = filamentObject.GetIDs();
        catch
            filamentIDs = [];
        end
    end

    exportedCount = 0;
    omittedCount = 0;
    for filamentIndex = 0:numberOfFilaments-1
        positions = filamentObject.GetPositions(filamentIndex);
        numberOfNodes = size(positions, 1);
        if numberOfNodes < minNodes
            omittedCount = omittedCount + 1;
            continue;
        end

        edges = filamentObject.GetEdges(filamentIndex) + 1;
        radii = filamentObject.GetRadii(filamentIndex);
        nodeTypes = filamentObject.GetTypes(filamentIndex);
        if isempty(radii) || numel(radii) ~= numberOfNodes
            radii = ones(numberOfNodes, 1);
        else
            radii = radii(:);
        end
        if isempty(nodeTypes) || numel(nodeTypes) ~= numberOfNodes
            nodeTypes = zeros(numberOfNodes, 1);
        else
            nodeTypes = nodeTypes(:);
        end

        root = filamentObject.GetBeginningVertexIndex(filamentIndex) + 1;
        root = max(1, min(numberOfNodes, root));
        parents = build_tree_sparse(edges, root, numberOfNodes);

        swc = zeros(numberOfNodes, 7);
        swc(:, 1) = (1:numberOfNodes).';
        swc(:, 2) = 3;
        swc(nodeTypes == 1, 2) = 4;
        swc(:, 3:5) = positions;
        swc(:, 6) = radii;
        swc(:, 7) = parents;
        swc(root, 2) = 1;
        swc(root, 7) = -1;
        swc(swc(:, 7) == 0, 7) = root;

        label = build_filament_label( ...
            objectName, filamentIndex, numberOfFilaments, filamentIDs);
        outputStem = sanitize_name(imsBase + "_" + label);
        outputPath = build_unique_output_path(outDir, outputStem, ".swc");
        writematrix(swc, char(outputPath), ...
            "Delimiter", " ", "FileType", "text");
        exportedCount = exportedCount + 1;
    end
end

function imsBase = resolve_ims_base_name(imsFileName)
    [~, baseName, ~] = fileparts(char(imsFileName));
    imsBase = string(baseName);
    if strlength(imsBase) == 0
        imsBase = "ImarisFile";
    end
end

function label = build_filament_label(baseName, index, count, filamentIDs)
    if count == 1
        label = baseName;
    elseif ~isempty(filamentIDs) && numel(filamentIDs) >= index + 1
        label = baseName + "_ID" + string(filamentIDs(index + 1));
    else
        label = baseName + "_f" + string(index);
    end
end

function outputPath = build_unique_output_path(outDir, fileStem, extension)
    outputPath = string(fullfile(char(outDir), char(fileStem + extension)));
    suffix = 1;
    while isfile(outputPath)
        candidate = fileStem + "_" + sprintf("%02d", suffix) + extension;
        outputPath = string(fullfile(char(outDir), char(candidate)));
        suffix = suffix + 1;
    end
end

function safeName = sanitize_name(rawName)
    safeName = strtrim(string(rawName));
    safeName = regexprep(safeName, "\s+", "_");
    safeName = regexprep(safeName, "[^a-zA-Z0-9_-]", "_");
    safeName = regexprep(safeName, "_+", "_");
    safeName = regexprep(safeName, "^_+|_+$", "");
end

function parent = build_tree_sparse(edges, root, numberOfNodes)
    if isempty(edges)
        parent = -ones(numberOfNodes, 1);
        parent(root) = -1;
        return;
    end

    edges = double(edges);
    edges(edges < 1) = 1;
    edges(edges > numberOfNodes) = numberOfNodes;
    graphMatrix = sparse( ...
        [edges(:, 1); edges(:, 2)], ...
        [edges(:, 2); edges(:, 1)], ...
        1, numberOfNodes, numberOfNodes);

    parent = -ones(numberOfNodes, 1);
    stack = root;
    while ~isempty(stack)
        current = stack(end);
        stack(end) = [];
        neighbours = find(graphMatrix(current, :));
        for neighbour = neighbours
            if parent(neighbour) == -1 && neighbour ~= root
                parent(neighbour) = current;
                stack(end + 1) = neighbour; %#ok<AGROW>
            end
        end
    end

    orphanMask = parent == -1 & (1:numberOfNodes).' ~= root;
    parent(orphanMask) = root;
    parent(root) = -1;
end
