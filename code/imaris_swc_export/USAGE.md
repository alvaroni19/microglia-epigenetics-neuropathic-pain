# Imaris filament to SWC export

Add this directory to the MATLAB path. The required input is a folder
containing Imaris `.ims` files; no study-specific folder structure or animal
identifier convention is assumed.

Export every `.ims` file in one folder:

```matlab
summary = batch_export_filaments( ...
    "D:\project\imaris", ...
    "D:\project\swc", ...
    MinNodes=1);
```

Search recursively:

```matlab
summary = batch_export_filaments( ...
    "D:\project\imaris", ...
    "D:\project\swc", ...
    MinNodes=1, ...
    Recursive=true);
```

Optional conversion with a separately installed copy of
NLMorphologyConverter:

```matlab
summary = batch_export_filaments( ...
    "D:\project\imaris", ...
    "D:\project\swc", ...
    RunNLM=true, ...
    NLMOutputDir="D:\project\swc_nlm", ...
    ConverterPath="C:\Program Files (x86)\Neuronland\NLMorphologyConverter\NLMorphologyConverter.exe");
```

The function returns counts of processed `.ims` files, exported filaments,
omitted filaments and failed input files. NLMorphologyConverter is not
downloaded or redistributed by these scripts.

