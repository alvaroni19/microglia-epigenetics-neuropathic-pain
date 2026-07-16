# Dependencies

## Required

- **MATLAB.** The scripts use argument validation, string arrays,
  `writematrix`, sparse matrices and HDF5 functionality. The exact MATLAB
  release used for the archived analysis was not recorded in the available
  provenance files.
- **ImarisReader**, Peter Beemiller, commit
  `9b71e5ca3b570c52698e9494ba339eb13c971de2`, MIT license:
  https://github.com/PeterBeemiller/ImarisReader

  Obtain the complete ImarisReader repository from its official source and
  add it to the MATLAB path. It is not selected for redistribution in this
  archive.
-  **NLMorphologyConverter.** Required only when `RunNLM=true` or when calling
    `batch_nlm_convert.m`. The local installer name suggests version 0.9.0;
    verify the installed version, authorship, official distribution source and
    installer integrity before execution. The executable is not
    redistributed or downloaded by these scripts.
