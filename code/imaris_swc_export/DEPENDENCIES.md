# Dependencies

## Required

- **MATLAB.** The scripts use argument validation, string arrays,
  `writematrix`, sparse matrices and HDF5 functionality.
  
- **ImarisReader**, Peter Beemiller, commit
  `9b71e5ca3b570c52698e9494ba339eb13c971de2`, MIT license:
  https://github.com/PeterBeemiller/ImarisReader

-  **NLMorphologyConverter.** Required only when `RunNLM=true` or when calling
    `batch_nlm_convert.m`. The executable is not
    redistributed or downloaded by these scripts.
