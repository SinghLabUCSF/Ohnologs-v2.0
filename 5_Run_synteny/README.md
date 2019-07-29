## Run macro-synteny comparisons

* [1_SyntenyScripts](1_SyntenyScripts): Scripts to run synteny for each outgroup and vertebrate genome, or with self for a defined window. 

* [2_ParameterFiles](2_ParameterFiles): Directory having parameter files for each organism -- for self and outgroup synteny -- for 2R and 3R (for fish).
    
* [3_SyntenyOutputFiles](3_SyntenyOutputFiles): This will have synteny output file for each organism, WGD and window.  
  
> **Run synteny scripts:**  When everything is set up. To run, cd into the [1_SyntenyScripts](1_SyntenyScripts) and run the OHNOLOG scripts. Two master scripts that needs to be run are:
  **`OHNOLOGS_SelfComparison_CL.pl`**   
  **`OHNOLOGS_OutgroupComparison_CL.pl`**   
  
> No warning or error should be printed.  

* [4_CombineWindows](4_CombineWindows) : Combine pairs from for different windows for self and outgroup comparisons. Also combine all pairs from different outgroups.

* [5_CombineSelfOutgroup](5_CombineSelfOutgroup): Combine pairs from self and outgroup comparisons. This is the final output for a single organism before input from other vertebrates.

