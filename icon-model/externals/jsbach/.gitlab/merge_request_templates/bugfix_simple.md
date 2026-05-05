### short description
* > information may be structured by bullet points and may be used as commit message
* brief description of the bug
* summary of core changes in science and technical aspects of the code



### detailed description
(add more details to the below sections, not right here)  
(these sections may not be included in the commit message)
#### bug description
* description of the bug (summary of the associated issue, see bug-report template)
* add here the merge request that introduced the bug (if known)

#### solution
* what has changed (science and technical aspects) in the code and where (process, module, routine/function)
* plus scientific justification, if applicable

#### ICON-Land testing, docu of changes in simulation results
* [ ] I tested the latest version of the branch
* [ ] simulation results do NOT change
* [ ] simulation results do change
  * provide a brief description of the changes
* compile with NAG / nvhpc at levante and run:
  * JSBACH:
    * [ ] land_jsbach_R2B4_test with `--enable-quincy` CPU (NAG)
    * [ ] land_jsbach_R2B4_test with `--disable-quincy` CPU (NAG)
  * QUINCY with JSBACH soil-physics enabled:
    * [ ] exp.land_quincy_canopy_R2B4_test (NAG)
    * [ ] exp.land_quincy_canopy_R2B4_test with gpu.nvhpc
    * [ ] exp.land_quincy_standalone_test_2y (NAG)

#### associated merge requests
* ICON-mpim MR: add link (full url) to the MR in the ICON-Land repository
* QS MR: add link (full url) to the MR in the ICON-Land repository
* iq-scripts MR: add link (full url) to the MR in the ICON-Land repository

closes issue [add link (full url) to issue]
