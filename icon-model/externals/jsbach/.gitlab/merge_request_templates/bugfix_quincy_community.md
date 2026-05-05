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
  * please generally avoid adding chunks of the code to this MR description
  * instead provide a brief description of your changes along with `module-name:routine-name` (and optionally the line number)
* plus scientific justification, if applicable

#### additional changes
* additional changes worth mentioning that do not directly belong to the issue / solution but are nevertheless connected to the issue in a way that they can be seen as part of the same feature

#### follow-up or remaining issues
* add here a list of new (follow-up) or remaining issues, if applicable
* link every item to an existing gitlab issue or create a new gitlab issue and provide the link to it

#### ICON-Land testing, docu of changes in simulation results
* [ ] I tested the latest version of the branch (using the NAG compiler, i.e., applying debug compile options)
* [ ] simulation results do NOT change
* [ ] simulation results do change
  * provide the links to the directory of the test simulations as a comment to this merge request, i.e., not here in the merge-request description
  * document for each test simulation whether the simulation results are changing
    * provide:
      * a brief description of the changes (and whether these are specific to any test simulation/usecase)
      * justification for how this is improving the model
    * docu any sites that die or not die anymore compared to the reference run
  * [ ] no changes in simulation results unrelated to this MR, i.e., unexpected side effects
    * if you find side effect please add an explanation and docu
* compile with NAG / nvhpc at levante and run:
  * JSBACH:
    * [ ] land_jsbach_R2B4_test with `--enable-quincy` CPU (NAG, nvcpu)
    * [ ] land_jsbach_R2B4_test with `--enable-quincy` GPU
    * [ ] land_jsbach_R2B4_test with `--disable-quincy` CPU (NAG)
    * [ ] atm_nwp_jsbach_test with `--disable-quincy` CPU (NAG) note: needs to set 80 GB of mem
    * [ ] aes_amip_test with `--disable-quincy` CPU (NAG)
  * QUINCY with JSBACH soil-physics enabled:
    * [ ] exp.land_quincy_canopy_R2B4_test (NAG)
    * [ ] exp.land_quincy_canopy_R2B4_test with gpu.nvhpc
    * [ ] exp.land_quincy_canopy_R2B4_test with cpu.nvhpc
      * [ ] identical results for cpu / gpu simulations
    * [ ] exp.land_quincy_standalone_test_1d (NAG)
    * [ ] exp.land_quincy_standalone_test_2m (NAG)
    * [ ] exp.land_quincy_standalone_test_2y (NAG)
    * [ ] results of a quincy land CN simulation with and without restart are identical (NAG)
      * any additional info
  * QUINCY with SPQ_ enabled:
    * [ ] exp.land_quincy_standalone_test_2m_spq (NAG)
* [ ] no water-balance issue
  * provide details about the test simulation

#### associated merge requests
* ICON-mpim MR: add link (full url) to the MR in the ICON-Land repository
* QS MR: add link (full url) to the MR in the ICON-Land repository
  * brief summary of IQ test simulations
* iq-scripts MR: add link (full url) to the MR in the ICON-Land repository

closes issue [add link (full url) to issue]
