# Steps to go through when expanding the model to account for a new data stream.

- add a flag to indicate the presence/absence of the data stream to the ``global_model_instance_parameters`` object in the ``RTM_StructDefs.h`` file.

- add a name for the variable to use in the mod_input.txt files to specify the value for this new flag. You will add this name into the specification of the global constant ``GLOBAL_MODEL_INSTANCE_MEMBERS`` in the ``RTM_StructAssign.h`` file. This variable is a `` : `` separated list.

- add a default value for this variable to the global constant ``GLOBAL_MODEL_INSTANCE_DEFAULTS``, this is a similarly delimited list and the default variable has to be in the same position within the list as its equivalent name in the above list.

- add a line to the function ``read_global_fixed_parameters`` in the ``RTM_Inputs.cc`` file that will read in the value of this flag. The position of this new line of code has to be located between the equivalent code for the variables defined before and after it in the ``GLOBAL_MODEL_INSTANCE_MEMBERS`` list.

- add a ``rtmData*`` object withing the ``Region`` structure defined in ``RTM_StructDefs.h``.

- copy and edit some allocation code for the new ``rtmData*`` object within the ``Region_alloc`` function in file ``RTM_StructAllocFree.cc``. This will use the initialising functions for the ``rtmData`` class. These are defined in ``rtm_data_class.cc`` and ``rtm_data_class.h``. Look at these functions to see how to most appropriately initiate this class in this new instance. Currently this will depend on whether or not there will be a likelihood attached to these data.

- copy and edit a line to free any memory allocated by the above function using the ``Region_free`` function in the same file.

- copy and edit two lines in the ``Region_memcpy`` function to be able to copy any memory allocated using the ``Region_alloc`` function, within the same file.

- add into the ``mod_inputs.Rmd`` template file a data flag for the new data. Must have the same name that was inserted into the ``GLOBAL MODEL_INSTANCE_MEMBERS`` string. Could be set equal to the output of some R code located in the ``config.R`` file.

- add a code chunk to handle the reading in of data from the ``mod_inputs.txt`` file, populating a variable of ``rtmData`` class. This is done using the ``read_metaregion_datatype`` function. The sixth and seventh arguments will typically give filenames within the same directory that contain data counts and a denominator respectively. The eighth parameter indicates any aggregation of age classes that is required. The last argument indicates whether or not the denominator is to be normalised using the population size.

- adapt the likelihood code accordingly.



TEST THE INITIALISATION OF THE VACCINATION DATACLASS OBJECT AND THE NDAYS MEMBER IN PARTICULAR.

