
# Steps to go through when adding a new parameter to the model.

- add a name for the parameter in the list of parameter names attached to the definition of ``GLOBAL_MODEL_PARAMETERS_MEMBERS`` in `RTM_StructAssign.h`. This will only determine the string that is searched for in the `mod_pars.txt` file for this parameter's info.

- add an index name to the enumerated typedef ``updateable_parameter_index``. The position in this list must correspond to the position of the parameter name in the above list.

- add a default parameter value to the string ``GLOBAL_MODEL_PARAMETERS_DEFAULT_FIXED_VALS``. The position in this list must correspond to the position of the parameter name in the ``GLOBAL_MODEL_PARAMETERS`` list.

- add a string of 0 and 1 flags to the string ``GLOBAL_MODEL_PARAMETERS_LIKELIHOOD_FLAGS``. These are indicators of quantities that will need to be recalculated if a new value for the parameter is proposed. In sequence, the flags correspond to the transmission mode, the reporting model, the GP_likelihood, the Hospitalisation likelihood, the serological likelihood, and the virological likelihood. The position in this list must correspond to the position of the parameter name in the ``GLOBAL_MODEL_PARAMETERS`` list.

- add a member names to the variable ``REGIONAL_MODEL_PARAMS_MEMBERS`` in file `RTM_StructDefs.h`. This will be the names of any new transformed parameters that will need to go into the parameter list sent to each region object.

- add this member to the ``regional_model_params`` structure in file ``RTM_StructDefs.h``.

- allocate some memory to the new component of the ``regional_model_params`` in the ``regional_model_params_alloc`` functions (this is an overloaded function name, so this needs to be done in two places). This isn't necessary for scalar quantities.

- free memory from the new component of the ``regional_model_params`` in the ``regional_model_params_free`` function. This isn't necessary for scalar quantities.

- similarly adapt the memory copy routine ``regional_model_params_memcpy``.

- in the function ``flagclass::UPItoRMP`` in the ``RTM_flagclass.cc`` files, add a case for the new parameter's index name. This should return a " : " delimited list of the regional parameters that will need modifying when you change the new parameter.

- add a line to the nameMap at the top of ``RTM_updParams.cc``, and choose a short constant name (in all caps).

- add the above constant name to the `paramIndex` enum (in the same order as in `nameMap`) further down ``RTM_updParams.cc``.

- if this parameter is a variance or overdispersion parameter that relates directly to data, re-size the relevant structures in the ``regional_model_params`` 

- add a call to the mapping function ``regional_matrix_parameter`` or ``regional_vector_parameter`` or modify exsiting ones as appropriate within the ``evaluate_regional_parameters``. This describes how the updateable model parameter is mapped into the regions.

- the above function is deprecated and is updated for consistency. Noting the change in format for the function arguments, enter similar blocks of code into the `block_regional_parameters` function.

- modify the likelihood functions to accommodate the new model parameter.
