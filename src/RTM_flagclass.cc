#include "RTM_flagclass.h"
#include "RTM_StructDefs.h"
#include "string_fns.h"
#include "RTM_StructAssign.h"

using namespace std;
using std::string;

/////// ///////
unsigned int flagclass::instances = 0;
unsigned int flagclass::iSize = 0;
string* flagclass::rmp_member_names = 0;
void flagclass::flagclass_static_new()
{
  if(instances == 0)
    {
      int int_delim_position = 0;
      string all_names = REGIONAL_MODEL_PARAMS_MEMBERS;
      iSize = count_delims_in_string(all_names, ",") + 1;
      rmp_member_names = new string[iSize];

      for(int int_i = 0; int_i < iSize; int_i++)
	{
	  rmp_member_names[int_i].assign(read_from_delim_string<string>(all_names,
									",;",
									int_delim_position));
	}
    }
}
void flagclass::new_and_set_all_flags(const bool& set_val)
{
  if(iSize > 0){
    regional_update_flags = new bool[iSize];
    for(int int_i = 0; int_i < iSize; int_i++)
      regional_update_flags[int_i] = set_val;
  }
}
flagclass::flagclass()
{
  flagclass_static_new();
  new_and_set_all_flags(true);
  ++instances;
}
string flagclass::UPItoRMP(const int& iUPI_index)
{
  string out_string;

  switch(iUPI_index){
  case GP_OVERDISP_INDEX :
    return "l_gp_negbin_overdispersion";
  case HOSP_OVERDISP_INDEX :
    return "l_hosp_negbin_overdispersion";
  case ALP_INDEX : 
    return "l_latent_period : l_R0_init : l_R0_Amplitude : l_I0";
  case AIP_INDEX :
    return "l_average_infectious_period : l_R0_init : l_I0 : l_R0_Amplitude";
  case AR1_INDEX :
    return "l_r1_period";
  case VAC1_DISEASE_INDEX :
    return "l_vacc1_disease";
  case VACN_DISEASE_INDEX :
    return "l_vaccn_disease";
  case VACB_DISEASE_INDEX :
    return "l_vaccb_disease";
  case VAC1_INFECT_INDEX :
    return "l_vacc1_infect";
  case VACN_INFECT_INDEX :
    return "l_vaccn_infect";
  case VACB_INFECT_INDEX :
    return "l_vaccb_infect";
  case REL_INFECT_INDEX :
    return "l_relative_infectiousness_I2_wrt_I1";
  case LBETA_RW_INDEX :
    return "l_lbeta_rw";
  case PROP_SYMP_INDEX :
    return "l_pr_symp";
  case CONTACT_INDEX :
    return "l_MIXMOD";
  case R0_AMP_INDEX :
    return "l_R0_Amplitude";
  case R0_PEAKDAY_INDEX :
    return "l_R0_peakday : d_R0_phase_differences : l_R0_Amplitude";
  case EGR_INDEX :
    return "l_EGR : l_R0_init : l_I0 : l_R0_Amplitude";
  case LPL0_INDEX :
    return "l_I0";
  case PROP_SUS_INDEX :
    return "l_init_prop_sus";
  case PROP_HI_GEQ_32_INDEX :
    return "l_init_prop_sus_HI_geq_32";
  case PROP_GP_INDEX :
    return "l_pr_onset_to_GP : l_I0";
  case PROP_HOSP_INDEX :
    return "l_pr_onset_to_Hosp";
  case PROP_DEATH_INDEX :
    return "l_pr_onset_to_Death";
  case IMPORTATION_INDEX :
    return "l_importation_rate";
  case BGR_INDEX :
    return "l_background_gps_counts";
  case SENS_INDEX :
    return "l_sensitivity";
  case SPEC_INDEX :
    return "l_specificity";
  case SSENS_INDEX :
    return "l_sero_sensitivity";
  case SSPEC_INDEX :
    return "l_sero_specificity";
  case IWAN_INDEX :
    return "l_waning_period";
  case DOW_EFFECTS_INDEX :
    return "l_day_of_week_effect";
  default :
    return REGIONAL_MODEL_PARAMS_MEMBERS;
  }
}
bool flagclass::getFlag(const string& reg_varname)
{
  for(int int_i = 0; int_i < iSize; int_i++)
    if(reg_varname.compare(rmp_member_names[int_i]) == 0)
      return regional_update_flags[int_i];
  return true;
}
unsigned int flagclass::getSize()
{
  return iSize;
}
void flagclass::switchFlag(const string& match_string)
{
  for(int int_i = 0; int_i < iSize; int_i++)
    {
      if(match_string.find(rmp_member_names[int_i]) != string::npos)
	{
	  // we have a match
	  regional_update_flags[int_i] = !regional_update_flags[int_i];
	}
    }
}
flagclass::flagclass(const int& iUPI_index)
{
  flagclass_static_new(); // allocates static memory and sets static variables if not already done so
  new_and_set_all_flags(false);
  ++instances;
  string member_names_to_update = UPItoRMP(iUPI_index);
  switchFlag(member_names_to_update);
}

flagclass::~flagclass()
{
  --instances;
  if(iSize > 0)
    {
      delete [] regional_update_flags;
      if(instances == 0)
	delete [] rmp_member_names;
    }
}
