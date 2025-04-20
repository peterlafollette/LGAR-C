#ifndef BMI_LGAR_CXX_INCLUDED
#define BMI_LGAR_CXX_INCLUDED


#include <stdio.h>
#include <string>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <iostream>
#include "../bmi/bmi.hxx"
#include "../include/bmi_lgar.hxx"
#include "../include/all.hxx"


// default verbosity is set to 'none' other option 'high' or 'low' needs to be specified in the config file
string verbosity="none";

/**
 * @brief Delete dynamic arrays allocated in Initialize() and held by this object
 * 
 */
BmiLGAR::~BmiLGAR(){
  delete [] giuh_ordinates;
  delete [] giuh_runoff_queue;
}

/* The `head` pointer stores the address in memory of the first member of the linked list containing
   all the wetting fronts. The contents of struct wetting_front are defined in "all.h" */

void BmiLGAR::
Initialize (std::string config_file)
{
  if (config_file.compare("") != 0 ) {
    this->state = new model_state;
    state->head = NULL;
    state->state_previous = NULL;
    state->head_frac = NULL;
    state->state_previous_frac = NULL;
    lgar_initialize(config_file, state);
  }

  num_giuh_ordinates = state->lgar_bmi_params.num_giuh_ordinates;

  /* giuh ordinates are static and read in the lgar.cxx, and we need to have a copy of it to pass to
     giuh.cxx, so allocating/copying here*/

  giuh_ordinates = new double[num_giuh_ordinates];
  giuh_runoff_queue = new double[num_giuh_ordinates+1];

  for (int i=0; i<num_giuh_ordinates;i++){
    giuh_ordinates[i] = state->lgar_bmi_params.giuh_ordinates[i+1]; // note lgar uses 1-indexing
  }

  for (int i=0; i<=num_giuh_ordinates;i++){
    giuh_runoff_queue[i] = 0.0;
  }

}

/**
 * @brief Allocate (or reallocate) storage for soil parameters
 * 
 */
void BmiLGAR::realloc_soil(){

  delete [] state->lgar_bmi_params.soil_depth_wetting_fronts;
  delete [] state->lgar_bmi_params.soil_moisture_wetting_fronts;

  state->lgar_bmi_params.soil_depth_wetting_fronts = new double[state->lgar_bmi_params.num_wetting_fronts];
  state->lgar_bmi_params.soil_moisture_wetting_fronts = new double[state->lgar_bmi_params.num_wetting_fronts];
}

/*
  This is the main function calling lgar subroutines for creating, moving, and merging wetting fronts.
  Calls to AET and mass balance module are also happening here
  If the model's timestep is smaller than the forcing's timestep then we take subtimesteps inside the subcycling loop
*/
void BmiLGAR::
Update()
{
  if (verbosity.compare("none") != 0) {
    std::cerr<<"---------------------------------------------------------\n";
    std::cerr<<"|****************** LASAM BMI Update... ******************|\n";
    std::cerr<<"---------------------------------------------------------\n";
  }


  double mm_to_cm = 0.1; // unit conversion
  double mm_to_m = 0.001;

  if (state->lgar_bmi_params.is_invalid_soil_type) {
    // add to mass balance accumulated variables
    state->lgar_mass_balance.volprecip_cm  += state->lgar_bmi_input_params->precipitation_mm_per_h * mm_to_cm;
    state->lgar_mass_balance.volin_cm       = 0.0;
    state->lgar_mass_balance.volin_cm_frac  = 0.0;
    state->lgar_mass_balance.volon_cm       = 0.0;
    state->lgar_mass_balance.volend_cm      = state->lgar_mass_balance.volstart_cm;
    state->lgar_mass_balance.volend_cm_frac = state->lgar_mass_balance.volstart_cm_frac;
    state->lgar_mass_balance.volAET_cm      = 0.0;
    state->lgar_mass_balance.volrech_cm     = 0.0;
    state->lgar_mass_balance.volrech_cm_frac= 0.0;
    state->lgar_mass_balance.volrunoff_cm  += state->lgar_bmi_input_params->precipitation_mm_per_h * mm_to_cm;
    state->lgar_mass_balance.volrunoff_cm_frac = 0.0;
    state->lgar_mass_balance.volQ_cm       += state->lgar_bmi_input_params->precipitation_mm_per_h * mm_to_cm;
    state->lgar_mass_balance.volQ_gw_cm     = 0.0;
    state->lgar_mass_balance.volPET_cm      = 0.0;
    state->lgar_mass_balance.volrunoff_giuh_cm  = 0.0;
    state->lgar_mass_balance.volchange_calib_cm = 0.0;
    state->lgar_mass_balance.mass_transfer_cm = 0.0;

    // converted values, a struct local to the BMI and has bmi output variables
    bmi_unit_conv.mass_balance_m        = 0.0;
    bmi_unit_conv.volprecip_timestep_m  = state->lgar_bmi_input_params->precipitation_mm_per_h * mm_to_m;
    bmi_unit_conv.volin_timestep_m      = 0.0;
    bmi_unit_conv.volin_timestep_m_frac = 0.0;
    bmi_unit_conv.volend_timestep_m     = 0.0;
    bmi_unit_conv.volend_timestep_m_frac= 0.0;
    bmi_unit_conv.volAET_timestep_m     = 0.0;
    bmi_unit_conv.volrech_timestep_m    = 0.0;
    bmi_unit_conv.volrech_timestep_m_frac= 0.0;
    bmi_unit_conv.volrunoff_timestep_m  = state->lgar_bmi_input_params->precipitation_mm_per_h * mm_to_m;
    bmi_unit_conv.volrunoff_timestep_m_frac  = 0.0;
    bmi_unit_conv.volQ_timestep_m       = state->lgar_bmi_input_params->precipitation_mm_per_h * mm_to_m;
    bmi_unit_conv.volQ_gw_timestep_m    = 0.0;
    bmi_unit_conv.volPET_timestep_m     = 0.0;
    bmi_unit_conv.volrunoff_giuh_timestep_m = 0.0;
    bmi_unit_conv.mass_transfer_timestep_m = 0.0;

    return;
  }

  // if lasam is coupled to soil freeze-thaw, frozen fraction module is called
  if (state->lgar_bmi_params.sft_coupled)
    frozen_factor_hydraulic_conductivity(state->lgar_bmi_params);

  double volchange_calib_cm = 0.0;

  if(state->lgar_bmi_params.calib_params_flag) {
    volchange_calib_cm = update_calibratable_parameters(state->lgar_bmi_params.dual_perm, state->lgar_bmi_params.ratio_fracture_vol_to_total_vol); // change in soil water volume due to calibratable parameters
    state->lgar_bmi_params.calib_params_flag = false;
  }

  // local variables for readibility
  int subcycles;
  int num_layers = state->lgar_bmi_params.num_layers;

  // local variables for a full timestep (i.e., timestep of the forcing data)
  // see 'struct lgar_mass_balance_variables' in all.hxx for full description of the variables
  double precip_timestep_cm = 0.0;
  double PET_timestep_cm    = 0.0;
  double AET_timestep_cm    = 0.0;
  double volend_timestep_cm = (1-state->lgar_bmi_params.ratio_fracture_vol_to_total_vol)*lgar_calc_mass_bal(state->lgar_bmi_params.cum_layer_thickness_cm, state->head); // this should not be reset to 0.0 in the for loop
  double volend_timestep_cm_frac = state->lgar_bmi_params.ratio_fracture_vol_to_total_vol*lgar_calc_mass_bal(state->lgar_bmi_params.cum_layer_thickness_cm, state->head_frac);
  double volin_timestep_cm  = 0.0;
  double volin_timestep_cm_frac  = 0.0;
  double volon_timestep_cm  = state->lgar_mass_balance.volon_timestep_cm;
  // double volon_timestep_cm_frac  = 0.0; volon is separated at the subtimestep level but aggregated for the time step level
  double volrunoff_timestep_cm      = 0.0;
  double volrunoff_timestep_cm_frac      = 0.0;
  double volrech_timestep_cm        = 0.0;
  double volrech_timestep_cm_frac        = 0.0;
  double surface_runoff_timestep_cm = 0.0; // direct surface runoff
  double volrunoff_giuh_timestep_cm = 0.0;
  double volQ_timestep_cm           = 0.0;
  double volQ_gw_timestep_cm        = 0.0;
  double QF_Q_timestep_cm           = 0.0;
  double mass_transfer_timestep_cm  = 0.0;
  
  // // local variables for a subtimestep (i.e., timestep of the model)
  // double precip_subtimestep_cm;
  // double precip_subtimestep_cm_per_h;
  // double PET_subtimestep_cm;
  // double PET_subtimestep_cm_per_h;
  // double ponded_depth_subtimestep_cm;
  // double ponded_depth_subtimestep_cm_frac;
  // double ponded_depth_subtimestep_cm_both_domains;
  // double AET_subtimestep_cm;
  // double volstart_subtimestep_cm;
  // double volstart_subtimestep_cm_frac;
  // double volend_subtimestep_cm = volend_timestep_cm; // this should not be reset to 0.0 in the for loop
  // double volend_subtimestep_cm_frac = volend_timestep_cm_frac; // this should not be reset to 0.0 in the for loop
  // double volin_subtimestep_cm;
  // double volin_subtimestep_cm_frac = 0.0;
  // double volon_subtimestep_cm;
  // double volon_subtimestep_cm_frac = 0.0;
  // double volrunoff_subtimestep_cm;
  // double volrunoff_subtimestep_cm_frac = 0.0;
  // double volrech_subtimestep_cm;
  // double volrech_subtimestep_cm_frac = 0.0;
  // double surface_runoff_subtimestep_cm; // direct surface runoff
  // double precip_previous_subtimestep_cm;
  // double volQ_gw_subtimestep_cm = 0.0; // fix it for non-zero values after adding groundwater reservoir

  double precip_subtimestep_cm = 0.0;
  double precip_subtimestep_cm_per_h = 0.0;
  double PET_subtimestep_cm = 0.0;
  double PET_subtimestep_cm_per_h = 0.0;
  double ponded_depth_subtimestep_cm = 0.0;
  double ponded_depth_subtimestep_cm_frac = 0.0;
  double ponded_depth_subtimestep_cm_both_domains = 0.0;
  double AET_subtimestep_cm = 0.0;
  double volstart_subtimestep_cm = 0.0;
  double volstart_subtimestep_cm_frac = 0.0;
  double volend_subtimestep_cm = volend_timestep_cm; // this should not be reset to 0.0 in the for loop
  double volend_subtimestep_cm_frac = volend_timestep_cm_frac; // this should not be reset to 0.0 in the for loop
  double volin_subtimestep_cm = 0.0;
  double volin_subtimestep_cm_frac = 0.0;
  double volon_subtimestep_cm = 0.0;
  double volon_subtimestep_cm_frac = 0.0;
  double volrunoff_subtimestep_cm = 0.0;
  double volrunoff_subtimestep_cm_frac = 0.0;
  double volrech_subtimestep_cm = 0.0;
  double volrech_subtimestep_cm_frac = 0.0;
  double surface_runoff_subtimestep_cm = 0.0; // direct surface runoff
  double precip_previous_subtimestep_cm = 0.0;
  double volQ_gw_subtimestep_cm = 0.0; // fix it for non-zero values after adding groundwater reservoir
  double mass_transfer_subtimestep_cm = 0.0;
  
  double subtimestep_h = state->lgar_bmi_params.timestep_h;
  int nint = state->lgar_bmi_params.nint;
  double wilting_point_psi_cm = state->lgar_bmi_params.wilting_point_psi_cm;
  double field_capacity_psi_cm = state->lgar_bmi_params.field_capacity_psi_cm;
  double a = state->lgar_bmi_params.a;
  double b = state->lgar_bmi_params.b;
  double frac_to_GW = state->lgar_bmi_params.frac_to_GW;
  double frac_to_GW_adjusted = 0.0;
  double spf_factor = state->lgar_bmi_params.spf_factor;
  double frac_to_pref = state->lgar_bmi_params.frac_to_pref;
  double ratio_fracture_vol_to_total_vol = state->lgar_bmi_params.ratio_fracture_vol_to_total_vol;
  double root_zone_depth_cm = state->lgar_bmi_params.root_zone_depth_cm;
  bool use_closed_form_G = state->lgar_bmi_params.use_closed_form_G; 
  bool adaptive_timestep = state->lgar_bmi_params.adaptive_timestep;
  bool PET_affects_precip = state->lgar_bmi_params.PET_affects_precip;
  double mbal_tol = state->lgar_bmi_params.mbal_tol;

  // constant value used in the AET function
  double AET_thresh_Theta = 0.85;    // scaled soil moisture (0-1) above which AET=PET (fix later!)
  double AET_expon        = 1.0;     // exponent that allows curvature of the rising portion of the Budyko curve (fix later!)

  double ponded_depth_max_cm = state->lgar_bmi_params.ponded_depth_max_cm;

  if (verbosity.compare("high") == 0) {
    std::cerr<<"Pr  [cm/h] (timestep) = "<<state->lgar_bmi_input_params->precipitation_mm_per_h * mm_to_cm <<"\n";
    std::cerr<<"PET [cm/h] (timestep) = "<<state->lgar_bmi_input_params->PET_mm_per_h * mm_to_cm <<"\n"; 
  }

  assert (state->lgar_bmi_input_params->precipitation_mm_per_h >= 0.0);
  assert(state->lgar_bmi_input_params->PET_mm_per_h >=0.0);

  // adaptive time step is set 
  if (adaptive_timestep) {
    subtimestep_h = state->lgar_bmi_params.forcing_resolution_h;
    if (state->lgar_bmi_input_params->precipitation_mm_per_h > 10.0) {
      subtimestep_h = state->lgar_bmi_params.minimum_timestep_h;  //case where precip > 1 cm/h
    }
    else if (state->lgar_bmi_input_params->precipitation_mm_per_h > 0.0) {
      subtimestep_h = state->lgar_bmi_params.minimum_timestep_h * 2.0;  //case where precip is less than 1 cm/h but greater than 0
    }
    subtimestep_h = fmin(subtimestep_h, state->lgar_bmi_params.forcing_resolution_h);  //just in case the user has specified a minimum time step that would make the subtimestep_h greater than the forcing resolution 
    state->lgar_bmi_params.timestep_h = subtimestep_h;
  }

  state->lgar_bmi_params.forcing_interval = int(state->lgar_bmi_params.forcing_resolution_h/state->lgar_bmi_params.timestep_h+1.0e-08); // add 1.0e-08 to prevent truncation error
  subcycles = state->lgar_bmi_params.forcing_interval;

  if (verbosity.compare("high") == 0) {
    printf("time step size in hours: %lf \n", state->lgar_bmi_params.timestep_h);
  }

  // ensure precip and PET are non-negative
  state->lgar_bmi_input_params->precipitation_mm_per_h = fmax(state->lgar_bmi_input_params->precipitation_mm_per_h, 0.0);
  state->lgar_bmi_input_params->PET_mm_per_h           = fmax(state->lgar_bmi_input_params->PET_mm_per_h, 0.0);

  if (PET_affects_precip){ // if the user wants PET subtracted from precip //has not beed updated with ratio_fracture_vol_to_total_vol yet
    if (state->lgar_bmi_input_params->precipitation_mm_per_h > state->lgar_bmi_input_params->PET_mm_per_h){
      state->lgar_bmi_input_params->precipitation_mm_per_h = state->lgar_bmi_input_params->precipitation_mm_per_h - state->lgar_bmi_input_params->PET_mm_per_h;
      state->lgar_bmi_input_params->PET_mm_per_h = 0.0;
    }
    else{
      state->lgar_bmi_input_params->PET_mm_per_h = state->lgar_bmi_input_params->PET_mm_per_h - state->lgar_bmi_input_params->precipitation_mm_per_h;
      state->lgar_bmi_input_params->precipitation_mm_per_h = 0.0;
    }
  }

  if ( (verbosity.compare("high") == 0) && (PET_affects_precip)) {
    std::cerr<<"Pr  [cm/h] (timestep), after PET is subtracted from precip = "<<state->lgar_bmi_input_params->precipitation_mm_per_h * mm_to_cm <<"\n";
    std::cerr<<"PET [cm/h] (timestep), after PET is subtracted from precip = "<<state->lgar_bmi_input_params->PET_mm_per_h * mm_to_cm <<"\n"; 
  }
  
  // subcycling loop (loop over model's timestep)
  for (int cycle=1; cycle <= subcycles; cycle++) {

    this->state->lgar_bmi_params.time_s    += subtimestep_h * state->units.hr_to_sec;
    this->state->lgar_bmi_params.timesteps ++;
    
    if (verbosity.compare("high") == 0 || verbosity.compare("low") == 0) {
      std::cerr<<"BMI Update |---------------------------------------------------------------|\n";
      std::cerr<<"BMI Update |Timesteps = "<< state->lgar_bmi_params.timesteps<<", Time [h] = "<<this->state->lgar_bmi_params.time_s / 3600.<<", Subcycle = "<< cycle <<" of "<<subcycles<<std::endl;
    }

    if( state->state_previous != NULL ){
      listDelete(state->state_previous);
      state->state_previous = NULL;
    }
    state->state_previous = listCopy(state->head);

    if( state->state_previous_frac != NULL ){
      listDelete(state->state_previous_frac);
      state->state_previous_frac = NULL;
    }
    state->state_previous_frac = listCopy(state->head_frac);

    // // ensure precip and PET are non-negative
    // state->lgar_bmi_input_params->precipitation_mm_per_h = fmax(state->lgar_bmi_input_params->precipitation_mm_per_h, 0.0);
    // state->lgar_bmi_input_params->PET_mm_per_h           = fmax(state->lgar_bmi_input_params->PET_mm_per_h, 0.0);

    /* Note unit conversion:
       Pr and PET are rates (fluxes) in mm/h
       Pr [mm/h] * 1h/3600sec = Pr [mm/3600sec]
       Model timestep (dt) = 300 sec (5 minutes for example)
       convert rate to amount
       Pr [mm/3600sec] * dt [300 sec] = Pr[mm] * 300/3600.
       in the code below, subtimestep_h is this 300/3600 factor (see initialize from config in lgar.cxx)
    */

    precip_subtimestep_cm_per_h = state->lgar_bmi_input_params->precipitation_mm_per_h * mm_to_cm; // rate [cm/hour]
    double precip_for_QF_subtimestep_cm_per_h = 0.0;

    if (state->lgar_bmi_params.runoff_in_prev_step){ // || state->lgar_mass_balance.QF_storage_cm > 0.001
      double precip_subtimestep_cm_per_h_total = precip_subtimestep_cm_per_h;
      if (frac_to_GW_adjusted<1.E-5){
        frac_to_GW_adjusted = 0.0;
      }
      frac_to_GW_adjusted = frac_to_GW_adjusted * frac_to_GW;
      precip_for_QF_subtimestep_cm_per_h = frac_to_GW_adjusted * precip_subtimestep_cm_per_h_total;
      precip_subtimestep_cm_per_h = (1.0 - frac_to_GW_adjusted) * precip_subtimestep_cm_per_h_total;
    }

    PET_subtimestep_cm_per_h = state->lgar_bmi_input_params->PET_mm_per_h * mm_to_cm;

    ponded_depth_subtimestep_cm = precip_subtimestep_cm_per_h * subtimestep_h; // the amount of water on the surface before any infiltration and runoff

    ponded_depth_subtimestep_cm += volon_timestep_cm; // add volume of water on the surface (from the last timestep) to ponded depth as well


    ponded_depth_subtimestep_cm_both_domains = ponded_depth_subtimestep_cm;
    ponded_depth_subtimestep_cm_frac = 0.0;
    if (frac_to_pref){ //if dual permeability is on, then separate ponded depth into the two model domains
      ponded_depth_subtimestep_cm = (1 - frac_to_pref) * ponded_depth_subtimestep_cm_both_domains;
      ponded_depth_subtimestep_cm_frac = frac_to_pref * ponded_depth_subtimestep_cm_both_domains;
      if (verbosity.compare("high") == 0){
        printf("frac_to_pref: %lf \n", frac_to_pref);
        printf("ponded_depth_subtimestep_cm: %lf \n", ponded_depth_subtimestep_cm);
        printf("ponded_depth_subtimestep_cm_frac: %lf \n", ponded_depth_subtimestep_cm_frac);
      }
    }

    if (ratio_fracture_vol_to_total_vol){
      ponded_depth_subtimestep_cm = ponded_depth_subtimestep_cm / (1-ratio_fracture_vol_to_total_vol);
      ponded_depth_subtimestep_cm_frac = ponded_depth_subtimestep_cm_frac / (ratio_fracture_vol_to_total_vol);
      if (verbosity.compare("high") == 0){
        printf("ratio_fracture_vol_to_total_vol: %lf \n", ratio_fracture_vol_to_total_vol);
        printf("ponded_depth_subtimestep_cm: %lf \n", ponded_depth_subtimestep_cm);
        printf("ponded_depth_subtimestep_cm_frac: %lf \n", ponded_depth_subtimestep_cm_frac);
      }
    }







    precip_subtimestep_cm = precip_subtimestep_cm_per_h * subtimestep_h; // rate x dt = amount (portion of the water on the suface for model's timestep [cm])
    PET_subtimestep_cm = PET_subtimestep_cm_per_h * subtimestep_h;      // potential ET for this subtimestep [cm]

    //using cerr instead of cout due to some cout buffering issues when running in the ngen framework, cerr doesn't buffer so it prints immediately to the sreeen.
    if (verbosity.compare("high") == 0 || verbosity.compare("low") == 0) {

      std::cerr<<"Pr [cm/h], Pr [cm] (subtimestep), subtimestep [h] = "<<state->lgar_bmi_input_params->precipitation_mm_per_h * mm_to_cm <<", "<< precip_subtimestep_cm <<", "<< subtimestep_h<<" ("<<subtimestep_h*3600<<" sec)"<<"\n";
      std::cerr<<"PET [cm/h], PET [cm] (subtimestep) = "<<state->lgar_bmi_input_params->PET_mm_per_h * mm_to_cm <<", "<< PET_subtimestep_cm<<"\n";
    }

    AET_subtimestep_cm            = 0.0;
    volstart_subtimestep_cm       = 0.0;
    volstart_subtimestep_cm_frac  = 0.0;
    volin_subtimestep_cm          = 0.0;
    volrunoff_subtimestep_cm      = 0.0;
    volrech_subtimestep_cm        = 0.0;
    surface_runoff_subtimestep_cm = 0.0;
    double free_drainage_subtimestep_cm = 0.0;
    double free_drainage_subtimestep_cm_frac = 0.0;
    if (!(state->lgar_bmi_params.TO_enabled) && (state->lgar_bmi_params.free_drainage_enabled)){
      free_drainage_subtimestep_cm = subtimestep_h*listFindFront(listLength(state->head), state->head, NULL)->K_cm_per_h;
      if (free_drainage_subtimestep_cm<1.E-7){
        free_drainage_subtimestep_cm = 0.0;
      }
      if (listFindFront(listLength(state->head), state->head, NULL)->psi_cm>1.E6){
        free_drainage_subtimestep_cm = 0.0;
      }

      if (state->lgar_bmi_params.dual_perm){
        free_drainage_subtimestep_cm_frac = subtimestep_h*listFindFront(listLength(state->head_frac), state->head_frac, NULL)->K_cm_per_h;
        if (free_drainage_subtimestep_cm_frac<1.E-7){
          free_drainage_subtimestep_cm_frac = 0.0;
        }
        if (listFindFront(listLength(state->head_frac), state->head_frac, NULL)->psi_cm>1.E6){
          free_drainage_subtimestep_cm_frac = 0.0;
        }
      }

    }
    if ((state->lgar_bmi_params.free_drainage_enabled) && verbosity.compare("high") == 0){
      printf("free_drainage_subtimestep_cm: %.10lf \n", free_drainage_subtimestep_cm);
      if (state->lgar_bmi_params.dual_perm){
        printf("free_drainage_subtimestep_cm_frac: %.10lf \n", free_drainage_subtimestep_cm_frac);
      }
    }
    
    precip_previous_subtimestep_cm = state->lgar_bmi_params.precip_previous_timestep_cm; // creation of a new wetting front depends on previous timestep's rainfall

    num_layers = state->lgar_bmi_params.num_layers;
    double delta_theta;   // the width of a front, such that its volume=depth*delta_theta
    double dry_depth;
    double surf_frac_rz = 0.0;
    double *surf_AET_vec = (double *) malloc(sizeof(double)*(listLength(state->head)+1)); //does it really have to be the whole list length? Can't it just be up until the end of surf WFs?
    double *surf_AET_vec_frac = (double *) malloc(sizeof(double)*(listLength(state->head_frac)+1));

    // Calculate AET from PET if PET is non-zero
    if (PET_subtimestep_cm_per_h > 0.0) {
      calc_aet(state->lgar_bmi_params.TO_enabled, PET_subtimestep_cm_per_h/(1 - ratio_fracture_vol_to_total_vol), subtimestep_h, wilting_point_psi_cm, field_capacity_psi_cm, root_zone_depth_cm, &surf_frac_rz, 
               state->lgar_bmi_params.layer_soil_type, AET_thresh_Theta, AET_expon, 
               state->head, state->soil_properties, surf_AET_vec);
      if (!(state->lgar_bmi_params.TO_enabled)){
        int wetting_front_free_drainage_temp = wetting_front_free_drainage(state->head);
        AET_subtimestep_cm = surf_AET_vec[wetting_front_free_drainage_temp];
      }
    }

    if (PET_subtimestep_cm_per_h == 0.0){
      for (int k=1; k<listLength(state->head)+1; k++) {
        surf_AET_vec[k] = 0.0;
      }
    }

    for (int k=1; k<listLength(state->head_frac)+1; k++) {
      surf_AET_vec_frac[k] = 0.0;
    }


    precip_timestep_cm += precip_subtimestep_cm + precip_for_QF_subtimestep_cm_per_h*subtimestep_h;
    PET_timestep_cm += fmax(PET_subtimestep_cm,0.0); // ensures non-negative PET

    volstart_subtimestep_cm = (1 - ratio_fracture_vol_to_total_vol)*lgar_calc_mass_bal(state->lgar_bmi_params.cum_layer_thickness_cm, state->head);
    if (state->lgar_bmi_params.dual_perm) {
      volstart_subtimestep_cm_frac = ratio_fracture_vol_to_total_vol*lgar_calc_mass_bal(state->lgar_bmi_params.cum_layer_thickness_cm, state->head_frac);
    }

    //addressed machine precision issues where volon_timestep_error could be for example -1E-17 or 1.E-20 or smaller
    volon_timestep_cm = fmax(volon_timestep_cm,0.0);
    volon_timestep_cm = volon_timestep_cm > 1.0E-12 ? volon_timestep_cm : 0.0;

    int wf_free_drainage_demand = wetting_front_free_drainage(state->head);


     /*----------------------------------------------------------------------*/
    // Should a new wetting front be created in the soil matrix?
    int soil_num = state->lgar_bmi_params.layer_soil_type[state->head->layer_num];
    double theta_e = state->soil_properties[soil_num].theta_e;
    bool is_top_wf_saturated = false;
    bool top_near_sat = false;
    double psi_below_which_precip_contribs_to_GW = 5000000000000.0; //if psi goes below this value and the GW reservoir is enabled, then some fraction of precipitation will be directed to GW.
    double factor_for_simple_pref_flow = spf_factor;
    double theta_above_which_precip_contribs_to_GW = theta_e * factor_for_simple_pref_flow;
    //The idea here is that it should be more necessary in humid environments with highly conductive soils, such that contribtions to GW are expected but runoff is rare. In most semi arid or arid environments, it's probably not necessary / can be set to a low value.
    if (!state->lgar_bmi_params.TO_enabled){
      is_top_wf_saturated = (state->head->theta+1.0E-12) >= theta_e ? true : false; //sometimes a machine precision error would erroneously create a new wetting front during saturated conditions. The + 1.0E-12 seems to prevent this.
      top_near_sat = state->head->psi_cm < psi_below_which_precip_contribs_to_GW ? true : false; //is the top WF near saturation
      frac_to_GW_adjusted = fmin(1, pow((state->head->theta / theta_e), spf_factor));
      // top_near_sat = state->head->theta > theta_above_which_precip_contribs_to_GW ? true : false; //is the top WF near saturation
    }
    else {
      if (listLength_surface(state->head)>0){
        struct wetting_front *top_most_surface_WF;
        top_most_surface_WF = state->head;

        while (top_most_surface_WF->is_WF_GW==1){
          top_most_surface_WF = top_most_surface_WF->next;
        }

        int soil_num_highest_surf = state->lgar_bmi_params.layer_soil_type[top_most_surface_WF->layer_num];
        double theta_e_highest_surf = state->soil_properties[soil_num_highest_surf].theta_e;

        is_top_wf_saturated = (top_most_surface_WF->theta+1.0E-12) >= theta_e_highest_surf ? true : false; 
        top_near_sat = top_most_surface_WF->psi_cm < psi_below_which_precip_contribs_to_GW ? true : false; //is the top WF near saturation
        frac_to_GW_adjusted = fmin(1, pow((top_most_surface_WF->theta / theta_e_highest_surf), spf_factor));
        // top_near_sat = top_most_surface_WF->theta > theta_above_which_precip_contribs_to_GW ? true : false; //is the top WF near saturation

      }
      else {
        is_top_wf_saturated = (state->head->theta+1.0E-12) >= theta_e ? true : false;
        top_near_sat = state->head->psi_cm < psi_below_which_precip_contribs_to_GW ? true : false; //is the top WF near saturation
        frac_to_GW_adjusted = fmin(1, pow((state->head->theta / theta_e), spf_factor));
        // top_near_sat = state->head->psi_cm > theta_above_which_precip_contribs_to_GW ? true : false; //is the top WF near saturation
      }
    }

    // // checks on creating a new surficial front
    // // 1. check current and previous timestep precipitation
    // bool create_surficial_front = (precip_previous_subtimestep_cm == 0.0 && precip_subtimestep_cm > 0.0);

    bool create_surficial_front = (precip_previous_subtimestep_cm == 0.0 && precip_subtimestep_cm > 0.0 && volon_timestep_cm == 0) || ( (precip_subtimestep_cm > 0.0 || volon_timestep_cm > 0) && (listLength(state->head)==num_layers) && !(state->lgar_bmi_params.TO_enabled) ); //the volon_timestep_cm == 0 condition is necessary; a new wetting front can't be created if there is already ponded head greater than 0, at least in LGAR mode
    if (is_top_wf_saturated)
      create_surficial_front = false;

    if ( (state->lgar_bmi_params.TO_enabled) && ( (precip_subtimestep_cm > 0.0) || (volon_timestep_cm > 0.0) ) && (listLength_surface(state->head)==0) && (state->head->theta<(theta_e-1.0E-12)) ){ 
      create_surficial_front = true; //however in LGARTO mode, a new WF can be created if there is nonzero ponded head in the event that a surface WF merged with a TO WF and became TO, and there is still remaining ponded head 
    }

    if (verbosity.compare("high") == 0 || verbosity.compare("low") == 0) {
      if (state->lgar_bmi_params.dual_perm){
        printf("conditions for soil matrix \n");
      }
      std::string flag        = (create_surficial_front && !is_top_wf_saturated) == true ? "Yes" : "No";
      std::string flag_top_wf = is_top_wf_saturated == true ? "Yes" : "No";
      std::cerr<<"Is top wetting front saturated? "<< flag_top_wf  << "\n";
      std::cerr<<"Create superficial wetting front? "<< flag << "\n";
    }

    double temp_excess_water = 0.0;

    /*----------------------------------------------------------------------*/
    /* create a new wetting front if the following is true. */
    // if(create_surficial_front && !is_top_wf_saturated) {
    if(create_surficial_front) {

      double temp_pd = 0.0; // necessary to assign zero precip due to the creation of new wetting front; AET will still be taken out of the layers
      double temp_rch = 0.0;
      wf_free_drainage_demand = wetting_front_free_drainage(state->head);

      // move the wetting fronts without adding any water; this is done to close the mass balance
      temp_rch = lgar_move_wetting_fronts(state->lgar_bmi_params.TO_enabled, subtimestep_h, &free_drainage_subtimestep_cm, PET_subtimestep_cm_per_h/(1 - ratio_fracture_vol_to_total_vol), wilting_point_psi_cm, field_capacity_psi_cm, root_zone_depth_cm, 
             &temp_pd, wf_free_drainage_demand, volend_subtimestep_cm/(1-ratio_fracture_vol_to_total_vol),
			       num_layers, surf_frac_rz, &AET_subtimestep_cm, state->lgar_bmi_params.cum_layer_thickness_cm,
			       state->lgar_bmi_params.layer_soil_type, state->lgar_bmi_params.frozen_factor,
			       &state->head, state->state_previous, state->soil_properties, surf_AET_vec);

      if (temp_pd != 0.0){ //if temp_pd != 0.0, that means that some water left the model through the lower model bdy. this has been refactored such that temp_rch handles this now.
        // volrech_subtimestep_cm = temp_pd;
        // volrech_timestep_cm += volrech_subtimestep_cm;
        // temp_pd = 0.0;
      }

      // depth of the surficial front to be created 
      dry_depth = lgar_calc_dry_depth(use_closed_form_G, false, state->lgar_bmi_params.TO_enabled, nint, subtimestep_h, &delta_theta, state->lgar_bmi_params.layer_soil_type,
				      state->lgar_bmi_params.cum_layer_thickness_cm, state->lgar_bmi_params.frozen_factor,
				      state->head, state->soil_properties);
      
      double theta_for_new_wf = 0.0;
      struct wetting_front *top_most_surface_WF;
      top_most_surface_WF = state->head;

      //in LGARTO, there can be GW WFs that technically have a depth of 0 and are to the left of surface WFs. So the top surface WF won't necessarily just be head 
      if (listLength_surface(state->head)>0){
        while (top_most_surface_WF->is_WF_GW==1){
          top_most_surface_WF = top_most_surface_WF->next;
        }
      }

      //in the case where there are no surface WFs, the shallowest GW WF with a nonzero depth is taken as top_most_surface_WF, which is simply used as the most superficial moisture for creating a new WF 
      //note that that in this case top_most_surface_WF is technically a GW WF and not a surface WF 
      if (top_most_surface_WF->depth_cm==0.0){
        while (top_most_surface_WF->depth_cm==0){
          top_most_surface_WF = top_most_surface_WF->next;
        }
      }

      theta_for_new_wf = top_most_surface_WF->theta;

      if (verbosity.compare("high") == 0) {
        printf("State before moving creating new WF...\n");
        listPrint(state->head);
      }

      temp_excess_water = lgar_create_surficial_front(state->lgar_bmi_params.TO_enabled, num_layers, &ponded_depth_subtimestep_cm, &volin_subtimestep_cm, dry_depth, theta_for_new_wf,
				  state->lgar_bmi_params.layer_soil_type, state->lgar_bmi_params.cum_layer_thickness_cm,
				  state->lgar_bmi_params.frozen_factor, &state->head, state->soil_properties);

      if (verbosity.compare("high") == 0) {
        printf("State after moving creating new WF...\n");
        listPrint(state->head);
      }

      temp_rch += temp_excess_water;

      if(state->state_previous != NULL ){
        listDelete(state->state_previous);
        state->state_previous = NULL;
      }
      state->state_previous = listCopy(state->head);

      // volin_timestep_cm += volin_subtimestep_cm;
      volrech_subtimestep_cm = temp_rch;
      // volrech_timestep_cm += volrech_subtimestep_cm;

      if (verbosity.compare("high") == 0) {
        std::cerr<<"New wetting front created...\n";
        if (state->lgar_bmi_params.dual_perm){
          printf("soil matrix \n");
        }
        listPrint(state->head);
      }
    }














    /*----------------------------------------------------------------------*/
    /* infiltrate water based on the infiltration capacity given no new wetting front
        is created and that there is water on the surface (or raining), in soil matrix. */

    double precip_subtimestep_cm_per_h_both_domains = precip_subtimestep_cm_per_h;
    double precip_subtimestep_cm_per_h_frac = 0.0;
    if (frac_to_pref){
      precip_subtimestep_cm_per_h = (1 - frac_to_pref) * precip_subtimestep_cm_per_h_both_domains;
      precip_subtimestep_cm_per_h_frac = frac_to_pref * precip_subtimestep_cm_per_h_both_domains;
    }

    if (ponded_depth_subtimestep_cm > 0.0 && !create_surficial_front) {

      wf_free_drainage_demand = wetting_front_free_drainage(state->head);
      volrunoff_subtimestep_cm = lgar_insert_water(use_closed_form_G, false, nint, subtimestep_h, &free_drainage_subtimestep_cm, &AET_subtimestep_cm, &ponded_depth_subtimestep_cm,
						   &volin_subtimestep_cm, precip_subtimestep_cm_per_h,
						   wf_free_drainage_demand, num_layers,
						   ponded_depth_max_cm, state->lgar_bmi_params.layer_soil_type,
						   state->lgar_bmi_params.cum_layer_thickness_cm,
						   state->lgar_bmi_params.frozen_factor, state->head,
						   state->soil_properties); 

      if (verbosity.compare("high") == 0) {
        printf("AET_subtimestep_cm after lgar_insert_water \n");
        printf("AET_subtimestep_cm: %lf \n", AET_subtimestep_cm);
      }

      // volin_timestep_cm += volin_subtimestep_cm;
      // volrunoff_timestep_cm += volrunoff_subtimestep_cm;
      volrech_subtimestep_cm = volin_subtimestep_cm;
      volon_subtimestep_cm = ponded_depth_subtimestep_cm;
      if (volrunoff_subtimestep_cm < 0.0) abort();
    }
    else {

      if (ponded_depth_subtimestep_cm < ponded_depth_max_cm) {
        volrunoff_timestep_cm += 0.0;
        volon_subtimestep_cm = ponded_depth_subtimestep_cm;
        ponded_depth_subtimestep_cm = 0.0;
        volrunoff_subtimestep_cm = 0.0;
      }
      else {
        volrunoff_subtimestep_cm = (ponded_depth_subtimestep_cm - ponded_depth_max_cm);
        // volrunoff_timestep_cm += (ponded_depth_subtimestep_cm - ponded_depth_max_cm);
        volon_subtimestep_cm = ponded_depth_max_cm;

        ponded_depth_subtimestep_cm = ponded_depth_max_cm;
      }
    }
    //finally, following example of HYDRUS, if there is runoff simulated from soil matrix, send this to fracture domain
    if (volrunoff_subtimestep_cm>0.0 && state->lgar_bmi_params.dual_perm){
      ponded_depth_subtimestep_cm_frac += volrunoff_subtimestep_cm*(1 - ratio_fracture_vol_to_total_vol)/ratio_fracture_vol_to_total_vol;
      volrunoff_subtimestep_cm = 0.0;
    }
    /*----------------------------------------------------------------------*/
























     /*----------------------------------------------------------------------*/
    // Should a new wetting front be created in the fracture domain?
    bool is_top_wf_saturated_frac = false;
    bool top_near_sat_frac = false;
    bool create_surficial_front_frac = false;
    if (state->lgar_bmi_params.dual_perm) {
      int soil_num_frac = state->lgar_bmi_params.layer_soil_type[state->head_frac->layer_num];
      double theta_e_frac = state->soil_properties_frac[soil_num_frac].theta_e;
      // double psi_below_which_precip_contribs_to_GW = 5.0; //if psi goes below this value and the GW reservoir is enabled, then some fraction of precipitation will be directed to GW.
      double theta_above_which_precip_contribs_to_GW_frac = theta_e_frac * factor_for_simple_pref_flow;
      // //The idea here is that it should be more necessary in humid environments with highly conductive soils, such that contribtions to GW are expected but runoff is rare. In most semi arid or arid environments, it's probably not necessary / can be set to a low value.
      if (!state->lgar_bmi_params.TO_enabled){
        is_top_wf_saturated_frac = (state->head_frac->theta+1.0E-12) >= theta_e_frac ? true : false; //sometimes a machine precision error would erroneously create a new wetting front during saturated conditions. The + 1.0E-12 seems to prevent this.
        top_near_sat_frac = state->head_frac->psi_cm < psi_below_which_precip_contribs_to_GW ? true : false; //is the top WF near saturation
        // top_near_sat_frac = state->head_frac->theta > theta_above_which_precip_contribs_to_GW_frac ? true : false; //is the top WF near saturation
      }
      else {
        if (listLength_surface(state->head_frac)>0){
          struct wetting_front *top_most_surface_WF;
          top_most_surface_WF = state->head_frac;

          while (top_most_surface_WF->is_WF_GW==1){
            top_most_surface_WF = top_most_surface_WF->next;
          }

          int soil_num_highest_surf = state->lgar_bmi_params.layer_soil_type[top_most_surface_WF->layer_num];
          double theta_e_highest_surf = state->soil_properties_frac[soil_num_highest_surf].theta_e;

          is_top_wf_saturated_frac = (top_most_surface_WF->theta+1.0E-12) >= theta_e_highest_surf ? true : false; 
          top_near_sat_frac = top_most_surface_WF->psi_cm < psi_below_which_precip_contribs_to_GW ? true : false; //is the top WF near saturation
          // top_near_sat_frac = top_most_surface_WF->theta > theta_above_which_precip_contribs_to_GW_frac ? true : false;

        }
        else {
          is_top_wf_saturated_frac = (state->head->theta+1.0E-12) >= theta_e_frac ? true : false;
          top_near_sat_frac = state->head->psi_cm < psi_below_which_precip_contribs_to_GW ? true : false; //is the top WF near saturation
          // top_near_sat_frac = state->head->theta > theta_above_which_precip_contribs_to_GW_frac ? true : false; //is the top WF near saturation
        }
      }

      // // checks on creating a new surficial front
      // // 1. check current and previous timestep precipitation
      // bool create_surficial_front = (precip_previous_subtimestep_cm == 0.0 && precip_subtimestep_cm > 0.0);

      create_surficial_front_frac = (precip_previous_subtimestep_cm == 0.0 && precip_subtimestep_cm > 0.0 && volon_timestep_cm == 0) || ( (precip_subtimestep_cm > 0.0 || volon_timestep_cm > 0) && (listLength(state->head_frac)==num_layers) && !(state->lgar_bmi_params.TO_enabled) ); //the volon_timestep_cm == 0 condition is necessary; a new wetting front can't be created if there is already ponded head greater than 0, at least in LGAR mode
      if (is_top_wf_saturated_frac)
        create_surficial_front_frac = false;

      if ( (state->lgar_bmi_params.TO_enabled) && ( (precip_subtimestep_cm > 0.0) || (volon_timestep_cm > 0.0) ) && (listLength_surface(state->head_frac)==0) && (state->head_frac->theta<(theta_e_frac-1.0E-12)) ){ 
        create_surficial_front_frac = true; //however in LGARTO mode, a new WF can be created if there is nonzero ponded head in the event that a surface WF merged with a TO WF and became TO, and there is still remaining ponded head 
      }

      if (verbosity.compare("high") == 0 || verbosity.compare("low") == 0) {
        printf("conditions for fracture domain \n");
        std::string flag        = (create_surficial_front_frac && !is_top_wf_saturated_frac) == true ? "Yes" : "No";
        std::string flag_top_wf = is_top_wf_saturated_frac == true ? "Yes" : "No";
        std::cerr<<"Is top wetting front saturated? "<< flag_top_wf  << "\n";
        std::cerr<<"Create superficial wetting front? "<< flag << "\n";
        printf("but here's the frac list: \n");
        listPrint(state->head_frac);
      }

      double temp_excess_water = 0.0;

      /*----------------------------------------------------------------------*/
      /* create a new wetting front if the following is true. */
      // if(create_surficial_front && !is_top_wf_saturated) {
      if(create_surficial_front_frac) {

        double temp_rch = 0.0;
        double temp_pd_dummy = 0.0;
        double AET_dummy = 0.0;
        int wf_free_drainage_demand_frac = wetting_front_free_drainage(state->head_frac);
        temp_rch += lgar_move_wetting_fronts(state->lgar_bmi_params.TO_enabled, subtimestep_h, &free_drainage_subtimestep_cm_frac, 0.0, wilting_point_psi_cm, field_capacity_psi_cm, root_zone_depth_cm, 
          &temp_pd_dummy, wf_free_drainage_demand_frac, volend_subtimestep_cm_frac/ratio_fracture_vol_to_total_vol,
          num_layers, surf_frac_rz, &AET_dummy, state->lgar_bmi_params.cum_layer_thickness_cm,
          state->lgar_bmi_params.layer_soil_type, state->lgar_bmi_params.frozen_factor,
          &state->head_frac, state->state_previous_frac, state->soil_properties_frac, surf_AET_vec_frac);
        // }

        if (temp_pd_dummy != 0.0){ //if temp_pd != 0.0, that means that some water left the model through the lower model bdy. this has been refactored such that temp_rch handles this now.
          // volrech_subtimestep_cm_frac = temp_pd_dummy;
          // volrech_timestep_cm_frac += volrech_subtimestep_cm_frac;
          // temp_pd_dummy = 0.0;
        }

        dry_depth = lgar_calc_dry_depth(use_closed_form_G, true, state->lgar_bmi_params.TO_enabled, nint, subtimestep_h, &delta_theta, state->lgar_bmi_params.layer_soil_type,
                state->lgar_bmi_params.cum_layer_thickness_cm, state->lgar_bmi_params.frozen_factor,
                state->head_frac, state->soil_properties_frac);
        
        double theta_for_new_wf = 0.0;
        struct wetting_front *top_most_surface_WF;
        top_most_surface_WF = state->head_frac;

        //in LGARTO, there can be GW WFs that technically have a depth of 0 and are to the left of surface WFs. So the top surface WF won't necessarily just be head 
        if (listLength_surface(state->head_frac)>0){
          while (top_most_surface_WF->is_WF_GW==1){
            top_most_surface_WF = top_most_surface_WF->next;
          }
        }

        //in the case where there are no surface WFs, the shallowest GW WF with a nonzero depth is taken as top_most_surface_WF, which is simply used as the most superficial moisture for creating a new WF 
        //note that that in this case top_most_surface_WF is technically a GW WF and not a surface WF 
        if (top_most_surface_WF->depth_cm==0.0){
          while (top_most_surface_WF->depth_cm==0){
            top_most_surface_WF = top_most_surface_WF->next;
          }
        }

        theta_for_new_wf = top_most_surface_WF->theta;

        if (verbosity.compare("high") == 0) {
          printf("State before moving creating new WF, fracture domain...\n");
          listPrint(state->head_frac);
        }
        
        temp_excess_water = lgar_create_surficial_front(state->lgar_bmi_params.TO_enabled, num_layers, &ponded_depth_subtimestep_cm_frac, &volin_subtimestep_cm_frac, dry_depth, theta_for_new_wf,
            state->lgar_bmi_params.layer_soil_type, state->lgar_bmi_params.cum_layer_thickness_cm,
            state->lgar_bmi_params.frozen_factor, &state->head_frac, state->soil_properties_frac);

        if (verbosity.compare("high") == 0) {
          printf("State after moving creating new WF, fracture domain...\n");
          listPrint(state->head_frac);
        }

        temp_rch += temp_excess_water;

        if (state->state_previous_frac != NULL ){
          listDelete(state->state_previous_frac);
          state->state_previous_frac = NULL;
        }
        state->state_previous_frac = listCopy(state->head_frac);


        // volin_timestep_cm_frac += volin_subtimestep_cm_frac;
        volrech_subtimestep_cm_frac = temp_rch;
        // volrech_timestep_cm_frac += volrech_subtimestep_cm_frac;

        if (verbosity.compare("high") == 0) {
          std::cerr<<"New wetting front created...\n";
          printf("fracture domain \n");
          listPrint(state->head_frac);
        }
      }
    }










// /////comment this out 
//     /*----------------------------------------------------------------------*/
//     /* infiltrate water based on the infiltration capacity given no new wetting front
//         is created and that there is water on the surface (or raining), in soil matrix. */

//     double precip_subtimestep_cm_per_h_both_domains = precip_subtimestep_cm_per_h;
//     double precip_subtimestep_cm_per_h_frac = 0.0;
//     if (frac_to_pref){
//       precip_subtimestep_cm_per_h = (1 - frac_to_pref) * precip_subtimestep_cm_per_h_both_domains;
//       precip_subtimestep_cm_per_h_frac = frac_to_pref * precip_subtimestep_cm_per_h_both_domains;
//     }

//     if (ponded_depth_subtimestep_cm > 0.0 && !create_surficial_front) {

//       wf_free_drainage_demand = wetting_front_free_drainage(state->head);
//       volrunoff_subtimestep_cm = lgar_insert_water(use_closed_form_G, nint, subtimestep_h, &free_drainage_subtimestep_cm, &AET_subtimestep_cm, &ponded_depth_subtimestep_cm,
// 						   &volin_subtimestep_cm, precip_subtimestep_cm_per_h,
// 						   wf_free_drainage_demand, num_layers,
// 						   ponded_depth_max_cm, state->lgar_bmi_params.layer_soil_type,
// 						   state->lgar_bmi_params.cum_layer_thickness_cm,
// 						   state->lgar_bmi_params.frozen_factor, state->head,
// 						   state->soil_properties); 

//       if (verbosity.compare("high") == 0) {
//         printf("AET_subtimestep_cm after lgar_insert_water \n");
//         printf("AET_subtimestep_cm: %lf \n", AET_subtimestep_cm);
//       }

//       // volin_timestep_cm += volin_subtimestep_cm;
//       // volrunoff_timestep_cm += volrunoff_subtimestep_cm;
//       volrech_subtimestep_cm = volin_subtimestep_cm;
//       volon_subtimestep_cm = ponded_depth_subtimestep_cm;
//       printf("segfault 10 \n");
//       if (volrunoff_subtimestep_cm < 0.0) abort();
//     }
//     else {

//       if (ponded_depth_subtimestep_cm < ponded_depth_max_cm) {
//         volrunoff_timestep_cm += 0.0;
//         volon_subtimestep_cm = ponded_depth_subtimestep_cm;
//         ponded_depth_subtimestep_cm = 0.0;
//         volrunoff_subtimestep_cm = 0.0;
//       }
//       else {
//         volrunoff_subtimestep_cm = (ponded_depth_subtimestep_cm - ponded_depth_max_cm);
//         // volrunoff_timestep_cm += (ponded_depth_subtimestep_cm - ponded_depth_max_cm);
//         volon_subtimestep_cm = ponded_depth_max_cm;

//         ponded_depth_subtimestep_cm = ponded_depth_max_cm;
//       }
//     }
//     //finally, following example of HYDRUS, if there is runoff simulated from soil matrix, send this to fracture domain
//     if (volrunoff_subtimestep_cm>0.0 && state->lgar_bmi_params.dual_perm){
//       ponded_depth_subtimestep_cm_frac += volrunoff_subtimestep_cm*(1 - ratio_fracture_vol_to_total_vol)/ratio_fracture_vol_to_total_vol;
//       volrunoff_subtimestep_cm = 0.0;
//     }
//     /*----------------------------------------------------------------------*/







    

    /*----------------------------------------------------------------------*/
    /* infiltrate water based on the infiltration capacity given no new wetting front
       is created and that there is water on the surface (or raining), in fracture domain. */

    //   //first, following example of HYDRUS, if there is runoff simulated from soil matrix, send this to fracture domain
    // if (volrunoff_subtimestep_cm>0.0 && state->lgar_bmi_params.dual_perm){
    //   ponded_depth_subtimestep_cm_frac += volrunoff_subtimestep_cm*(1 - ratio_fracture_vol_to_total_vol)/ratio_fracture_vol_to_total_vol;
    //   volrunoff_subtimestep_cm = 0.0;
    // }
    if (state->lgar_bmi_params.dual_perm){
      if (ponded_depth_subtimestep_cm_frac > 0.0 && !create_surficial_front_frac) {
        double AET_dummy = 0.0;
        int wf_free_drainage_demand_frac = wetting_front_free_drainage(state->head_frac);
        volrunoff_subtimestep_cm_frac = lgar_insert_water(use_closed_form_G, true, nint, subtimestep_h, &free_drainage_subtimestep_cm_frac, &AET_dummy, &ponded_depth_subtimestep_cm_frac,
                &volin_subtimestep_cm_frac, precip_subtimestep_cm_per_h_frac,
                wf_free_drainage_demand_frac, num_layers,
                ponded_depth_max_cm, state->lgar_bmi_params.layer_soil_type,
                state->lgar_bmi_params.cum_layer_thickness_cm,
                state->lgar_bmi_params.frozen_factor, state->head_frac,
                state->soil_properties_frac); 
        if (verbosity.compare("high") == 0) {
          printf("AET_subtimestep_cm after lgar_insert_water for frac domain. Should be unaffected because frac domain does not contribute to AET \n");
          printf("AET_subtimestep_cm: %lf \n", AET_subtimestep_cm);
        }

        // volin_timestep_cm_frac += volin_subtimestep_cm_frac;
        // volrunoff_timestep_cm_frac += volrunoff_subtimestep_cm_frac;
        volrech_subtimestep_cm_frac = volin_subtimestep_cm_frac;

        volon_subtimestep_cm_frac = ponded_depth_subtimestep_cm_frac;
        if (volrunoff_subtimestep_cm_frac < 0.0) abort();
      }
      else {
        

        if (ponded_depth_subtimestep_cm_frac < ponded_depth_max_cm) {
          volrunoff_timestep_cm_frac += 0.0;
          volon_subtimestep_cm_frac = ponded_depth_subtimestep_cm_frac;
          ponded_depth_subtimestep_cm_frac = 0.0;
          volrunoff_subtimestep_cm_frac = 0.0;
        }
        else {
          volrunoff_subtimestep_cm_frac = (ponded_depth_subtimestep_cm_frac - ponded_depth_max_cm);
          // volrunoff_timestep_cm_frac += (ponded_depth_subtimestep_cm_frac - ponded_depth_max_cm);
          volon_subtimestep_cm_frac = ponded_depth_max_cm;
          ponded_depth_subtimestep_cm_frac = ponded_depth_max_cm;
        }
      }
    }
    /*----------------------------------------------------------------------*/


    /* move wetting fronts if no new wetting front is created in soil matrix. Otherwise, movement
       of wetting fronts has already happened at the time of creating surficial front,
       so no need to move them here. */
    if (!create_surficial_front) {
      double temp_rch = 0.0;
      double volin_subtimestep_cm_temp = volin_subtimestep_cm;  /* passing this for mass balance only, the method modifies it
								   and returns percolated value, so we need to keep its original
								   value stored to copy it back*/
      wf_free_drainage_demand = wetting_front_free_drainage(state->head);
      temp_rch = lgar_move_wetting_fronts(state->lgar_bmi_params.TO_enabled, subtimestep_h, &free_drainage_subtimestep_cm, PET_subtimestep_cm_per_h/(1 - ratio_fracture_vol_to_total_vol), wilting_point_psi_cm, field_capacity_psi_cm, root_zone_depth_cm, 
             &volin_subtimestep_cm, wf_free_drainage_demand, volend_subtimestep_cm/(1-ratio_fracture_vol_to_total_vol),
			       num_layers, surf_frac_rz, &AET_subtimestep_cm, state->lgar_bmi_params.cum_layer_thickness_cm,
			       state->lgar_bmi_params.layer_soil_type, state->lgar_bmi_params.frozen_factor,
			       &state->head, state->state_previous, state->soil_properties, surf_AET_vec);

      volrech_subtimestep_cm = temp_rch;
      // volrech_timestep_cm += volrech_subtimestep_cm;

      volin_subtimestep_cm = volin_subtimestep_cm_temp;
    }


    /* move wetting fronts if no new wetting front is created in fracture domain. Otherwise, movement
       of wetting fronts has already happened at the time of creating surficial front,
       so no need to move them here. */
    if (state->lgar_bmi_params.dual_perm){
      if (!create_surficial_front_frac) {
        double temp_rch = 0.0;
        double volin_subtimestep_cm_temp = volin_subtimestep_cm_frac;  /* passing this for mass balance only, the method modifies it
                    and returns percolated value, so we need to keep its original
                    value stored to copy it back*/
        double AET_dummy = 0.0;
        int wf_free_drainage_demand_frac = wetting_front_free_drainage(state->head_frac);
        temp_rch = lgar_move_wetting_fronts(state->lgar_bmi_params.TO_enabled, subtimestep_h, &free_drainage_subtimestep_cm_frac, 0.0, wilting_point_psi_cm, field_capacity_psi_cm, root_zone_depth_cm, 
              &volin_subtimestep_cm_frac, wf_free_drainage_demand_frac, volend_subtimestep_cm_frac/ratio_fracture_vol_to_total_vol,
              num_layers, surf_frac_rz, &AET_dummy, state->lgar_bmi_params.cum_layer_thickness_cm,
              state->lgar_bmi_params.layer_soil_type, state->lgar_bmi_params.frozen_factor,
              &state->head_frac, state->state_previous_frac, state->soil_properties_frac, surf_AET_vec_frac);
        if (AET_dummy!=0.0){
          printf("AET_dummy: %lf \n", AET_dummy);
          // abort();
        }
        volrech_subtimestep_cm_frac = temp_rch;
        // volrech_timestep_cm_frac += volrech_subtimestep_cm_frac;
        volin_subtimestep_cm_frac = volin_subtimestep_cm_temp;
      }
    }


    lgar_clean_redundant_fronts(&state->head, state->lgar_bmi_params.layer_soil_type, state->soil_properties); //in the event that the soil profile was completely saturated and there was a new WF that was intense enough to maintain complete saturation, a redundant WF was created and must be deleted
    //not sure if this is still necessary however
    if (state->lgar_bmi_params.dual_perm){
      lgar_clean_redundant_fronts(&state->head_frac, state->lgar_bmi_params.layer_soil_type, state->soil_properties_frac);
    }

    volrech_subtimestep_cm += free_drainage_subtimestep_cm;
    // volrech_timestep_cm += free_drainage_subtimestep_cm;

    volrech_subtimestep_cm_frac += free_drainage_subtimestep_cm_frac;
    // volrech_timestep_cm_frac += free_drainage_subtimestep_cm_frac;

    if (ratio_fracture_vol_to_total_vol){
      volon_subtimestep_cm = (1 - ratio_fracture_vol_to_total_vol)*volon_subtimestep_cm;
      volon_subtimestep_cm_frac = ratio_fracture_vol_to_total_vol*volon_subtimestep_cm_frac;
      volon_subtimestep_cm = volon_subtimestep_cm + volon_subtimestep_cm_frac; //aggregate all surface water for next time step and for mass balance
    }



    // volon_subtimestep_cm += volon_subtimestep_cm_frac; //aggregate all surface water


    // if (state->lgar_bmi_params.dual_perm && (listLength(state->head_frac)!=listLength(state->head))){
    if (state->lgar_bmi_params.dual_perm){
      double inc_matrix_rch = 0.0;
      double inc_fracture_rch = 0.0;
      mass_transfer_subtimestep_cm = lgar_dual_perm_mass_transfer(state->lgar_bmi_params.TO_enabled, state->lgar_bmi_params.num_layers, ratio_fracture_vol_to_total_vol, subtimestep_h, &inc_matrix_rch, &inc_fracture_rch, &state->head, &state->head_frac, state->lgar_bmi_params.cum_layer_thickness_cm,
                                   state->lgar_bmi_params.frozen_factor, state->lgar_bmi_params.layer_soil_type, state->soil_properties, state->soil_properties_frac, state->mass_transfer_soil_properties);
      volrech_subtimestep_cm += inc_matrix_rch*(1 - ratio_fracture_vol_to_total_vol);
      // volrech_timestep_cm += inc_matrix_rch;

      volrech_subtimestep_cm_frac += inc_fracture_rch*ratio_fracture_vol_to_total_vol;
      // volrech_timestep_cm_frac += inc_fracture_rch;
      // abort();
    }


    /*----------------------------------------------------------------------*/
    // calculate derivative (dz/dt) for all wetting fronts
    lgar_dzdt_calc(use_closed_form_G, nint, subtimestep_h, ponded_depth_subtimestep_cm, state->lgar_bmi_params.layer_soil_type,
		   state->lgar_bmi_params.cum_layer_thickness_cm, state->lgar_bmi_params.frozen_factor,
		   state->head, state->soil_properties, state->lgar_bmi_params.num_layers, false);

    if (state->lgar_bmi_params.dual_perm){
      lgar_dzdt_calc(use_closed_form_G, nint, subtimestep_h, ponded_depth_subtimestep_cm, state->lgar_bmi_params.layer_soil_type, //ponded_depth_subtimestep_cm_frac instead?
        state->lgar_bmi_params.cum_layer_thickness_cm, state->lgar_bmi_params.frozen_factor,
        state->head_frac, state->soil_properties_frac, state->lgar_bmi_params.num_layers, true);
    }

    volend_subtimestep_cm = (1-ratio_fracture_vol_to_total_vol)*lgar_calc_mass_bal(state->lgar_bmi_params.cum_layer_thickness_cm, state->head);
    volend_timestep_cm = volend_subtimestep_cm;
    state->lgar_bmi_params.precip_previous_timestep_cm = precip_subtimestep_cm;

    volend_subtimestep_cm_frac = ratio_fracture_vol_to_total_vol*lgar_calc_mass_bal(state->lgar_bmi_params.cum_layer_thickness_cm, state->head_frac);
    volend_timestep_cm_frac = volend_subtimestep_cm_frac;

    double QF_storage_start_cm = state->lgar_mass_balance.QF_storage_cm;
    double QF_Q_subtimestep_cm = calc_QF_Q(subtimestep_h, a, b, precip_for_QF_subtimestep_cm_per_h, &state->lgar_mass_balance.QF_storage_cm);
    state->lgar_mass_balance.volrunoff_QF_cm += QF_Q_subtimestep_cm;
    QF_Q_timestep_cm += QF_Q_subtimestep_cm;

    // set runoff_in_prev_step for next step
    if ((volrunoff_subtimestep_cm > 1.E-9) || (top_near_sat) || (top_near_sat_frac)){
      state->lgar_bmi_params.runoff_in_prev_step = true;
    }
    else {
      state->lgar_bmi_params.runoff_in_prev_step = false;
    }

    // add QF_Q_subtimestep_cm into volrunoff for the time step for routing through the GIUH
    // volrunoff_timestep_cm += QF_Q_subtimestep_cm;
    //add precip_for_QF_subtimestep_cm_per_h back into precip for mass balance
    precip_subtimestep_cm += precip_for_QF_subtimestep_cm_per_h * subtimestep_h;

    /*----------------------------------------------------------------------*/
    // mass balance at the subtimestep (local mass balance)

    //normalize fluxes considering fraction of total volume that is fracture domain and volume that is soil matrix
    volrunoff_subtimestep_cm = (1 - ratio_fracture_vol_to_total_vol)*volrunoff_subtimestep_cm;
    volrunoff_subtimestep_cm_frac = ratio_fracture_vol_to_total_vol*volrunoff_subtimestep_cm_frac;
    AET_subtimestep_cm = (1 - ratio_fracture_vol_to_total_vol)*AET_subtimestep_cm;
    volrech_subtimestep_cm = (1 - ratio_fracture_vol_to_total_vol)*volrech_subtimestep_cm;
    volrech_subtimestep_cm_frac = ratio_fracture_vol_to_total_vol*volrech_subtimestep_cm_frac;
    volin_subtimestep_cm = (1 - ratio_fracture_vol_to_total_vol)*volin_subtimestep_cm;
    volin_subtimestep_cm_frac = ratio_fracture_vol_to_total_vol*volin_subtimestep_cm_frac;


    double local_mb = volstart_subtimestep_cm + volstart_subtimestep_cm_frac + precip_subtimestep_cm + volon_timestep_cm - volrunoff_subtimestep_cm - volrunoff_subtimestep_cm_frac - QF_Q_subtimestep_cm - state->lgar_mass_balance.QF_storage_cm + QF_storage_start_cm 
                      - AET_subtimestep_cm - volon_subtimestep_cm - volrech_subtimestep_cm - volrech_subtimestep_cm_frac - volend_subtimestep_cm - volend_subtimestep_cm_frac;


    // increment runoff for the subtimestep
    surface_runoff_subtimestep_cm = volrunoff_subtimestep_cm + volrunoff_subtimestep_cm_frac;
    surface_runoff_timestep_cm += surface_runoff_subtimestep_cm ;

    // adding groundwater flux to stream channel (note: this will be updated/corrected after adding the groundwater reservoir)
    // note that after adding free drainage, these fluxes will go to a nonlinear GW reservoir that will contribute to streamflow in LGARTO mode, but are assumed to be lost to deep GW in LGAR mode
    volQ_gw_timestep_cm += volQ_gw_subtimestep_cm;














    //separating code such that most non substep vars (so xxx_timestep and not xxx_subtimestep) are updated in just one place. not all, because some must be set before substepping.
    volin_timestep_cm += volin_subtimestep_cm;
    volrech_timestep_cm += volrech_subtimestep_cm;
    volin_timestep_cm_frac += volin_subtimestep_cm_frac;
    volrech_timestep_cm_frac += volrech_subtimestep_cm_frac;


    volrunoff_timestep_cm += volrunoff_subtimestep_cm;
    volrunoff_timestep_cm_frac += volrunoff_subtimestep_cm_frac;

    AET_timestep_cm += AET_subtimestep_cm;
    volon_timestep_cm = volon_subtimestep_cm; // surface ponded water at the end of the timestep 

    mass_transfer_timestep_cm += mass_transfer_subtimestep_cm;

















    if (verbosity.compare("high") == 0 || verbosity.compare("low") == 0) {
      printf("Printing wetting fronts at this subtimestep... \n");
      printf("soil matrix: \n ");
      listPrint(state->head);
      if (state->lgar_bmi_params.dual_perm){
        printf("fracture domain: \n ");
        listPrint(state->head_frac);
      }
    }

    bool unexpected_local_error = fabs(local_mb) > mbal_tol ? true : false; //1.0E-4 was the default for LASAM stability testing 
    if (isinf(local_mb)){
      unexpected_local_error = true;
    }

    // if (volrech_subtimestep_cm<-0.001){
    //   unexpected_local_error = true; //testing some weird recharge behavior in dual perm mode
    // }

    
    if (verbosity.compare("high") == 0 || verbosity.compare("low") == 0 || unexpected_local_error) {

      if (state->lgar_bmi_params.dual_perm) {
        printf("\nLocal mass balance at this timestep... \n\
        Error              = %14.10f \n\
        Initial water      = %14.10f \n\
        Initial water frac = %14.10f \n\
        Water added        = %14.10f \n\
        Ponded water       = %14.10f \n\
        Infiltration       = %14.10f \n\
        Infiltration frac  = %14.10f \n\
        Runoff             = %14.10f \n\
        Runoff frac        = %14.10f \n\
        AET                = %14.10f \n\
        Percolation        = %14.10f \n\
        Percolation frac   = %14.10f \n\
        Mass transfer      = %14.10f \n\
        Final water        = %14.10f \n\
        Final water frac   = %14.10f \n", local_mb, volstart_subtimestep_cm, volstart_subtimestep_cm_frac, precip_subtimestep_cm, volon_subtimestep_cm,
        volin_subtimestep_cm, volin_subtimestep_cm_frac, volrunoff_subtimestep_cm, volrunoff_subtimestep_cm_frac, AET_subtimestep_cm, volrech_subtimestep_cm,
        volrech_subtimestep_cm_frac, mass_transfer_subtimestep_cm, volend_subtimestep_cm, volend_subtimestep_cm_frac);
      }
      else {
        printf("\nLocal mass balance at this timestep... \n\
        Error         = %14.10f \n\
        Initial water = %14.10f \n\
        Water added   = %14.10f \n\
        Ponded water  = %14.10f \n\
        Infiltration  = %14.10f \n\
        Runoff        = %14.10f \n\
        AET           = %14.10f \n\
        Percolation   = %14.10f \n\
        Final water   = %14.10f \n", local_mb, volstart_subtimestep_cm, precip_subtimestep_cm, volon_subtimestep_cm,
        volin_subtimestep_cm, volrunoff_subtimestep_cm, AET_subtimestep_cm, volrech_subtimestep_cm,
        volend_subtimestep_cm);
      }

      if (unexpected_local_error) {
	printf("Local mass balance (in this timestep) is %14.10f, larger than expected, needs some debugging...\n ",local_mb);
	abort();
      }

    }

    if (isnan(volrech_subtimestep_cm)){
      listPrint(state->head);
      printf("volrech_subtimestep_cm is not a number \n");
      abort();
    }

    free(surf_AET_vec);
    free(surf_AET_vec_frac);

    // store local mass balance error to the struct
    state->lgar_mass_balance.local_mass_balance = local_mb;
    state->lgar_mass_balance.cumulative_mbal += local_mb;

    assert (state->head->depth_cm > -1e-7); // check on negative layer depth --> move this to somewhere else AJ (later)
    if (state->lgar_bmi_params.dual_perm){
      assert (state->head_frac->depth_cm > -1e-7); // check on negative layer depth --> move this to somewhere else AJ (later)
    }

    bool lasam_standalone = true;
#ifdef NGEN
    lasam_standalone = false;
#endif
    // simuation time can't exceed the endtime when running standalone
    if ( (this->state->lgar_bmi_params.time_s >= this->state->lgar_bmi_params.endtime_s) && lasam_standalone)
      break;

  } // end of subcycling

  //update giuh at the time step level (was previously updated at the sub time step level)
  volrunoff_giuh_timestep_cm = giuh_convolution_integral(volrunoff_timestep_cm + QF_Q_timestep_cm + volrunoff_timestep_cm_frac, num_giuh_ordinates, giuh_ordinates, giuh_runoff_queue);

  // total mass of water leaving the system, at this time it is the giuh-only, but later will add groundwater component as well.
  // for LGARTO, I believe that volQ will be streamflow via GIUH, and indeed a streamflow component due to GW will be added once we add a GW reservoir
  volQ_timestep_cm = volrunoff_giuh_timestep_cm;



  /*----------------------------------------------------------------------*/
  // Everything related to lgar state is done at this point, now time to update some dynamic variables

  // update number of wetting fronts
  state->lgar_bmi_params.num_wetting_fronts = listLength(state->head);

  // allocate new memory based on updated wetting fronts; we could make it conditional i.e. create only if no. of wf are changed
  realloc_soil();

  // update thickness/depth and soil moisture of wetting fronts (used for state coupling)
  struct wetting_front *current = state->head;
  for (int i=0; i<state->lgar_bmi_params.num_wetting_fronts; i++) {
    assert (current != NULL);
    state->lgar_bmi_params.soil_moisture_wetting_fronts[i] = current->theta;
    state->lgar_bmi_params.soil_depth_wetting_fronts[i] = current->depth_cm * state->units.cm_to_m;
    current = current->next;
    if (verbosity.compare("high") == 0)
      std::cerr<<"Wetting fronts (bmi outputs) (depth in meters, theta)= "
	       <<state->lgar_bmi_params.soil_depth_wetting_fronts[i]
	       <<" "<<state->lgar_bmi_params.soil_moisture_wetting_fronts[i]<<"\n";
  }
  
  // add to mass balance timestep variables
  state->lgar_mass_balance.volprecip_timestep_cm      = precip_timestep_cm;
  state->lgar_mass_balance.volin_timestep_cm          = volin_timestep_cm;
  state->lgar_mass_balance.volin_timestep_cm_frac     = volin_timestep_cm_frac;
  state->lgar_mass_balance.volon_timestep_cm          = volon_timestep_cm;
  state->lgar_mass_balance.volend_timestep_cm         = volend_timestep_cm;
  state->lgar_mass_balance.volend_timestep_cm_frac    = volend_timestep_cm_frac;
  state->lgar_mass_balance.volAET_timestep_cm         = AET_timestep_cm;
  state->lgar_mass_balance.volrech_timestep_cm        = volrech_timestep_cm;
  state->lgar_mass_balance.volrech_timestep_cm_frac   = volrech_timestep_cm_frac;
  state->lgar_mass_balance.volrunoff_timestep_cm      = volrunoff_timestep_cm;
  state->lgar_mass_balance.volrunoff_timestep_cm_frac = volrunoff_timestep_cm_frac;
  state->lgar_mass_balance.volQ_timestep_cm           = volQ_timestep_cm;
  state->lgar_mass_balance.volQ_gw_timestep_cm        = volQ_gw_timestep_cm;
  state->lgar_mass_balance.volPET_timestep_cm         = PET_timestep_cm;
  state->lgar_mass_balance.volrunoff_giuh_timestep_cm = volrunoff_giuh_timestep_cm;
  state->lgar_mass_balance.mass_transfer_timestep_cm  = mass_transfer_timestep_cm;

  // add to mass balance accumulated variables
  state->lgar_mass_balance.volprecip_cm       += precip_timestep_cm;
  state->lgar_mass_balance.volin_cm           += volin_timestep_cm;
  state->lgar_mass_balance.volin_cm_frac      += volin_timestep_cm_frac;
  state->lgar_mass_balance.volon_cm            = volon_timestep_cm;
  state->lgar_mass_balance.volend_cm           = volend_timestep_cm;
  state->lgar_mass_balance.volend_cm_frac      = volend_timestep_cm_frac;
  state->lgar_mass_balance.volAET_cm          += AET_timestep_cm;
  state->lgar_mass_balance.volrech_cm         += volrech_timestep_cm;
  state->lgar_mass_balance.volrech_cm_frac    += volrech_timestep_cm_frac;
  state->lgar_mass_balance.volrunoff_cm       += volrunoff_timestep_cm;
  state->lgar_mass_balance.volrunoff_cm_frac  += volrunoff_timestep_cm_frac;
  state->lgar_mass_balance.volQ_cm            += volQ_timestep_cm;
  state->lgar_mass_balance.volQ_gw_cm         += volQ_gw_timestep_cm;
  state->lgar_mass_balance.volPET_cm          += PET_timestep_cm;
  state->lgar_mass_balance.volrunoff_giuh_cm  += volrunoff_giuh_timestep_cm;
  state->lgar_mass_balance.volchange_calib_cm += volchange_calib_cm ;
  state->lgar_mass_balance.mass_transfer_cm   += mass_transfer_timestep_cm;
 
  // converted values, a struct local to the BMI and has bmi output variables
  //needs volon?
  bmi_unit_conv.mass_balance_m            = state->lgar_mass_balance.local_mass_balance * state->units.cm_to_m;
  bmi_unit_conv.volprecip_timestep_m      = precip_timestep_cm * state->units.cm_to_m;
  bmi_unit_conv.volin_timestep_m          = volin_timestep_cm * state->units.cm_to_m;
  bmi_unit_conv.volin_timestep_m_frac     = volin_timestep_cm_frac * state->units.cm_to_m;
  bmi_unit_conv.volend_timestep_m         = volend_timestep_cm * state->units.cm_to_m;
  bmi_unit_conv.volend_timestep_m_frac    = volend_timestep_cm_frac * state->units.cm_to_m;
  bmi_unit_conv.volAET_timestep_m         = AET_timestep_cm * state->units.cm_to_m;
  bmi_unit_conv.volrech_timestep_m        = volrech_timestep_cm * state->units.cm_to_m;
  bmi_unit_conv.volrech_timestep_m_frac   = volrech_timestep_cm_frac * state->units.cm_to_m;
  bmi_unit_conv.volrunoff_timestep_m      = volrunoff_timestep_cm * state->units.cm_to_m;
  bmi_unit_conv.volrunoff_timestep_m_frac = volrunoff_timestep_cm_frac * state->units.cm_to_m;
  bmi_unit_conv.volQ_timestep_m           = volQ_timestep_cm * state->units.cm_to_m;
  bmi_unit_conv.volQ_gw_timestep_m        = volQ_gw_timestep_cm * state->units.cm_to_m;
  bmi_unit_conv.volPET_timestep_m         = PET_timestep_cm * state->units.cm_to_m;
  bmi_unit_conv.volrunoff_giuh_timestep_m = volrunoff_giuh_timestep_cm * state->units.cm_to_m;
  bmi_unit_conv.mass_transfer_timestep_m  = mass_transfer_timestep_cm * state->units.cm_to_m;

}


void BmiLGAR::
UpdateUntil(double t)
{
  assert (t > 0.0);
  this->Update();
}

struct model_state* BmiLGAR::get_model()
{
  return state;
}

void BmiLGAR::
global_mass_balance()
{
  lgar_global_mass_balance(this->state, giuh_runoff_queue);
}

double BmiLGAR::
update_calibratable_parameters(bool dual_perm, double ratio_fracture_vol_to_total_vol)
{
  int soil, layer_num;
  struct wetting_front *current = state->head;

  if (verbosity.compare("high") == 0)
    listPrint(state->head);
  
  double volstart_before = (1 - ratio_fracture_vol_to_total_vol)*lgar_calc_mass_bal(state->lgar_bmi_params.cum_layer_thickness_cm, state->head);

  if (dual_perm){
    volstart_before += ratio_fracture_vol_to_total_vol*lgar_calc_mass_bal(state->lgar_bmi_params.cum_layer_thickness_cm, state->head_frac);
  }

  // for (int i=0; i<state->lgar_bmi_params.num_wetting_fronts; i++) {//first we update the parameters that depend on soil layer, for each layer
  for (int i=1; i<state->lgar_bmi_params.num_layers + 1; i++) {//first we update the parameters that depend on soil layer, for each layer (not each WF)
    // layer_num  = current->layer_num;
    layer_num = i;
    soil = state->lgar_bmi_params.layer_soil_type[layer_num];
    
    assert (current != NULL);

    if (verbosity.compare("high") == 0 || verbosity.compare("low") == 0) {
      std::cerr<<"----------- Calibratable parameters depending on soil layer (initial values) ----------- \n";
      std::cerr<<"| soil_type = "<< soil <<", layer = "<<layer_num
	       <<", smcmax = "   << state->soil_properties[soil].theta_e
	       <<", smcmin = "   << state->soil_properties[soil].theta_r
	       <<", vg_n = "     << state->soil_properties[soil].vg_n
	       <<", vg_alpha = " << state->soil_properties[soil].vg_alpha_per_cm
	       <<", Ksat = "     << state->soil_properties[soil].Ksat_cm_per_h <<"\n";
	      //  <<", theta = "    << current->theta <<"\n"; //theta is not a calibratable parameter

      if (dual_perm){
        std::cerr<<"| soil_type = "<< soil <<", layer = "<<layer_num
          <<", smcmax = "   << state->soil_properties_frac[soil].theta_e
          <<", smcmin = "   << state->soil_properties_frac[soil].theta_r
          <<", vg_n = "     << state->soil_properties_frac[soil].vg_n
          <<", vg_alpha = " << state->soil_properties_frac[soil].vg_alpha_per_cm
          <<", Ksat = "     << state->soil_properties_frac[soil].Ksat_cm_per_h 

          <<", a_f = "      << state->mass_transfer_soil_properties[soil].a_f
          <<", beta_f = "   << state->mass_transfer_soil_properties[soil].beta_f
          <<", K_sa_f = "   << state->mass_transfer_soil_properties[soil].K_sa_f
          <<", gamma_f = "  << state->mass_transfer_soil_properties[soil].gamma_f
          
          <<"\n";
      }
    }
    
    state->soil_properties[soil].theta_e = state->lgar_calib_params.theta_e[layer_num-1];
    state->soil_properties[soil].theta_r = state->lgar_calib_params.theta_r[layer_num-1];
    state->soil_properties[soil].vg_n    = state->lgar_calib_params.vg_n[layer_num-1];
    state->soil_properties[soil].vg_m    = 1.0 - 1.0/state->soil_properties[soil].vg_n;
    state->soil_properties[soil].vg_alpha_per_cm = state->lgar_calib_params.vg_alpha[layer_num-1];
    state->soil_properties[soil].Ksat_cm_per_h   = state->lgar_calib_params.Ksat[layer_num-1];
    
    // current->theta = calc_theta_from_h(current->psi_cm, state->soil_properties[soil].vg_alpha_per_cm,
		// 		       state->soil_properties[soil].vg_m, state->soil_properties[soil].vg_n,
		// 		       state->soil_properties[soil].theta_e, state->soil_properties[soil].theta_r);

    if (dual_perm){
      state->soil_properties_frac[soil].theta_e = state->lgar_calib_params.theta_e_f[layer_num-1];
      state->soil_properties_frac[soil].theta_r = state->lgar_calib_params.theta_r_f[layer_num-1];
      state->soil_properties_frac[soil].vg_n    = state->lgar_calib_params.vg_n_f[layer_num-1];
      state->soil_properties_frac[soil].vg_m    = 1.0 - 1.0/state->soil_properties_frac[soil].vg_n;
      state->soil_properties_frac[soil].vg_alpha_per_cm = state->lgar_calib_params.vg_alpha_f[layer_num-1];
      state->soil_properties_frac[soil].Ksat_cm_per_h   = state->lgar_calib_params.Ksat_f[layer_num-1];

      // current->theta = calc_theta_from_h(current->psi_cm, state->soil_properties[soil].vg_alpha_per_cm,
			// 	       state->soil_properties[soil].vg_m, state->soil_properties[soil].vg_n,
			// 	       state->soil_properties[soil].theta_e, state->soil_properties[soil].theta_r_f);

      state->mass_transfer_soil_properties[soil].a_f     = state->lgar_calib_params.a_f[layer_num-1];
      state->mass_transfer_soil_properties[soil].beta_f  = state->lgar_calib_params.beta_f[layer_num-1];
      state->mass_transfer_soil_properties[soil].K_sa_f  = state->lgar_calib_params.K_sa_f[layer_num-1];
      state->mass_transfer_soil_properties[soil].gamma_f = state->lgar_calib_params.gamma_f[layer_num-1];

    }

    if (verbosity.compare("high") == 0 || verbosity.compare("low") == 0) {
      std::cerr<<"----------- Calibratable parameters depending on soil layer (updated values) ----------- \n";
      std::cerr<<"| soil_type = "<< soil <<", layer = "<<layer_num
	       <<", smcmax = "   << state->soil_properties[soil].theta_e
	       <<", smcmin = "   << state->soil_properties[soil].theta_r
	       <<", vg_n = "     << state->soil_properties[soil].vg_n
	       <<", vg_alpha = " << state->soil_properties[soil].vg_alpha_per_cm
	       <<", Ksat = "     << state->soil_properties[soil].Ksat_cm_per_h <<"\n";
	      //  <<", theta = "    << current->theta <<"\n";

      if (dual_perm){
        std::cerr<<"| soil_type = "<< soil <<", layer = "<<layer_num
          <<", smcmax = "   << state->soil_properties_frac[soil].theta_e
          <<", smcmin = "   << state->soil_properties_frac[soil].theta_r
          <<", vg_n = "     << state->soil_properties_frac[soil].vg_n
          <<", vg_alpha = " << state->soil_properties_frac[soil].vg_alpha_per_cm
          <<", Ksat = "     << state->soil_properties_frac[soil].Ksat_cm_per_h 

          <<", a_f = "      << state->mass_transfer_soil_properties[soil].a_f
          <<", beta_f = "   << state->mass_transfer_soil_properties[soil].beta_f
          <<", K_sa_f = "   << state->mass_transfer_soil_properties[soil].K_sa_f
          <<", gamma_f = "  << state->mass_transfer_soil_properties[soil].gamma_f
          
          <<"\n";
      }
    }
    
    // current = current->next;
  }

  //next we update the parameters that apply to the whole model domain and do not depend on soil layer
  if (verbosity.compare("high") == 0 || verbosity.compare("low") == 0) {
    std::cerr<<"----------- Calibratable parameters independent of soil layer (initial values) ----------- \n";
    std::cerr<<"field_capacity_psi = "   << state->lgar_bmi_params.field_capacity_psi_cm
      <<", ponded_depth_max = "     << state->lgar_bmi_params.ponded_depth_max_cm
      <<", a = "     << state->lgar_bmi_params.a
      <<", b = "     << state->lgar_bmi_params.b
      <<", frac_to_GW = "     << state->lgar_bmi_params.frac_to_GW
      <<", frac_to_pref = "     << state->lgar_bmi_params.frac_to_pref
      <<", spf_factor = "     << state->lgar_bmi_params.spf_factor
      <<", ratio_fracture_vol_to_total_vol = "     << state->lgar_bmi_params.ratio_fracture_vol_to_total_vol <<
      "\n";
  }

  state->lgar_bmi_params.field_capacity_psi_cm = state->lgar_calib_params.field_capacity_psi;
  state->lgar_bmi_params.ponded_depth_max_cm   = state->lgar_calib_params.ponded_depth_max;
  state->lgar_bmi_params.a                     = state->lgar_calib_params.a;
  state->lgar_bmi_params.b                     = state->lgar_calib_params.b;
  state->lgar_bmi_params.frac_to_GW            = state->lgar_calib_params.frac_to_GW;
  state->lgar_bmi_params.spf_factor            = state->lgar_calib_params.spf_factor;
  state->lgar_bmi_params.frac_to_pref          = state->lgar_calib_params.frac_to_pref;
  state->lgar_bmi_params.ratio_fracture_vol_to_total_vol          = state->lgar_calib_params.ratio_fracture_vol_to_total_vol;

  if (verbosity.compare("high") == 0 || verbosity.compare("low") == 0) {
    std::cerr<<"----------- Calibratable parameters independent of soil layer (updated values) ----------- \n";
    std::cerr<<"field_capacity_psi = "   << state->lgar_bmi_params.field_capacity_psi_cm
      <<", ponded_depth_max = "     << state->lgar_bmi_params.ponded_depth_max_cm
      <<", a = "     << state->lgar_bmi_params.a
      <<", b = "     << state->lgar_bmi_params.b
      <<", frac_to_GW = "     << state->lgar_bmi_params.frac_to_GW
      <<", frac_to_pref = "     << state->lgar_bmi_params.frac_to_pref
      <<", spf_factor = "     << state->lgar_bmi_params.spf_factor
      <<", ratio_fracture_vol_to_total_vol = "     << state->lgar_bmi_params.ratio_fracture_vol_to_total_vol <<
      "\n";
  }
  
  if (verbosity.compare("high") == 0)
    listPrint(state->head);
  
  double volstart_after = (1 - ratio_fracture_vol_to_total_vol)*lgar_calc_mass_bal(state->lgar_bmi_params.cum_layer_thickness_cm, state->head);

  if (dual_perm){
    volstart_after += ratio_fracture_vol_to_total_vol*lgar_calc_mass_bal(state->lgar_bmi_params.cum_layer_thickness_cm, state->head_frac);
  }

  if (verbosity.compare("high") == 0 || verbosity.compare("low") == 0)
    std::cerr<<"Mass of water (before and after) = "<< volstart_before<<", "<< volstart_after <<"\n";
  
  return volstart_after - volstart_before;
}

void BmiLGAR::
Finalize()
{
  global_mass_balance();
  listDelete(state->head);
  listDelete(state->state_previous);

  listDelete(state->head_frac);
  listDelete(state->state_previous_frac);

  delete [] state->soil_properties;
  delete [] state->soil_properties_frac;
  delete [] state->mass_transfer_soil_properties;

  delete [] state->lgar_bmi_params.soil_depth_wetting_fronts;
  delete [] state->lgar_bmi_params.soil_moisture_wetting_fronts;

  delete [] state->lgar_bmi_params.soil_temperature;
  delete [] state->lgar_bmi_params.soil_temperature_z;
  delete [] state->lgar_bmi_params.layer_soil_type;

  delete [] state->lgar_calib_params.theta_e;
  delete [] state->lgar_calib_params.theta_r;
  delete [] state->lgar_calib_params.vg_n;
  delete [] state->lgar_calib_params.vg_alpha;
  delete [] state->lgar_calib_params.Ksat;

  delete [] state->lgar_calib_params.theta_e_f;
  delete [] state->lgar_calib_params.theta_r_f;
  delete [] state->lgar_calib_params.vg_n_f;
  delete [] state->lgar_calib_params.vg_alpha_f;
  delete [] state->lgar_calib_params.Ksat_f;

  delete [] state->lgar_bmi_params.layer_thickness_cm;
  delete [] state->lgar_bmi_params.cum_layer_thickness_cm;
  delete [] state->lgar_bmi_params.giuh_ordinates;
  delete [] state->lgar_bmi_params.frozen_factor;
  delete state->lgar_bmi_input_params;
  delete state;
}


int BmiLGAR::
GetVarGrid(std::string name)
{
  if (name.compare("soil_storage_model") == 0 || name.compare("soil_num_wetting_fronts") == 0)   // int
    return 0;
  else if (name.compare("precipitation_rate") == 0 || name.compare("precipitation") == 0)
    return 1;
  else if (name.compare("potential_evapotranspiration_rate") == 0
	   || name.compare("potential_evapotranspiration") == 0
	   || name.compare("actual_evapotranspiration") == 0) // double
    return 1;
  else if (name.compare("surface_runoff") == 0 || name.compare("giuh_runoff") == 0 || name.compare("surface_runoff_frac") == 0 || name.compare("ratio_fracture_vol_to_total_vol") == 0
	   || name.compare("soil_storage") == 0 || name.compare("soil_storage_frac") == 0 || name.compare("field_capacity") == 0 || name.compare("a") == 0 || name.compare("b") == 0 || name.compare("frac_to_GW") == 0 || name.compare("spf_factor") == 0 || name.compare("frac_to_pref") == 0 || name.compare("ponded_depth_max") == 0)// double
    return 1;
  else if (name.compare("total_discharge") == 0 || name.compare("infiltration") == 0 || name.compare("infiltration_frac") == 0 || name.compare("initial_psi_cm") == 0
	   || name.compare("percolation") == 0 || name.compare("percolation_frac") == 0 || name.compare("groundwater_to_stream_recharge") == 0) // double
    return 1;
  else if (name.compare("mass_balance") == 0 || name.compare("mass_transfer") == 0)
    return 1;
  else if (name.compare("soil_depth_layers") == 0  || name.compare("smcmax") == 0 || name.compare("smcmax_f") == 0 || name.compare("smcmin") == 0 || name.compare("smcmin_f") == 0
	   || name.compare("van_genuchten_m") == 0 || name.compare("van_genuchten_m_f") == 0 || name.compare("van_genuchten_alpha") == 0 || name.compare("van_genuchten_alpha_f") == 0 || name.compare("van_genuchten_n") == 0 || name.compare("van_genuchten_n_f") == 0 
	   || name.compare("hydraulic_conductivity") == 0 || name.compare("hydraulic_conductivity_f") == 0 
     || name.compare("a_f") == 0 || name.compare("beta_f") == 0 || name.compare("K_sa_f") == 0 || name.compare("gamma_f") == 0)// array of doubles (fixed length) 
    return 2;
  else if (name.compare("soil_moisture_wetting_fronts") == 0 || name.compare("soil_depth_wetting_fronts") == 0) // array of doubles (dynamic length)
    return 3;
  else if (name.compare("soil_temperature_profile") == 0) // array of doubles (fixed and of the size of soil temperature profile)
    return 4;
  else
    return -1;
}


std::string BmiLGAR::
GetVarType(std::string name)
{
  int var_grid = GetVarGrid(name);

  if (var_grid == 0)
    return "int";
  else if (var_grid == 1 || var_grid == 2 || var_grid == 3 || var_grid == 4)
    return "double";
  else
    return "none";
}


int BmiLGAR::
GetVarItemsize(std::string name)
{
  int var_grid = GetVarGrid(name);

   if (var_grid == 0)
    return sizeof(int);
  else if (var_grid == 1 || var_grid == 2 || var_grid == 3 || var_grid == 4)
    return sizeof(double);
  else
    return 0;
}


std::string BmiLGAR::
GetVarUnits(std::string name)
{
  if (name.compare("precipitation_rate") == 0 || name.compare("potential_evapotranspiration_rate") == 0)
    return "mm h^-1";
  else if (name.compare("precipitation") == 0 || name.compare("potential_evapotranspiration") == 0
	   || name.compare("actual_evapotranspiration") == 0) // double
    return "m";
  else if (name.compare("surface_runoff") == 0 || name.compare("giuh_runoff") == 0 || name.compare("surface_runoff_frac") == 0
	   || name.compare("soil_storage") == 0 || name.compare("soil_storage_frac") == 0) // double
    return "m";
  else if (name.compare("total_discharge") == 0 || name.compare("infiltration") == 0 || name.compare("infiltration_frac") == 0
	   || name.compare("percolation") == 0 || name.compare("percolation_frac") == 0) // double
    return "m";
  else if (name.compare("mass_balance") == 0 || name.compare("groundwater_to_stream_recharge") == 0 || name.compare("mass_transfer") == 0)
    return "m";
  else if (name.compare("soil_moisture_wetting_fronts") == 0) // array of doubles
    return "none";
  else if (name.compare("soil_depth_layers") == 0 || name.compare("soil_depth_wetting_fronts") == 0) // array of doubles
    return "m";
  else if (name.compare("soil_temperature_profile") == 0)
    return "K";
  else
    return "none";

}


int BmiLGAR::
GetVarNbytes(std::string name)
{
  int itemsize;
  int gridsize;

  itemsize = this->GetVarItemsize(name);
  gridsize = this->GetGridSize(this->GetVarGrid(name));
  return itemsize * gridsize;
}


std::string BmiLGAR::
GetVarLocation(std::string name)
{
  if (name.compare("precipitation_rate") == 0 || name.compare("precipitation") == 0 || name.compare("initial_psi_cm") == 0 ||
      name.compare("potential_evapotranspiration") == 0 || name.compare("potential_evapotranspiration_rate") == 0
      || name.compare("actual_evapotranspiration") == 0) // double
    return "node";
  else if (name.compare("surface_runoff") == 0 || name.compare("surface_runoff_frac") == 0 || name.compare("giuh_runoff") == 0
	   || name.compare("soil_storage") == 0 || name.compare("soil_storage_frac") == 0)  // double
    return "node";
   else if (name.compare("total_discharge") == 0 || name.compare("infiltration") == 0 || name.compare("infiltration_frac") == 0
	    || name.compare("percolation") == 0 || name.compare("percolation_frac") == 0 || name.compare("groundwater_to_stream_recharge") == 0) // double
    return "node";
  else if (name.compare("soil_moisture_wetting_fronts") == 0) // array of doubles
    return "node";
  else if (name.compare("mass_balance") == 0 || name.compare("mass_transfer") == 0)
    return "node";
  else if (name.compare("soil_depth_layers") == 0 || name.compare("soil_depth_wetting_fronts") == 0
	   || name.compare("soil_num_wetting_fronts") == 0) // array of doubles
    return "node";
  else if (name.compare("soil_temperature_profile") == 0)
    return "node";
  else
    return "none";
}


void BmiLGAR::
GetGridShape(const int grid, int *shape)
{
  if (grid == 2)
    shape[0] = this->state->lgar_bmi_params.num_layers;
  else if (grid == 3) // number of wetting fronts (dynamic)
    shape[1] = this->state->lgar_bmi_params.num_wetting_fronts;
}


void BmiLGAR::
GetGridSpacing (const int grid, double * spacing)
{
  if (grid == 0) {
    spacing[0] = this->state->lgar_bmi_params.spacing[0];
  }
}


void BmiLGAR::
GetGridOrigin (const int grid, double *origin)
{
  if (grid == 0) {
    origin[0] = this->state->lgar_bmi_params.origin[0];
  }
}


int BmiLGAR::
GetGridRank(const int grid)
{
  if (grid == 0 || grid == 1 || grid == 2 || grid == 3 || grid == 4)
    return 1;
  else
    return -1;
}


int BmiLGAR::
GetGridSize(const int grid)
{
  if (grid == 0 || grid == 1)
    return 1;
  else if (grid == 2) // number of layers (fixed)
    return this->state->lgar_bmi_params.num_layers;
  else if (grid == 3) // number of wetting fronts (dynamic)
    return this->state->lgar_bmi_params.num_wetting_fronts;
  else if (grid == 4) // number of cells (discretized temperature profile, input from SFT)
    return this->state->lgar_bmi_params.num_cells_temp;
  else
    return -1;
}



void BmiLGAR::
GetValue (std::string name, void *dest)
{
  void * src = NULL;
  int nbytes = 0;

  src = this->GetValuePtr(name);
  nbytes = this->GetVarNbytes(name);
  memcpy (dest, src, nbytes);
}


void *BmiLGAR::
GetValuePtr (std::string name)
{
  if (name.compare("precipitation_rate") == 0)
    return (void*)(&this->state->lgar_bmi_input_params->precipitation_mm_per_h);
  else if (name.compare("precipitation") == 0)
    return (void*)(&bmi_unit_conv.volprecip_timestep_m);
  else if (name.compare("potential_evapotranspiration_rate") == 0)
    return (void*)(&this->state->lgar_bmi_input_params->PET_mm_per_h);
  else if (name.compare("potential_evapotranspiration") == 0)
    return (void*)(&bmi_unit_conv.volPET_timestep_m);
  else if (name.compare("actual_evapotranspiration") == 0)
    return (void*)(&bmi_unit_conv.volAET_timestep_m);
  else if (name.compare("surface_runoff") == 0)
    return (void*)(&bmi_unit_conv.volrunoff_timestep_m);
  else if (name.compare("surface_runoff_frac") == 0)
    return (void*)(&bmi_unit_conv.volrunoff_timestep_m_frac);
  else if (name.compare("giuh_runoff") == 0)
    return (void*)(&bmi_unit_conv.volrunoff_giuh_timestep_m);
  else if (name.compare("soil_storage") == 0)
    return (void*)(&bmi_unit_conv.volend_timestep_m);
  else if (name.compare("soil_storage_frac") == 0)
    return (void*)(&bmi_unit_conv.volend_timestep_m_frac);
  else if (name.compare("total_discharge") == 0)
    return (void*)(&bmi_unit_conv.volQ_timestep_m);
  else if (name.compare("infiltration") == 0)
    return (void*)(&bmi_unit_conv.volin_timestep_m);
  else if (name.compare("infiltration_frac") == 0)
    return (void*)(&bmi_unit_conv.volin_timestep_m_frac);
  else if (name.compare("percolation") == 0)
    return (void*)(&bmi_unit_conv.volrech_timestep_m);
  else if (name.compare("percolation_frac") == 0)
    return (void*)(&bmi_unit_conv.volrech_timestep_m_frac);
  else if (name.compare("groundwater_to_stream_recharge") == 0)
    return (void*)(&bmi_unit_conv.volQ_gw_timestep_m);
  else if (name.compare("mass_balance") == 0)
    return (void*)(&bmi_unit_conv.mass_balance_m);
  else if (name.compare("mass_transfer") == 0)
    return (void*)(&bmi_unit_conv.mass_transfer_timestep_m);
  else if (name.compare("soil_depth_layers") == 0)
    return (void*)this->state->lgar_bmi_params.cum_layer_thickness_cm;  // this too and, if needed, change soil_moisture_layers to soil_thickness_layers
  else if (name.compare("soil_moisture_wetting_fronts") == 0)
    return (void*)this->state->lgar_bmi_params.soil_moisture_wetting_fronts;
  else if (name.compare("soil_depth_wetting_fronts") == 0)
    return (void*)this->state->lgar_bmi_params.soil_depth_wetting_fronts;
  else if (name.compare("soil_num_wetting_fronts") == 0)
    return (void*)(&state->lgar_bmi_params.num_wetting_fronts);
  else if (name.compare("initial_psi_cm") == 0)
    return (void*)(&state->lgar_bmi_params.initial_psi_cm);
  else if (name.compare("soil_temperature_profile") == 0)
    return (void*)this->state->lgar_bmi_params.soil_temperature;
  else if (name.compare("smcmax") == 0)
    return (void*)this->state->lgar_calib_params.theta_e;
  else if (name.compare("smcmin") == 0)
    return (void*)this->state->lgar_calib_params.theta_r;
  else if (name.compare("van_genuchten_n") == 0)
    return (void*)this->state->lgar_calib_params.vg_n;
  else if (name.compare("van_genuchten_alpha") == 0)
    return (void*)this->state->lgar_calib_params.vg_alpha;
  else if (name.compare("hydraulic_conductivity") == 0)
    return (void*)this->state->lgar_calib_params.Ksat;
  else if (name.compare("smcmax_f") == 0)
    return (void*)this->state->lgar_calib_params.theta_e_f;
  else if (name.compare("smcmin_f") == 0)
    return (void*)this->state->lgar_calib_params.theta_r_f;
  else if (name.compare("van_genuchten_n_f") == 0)
    return (void*)this->state->lgar_calib_params.vg_n_f;
  else if (name.compare("van_genuchten_alpha_f") == 0)
    return (void*)this->state->lgar_calib_params.vg_alpha_f;
  else if (name.compare("hydraulic_conductivity_f") == 0)
      return (void*)this->state->lgar_calib_params.Ksat_f;
  else if (name.compare("a_f") == 0)
    return (void*)this->state->lgar_calib_params.a_f;
  else if (name.compare("beta_f") == 0)
    return (void*)this->state->lgar_calib_params.beta_f;
  else if (name.compare("K_sa_f") == 0)
    return (void*)this->state->lgar_calib_params.K_sa_f;
  else if (name.compare("gamma_f") == 0)
    return (void*)this->state->lgar_calib_params.gamma_f;
  else if (name.compare("ponded_depth_max") == 0)
    return (void*)&this->state->lgar_calib_params.ponded_depth_max;
  else if (name.compare("field_capacity") == 0)
    return (void*)&this->state->lgar_calib_params.field_capacity_psi;
  else if (name.compare("a") == 0)
    return (void*)&this->state->lgar_calib_params.a;
  else if (name.compare("b") == 0)
    return (void*)&this->state->lgar_calib_params.b;
  else if (name.compare("frac_to_GW") == 0)
    return (void*)&this->state->lgar_calib_params.frac_to_GW;
    else if (name.compare("spf_factor") == 0)
    return (void*)&this->state->lgar_calib_params.spf_factor;
  else if (name.compare("frac_to_pref") == 0)
    return (void*)&this->state->lgar_calib_params.frac_to_pref;
  else if (name.compare("ratio_fracture_vol_to_total_vol") == 0)
    return (void*)&this->state->lgar_calib_params.ratio_fracture_vol_to_total_vol;
  else {
    std::stringstream errMsg;
    errMsg << "variable "<< name << " does not exist";
    throw std::runtime_error(errMsg.str());
    return NULL;
  }
  
  // delete it later
  return NULL;
}

void BmiLGAR::
GetValueAtIndices (std::string name, void *dest, int *inds, int len)
{
  void * src = NULL;

  src = this->GetValuePtr(name);

  if (src) {
    int i;
    int itemsize = 0;
    int offset;
    char *ptr;

    itemsize = this->GetVarItemsize(name);

    for (i=0, ptr=(char *)dest; i<len; i++, ptr+=itemsize) {
      offset = inds[i] * itemsize;
      memcpy(ptr, (char *)src + offset, itemsize);
    }
  }
}


void BmiLGAR::
SetValue (std::string name, void *src)
{
  void * dest = NULL;
  dest = this->GetValuePtr(name);

  if (dest) {
    int nbytes = 0;
    nbytes = this->GetVarNbytes(name);
    memcpy(dest, src, nbytes);
  }

}


void BmiLGAR::
SetValueAtIndices (std::string name, int * inds, int len, void *src)
{
  void * dest = NULL;

  dest = this->GetValuePtr(name);

  if (dest) {
    int i;
    int itemsize = 0;
    int offset;
    char *ptr;

    itemsize = this->GetVarItemsize(name);

    for (i=0, ptr=(char *)src; i<len; i++, ptr+=itemsize) {
      offset = inds[i] * itemsize;
      memcpy((char *)dest + offset, ptr, itemsize);
    }
  }
}


std::string BmiLGAR::
GetComponentName()
{
  return "LASAM (Lumped Arid/Semi-arid Model)";
}


int BmiLGAR::
GetInputItemCount()
{
  return this->input_var_name_count;
}


int BmiLGAR::
GetOutputItemCount()
{
  return this->output_var_name_count;
}


std::vector<std::string> BmiLGAR::
GetInputVarNames()
{
  std::vector<std::string> names;

  for (int i=0; i<this->input_var_name_count; i++)
    names.push_back(this->input_var_names[i]);

  return names;
}


std::vector<std::string> BmiLGAR::
GetOutputVarNames()
{
  std::vector<std::string> names;

  for (int i=0; i<this->output_var_name_count; i++)
    names.push_back(this->output_var_names[i]);

  return names;
}


double BmiLGAR::
GetStartTime () {
  return 0.0;
}


double BmiLGAR::
GetEndTime () {
  return this->state->lgar_bmi_params.endtime_s;
}


double BmiLGAR::
GetCurrentTime () {
  return this->state->lgar_bmi_params.time_s;
}


std::string BmiLGAR::
GetTimeUnits() {
  return "s";
}


double BmiLGAR::
GetTimeStep () {
  return this->state->lgar_bmi_params.forcing_resolution_h * 3600.; // convert hours to seconds
}

std::string BmiLGAR::
GetGridType(const int grid)
{
  if (grid == 0)
    return "uniform_rectilinear";
  else
    return "";
}


void BmiLGAR::
GetGridX(const int grid, double *x)
{
  // this is not needed but printing here to avoid compiler warnings
  std::cerr<<"GetGridX: "<<grid<<" "<<x[0]<<"\n";
  throw bmi_lgar::NotImplemented();
}


void BmiLGAR::
GetGridY(const int grid, double *y)
{
  // this is not needed but printing here to avoid compiler warnings
  std::cerr<<"GetGridY: "<<grid<<" "<<y[0]<<"\n";
  throw bmi_lgar::NotImplemented();
}


void BmiLGAR::
GetGridZ(const int grid, double *z)
{
  // this is not needed but printing here to avoid compiler warnings
  std::cerr<<"GetGridZ: "<<grid<<" "<<z[0]<<"\n";
  throw bmi_lgar::NotImplemented();
}


int BmiLGAR::
GetGridNodeCount(const int grid)
{
  // this is not needed but printing here to avoid compiler warnings
  std::cerr<<"GetGridNodeCount: "<<grid<<"\n";
  throw bmi_lgar::NotImplemented();
}


int BmiLGAR::
GetGridEdgeCount(const int grid)
{
  // this is not needed but printing here to avoid compiler warnings
  std::cerr<<"GetGridEdgeCount: "<<grid<<"\n";
  throw bmi_lgar::NotImplemented();
}


int BmiLGAR::
GetGridFaceCount(const int grid)
{
  // this is not needed but printing here to avoid compiler warnings
  std::cerr<<"GetGridFaceCount: "<<grid<<"\n";
  throw bmi_lgar::NotImplemented();
}


void BmiLGAR::
GetGridEdgeNodes(const int grid, int *edge_nodes)
{
  // this is not needed but printing here to avoid compiler warnings
  std::cerr<<"GetGridEdgeNodes: "<<grid<<" "<<edge_nodes[0]<<"\n";
  throw bmi_lgar::NotImplemented();
}


void BmiLGAR::
GetGridFaceEdges(const int grid, int *face_edges)
{
  // this is not needed but printing here to avoid compiler warnings
  std::cerr<<"GetGridFaceNodes: "<<grid<<" "<<face_edges[0]<<"\n";
  throw bmi_lgar::NotImplemented();
}


void BmiLGAR::
GetGridFaceNodes(const int grid, int *face_nodes)
{
  // this is not needed but printing here to avoid compiler warnings
  std::cerr<<"GetGridFaceNodes: "<<grid<<" "<<face_nodes[0]<<"\n";
  throw bmi_lgar::NotImplemented();
}


void BmiLGAR::
GetGridNodesPerFace(const int grid, int *nodes_per_face)
{
  // this is not needed but printing here to avoid compiler warnings
  std::cerr<<"GetGridNodesPerFace: "<<grid<<" "<<nodes_per_face[0]<<"\n";
  throw bmi_lgar::NotImplemented();
}

#endif
