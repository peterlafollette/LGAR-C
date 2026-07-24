#ifndef BMI_LGAR_CXX_INCLUDED
#define BMI_LGAR_CXX_INCLUDED


#include <stdio.h>
#include <string>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <iostream>
#include "../bmi/bmi.hxx"
#include "../include/bmi_lgar.hxx"
#include "../include/all.hxx"


// default verbosity is set to 'none' other option 'high' or 'low' needs to be specified in the config file
string verbosity="none";

// some macros for flux caching, these are not relevant if flux caching is disabled in the config file.
// precip intensity above which flux caching does not occur
#ifndef PRECIP_THRESHOLD_MM_PER_H
#define PRECIP_THRESHOLD_MM_PER_H ((double)1.0e-6)
#endif

// Small epsilon added to PET denominator to prevent division by 0, and in another case to check if there is significant PET
#ifndef PET_EPSILON
#define PET_EPSILON ((double)1.0e-8)
#endif

// Threshold for AET/PET ratio to indicate low enough AET for flux caching to start
#ifndef AET_PET_RATIO_THRESHOLD
#define AET_PET_RATIO_THRESHOLD ((double)0.75)
#endif

// while ponded head is greater than this, flux caching will not occur
#ifndef VOLON_TIMESTEP_THRESHOLD_CM
#define VOLON_TIMESTEP_THRESHOLD_CM ((double)1.0e-16)
#endif

// While the vadose zone lower boundary flux is greater than this, flux caching will not happen
#ifndef BOTTOM_BDY_FLUX_THRESHOLD_CM
#define BOTTOM_BDY_FLUX_THRESHOLD_CM ((double)1.0e-2)
#endif

// While the top most wetting front has a dzdt value greater than this, flux caching will not happen
#ifndef THRESHOLD_DZDT_CM_PER_H
#define THRESHOLD_DZDT_CM_PER_H ((double)1.0e0)
#endif

// Number of timesteps before flux cache reset
#ifndef NUM_TIMESTEPS_BEFORE_RESET_CACHE
#define NUM_TIMESTEPS_BEFORE_RESET_CACHE 24
#endif

// Fraction of the local mass-balance tolerance that cached lower-boundary
// correction is allowed to accumulate before forcing a recomputed step.
#ifndef CACHE_LOWER_BOUNDARY_MBAL_FRACTION
#define CACHE_LOWER_BOUNDARY_MBAL_FRACTION ((double)0.5)
#endif

// small epsillon that is used to determine if the difference between two quantities is 0 while avoiding machine precision errors
#define SMALL_EPS 1.E-12

#ifndef LGARTO_SATURATED_CREATION_PRESSURE_GRADIENT_CAP
#define LGARTO_SATURATED_CREATION_PRESSURE_GRADIENT_CAP ((double)1.0)
#endif

#ifndef LGARTO_SATURATED_CREATION_PSI_TOL_CM
#define LGARTO_SATURATED_CREATION_PSI_TOL_CM ((double)1.0e-6)
#endif

#ifndef LGARTO_SATURATED_CREATION_THETA_TOL
#define LGARTO_SATURATED_CREATION_THETA_TOL ((double)1.0e-8)
#endif

#ifndef LGARTO_EVENT_SPLIT_BOUNDARY_EPS_CM
#define LGARTO_EVENT_SPLIT_BOUNDARY_EPS_CM ((double)1.0e-9)
#endif

#ifndef LGARTO_EVENT_SPLIT_MAX_PER_FORCING
#define LGARTO_EVENT_SPLIT_MAX_PER_FORCING 1000
#endif

#ifndef LGARTO_TO_PACKET_EVENT_SPLIT_MIN_SPACING_CM
#define LGARTO_TO_PACKET_EVENT_SPLIT_MIN_SPACING_CM ((double)1.0e-5)
#endif

#ifndef LGARTO_TO_PACKET_EVENT_SPLIT_MERGE_READY_SPACING_CM
// Hand near-contact mobile TO/GW packets to the normal merge/correction pass.
#define LGARTO_TO_PACKET_EVENT_SPLIT_MERGE_READY_SPACING_CM ((double)1.0e-3)
#endif

#ifndef LGARTO_TO_PACKET_EVENT_SPLIT_REPEAT_HANDOFF_COUNT
#define LGARTO_TO_PACKET_EVENT_SPLIT_REPEAT_HANDOFF_COUNT 8
#endif

#ifndef LGARTO_TO_PACKET_EVENT_SPLIT_REPEAT_HANDOFF_MAX_GAP_CM
#define LGARTO_TO_PACKET_EVENT_SPLIT_REPEAT_HANDOFF_MAX_GAP_CM ((double)2.5e-1)
#endif

#ifndef LGARTO_TO_PACKET_EVENT_SPLIT_HANDOFF_OVERTAKE_CM
#define LGARTO_TO_PACKET_EVENT_SPLIT_HANDOFF_OVERTAKE_CM ((double)2.5e-1)
#endif

#ifndef LGARTO_TO_PACKET_EVENT_SPLIT_MIN_DT_H
#define LGARTO_TO_PACKET_EVENT_SPLIT_MIN_DT_H ((double)1.0e-8)
#endif

#ifndef LGARTO_TO_PACKET_INTERFLOW_EVENT_BISECTION_ITERATIONS
#define LGARTO_TO_PACKET_INTERFLOW_EVENT_BISECTION_ITERATIONS 60
#endif

#ifndef LGARTO_MOBILE_GROUNDWATER_STORAGE_COEFF_MIN
#define LGARTO_MOBILE_GROUNDWATER_STORAGE_COEFF_MIN ((double)0.03)
#endif

#ifndef LGARTO_MOBILE_GROUNDWATER_STORAGE_COEFF_FALLBACK
#define LGARTO_MOBILE_GROUNDWATER_STORAGE_COEFF_FALLBACK ((double)1.0)
#endif

#ifndef LGARTO_MOBILE_GW_RECESSION_REWET_START_MAX_DEPTH_CM
#define LGARTO_MOBILE_GW_RECESSION_REWET_START_MAX_DEPTH_CM ((double)5.0)
#endif

static double lgar_fixed_soil_depth_cm(const lgar_bmi_parameters *params)
{
  if (params == NULL) {
    return 0.0;
  }

  if (params->num_layers > 0 && params->cum_layer_thickness_cm != NULL) {
    return params->cum_layer_thickness_cm[params->num_layers];
  }

  return params->soil_depth_cm;
}

static double lgar_effective_groundwater_depth_cm(const lgar_bmi_parameters *params)
{
  const double fixed_depth_cm = lgar_fixed_soil_depth_cm(params);
  if (params == NULL || !params->mobile_groundwater_level) {
    return fixed_depth_cm;
  }

  if (!std::isfinite(params->groundwater_depth_cm) || params->groundwater_depth_cm < 0.0) {
    return fixed_depth_cm;
  }

  return params->groundwater_depth_cm;
}

static int lgar_mobile_groundwater_layer_num(const lgar_bmi_parameters *params,
                                             double groundwater_depth_cm)
{
  if (params == NULL || params->num_layers <= 0 ||
      params->cum_layer_thickness_cm == NULL ||
      !std::isfinite(groundwater_depth_cm)) {
    return 0;
  }

  for (int layer_num = 1; layer_num <= params->num_layers; layer_num++) {
    if (groundwater_depth_cm <=
        params->cum_layer_thickness_cm[layer_num] + LGARTO_EVENT_SPLIT_BOUNDARY_EPS_CM) {
      return layer_num;
    }
  }

  return params->num_layers;
}

static double lgar_mobile_groundwater_storage_coefficient(const model_state *state,
                                                          double groundwater_depth_cm)
{
  if (state == NULL || state->lgar_bmi_params.layer_soil_type == NULL ||
      state->lgar_bmi_params.cum_layer_thickness_cm == NULL ||
      state->soil_properties == NULL) {
    return LGARTO_MOBILE_GROUNDWATER_STORAGE_COEFF_FALLBACK;
  }

  const lgar_bmi_parameters *params = &state->lgar_bmi_params;
  const int groundwater_layer =
    lgar_mobile_groundwater_layer_num(params, groundwater_depth_cm);
  if (groundwater_layer < 1 || groundwater_layer > params->num_layers) {
    return LGARTO_MOBILE_GROUNDWATER_STORAGE_COEFF_FALLBACK;
  }

  const int soil_num = params->layer_soil_type[groundwater_layer];
  const double theta_e = state->soil_properties[soil_num].theta_e;
  if (!std::isfinite(theta_e) || theta_e <= 0.0) {
    return LGARTO_MOBILE_GROUNDWATER_STORAGE_COEFF_FALLBACK;
  }

  const double fixed_column_depth_cm = lgar_fixed_soil_depth_cm(params);
  const double search_depth_cm =
    std::isfinite(fixed_column_depth_cm) && fixed_column_depth_cm > 0.0
      ? fmin(groundwater_depth_cm, fixed_column_depth_cm)
      : groundwater_depth_cm;
  const double depth_tol_cm =
    fmax(1.0e-8, 1.0e-10 * fmax(1.0, fmax(search_depth_cm, fixed_column_depth_cm)));

  bool found_unsaturated_theta = false;
  double deepest_unsaturated_depth_cm = -HUGE_VAL;
  double deepest_unsaturated_theta = theta_e;
  for (const wetting_front *current = state->head; current != NULL; current = current->next) {
    if (current->layer_num != groundwater_layer ||
        !std::isfinite(current->depth_cm) ||
        current->depth_cm > search_depth_cm + depth_tol_cm ||
        !std::isfinite(current->theta)) {
      continue;
    }

    if (current->theta >= theta_e - LGARTO_SATURATED_CREATION_THETA_TOL) {
      continue;
    }

    if (!found_unsaturated_theta ||
        current->depth_cm > deepest_unsaturated_depth_cm) {
      found_unsaturated_theta = true;
      deepest_unsaturated_depth_cm = current->depth_cm;
      deepest_unsaturated_theta = current->theta;
    }
  }

  const double theta_above_groundwater =
    found_unsaturated_theta ? deepest_unsaturated_theta : theta_e;
  double storage_coefficient = theta_e - theta_above_groundwater;
  const double storage_coefficient_min =
    fmax(1.0e-12, LGARTO_MOBILE_GROUNDWATER_STORAGE_COEFF_MIN);
  const double storage_coefficient_max = fmax(storage_coefficient_min, theta_e);
  if (!std::isfinite(storage_coefficient)) {
    storage_coefficient = LGARTO_MOBILE_GROUNDWATER_STORAGE_COEFF_FALLBACK;
  }

  return fmax(storage_coefficient_min,
              fmin(storage_coefficient, storage_coefficient_max));
}

struct lgar_mobile_groundwater_profile_point {
  double depth_cm;
  double theta;
  int order;
};

static bool lgar_mobile_groundwater_depth_before(
    const lgar_mobile_groundwater_profile_point &a,
    const lgar_mobile_groundwater_profile_point &b)
{
  if (a.depth_cm < b.depth_cm) {
    return true;
  }
  if (a.depth_cm > b.depth_cm) {
    return false;
  }
  return a.order < b.order;
}

static double lgar_profile_theta_at_depth_in_layer(const model_state *state,
                                                   int layer_num,
                                                   double depth_cm)
{
  if (state == NULL || state->head == NULL || layer_num < 1 ||
      layer_num > state->lgar_bmi_params.num_layers ||
      !std::isfinite(depth_cm)) {
    return NAN;
  }

  std::vector<lgar_mobile_groundwater_profile_point> points;
  int order = 0;
  for (const wetting_front *current = state->head; current != NULL;
       current = current->next, order++) {
    if (current->layer_num != layer_num ||
        !std::isfinite(current->depth_cm) ||
        !std::isfinite(current->theta)) {
      continue;
    }
    points.push_back({current->depth_cm, current->theta, order});
  }

  if (points.empty()) {
    return NAN;
  }

  std::stable_sort(points.begin(), points.end(),
                   lgar_mobile_groundwater_depth_before);

  const double depth_tol_cm =
    fmax(1.0e-8, 1.0e-10 * fmax(1.0, fabs(depth_cm)));
  for (const auto &point : points) {
    if (depth_cm <= point.depth_cm + depth_tol_cm) {
      return point.theta;
    }
  }

  // Below the deepest represented front, continue the deepest profile segment.
  return points.back().theta;
}

static double lgar_mobile_groundwater_interval_storage_coefficient(
    const model_state *state,
    double interval_midpoint_depth_cm)
{
  if (state == NULL || state->lgar_bmi_params.layer_soil_type == NULL ||
      state->soil_properties == NULL ||
      !std::isfinite(interval_midpoint_depth_cm)) {
    return LGARTO_MOBILE_GROUNDWATER_STORAGE_COEFF_FALLBACK;
  }

  const lgar_bmi_parameters *params = &state->lgar_bmi_params;
  const int layer_num =
    lgar_mobile_groundwater_layer_num(params,
                                      fmax(0.0, interval_midpoint_depth_cm));
  if (layer_num < 1 || layer_num > params->num_layers) {
    return LGARTO_MOBILE_GROUNDWATER_STORAGE_COEFF_FALLBACK;
  }

  const int soil_num = params->layer_soil_type[layer_num];
  const double theta_e = state->soil_properties[soil_num].theta_e;
  if (!std::isfinite(theta_e) || theta_e <= 0.0) {
    return LGARTO_MOBILE_GROUNDWATER_STORAGE_COEFF_FALLBACK;
  }

  const double profile_theta =
    lgar_profile_theta_at_depth_in_layer(state, layer_num,
                                         interval_midpoint_depth_cm);
  if (!std::isfinite(profile_theta)) {
    return lgar_mobile_groundwater_storage_coefficient(
      state, interval_midpoint_depth_cm);
  }

  double storage_coefficient = theta_e - profile_theta;
  const double storage_coefficient_min =
    fmax(1.0e-12, LGARTO_MOBILE_GROUNDWATER_STORAGE_COEFF_MIN);
  const double storage_coefficient_max = fmax(storage_coefficient_min, theta_e);
  if (!std::isfinite(storage_coefficient)) {
    storage_coefficient = LGARTO_MOBILE_GROUNDWATER_STORAGE_COEFF_FALLBACK;
  }

  return fmax(storage_coefficient_min,
              fmin(storage_coefficient, storage_coefficient_max));
}

static bool lgar_mobile_groundwater_boundary_exists(std::vector<double> &boundaries,
                                                    double candidate_cm)
{
  if (!std::isfinite(candidate_cm) || candidate_cm < 0.0) {
    return true;
  }

  const double tol_cm = fmax(1.0e-8, 1.0e-10 * fmax(1.0, candidate_cm));
  for (double boundary_cm : boundaries) {
    if (fabs(boundary_cm - candidate_cm) <= tol_cm) {
      return true;
    }
  }
  return false;
}

static std::vector<double> lgar_mobile_groundwater_storage_boundaries(
    const model_state *state)
{
  std::vector<double> boundaries;
  if (state == NULL || state->lgar_bmi_params.cum_layer_thickness_cm == NULL) {
    boundaries.push_back(0.0);
    return boundaries;
  }

  const lgar_bmi_parameters *params = &state->lgar_bmi_params;
  for (int layer = 0; layer <= params->num_layers; layer++) {
    const double boundary_cm = params->cum_layer_thickness_cm[layer];
    if (!lgar_mobile_groundwater_boundary_exists(boundaries, boundary_cm)) {
      boundaries.push_back(boundary_cm);
    }
  }

  for (const wetting_front *current = state->head; current != NULL;
       current = current->next) {
    if (!std::isfinite(current->depth_cm) || current->depth_cm < 0.0) {
      continue;
    }
    if (!lgar_mobile_groundwater_boundary_exists(boundaries,
                                                 current->depth_cm)) {
      boundaries.push_back(current->depth_cm);
    }
  }

  std::sort(boundaries.begin(), boundaries.end());
  return boundaries;
}

static double lgar_next_mobile_groundwater_storage_boundary(
    const std::vector<double> &boundaries,
    double current_depth_cm,
    bool moving_up)
{
  const double tol_cm =
    fmax(1.0e-8, 1.0e-10 * fmax(1.0, fabs(current_depth_cm)));

  if (moving_up) {
    double next_boundary_cm = 0.0;
    for (double boundary_cm : boundaries) {
      if (boundary_cm < current_depth_cm - tol_cm) {
        next_boundary_cm = fmax(next_boundary_cm, boundary_cm);
      }
      else {
        break;
      }
    }
    return fmax(0.0, next_boundary_cm);
  }

  for (double boundary_cm : boundaries) {
    if (boundary_cm > current_depth_cm + tol_cm) {
      return boundary_cm;
    }
  }

  return HUGE_VAL;
}

static double lgar_solve_mobile_groundwater_depth_from_storage_exchange(
    const model_state *state,
    double current_depth_cm,
    double net_groundwater_storage_gain_cm)
{
  if (state == NULL || !std::isfinite(current_depth_cm) ||
      !std::isfinite(net_groundwater_storage_gain_cm) ||
      fabs(net_groundwater_storage_gain_cm) <= SMALL_EPS) {
    return current_depth_cm;
  }

  double depth_cm = fmax(0.0, current_depth_cm);
  double remaining_storage_cm = fabs(net_groundwater_storage_gain_cm);
  const bool moving_up = net_groundwater_storage_gain_cm > 0.0;
  const std::vector<double> boundaries =
    lgar_mobile_groundwater_storage_boundaries(state);

  for (int iter = 0;
       iter < MAX_NUM_WETTING_FRONTS * 2 && remaining_storage_cm > SMALL_EPS;
       iter++) {
    if (moving_up && depth_cm <= 0.0) {
      depth_cm = 0.0;
      break;
    }

    const double boundary_cm =
      lgar_next_mobile_groundwater_storage_boundary(boundaries, depth_cm,
                                                    moving_up);
    double interval_thickness_cm = 0.0;
    if (moving_up) {
      interval_thickness_cm = depth_cm - boundary_cm;
    }
    else if (std::isfinite(boundary_cm)) {
      interval_thickness_cm = boundary_cm - depth_cm;
    }

    if (!moving_up && !std::isfinite(boundary_cm)) {
      const double storage_coefficient =
        lgar_mobile_groundwater_interval_storage_coefficient(state,
                                                            depth_cm + 0.5);
      if (!std::isfinite(storage_coefficient) || storage_coefficient <= 0.0) {
        break;
      }
      depth_cm += remaining_storage_cm / storage_coefficient;
      remaining_storage_cm = 0.0;
      break;
    }

    if (interval_thickness_cm <= SMALL_EPS) {
      depth_cm = moving_up ? fmax(0.0, depth_cm - 1.0e-8)
                           : depth_cm + 1.0e-8;
      continue;
    }

    const double midpoint_depth_cm =
      moving_up ? depth_cm - 0.5 * interval_thickness_cm
                : depth_cm + 0.5 * interval_thickness_cm;
    const double storage_coefficient =
      lgar_mobile_groundwater_interval_storage_coefficient(state,
                                                          midpoint_depth_cm);
    if (!std::isfinite(storage_coefficient) || storage_coefficient <= 0.0) {
      break;
    }

    const double interval_storage_capacity_cm =
      storage_coefficient * interval_thickness_cm;
    if (remaining_storage_cm <= interval_storage_capacity_cm + SMALL_EPS) {
      const double depth_change_cm = remaining_storage_cm / storage_coefficient;
      depth_cm = moving_up ? fmax(0.0, depth_cm - depth_change_cm)
                           : depth_cm + depth_change_cm;
      remaining_storage_cm = 0.0;
      break;
    }

    remaining_storage_cm -= interval_storage_capacity_cm;
    depth_cm = moving_up ? fmax(0.0, boundary_cm) : boundary_cm;
  }

  if (moving_up) {
    depth_cm = fmax(0.0, depth_cm);
  }
  return depth_cm;
}

static double lgar_total_CR_storage_cm(const model_state *state)
{
  if (state == NULL) {
    return 0.0;
  }

  return state->lgar_mass_balance.CR_fast_storage_cm +
         state->lgar_mass_balance.CR_slow_storage_cm;
}

static double lgar_extract_from_CR_storage_fast_then_slow(double demand_cm,
                                                          model_state *state)
{
  if (state == NULL || demand_cm <= 0.0) {
    return 0.0;
  }

  double remaining_demand_cm = demand_cm;
  const double fast_extraction_cm =
    fmin(fmax(state->lgar_mass_balance.CR_fast_storage_cm, 0.0),
         remaining_demand_cm);
  state->lgar_mass_balance.CR_fast_storage_cm -= fast_extraction_cm;
  remaining_demand_cm -= fast_extraction_cm;

  const double slow_extraction_cm =
    fmin(fmax(state->lgar_mass_balance.CR_slow_storage_cm, 0.0),
         remaining_demand_cm);
  state->lgar_mass_balance.CR_slow_storage_cm -= slow_extraction_cm;

  if (state->lgar_mass_balance.CR_fast_storage_cm < 1.0e-12) {
    state->lgar_mass_balance.CR_fast_storage_cm = 0.0;
  }
  if (state->lgar_mass_balance.CR_slow_storage_cm < 1.0e-12) {
    state->lgar_mass_balance.CR_slow_storage_cm = 0.0;
  }

  return fast_extraction_cm + slow_extraction_cm;
}

static double lgar_explicit_aet_lower_boundary_depth_cm(
    const lgar_bmi_parameters *params)
{
  if (params == NULL || params->root_zone_depth_cm <= 0.0) {
    return 0.0;
  }

  if (!params->TO_enabled || !params->mobile_groundwater_level) {
    return params->root_zone_depth_cm;
  }

  const double groundwater_depth_cm = lgar_effective_groundwater_depth_cm(params);
  if (!std::isfinite(groundwater_depth_cm) || groundwater_depth_cm < 0.0) {
    return params->root_zone_depth_cm;
  }

  return fmax(0.0, fmin(params->root_zone_depth_cm, groundwater_depth_cm));
}

static double lgar_groundwater_supported_aet_potential_cm(
    const lgar_bmi_parameters *params,
    double PET_timestep_cm)
{
  if (params == NULL || PET_timestep_cm <= 0.0 ||
      params->root_zone_depth_cm <= 0.0 ||
      !params->TO_enabled || !params->mobile_groundwater_level) {
    return 0.0;
  }

  const double explicit_aet_depth_cm =
    lgar_explicit_aet_lower_boundary_depth_cm(params);
  const double groundwater_supported_root_depth_cm =
    fmax(0.0, params->root_zone_depth_cm - explicit_aet_depth_cm);
  const double groundwater_supported_fraction =
    fmin(1.0, groundwater_supported_root_depth_cm / params->root_zone_depth_cm);

  return PET_timestep_cm * groundwater_supported_fraction;
}

static double lgar_mobile_groundwater_depth_from_current_CR_storage(
    const model_state *state)
{
  if (state == NULL) {
    return NAN;
  }

  double reference_depth_cm =
    state->lgar_bmi_params.mobile_groundwater_reference_depth_cm;
  if (!std::isfinite(reference_depth_cm) || reference_depth_cm < 0.0) {
    reference_depth_cm = lgar_fixed_soil_depth_cm(&state->lgar_bmi_params);
  }
  if (!std::isfinite(reference_depth_cm)) {
    reference_depth_cm =
      lgar_effective_groundwater_depth_cm(&state->lgar_bmi_params);
  }

  double reference_CR_storage_cm =
    state->lgar_bmi_params.mobile_groundwater_reference_CR_storage_cm;
  if (!std::isfinite(reference_CR_storage_cm)) {
    reference_CR_storage_cm = 0.0;
  }

  return lgar_solve_mobile_groundwater_depth_from_storage_exchange(
    state,
    reference_depth_cm,
    lgar_total_CR_storage_cm(state) - reference_CR_storage_cm);
}

static double lgar_mobile_groundwater_storage_capacity_to_surface_cm(
    const model_state *state,
    double groundwater_depth_cm)
{
  if (state == NULL || !std::isfinite(groundwater_depth_cm) ||
      groundwater_depth_cm <= 0.0) {
    return 0.0;
  }

  const std::vector<double> boundaries =
    lgar_mobile_groundwater_storage_boundaries(state);
  double depth_cm = groundwater_depth_cm;
  double capacity_cm = 0.0;

  for (int iter = 0;
       iter < MAX_NUM_WETTING_FRONTS * 2 && depth_cm > SMALL_EPS;
       iter++) {
    const double next_boundary_cm =
      lgar_next_mobile_groundwater_storage_boundary(boundaries,
                                                    depth_cm,
                                                    true);
    const double interval_thickness_cm =
      fmax(0.0, depth_cm - next_boundary_cm);
    if (interval_thickness_cm <= SMALL_EPS) {
      depth_cm = fmax(0.0, depth_cm - 1.0e-8);
      continue;
    }

    const double midpoint_depth_cm =
      depth_cm - 0.5 * interval_thickness_cm;
    const double storage_coefficient =
      lgar_mobile_groundwater_interval_storage_coefficient(state,
                                                          midpoint_depth_cm);
    if (!std::isfinite(storage_coefficient) || storage_coefficient <= 0.0) {
      break;
    }

    capacity_cm += storage_coefficient * interval_thickness_cm;
    depth_cm = fmax(0.0, next_boundary_cm);
  }

  return fmax(0.0, capacity_cm);
}

static double lgar_predict_CR_discharge_cm(
    const model_state *state,
    double subtimestep_h,
    double a_fast,
    double a_slow,
    double b_fast,
    double b_slow,
    double fast_discharge_threshold_cm,
    double slow_discharge_threshold_cm,
    double frac_slow,
    double CR_input_rate_cm_per_h)
{
  if (state == NULL || subtimestep_h <= 0.0 ||
      !std::isfinite(CR_input_rate_cm_per_h)) {
    return 0.0;
  }

  double fast_storage_cm = state->lgar_mass_balance.CR_fast_storage_cm;
  double slow_storage_cm = state->lgar_mass_balance.CR_slow_storage_cm;
  return calc_CR_Q(subtimestep_h,
                   a_fast, a_slow,
                   b_fast, b_slow,
                   fast_discharge_threshold_cm,
                   slow_discharge_threshold_cm,
                   frac_slow,
                   CR_input_rate_cm_per_h,
                   &fast_storage_cm,
                   &slow_storage_cm);
}

static double lgar_cap_positive_CR_inputs_by_mobile_groundwater_surface(
    model_state *state,
    double same_substep_drainage_allowance_cm,
    double *preferential_CR_input_cm,
    double *lgarto_domain_CR_input_cm,
    const char *context_label)
{
  if (state == NULL ||
      preferential_CR_input_cm == NULL ||
      lgarto_domain_CR_input_cm == NULL ||
      !state->lgar_bmi_params.mobile_groundwater_level ||
      !state->lgar_bmi_params.lower_bdy_flux_to_CR) {
    return 0.0;
  }

  const double raw_preferential_cm =
    fmax(0.0, *preferential_CR_input_cm);
  const double raw_domain_cm =
    fmax(0.0, *lgarto_domain_CR_input_cm);
  const double raw_CR_input_cm = raw_preferential_cm + raw_domain_cm;
  if (raw_CR_input_cm <= SMALL_EPS) {
    return 0.0;
  }

  const double cap_depth_cm =
    lgar_mobile_groundwater_depth_from_current_CR_storage(state);
  const double available_storage_cm =
    lgar_mobile_groundwater_storage_capacity_to_surface_cm(state,
                                                          cap_depth_cm);

  double accepted_CR_input_cm = raw_CR_input_cm;
  const double rejected_CR_input_cm =
    lgar_cap_CR_input_by_available_storage(
      raw_CR_input_cm,
      available_storage_cm,
      same_substep_drainage_allowance_cm,
      &accepted_CR_input_cm);

  const double accepted_fraction =
    raw_CR_input_cm > SMALL_EPS
      ? fmax(0.0, fmin(1.0, accepted_CR_input_cm / raw_CR_input_cm))
      : 1.0;
  *preferential_CR_input_cm = raw_preferential_cm * accepted_fraction;
  *lgarto_domain_CR_input_cm = raw_domain_cm * accepted_fraction;

  if (verbosity.compare("high") == 0 && rejected_CR_input_cm > SMALL_EPS) {
    printf("Mobile groundwater surface cap rejected %.17lf cm of CR-bound "
           "water as saturation-excess runoff (%s): raw_CR_input_cm=%.17lf "
           "accepted_CR_input_cm=%.17lf available_storage_cm=%.17lf "
           "same_substep_drainage_allowance_cm=%.17lf cap_depth_cm=%.17lf.\n",
           rejected_CR_input_cm,
           context_label != NULL ? context_label : "CR input",
           raw_CR_input_cm,
           accepted_CR_input_cm,
           available_storage_cm,
           same_substep_drainage_allowance_cm,
           cap_depth_cm);
  }

  return rejected_CR_input_cm;
}

static void lgar_set_mobile_groundwater_depth_from_storage_exchange(
    model_state *state,
    double reference_depth_cm,
    double net_groundwater_storage_gain_cm,
    const char *update_label)
{
  lgar_bmi_parameters *params = state != NULL ? &state->lgar_bmi_params : NULL;
  if (params == NULL || !params->mobile_groundwater_level) {
    return;
  }

  if (!std::isfinite(reference_depth_cm)) {
    reference_depth_cm = lgar_effective_groundwater_depth_cm(params);
  }
  if (!std::isfinite(reference_depth_cm)) {
    reference_depth_cm = lgar_fixed_soil_depth_cm(params);
  }

  if (!std::isfinite(net_groundwater_storage_gain_cm)) {
    net_groundwater_storage_gain_cm = 0.0;
  }

  // Depth is a diagnostic/proxy for existing groundwater-reservoir storage.
  // It is not added as a separate mass-balance storage term.
  const double updated_depth_cm =
    lgar_solve_mobile_groundwater_depth_from_storage_exchange(
      state, reference_depth_cm, net_groundwater_storage_gain_cm);
  const double effective_storage_coefficient =
    (std::isfinite(updated_depth_cm) &&
     fabs(reference_depth_cm - updated_depth_cm) > SMALL_EPS)
        ? net_groundwater_storage_gain_cm / (reference_depth_cm - updated_depth_cm)
        : 0.0;

  if (verbosity.compare("high") == 0 &&
      std::fabs(net_groundwater_storage_gain_cm) > SMALL_EPS) {
    printf("Mobile groundwater depth update (%s): reference_depth_cm=%.17lf "
           "net_storage_gain_cm=%.17lf effective_storage_coefficient=%.17lf "
           "new_depth_cm=%.17lf.\n",
           update_label != NULL ? update_label : "incremental",
           reference_depth_cm,
           net_groundwater_storage_gain_cm,
           effective_storage_coefficient,
           updated_depth_cm);
  }

  if (std::isfinite(updated_depth_cm)) {
    params->groundwater_depth_cm = fmax(0.0, updated_depth_cm);
  }
}

static void lgar_update_mobile_groundwater_depth(model_state *state,
                                                 double groundwater_storage_exchange_cm,
                                                 double reservoir_discharge_cm)
{
  lgar_bmi_parameters *params = state != NULL ? &state->lgar_bmi_params : NULL;
  if (params == NULL || !params->mobile_groundwater_level) {
    return;
  }

  double current_depth_cm = lgar_effective_groundwater_depth_cm(params);
  if (!std::isfinite(current_depth_cm)) {
    current_depth_cm = lgar_fixed_soil_depth_cm(params);
  }

  if (!std::isfinite(groundwater_storage_exchange_cm)) {
    groundwater_storage_exchange_cm = 0.0;
  }
  if (!std::isfinite(reservoir_discharge_cm)) {
    reservoir_discharge_cm = 0.0;
  }

  lgar_set_mobile_groundwater_depth_from_storage_exchange(
    state,
    current_depth_cm,
    groundwater_storage_exchange_cm - reservoir_discharge_cm,
    "incremental");
}

static void lgar_update_mobile_groundwater_depth_from_CR_storage(model_state *state)
{
  if (state == NULL ||
      !state->lgar_bmi_params.mobile_groundwater_level ||
      !state->lgar_bmi_params.lower_bdy_flux_to_CR) {
    return;
  }

  double reference_depth_cm =
    state->lgar_bmi_params.mobile_groundwater_reference_depth_cm;
  if (!std::isfinite(reference_depth_cm) || reference_depth_cm < 0.0) {
    reference_depth_cm = lgar_fixed_soil_depth_cm(&state->lgar_bmi_params);
    state->lgar_bmi_params.mobile_groundwater_reference_depth_cm =
      reference_depth_cm;
  }

  double reference_CR_storage_cm =
    state->lgar_bmi_params.mobile_groundwater_reference_CR_storage_cm;
  if (!std::isfinite(reference_CR_storage_cm)) {
    reference_CR_storage_cm = 0.0;
    state->lgar_bmi_params.mobile_groundwater_reference_CR_storage_cm =
      reference_CR_storage_cm;
  }

  const double net_CR_storage_gain_cm =
    lgar_total_CR_storage_cm(state) - reference_CR_storage_cm;
  lgar_set_mobile_groundwater_depth_from_storage_exchange(
    state,
    reference_depth_cm,
    net_CR_storage_gain_cm,
    "CR-storage");
}

static double lgar_sync_mobile_groundwater_chain_from_CR_storage(model_state *state)
{
  if (state == NULL ||
      !state->lgar_bmi_params.mobile_groundwater_level ||
      !state->lgar_bmi_params.lower_bdy_flux_to_CR) {
    return 0.0;
  }

  if (!state->lgar_bmi_params.TO_enabled) {
    lgar_update_mobile_groundwater_depth_from_CR_storage(state);
    return 0.0;
  }

  const int num_layers = state->lgar_bmi_params.num_layers;
  double *cum_layer_thickness_cm =
    state->lgar_bmi_params.cum_layer_thickness_cm;
  int *soil_type = state->lgar_bmi_params.layer_soil_type;
  double *frozen_factor = state->lgar_bmi_params.frozen_factor;
  const double fixed_column_depth_cm = lgar_fixed_soil_depth_cm(&state->lgar_bmi_params);
  const double target_CR_storage_cm = fmax(0.0, lgar_total_CR_storage_cm(state));
  double current_groundwater_depth_cm =
    lgar_effective_groundwater_depth_cm(&state->lgar_bmi_params);
  if (std::isfinite(fixed_column_depth_cm) && fixed_column_depth_cm > 0.0) {
    current_groundwater_depth_cm =
      fmin(current_groundwater_depth_cm, fixed_column_depth_cm);
  }
  const double storage_tol_cm =
    fmax(1.0e-4, 1.0e-4 * fmax(1.0, target_CR_storage_cm));
  const double mass_before_cm =
    lgar_calc_mass_bal(cum_layer_thickness_cm, state->head);
  auto sync_trial_is_usable = [&](double chain_storage_after_cm,
                                  double actual_chain_storage_change_cm,
                                  double explicit_mass_change_cm,
                                  double trial_mass_after_cm) -> bool {
    const double explicit_chain_mismatch_tol_cm =
      fmax(1.0e-6, 1.0e-4 * fmax(1.0, target_CR_storage_cm));
    const bool near_zero_chain_change =
      fabs(actual_chain_storage_change_cm) <= storage_tol_cm;
    const bool suspicious_explicit_jump =
      near_zero_chain_change &&
      fabs(explicit_mass_change_cm - actual_chain_storage_change_cm) >
        explicit_chain_mismatch_tol_cm;
    return std::isfinite(actual_chain_storage_change_cm) &&
           std::isfinite(trial_mass_after_cm) &&
           fabs(chain_storage_after_cm - target_CR_storage_cm) <=
             storage_tol_cm &&
           !suspicious_explicit_jump;
  };

  wetting_front *trial_head = listCopy(state->head, NULL);
  if (trial_head == NULL) {
    lgar_update_mobile_groundwater_depth_from_CR_storage(state);
    return 0.0;
  }

  double chain_storage_before_cm =
    lgarto_mobile_groundwater_CR_storage_cm(num_layers, cum_layer_thickness_cm,
                                            soil_type, trial_head,
                                            state->soil_properties);
  double updated_groundwater_depth_cm = current_groundwater_depth_cm;
  double trial_explicit_mass_change_cm = 0.0;
  double actual_chain_storage_change_cm =
    lgarto_sync_mobile_groundwater_support_to_CR_storage(
      target_CR_storage_cm,
      current_groundwater_depth_cm,
      num_layers,
      cum_layer_thickness_cm,
      soil_type,
      frozen_factor,
      &trial_head,
      state->soil_properties,
      &updated_groundwater_depth_cm,
      &trial_explicit_mass_change_cm);
  double chain_storage_after_cm =
    lgarto_mobile_groundwater_CR_storage_cm(num_layers, cum_layer_thickness_cm,
                                            soil_type, trial_head,
                                            state->soil_properties);
  double mass_after_cm = lgar_calc_mass_bal(cum_layer_thickness_cm, trial_head);
  bool sync_usable =
    sync_trial_is_usable(chain_storage_after_cm,
                         actual_chain_storage_change_cm,
                         trial_explicit_mass_change_cm,
                         mass_after_cm);

  if (!sync_usable) {
    listDelete(trial_head);
    trial_head = listCopy(state->head, NULL);
    if (trial_head == NULL) {
      lgar_update_mobile_groundwater_depth_from_CR_storage(state);
      return 0.0;
    }

    (void)
      lgarto_maintain_mobile_groundwater_support(
        current_groundwater_depth_cm,
        num_layers,
        cum_layer_thickness_cm,
        soil_type,
        frozen_factor,
        &trial_head,
        state->soil_properties);
    chain_storage_before_cm =
      lgarto_mobile_groundwater_CR_storage_cm(num_layers, cum_layer_thickness_cm,
                                              soil_type, trial_head,
                                              state->soil_properties);
    updated_groundwater_depth_cm = current_groundwater_depth_cm;
    trial_explicit_mass_change_cm = 0.0;
    actual_chain_storage_change_cm =
      lgarto_sync_mobile_groundwater_support_to_CR_storage(
        target_CR_storage_cm,
        current_groundwater_depth_cm,
        num_layers,
        cum_layer_thickness_cm,
        soil_type,
        frozen_factor,
        &trial_head,
        state->soil_properties,
        &updated_groundwater_depth_cm,
        &trial_explicit_mass_change_cm);
    chain_storage_after_cm =
      lgarto_mobile_groundwater_CR_storage_cm(num_layers, cum_layer_thickness_cm,
                                              soil_type, trial_head,
                                              state->soil_properties);
    mass_after_cm = lgar_calc_mass_bal(cum_layer_thickness_cm, trial_head);
    sync_usable =
      sync_trial_is_usable(chain_storage_after_cm,
                           actual_chain_storage_change_cm,
                           trial_explicit_mass_change_cm,
                           mass_after_cm);
  }

  if (!sync_usable) {
    const double live_chain_storage_cm =
      lgarto_mobile_groundwater_CR_storage_cm(num_layers, cum_layer_thickness_cm,
                                              soil_type, state->head,
                                              state->soil_properties);
    const double maintained_trial_tol_cm =
      fmax(1.0e-2, 1.0e-3 * fmax(1.0, target_CR_storage_cm));
    if (live_chain_storage_cm <= storage_tol_cm &&
        chain_storage_before_cm > storage_tol_cm &&
        std::isfinite(chain_storage_after_cm) &&
        std::isfinite(mass_after_cm) &&
        fabs(chain_storage_after_cm - target_CR_storage_cm) <=
          maintained_trial_tol_cm) {
      sync_usable = true;
    }
  }

  if (!sync_usable) {
    if (verbosity.compare("high") == 0) {
      printf("Mobile groundwater CR-chain sync skipped: target_CR_storage_cm=%.17lf "
             "chain_before_cm=%.17lf chain_after_cm=%.17lf depth_cm=%.17lf.\n",
             target_CR_storage_cm,
             chain_storage_before_cm,
             chain_storage_after_cm,
             current_groundwater_depth_cm);
    }
    listDelete(trial_head);
	    if (chain_storage_before_cm > storage_tol_cm) {
	      if (verbosity.compare("high") == 0) {
	        printf("Mobile groundwater CR-chain sync preserving existing support chain "
	               "instead of applying proxy depth fallback.\n");
	      }
	      return lgarto_cleanup_redundant_colocated_psi0_support_stacks(
	        cum_layer_thickness_cm, soil_type, &state->head,
	        state->soil_properties);
	    }
	    lgar_update_mobile_groundwater_depth_from_CR_storage(state);
	    return lgarto_cleanup_redundant_colocated_psi0_support_stacks(
	      cum_layer_thickness_cm, soil_type, &state->head,
	      state->soil_properties);
	  }

	  double explicit_mass_change_cm = mass_after_cm - mass_before_cm;
	  listDelete(state->head);
	  state->head = trial_head;

  if (std::isfinite(updated_groundwater_depth_cm)) {
	    state->lgar_bmi_params.groundwater_depth_cm =
	      fmax(0.0, updated_groundwater_depth_cm);
	  }
	  explicit_mass_change_cm +=
	    lgarto_cleanup_redundant_colocated_psi0_support_stacks(
	      cum_layer_thickness_cm, soil_type, &state->head,
	      state->soil_properties);

	  if (verbosity.compare("high") == 0 &&
      (fabs(actual_chain_storage_change_cm) > SMALL_EPS ||
       fabs(explicit_mass_change_cm) > SMALL_EPS ||
       fabs(chain_storage_after_cm - target_CR_storage_cm) > storage_tol_cm)) {
    printf("Mobile groundwater CR-chain sync: target_CR_storage_cm=%.17lf "
           "chain_before_cm=%.17lf chain_after_cm=%.17lf "
           "actual_chain_storage_change_cm=%.17lf "
           "explicit_mass_change_cm=%.17lf depth_cm=%.17lf.\n",
           target_CR_storage_cm,
           chain_storage_before_cm,
           chain_storage_after_cm,
           actual_chain_storage_change_cm,
           explicit_mass_change_cm,
           state->lgar_bmi_params.groundwater_depth_cm);
  }

  return explicit_mass_change_cm;
}

static double lgar_CR_capillary_supply_scale(const model_state *state,
                                             double upward_demand_cm)
{
  if (state == NULL || upward_demand_cm <= SMALL_EPS ||
      !state->lgar_bmi_params.TO_enabled ||
      !state->lgar_bmi_params.mobile_groundwater_level ||
      !state->lgar_bmi_params.lower_bdy_flux_to_CR) {
    return 1.0;
  }

  const double threshold_cm =
    state->lgar_bmi_params.CR_capillary_supply_threshold_cm;
  if (!std::isfinite(threshold_cm) || threshold_cm <= 0.0) {
    return 1.0;
  }

  const double CR_storage_cm =
    fmax(state->lgar_mass_balance.CR_fast_storage_cm, 0.0) +
    fmax(state->lgar_mass_balance.CR_slow_storage_cm, 0.0);
  if (CR_storage_cm <= SMALL_EPS) {
    return 0.0;
  }

  double scale = fmin(1.0, CR_storage_cm / threshold_cm);
  scale = fmin(scale, CR_storage_cm / upward_demand_cm);
  return fmax(0.0, fmin(1.0, scale));
}

static double lgar_project_upward_TO_flux_demand_cm(const model_state *state,
                                                    double subtimestep_h)
{
  if (state == NULL || state->head == NULL || subtimestep_h <= 0.0) {
    return 0.0;
  }

  double upward_demand_cm = 0.0;
  const int wetting_front_count = listLength(state->head);
  for (int wf = wetting_front_count - 1; wf >= 1; wf--) {
    const wetting_front *current = listFindFront(wf, state->head, NULL);
    if (current == NULL || !current->is_WF_GW || current->to_bottom ||
        current->dzdt_cm_per_h >= 0.0) {
      continue;
    }

    const wetting_front *next_to_use = current->next;
    while (next_to_use != NULL && !next_to_use->is_WF_GW) {
      next_to_use = next_to_use->next;
    }

    if (next_to_use == NULL ||
        current->layer_num < 1 ||
        current->layer_num > state->lgar_bmi_params.num_layers) {
      continue;
    }

    double delta_depth_cm = current->dzdt_cm_per_h * subtimestep_h;
    if (current->depth_cm + delta_depth_cm < 0.0) {
      delta_depth_cm = -current->depth_cm;
    }
    if (delta_depth_cm >= 0.0) {
      continue;
    }

    double delta_theta = 0.0;
    if (current->layer_num == next_to_use->layer_num) {
      delta_theta = next_to_use->theta - current->theta;
    }
    else {
      const int soil_num_current =
        state->lgar_bmi_params.layer_soil_type[current->layer_num];
      const double theta_e_current =
        state->soil_properties[soil_num_current].theta_e;
      const double theta_r_current =
        state->soil_properties[soil_num_current].theta_r;
      const double vg_a_current =
        state->soil_properties[soil_num_current].vg_alpha_per_cm;
      const double vg_m_current =
        state->soil_properties[soil_num_current].vg_m;
      const double vg_n_current =
        state->soil_properties[soil_num_current].vg_n;
      const double equiv_next_theta =
        calc_theta_from_h(next_to_use->psi_cm, vg_a_current, vg_m_current,
                          vg_n_current, theta_e_current, theta_r_current);
      delta_theta = equiv_next_theta - current->theta;
    }

    const double projected_flux_cm = delta_depth_cm * delta_theta;
    if (std::isfinite(projected_flux_cm) && projected_flux_cm < 0.0) {
      upward_demand_cm += -projected_flux_cm;
    }
  }

  return upward_demand_cm;
}

static void lgar_limit_upward_TO_dzdt_by_CR_supply(model_state *state,
                                                   double subtimestep_h)
{
  const double upward_demand_cm =
    lgar_project_upward_TO_flux_demand_cm(state, subtimestep_h);
  const double supply_scale =
    lgar_CR_capillary_supply_scale(state, upward_demand_cm);

  if (supply_scale >= 1.0 - SMALL_EPS) {
    return;
  }

  for (wetting_front *current = state->head; current != NULL; current = current->next) {
    if (current->is_WF_GW && !current->to_bottom &&
        current->dzdt_cm_per_h < 0.0) {
      current->dzdt_cm_per_h =
        (supply_scale <= SMALL_EPS) ? 0.0 : current->dzdt_cm_per_h * supply_scale;
    }
  }

  if (verbosity.compare("high") == 0) {
    printf("Limited upward TO dzdt by CR capillary supply: "
           "upward_demand_cm=%.12e scale=%.12e CR_storage_cm=%.12e "
           "threshold_cm=%.12e.\n",
           upward_demand_cm,
           supply_scale,
           fmax(state->lgar_mass_balance.CR_fast_storage_cm, 0.0) +
             fmax(state->lgar_mass_balance.CR_slow_storage_cm, 0.0),
           state->lgar_bmi_params.CR_capillary_supply_threshold_cm);
  }
}

static double lgar_limit_cached_negative_lower_boundary_flux_by_CR_supply(
    const model_state *state,
    double lower_boundary_flux_cm)
{
  if (lower_boundary_flux_cm >= 0.0) {
    return lower_boundary_flux_cm;
  }

  const double upward_demand_cm = -lower_boundary_flux_cm;
  const double supply_scale =
    lgar_CR_capillary_supply_scale(state, upward_demand_cm);
  if (supply_scale >= 1.0 - SMALL_EPS) {
    return lower_boundary_flux_cm;
  }

  return (supply_scale <= SMALL_EPS) ? 0.0 : lower_boundary_flux_cm * supply_scale;
}

static double lgar_max_abs_mobile_dzdt_for_flux_cache(const wetting_front *head)
{
  double max_abs_dzdt = 0.0;

  for (const wetting_front *current = head; current != NULL; current = current->next) {
    if (current->to_bottom) {
      continue;
    }

    max_abs_dzdt = fmax(max_abs_dzdt, fabs(current->dzdt_cm_per_h));
  }

  return max_abs_dzdt;
}

static double lgarto_remaining_lower_boundary_capacity_cm(double subtimestep_h,
                                                          int num_layers,
                                                          const double *cum_layer_thickness_cm,
                                                          const wetting_front *head,
                                                          double already_booked_flux_cm)
{
  if (subtimestep_h <= 0.0 || num_layers <= 0 || cum_layer_thickness_cm == NULL) {
    return 0.0;
  }

  const double domain_depth_cm = cum_layer_thickness_cm[num_layers];
  const double depth_tol_cm = fmax(1.0e-8, 1.0e-10 * fmax(1.0, domain_depth_cm));
  const wetting_front *bottom_front = NULL;
  for (const wetting_front *current = head; current != NULL; current = current->next) {
    if (current->to_bottom && current->layer_num == num_layers &&
        std::fabs(current->depth_cm - domain_depth_cm) <= depth_tol_cm) {
      bottom_front = current;
    }
  }

  if (bottom_front == NULL ||
      !std::isfinite(bottom_front->K_cm_per_h) ||
      bottom_front->K_cm_per_h <= 0.0) {
    return 0.0;
  }

  const double hydraulic_capacity_cm = bottom_front->K_cm_per_h * subtimestep_h;
  return fmax(0.0, hydraulic_capacity_cm - fmax(0.0, already_booked_flux_cm));
}

static double lgarto_saturated_creation_lower_boundary_capacity_cm(double subtimestep_h,
                                                                   int num_layers,
                                                                   const double *cum_layer_thickness_cm,
                                                                   int *soil_type,
                                                                   double *frozen_factor,
                                                                   const wetting_front *head,
                                                                   struct soil_properties_ *soil_properties,
                                                                   double already_booked_flux_cm)
{
  if (subtimestep_h <= 0.0 || num_layers <= 0 || cum_layer_thickness_cm == NULL ||
      soil_type == NULL || frozen_factor == NULL || head == NULL || soil_properties == NULL ||
      LGARTO_SATURATED_CREATION_PRESSURE_GRADIENT_CAP <= 0.0) {
    return 0.0;
  }

  const double domain_depth_cm = cum_layer_thickness_cm[num_layers];
  if (domain_depth_cm <= 0.0) {
    return 0.0;
  }

  const double depth_tol_cm = fmax(1.0e-8, 1.0e-10 * fmax(1.0, domain_depth_cm));
  std::vector<int> saturated_to_bottom_layers(num_layers + 1, 0);

  for (const wetting_front *current = head; current != NULL; current = current->next) {
    if (!current->to_bottom) {
      continue;
    }

    const int layer_num = current->layer_num;
    if (layer_num < 1 || layer_num > num_layers) {
      return 0.0;
    }

    const double expected_depth_cm = cum_layer_thickness_cm[layer_num];
    if (std::fabs(current->depth_cm - expected_depth_cm) > depth_tol_cm) {
      return 0.0;
    }

    const int soil_num = soil_type[layer_num];
    const double theta_e = soil_properties[soil_num].theta_e;
    const bool near_theta_e = current->theta >= theta_e - LGARTO_SATURATED_CREATION_THETA_TOL;
    const bool near_zero_psi =
      std::isfinite(current->psi_cm) &&
      current->psi_cm <= LGARTO_SATURATED_CREATION_PSI_TOL_CM;
    if (!(near_theta_e || near_zero_psi)) {
      return 0.0;
    }

    saturated_to_bottom_layers[layer_num] = 1;
  }

  for (int layer_num = 1; layer_num <= num_layers; layer_num++) {
    if (!saturated_to_bottom_layers[layer_num]) {
      return 0.0;
    }
  }

  double saturated_resistance_h = 0.0;
  for (int layer_num = 1; layer_num <= num_layers; layer_num++) {
    const double layer_top_cm = cum_layer_thickness_cm[layer_num - 1];
    const double layer_bottom_cm = cum_layer_thickness_cm[layer_num];
    const double layer_thickness_cm = layer_bottom_cm - layer_top_cm;
    const int soil_num = soil_type[layer_num];
    const double Ksat_cm_per_h =
      soil_properties[soil_num].Ksat_cm_per_h * frozen_factor[layer_num];

    if (layer_thickness_cm <= 0.0 ||
        !std::isfinite(Ksat_cm_per_h) ||
        Ksat_cm_per_h <= 0.0) {
      return 0.0;
    }

    saturated_resistance_h += layer_thickness_cm / Ksat_cm_per_h;
  }

  if (saturated_resistance_h <= 0.0 || !std::isfinite(saturated_resistance_h)) {
    return 0.0;
  }

  /* Surface-creation repair can leave excess water only because the
     boundary-pinned TO scaffold is effectively saturated. In that narrow
     regime, treat the lower boundary as a saturated layered column and allow
     a bounded positive-pressure gradient to augment unit-gradient drainage.
     The returned value is remaining capacity after fluxes already booked in
     this substep, so this does not stack a second drainage term on top of the
     natural TO/GW flux bookkeeping. */
  const double effective_Ksat_cm_per_h = domain_depth_cm / saturated_resistance_h;
  const double saturated_column_flux_capacity_cm =
    effective_Ksat_cm_per_h *
    (1.0 + LGARTO_SATURATED_CREATION_PRESSURE_GRADIENT_CAP) *
    subtimestep_h;

  return fmax(0.0, saturated_column_flux_capacity_cm -
                   fmax(0.0, already_booked_flux_cm));
}

static bool lgarto_has_TO_fronts(const wetting_front *head)
{
  for (const wetting_front *current = head; current != NULL; current = current->next) {
    if (current->is_WF_GW) {
      return true;
    }
  }

  return false;
}

struct lgarto_infiltration_limit_trace_row
{
  int timestep;
  double time_h;
  int cycle;
  int subcycles;
  double subtimestep_h;
  const char *branch;
  int cache_fluxes;
  int create_surficial_front;
  int top_wf_saturated;
  int fronts_start;
  int surface_fronts_start;
  int fronts_end;
  int surface_fronts_end;
  double precip_cm;
  double PET_cm;
  double ponded_start_cm;
  double ponded_end_cm;
  double volin_cm;
  double runoff_before_creation_excess_cm;
  double creation_excess_runoff_cm;
  double runoff_after_creation_excess_cm;
  double creation_excess_gw_flux_cm;
  double surface_creation_gw_capacity_cm;
  double surface_creation_post_creation_TO_release_cm;
  double temp_rch_cm;
  double free_drainage_cm;
  double projected_TO_storage_release_cm;
  double insertion_storage_release_cm;
  double insert_raw_fp_cm_per_h;
  double insert_storage_limit_fp_cm_per_h;
  double insert_capped_fp_cm_per_h;
  int insert_storage_cap_active;
  double lower_boundary_flux_cm;
};

static FILE *lgarto_infiltration_limit_trace_file()
{
  const char *path = getenv("LGARTO_INF_LIMIT_TRACE");
  if (path == NULL || path[0] == '\0') {
    return NULL;
  }

  static FILE *trace_file = NULL;
  static bool failed_open = false;
  if (trace_file != NULL || failed_open) {
    return trace_file;
  }

  trace_file = fopen(path, "w");
  if (trace_file == NULL) {
    failed_open = true;
    fprintf(stderr, "Warning: could not open LGARTO_INF_LIMIT_TRACE file '%s'.\n", path);
    return NULL;
  }

  fprintf(trace_file,
          "timestep,time_h,cycle,subcycles,subtimestep_h,branch,cache_fluxes,"
          "create_surficial_front,top_wf_saturated,fronts_start,surface_fronts_start,"
          "fronts_end,surface_fronts_end,precip_cm,PET_cm,ponded_start_cm,ponded_end_cm,"
          "volin_cm,runoff_before_creation_excess_cm,creation_excess_runoff_cm,"
          "runoff_after_creation_excess_cm,creation_excess_gw_flux_cm,"
          "surface_creation_gw_capacity_cm,surface_creation_post_creation_TO_release_cm,temp_rch_cm,"
          "free_drainage_cm,projected_TO_storage_release_cm,insertion_storage_release_cm,"
          "insert_raw_fp_cm_per_h,insert_storage_limit_fp_cm_per_h,insert_capped_fp_cm_per_h,"
          "insert_storage_cap_active,lower_boundary_flux_cm\n");
  fflush(trace_file);

  return trace_file;
}

static void lgarto_write_infiltration_limit_trace(
  const lgarto_infiltration_limit_trace_row &row)
{
  FILE *trace_file = lgarto_infiltration_limit_trace_file();
  if (trace_file == NULL) {
    return;
  }

  fprintf(trace_file,
          "%d,%.17g,%d,%d,%.17g,%s,%d,%d,%d,%d,%d,%d,%d,"
          "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,"
          "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d,%.17g\n",
          row.timestep,
          row.time_h,
          row.cycle,
          row.subcycles,
          row.subtimestep_h,
          row.branch,
          row.cache_fluxes,
          row.create_surficial_front,
          row.top_wf_saturated,
          row.fronts_start,
          row.surface_fronts_start,
          row.fronts_end,
          row.surface_fronts_end,
          row.precip_cm,
          row.PET_cm,
          row.ponded_start_cm,
          row.ponded_end_cm,
          row.volin_cm,
          row.runoff_before_creation_excess_cm,
          row.creation_excess_runoff_cm,
          row.runoff_after_creation_excess_cm,
          row.creation_excess_gw_flux_cm,
          row.surface_creation_gw_capacity_cm,
          row.surface_creation_post_creation_TO_release_cm,
          row.temp_rch_cm,
          row.free_drainage_cm,
          row.projected_TO_storage_release_cm,
          row.insertion_storage_release_cm,
          row.insert_raw_fp_cm_per_h,
          row.insert_storage_limit_fp_cm_per_h,
          row.insert_capped_fp_cm_per_h,
          row.insert_storage_cap_active,
          row.lower_boundary_flux_cm);
  fflush(trace_file);
}

static double lgarto_limit_subtimestep_for_surface_boundary_crossing(
  double proposed_subtimestep_h,
  double groundwater_depth_cm,
  int num_layers,
  const double *cum_layer_thickness_cm,
  const wetting_front *head)
{
  if (proposed_subtimestep_h <= 0.0 || head == NULL || num_layers < 2 ||
      !lgarto_has_TO_fronts(head)) {
    return proposed_subtimestep_h;
  }

  double limited_subtimestep_h = proposed_subtimestep_h;

  for (const wetting_front *current = head; current != NULL; current = current->next) {
    const wetting_front *next = current->next;
    if (next == NULL || next->next == NULL) {
      continue;
    }

    const int layer_num = current->layer_num;
    if (current->is_WF_GW || current->to_bottom || current->dzdt_cm_per_h <= 0.0 ||
        layer_num < 1 || layer_num >= num_layers) {
      continue;
    }

    const double boundary_depth_cm = cum_layer_thickness_cm[layer_num];
    const double boundary_tol_cm =
      1.0e-10 * fmax(1.0, std::fabs(boundary_depth_cm));
    if (std::isfinite(groundwater_depth_cm) && groundwater_depth_cm > 0.0 &&
        groundwater_depth_cm < boundary_depth_cm - boundary_tol_cm) {
      continue;
    }

    const double projected_depth_cm =
      current->depth_cm + current->dzdt_cm_per_h * proposed_subtimestep_h;
    if (current->depth_cm >= boundary_depth_cm || projected_depth_cm <= boundary_depth_cm ||
        !next->to_bottom || next->depth_cm != boundary_depth_cm) {
      continue;
    }

    // Stop at the geometric event; the mover remains the sole owner of the
    // hydraulic remap performed immediately after the boundary is reached.
    const double event_subtimestep_h =
      (boundary_depth_cm + LGARTO_EVENT_SPLIT_BOUNDARY_EPS_CM - current->depth_cm) /
      current->dzdt_cm_per_h;
    if (event_subtimestep_h <= 0.0 || event_subtimestep_h >= limited_subtimestep_h) {
      continue;
    }

    limited_subtimestep_h = event_subtimestep_h;

    if (verbosity.compare("high") == 0) {
      printf("Adaptive LGARTO event split at surface layer-boundary crossing: "
             "front=%d layer=%d old_dt_h=%.17lf new_dt_h=%.17lf "
             "depth_cm=%.17lf projected_depth_cm=%.17lf boundary_cm=%.17lf\n",
             current->front_num,
             layer_num,
             proposed_subtimestep_h,
             limited_subtimestep_h,
             current->depth_cm,
             projected_depth_cm,
             boundary_depth_cm);
    }
  }

  return limited_subtimestep_h;
}

/*
  Project the same interflow-before-dzdt, deepest-to-shallowest depth operator
  used for movable TO/GW fronts in lgar_move_wetting_fronts().  This is only an
  event-time predictor: it does not change theta, psi, K, mass, or the live
  wetting-front list.

  Return false when the requested pair would require one of the more complex
  surface/layer-boundary caps in the mover.  The targeted TO/to_bottom pair is
  the exception: its uncapped projection is needed to locate boundary contact.
*/
static bool lgarto_project_mobile_TO_packet_gap_after_interflow(
  double subtimestep_h,
  const model_state *state,
  int current_front_num,
  int next_front_num,
  double *projected_current_depth_cm,
  double *projected_next_depth_cm,
  double *projected_gap_cm)
{
  if (state == NULL || state->head == NULL || subtimestep_h < 0.0 ||
      current_front_num < 1 || next_front_num < 1 ||
      projected_current_depth_cm == NULL ||
      projected_next_depth_cm == NULL ||
      projected_gap_cm == NULL) {
    return false;
  }

  const lgar_bmi_parameters& params = state->lgar_bmi_params;
  if (params.num_layers < 1 ||
      params.cum_layer_thickness_cm == NULL) {
    return false;
  }

  const int wetting_front_count = listLength(state->head);
  if (wetting_front_count < 2 ||
      current_front_num > wetting_front_count ||
      next_front_num > wetting_front_count) {
    return false;
  }

  std::vector<wetting_front *> front_by_num(wetting_front_count + 1, NULL);
  std::vector<double> projected_depth_by_num(
    wetting_front_count + 1, std::numeric_limits<double>::quiet_NaN());
  for (wetting_front *front = state->head;
       front != NULL;
       front = front->next) {
    if (front->front_num < 1 || front->front_num > wetting_front_count ||
        !std::isfinite(front->depth_cm)) {
      return false;
    }
    front_by_num[front->front_num] = front;
    projected_depth_by_num[front->front_num] = front->depth_cm;
  }

  const double column_depth_cm =
    params.cum_layer_thickness_cm[params.num_layers];
  if (!std::isfinite(column_depth_cm) || column_depth_cm <= 0.0) {
    return false;
  }

  const double upward_TO_supply_scale =
    lgar_CR_capillary_supply_scale(
      state, lgar_project_upward_TO_flux_demand_cm(state, subtimestep_h));

  for (int wf = wetting_front_count - 1; wf >= 1; --wf) {
    wetting_front *current = front_by_num[wf];
    if (current == NULL || !current->is_WF_GW || current->to_bottom) {
      continue;
    }
    if (current->layer_num < 1 || current->layer_num > params.num_layers ||
        current->next == NULL ||
        current->next->front_num < 1 ||
        current->next->front_num > wetting_front_count) {
      return false;
    }

    const int layer_num = current->layer_num;
    const double layer_top_cm = params.cum_layer_thickness_cm[layer_num - 1];
    const double layer_bottom_cm = params.cum_layer_thickness_cm[layer_num];
    double current_depth_cm = projected_depth_by_num[wf];

    if (params.lateral_flow_enabled && params.lateral_flow_factor > 0.0 &&
        current->next->layer_num == layer_num) {
      const double requested_lateral_flux_cm =
        lgarto_lateral_flux_candidate_cm(
          subtimestep_h, params.num_layers,
          params.lateral_flow_psi_threshold_cm, params.lateral_flow_factor,
          params.cum_layer_thickness_cm, state->head,
          wf > 1 ? front_by_num[wf - 1] : NULL, current,
          params.mobile_groundwater_level, params.layer_soil_type,
          state->soil_properties);
      current_depth_cm =
        lgarto_project_TO_interflow_depth_cm(
          requested_lateral_flux_cm, current_depth_cm, layer_bottom_cm,
          projected_depth_by_num[current->next->front_num],
          current->theta, current->next->theta, NULL);
    }

    double dzdt_cm_per_h = current->dzdt_cm_per_h;
    if (!std::isfinite(dzdt_cm_per_h)) {
      return false;
    }
    if (dzdt_cm_per_h < 0.0 &&
        upward_TO_supply_scale < 1.0 - SMALL_EPS) {
      dzdt_cm_per_h *= upward_TO_supply_scale;
    }
    const double vertically_projected_depth_cm =
      current_depth_cm + dzdt_cm_per_h * subtimestep_h;
    const bool targeted_to_bottom_boundary =
      current->front_num == current_front_num &&
      current->next->front_num == next_front_num &&
      current->next->to_bottom;

    /*
      The failing chatter occurs wholly within a layer.  Do not substitute this
      lightweight predictor where the mover would invoke a surface or
      unrelated layer-boundary cap/canonicalization.
    */
    if (vertically_projected_depth_cm < layer_top_cm ||
        (vertically_projected_depth_cm > layer_bottom_cm &&
         !targeted_to_bottom_boundary)) {
      return false;
    }
    if (dzdt_cm_per_h < 0.0) {
      for (const wetting_front *candidate = state->head;
           candidate != NULL;
           candidate = candidate->next) {
        if (candidate->is_WF_GW || candidate->to_bottom ||
            candidate->layer_num != layer_num) {
          continue;
        }
        if (candidate->depth_cm <= current_depth_cm &&
            candidate->depth_cm > vertically_projected_depth_cm) {
          return false;
        }
      }
    }

    projected_depth_by_num[wf] = vertically_projected_depth_cm;
    if (targeted_to_bottom_boundary) {
      *projected_current_depth_cm = vertically_projected_depth_cm;
      *projected_next_depth_cm =
        projected_depth_by_num[next_front_num];
      *projected_gap_cm =
        *projected_next_depth_cm - *projected_current_depth_cm;
      return std::isfinite(*projected_gap_cm);
    }
  }

  const double current_depth_cm =
    projected_depth_by_num[current_front_num];
  const double next_depth_cm =
    projected_depth_by_num[next_front_num];
  if (!std::isfinite(current_depth_cm) || !std::isfinite(next_depth_cm)) {
    return false;
  }

  *projected_current_depth_cm = current_depth_cm;
  *projected_next_depth_cm = next_depth_cm;
  *projected_gap_cm = next_depth_cm - current_depth_cm;
  return std::isfinite(*projected_gap_cm);
}

extern double lgarto_limit_subtimestep_for_mobile_TO_packet_overtake(
  double proposed_subtimestep_h,
  const model_state *state,
  int *repeat_front_num,
  int *repeat_next_front_num,
  int *repeat_layer_num,
  int *repeat_count)
{
  if (state == NULL || proposed_subtimestep_h <= 0.0 ||
      state->head == NULL ||
      !state->lgar_bmi_params.TO_enabled ||
      !state->lgar_bmi_params.mobile_groundwater_level ||
      !lgarto_has_TO_fronts(state->head)) {
    return proposed_subtimestep_h;
  }

  // The normal TO correction resolves interflow-assisted to_bottom crossings.
  double limited_subtimestep_h = proposed_subtimestep_h;
  if (state->lgar_bmi_params.lower_bdy_flux_to_CR &&
      lgar_total_CR_storage_cm(state) <= SMALL_EPS) {
    return limited_subtimestep_h;
  }

  const double min_spacing_cm =
    fmax(LGARTO_TO_PACKET_EVENT_SPLIT_MIN_SPACING_CM,
         10.0 * LGARTO_EVENT_SPLIT_BOUNDARY_EPS_CM);
  const wetting_front *limiting_current = NULL;
  const wetting_front *limiting_next = NULL;
  double limiting_current_gap_cm = 0.0;
  double limiting_closing_speed_cm_per_h = 0.0;
  double limiting_projected_current_depth_cm = 0.0;
  double limiting_projected_next_depth_cm = 0.0;
  double limiting_current_dzdt_cm_per_h = 0.0;
  double limiting_next_dzdt_cm_per_h = 0.0;
  bool limiting_event_is_interflow_aware = false;

  for (const wetting_front *current = state->head;
       current != NULL && current->next != NULL;
       current = current->next) {
    const wetting_front *next = current->next;
    if (!current->is_WF_GW || !next->is_WF_GW ||
        current->to_bottom || next->to_bottom ||
        current->layer_num != next->layer_num ||
        !std::isfinite(current->depth_cm) ||
        !std::isfinite(next->depth_cm) ||
        !std::isfinite(current->dzdt_cm_per_h) ||
        !std::isfinite(next->dzdt_cm_per_h)) {
      continue;
    }

    const double current_gap_cm = next->depth_cm - current->depth_cm;
    if (current_gap_cm <= LGARTO_TO_PACKET_EVENT_SPLIT_MERGE_READY_SPACING_CM) {
      continue;
    }
    if (current_gap_cm <= 2.0 * min_spacing_cm) {
      continue;
    }

    double ordered_projected_current_depth_cm = 0.0;
    double ordered_projected_next_depth_cm = 0.0;
    double ordered_projected_gap_cm = 0.0;
    const bool interflow_enabled =
      state->lgar_bmi_params.lateral_flow_enabled &&
      state->lgar_bmi_params.lateral_flow_factor > 0.0;
    const bool ordered_projection_valid =
      interflow_enabled &&
      lgarto_project_mobile_TO_packet_gap_after_interflow(
        proposed_subtimestep_h, state,
        current->front_num, next->front_num,
        &ordered_projected_current_depth_cm,
        &ordered_projected_next_depth_cm,
        &ordered_projected_gap_cm);

    // Match the event prediction to the CR-supply cap applied before TO motion.
    const double upward_TO_supply_scale =
      lgar_CR_capillary_supply_scale(
        state, lgar_project_upward_TO_flux_demand_cm(
                 state, proposed_subtimestep_h));
    double current_dzdt_cm_per_h = current->dzdt_cm_per_h;
    double next_dzdt_cm_per_h = next->dzdt_cm_per_h;
    if (upward_TO_supply_scale < 1.0 - SMALL_EPS) {
      if (current_dzdt_cm_per_h < 0.0) {
        current_dzdt_cm_per_h *= upward_TO_supply_scale;
      }
      if (next_dzdt_cm_per_h < 0.0) {
        next_dzdt_cm_per_h *= upward_TO_supply_scale;
      }
    }

    const double projected_current_depth_cm =
      current->depth_cm + current_dzdt_cm_per_h * proposed_subtimestep_h;
    const double projected_next_depth_cm =
      next->depth_cm + next_dzdt_cm_per_h * proposed_subtimestep_h;
    const double projected_gap_cm =
      projected_next_depth_cm - projected_current_depth_cm;
    if ((ordered_projection_valid ? ordered_projected_gap_cm :
         projected_gap_cm) >= min_spacing_cm) {
      continue;
    }

    const double closing_speed_cm_per_h =
      ordered_projection_valid
        ? (current_gap_cm - ordered_projected_gap_cm) /
            proposed_subtimestep_h
        : current_dzdt_cm_per_h - next_dzdt_cm_per_h;
    if (closing_speed_cm_per_h <= 0.0) {
      continue;
    }

    double event_subtimestep_h =
      (current_gap_cm - min_spacing_cm) / closing_speed_cm_per_h;
    double event_projected_current_depth_cm =
      ordered_projection_valid ? ordered_projected_current_depth_cm :
                                 projected_current_depth_cm;
    double event_projected_next_depth_cm =
      ordered_projection_valid ? ordered_projected_next_depth_cm :
                                 projected_next_depth_cm;
    bool event_is_interflow_aware = false;

    if (interflow_enabled) {
      double interflow_projected_current_depth_cm =
        ordered_projected_current_depth_cm;
      double interflow_projected_next_depth_cm =
        ordered_projected_next_depth_cm;
      double interflow_projected_gap_cm = ordered_projected_gap_cm;
      const bool projected = ordered_projection_valid;

      if (projected && interflow_projected_gap_cm >= min_spacing_cm) {
        // The actual ordered motion keeps this pair separated for the full step.
        continue;
      }

      if (!projected) {
        /*
          A full-step projection can leave the current layer before the packet
          would meet.  Recheck the earlier dzdt-only contact time with the actual
          interflow-before-dzdt operator.  If the pair is still separated there,
          this is not a packet event; the normal mover/correction pass will handle
          any later layer-boundary crossing.
        */
        const bool candidate_projected =
          lgarto_project_mobile_TO_packet_gap_after_interflow(
            event_subtimestep_h, state,
            current->front_num, next->front_num,
            &interflow_projected_current_depth_cm,
            &interflow_projected_next_depth_cm,
            &interflow_projected_gap_cm);
        if (candidate_projected &&
            interflow_projected_gap_cm >= min_spacing_cm) {
          continue;
        }
      }

      if (projected) {
        double safe_dt_h = 0.0;
        double crossing_dt_h = proposed_subtimestep_h;
        bool projection_valid = true;
        for (int iteration = 0;
             iteration < LGARTO_TO_PACKET_INTERFLOW_EVENT_BISECTION_ITERATIONS;
             ++iteration) {
          const double trial_dt_h = 0.5 * (safe_dt_h + crossing_dt_h);
          double trial_current_depth_cm = 0.0;
          double trial_next_depth_cm = 0.0;
          double trial_gap_cm = 0.0;
          if (!lgarto_project_mobile_TO_packet_gap_after_interflow(
                trial_dt_h, state,
                current->front_num, next->front_num,
                &trial_current_depth_cm,
                &trial_next_depth_cm,
                &trial_gap_cm)) {
            projection_valid = false;
            break;
          }

          if (trial_gap_cm >= min_spacing_cm) {
            safe_dt_h = trial_dt_h;
            event_projected_current_depth_cm = trial_current_depth_cm;
            event_projected_next_depth_cm = trial_next_depth_cm;
          }
          else {
            crossing_dt_h = trial_dt_h;
          }
        }

        if (projection_valid && safe_dt_h > 0.0) {
          event_subtimestep_h = safe_dt_h;
          event_is_interflow_aware = true;
        }
      }
    }

    if (event_subtimestep_h <= LGARTO_TO_PACKET_EVENT_SPLIT_MIN_DT_H ||
        event_subtimestep_h >= limited_subtimestep_h) {
      continue;
    }

    limited_subtimestep_h = event_subtimestep_h;
    limiting_current = current;
    limiting_next = next;
    limiting_current_gap_cm = current_gap_cm;
    limiting_closing_speed_cm_per_h = closing_speed_cm_per_h;
    limiting_projected_current_depth_cm = event_projected_current_depth_cm;
    limiting_projected_next_depth_cm = event_projected_next_depth_cm;
    limiting_current_dzdt_cm_per_h = current_dzdt_cm_per_h;
    limiting_next_dzdt_cm_per_h = next_dzdt_cm_per_h;
    limiting_event_is_interflow_aware = event_is_interflow_aware;
  }

  if (limiting_current == NULL || limiting_next == NULL) {
    if (repeat_front_num != NULL && repeat_next_front_num != NULL &&
        repeat_layer_num != NULL &&
        repeat_count != NULL) {
      *repeat_front_num = -1;
      *repeat_next_front_num = -1;
      *repeat_layer_num = -1;
      *repeat_count = 0;
    }
    return limited_subtimestep_h;
  }

  if (repeat_front_num != NULL && repeat_next_front_num != NULL &&
      repeat_layer_num != NULL && repeat_count != NULL) {
    const bool same_pair =
      *repeat_front_num == limiting_current->front_num &&
      *repeat_next_front_num == limiting_next->front_num;
    const bool same_nearby_layer_cluster =
      *repeat_layer_num == limiting_current->layer_num &&
      limiting_current_gap_cm <= LGARTO_TO_PACKET_EVENT_SPLIT_REPEAT_HANDOFF_MAX_GAP_CM;
    if (same_pair || same_nearby_layer_cluster) {
      (*repeat_count)++;
    }
    else {
      *repeat_count = 1;
    }
    *repeat_front_num = limiting_current->front_num;
    *repeat_next_front_num = limiting_next->front_num;
    *repeat_layer_num = limiting_current->layer_num;

    if (*repeat_count >= LGARTO_TO_PACKET_EVENT_SPLIT_REPEAT_HANDOFF_COUNT &&
        limiting_current_gap_cm <= LGARTO_TO_PACKET_EVENT_SPLIT_REPEAT_HANDOFF_MAX_GAP_CM) {
      // A repeated near-contact split is a stall. Step just past contact so
      // the existing TO/GW merge/correction logic resolves the local overtake.
      limited_subtimestep_h =
        fmin(proposed_subtimestep_h,
             (limiting_current_gap_cm + LGARTO_TO_PACKET_EVENT_SPLIT_HANDOFF_OVERTAKE_CM) /
             limiting_closing_speed_cm_per_h);
      if (verbosity.compare("high") == 0) {
        printf("Adaptive LGARTO event split handoff for repeated mobile TO/GW "
               "packet near-contact: front=%d next_front=%d layer=%d repeats=%d "
               "gap_cm=%.17lf proposed_dt_h=%.17lf handoff_dt_h=%.17lf\n",
               limiting_current->front_num,
               limiting_next->front_num,
               limiting_current->layer_num,
               *repeat_count,
               limiting_current_gap_cm,
               proposed_subtimestep_h,
               limited_subtimestep_h);
      }
      return limited_subtimestep_h;
    }
  }

  if (verbosity.compare("high") == 0) {
    printf("Adaptive LGARTO event split before mobile TO/GW packet overtake: "
           "front=%d next_front=%d layer=%d old_dt_h=%.17lf new_dt_h=%.17lf "
           "depth_cm=%.17lf next_depth_cm=%.17lf projected_depth_cm=%.17lf "
           "projected_next_depth_cm=%.17lf dzdt_cm_per_h=%.17lf "
           "next_dzdt_cm_per_h=%.17lf min_spacing_cm=%.17lf "
           "interflow_aware=%d\n",
           limiting_current->front_num,
           limiting_next->front_num,
           limiting_current->layer_num,
           proposed_subtimestep_h,
           limited_subtimestep_h,
           limiting_current->depth_cm,
           limiting_next->depth_cm,
           limiting_projected_current_depth_cm,
           limiting_projected_next_depth_cm,
           limiting_current_dzdt_cm_per_h,
           limiting_next_dzdt_cm_per_h,
           min_spacing_cm,
           limiting_event_is_interflow_aware ? 1 : 0);
  }

  return limited_subtimestep_h;
}

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

  if (!state->lgar_bmi_params.init_giuh_state_path.empty()) {
    InitializeGIUHRunoffQueueFromCSV(
        state->lgar_bmi_params.init_giuh_state_path.c_str(),
        giuh_runoff_queue,
        num_giuh_ordinates);
  }
  else {
    for (int i=0; i<=num_giuh_ordinates; i++) {
      giuh_runoff_queue[i] = 0.0;
    }
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
    state->lgar_mass_balance.volon_cm       = 0.0;
    state->lgar_mass_balance.volend_cm      = state->lgar_mass_balance.volstart_cm;
    state->lgar_mass_balance.volCRend_cm    = state->lgar_mass_balance.volCRstart_cm;
    state->lgar_mass_balance.volAET_cm      = 0.0;
    state->lgar_mass_balance.volrech_cm     = 0.0;
    state->lgar_mass_balance.volrunoff_cm  += state->lgar_bmi_input_params->precipitation_mm_per_h * mm_to_cm;
    state->lgar_mass_balance.volQ_cm       += state->lgar_bmi_input_params->precipitation_mm_per_h * mm_to_cm;
    state->lgar_mass_balance.volQ_CR_cm     = 0.0;
    state->lgar_mass_balance.volpref_flow_to_CR_cm = 0.0;
    state->lgar_mass_balance.vollgarto_domain_to_CR_cm = 0.0;
    state->lgar_mass_balance.vollateral_flow_cm = 0.0;
    state->lgar_mass_balance.vollateral_flow_timestep_cm = 0.0;
    state->lgar_mass_balance.volPET_cm      = 0.0;
    state->lgar_mass_balance.volrunoff_giuh_cm  = 0.0;
    state->lgar_mass_balance.volchange_calib_cm = 0.0;

    // converted values, a struct local to the BMI and has bmi output variables
    bmi_unit_conv.mass_balance_m        = 0.0;
    bmi_unit_conv.volprecip_timestep_m  = state->lgar_bmi_input_params->precipitation_mm_per_h * mm_to_m;
    bmi_unit_conv.volin_timestep_m      = 0.0;
    bmi_unit_conv.volend_timestep_m     = 0.0;
    bmi_unit_conv.volCRend_timestep_m   = 0.0;
    bmi_unit_conv.volAET_timestep_m     = 0.0;
    bmi_unit_conv.volrech_timestep_m    = 0.0;
    bmi_unit_conv.volrunoff_timestep_m  = state->lgar_bmi_input_params->precipitation_mm_per_h * mm_to_m;
    bmi_unit_conv.volQ_timestep_m       = state->lgar_bmi_input_params->precipitation_mm_per_h * mm_to_m;
    bmi_unit_conv.volQ_CR_timestep_m    = 0.0;
    bmi_unit_conv.volpref_flow_to_CR_timestep_m = 0.0;
    bmi_unit_conv.vollgarto_domain_to_CR_timestep_m = 0.0;
    bmi_unit_conv.vollateral_flow_timestep_m = 0.0;
    bmi_unit_conv.volPET_timestep_m     = 0.0;
    bmi_unit_conv.volrunoff_giuh_timestep_m = 0.0;

    return;
  }
  
  // if lasam is coupled to soil freeze-thaw, frozen fraction module is called
  if (state->lgar_bmi_params.sft_coupled)
    frozen_factor_hydraulic_conductivity(state->lgar_bmi_params);

  double volchange_calib_cm = 0.0;

  if(state->lgar_bmi_params.calib_params_flag) {
    volchange_calib_cm = update_calibratable_parameters(); // change in soil water volume due to calibratable parameters
    state->lgar_bmi_params.calib_params_flag = false;
  }


  // local variables for readibility
  int subcycles;
  int num_layers = state->lgar_bmi_params.num_layers;

  // local variables for a full timestep (i.e., timestep of the forcing data)
  // see 'struct lgar_mass_balance_variables' in all.hxx for full description of the variables
  double precip_timestep_cm   = 0.0;
  double PET_timestep_cm      = 0.0;
  double AET_timestep_cm      = 0.0;
  double volend_timestep_cm   = lgar_calc_mass_bal(state->lgar_bmi_params.cum_layer_thickness_cm, state->head); // this should not be reset to 0.0 in the for loop
  double volCRend_timestep_cm = state->lgar_mass_balance.CR_fast_storage_cm + state->lgar_mass_balance.CR_slow_storage_cm;
  double volin_timestep_cm    = 0.0;
  double volon_timestep_cm    = state->lgar_mass_balance.volon_timestep_cm;
  double volrunoff_timestep_cm      = 0.0;
  double volrech_timestep_cm        = 0.0;
  double volrunoff_giuh_timestep_cm = 0.0;
  double volQ_timestep_cm           = 0.0;
  double volQ_CR_timestep_cm        = 0.0;
  double volpref_flow_to_CR_timestep_cm = 0.0;
  double vollgarto_domain_to_CR_timestep_cm = 0.0;
  double vollateral_flow_timestep_cm = 0.0;
  double mobile_groundwater_explicit_mass_change_timestep_cm = 0.0;
  
  // local variables for a subtimestep (i.e., timestep of the model)
  double precip_subtimestep_cm;
  double precip_subtimestep_cm_per_h;
  double PET_subtimestep_cm;
  double PET_subtimestep_cm_per_h;
  double ponded_depth_subtimestep_cm;
  double AET_subtimestep_cm;
  double volstart_subtimestep_cm;
  double volend_subtimestep_cm = volend_timestep_cm; // this should not be reset to 0.0 in the for loop
  double volin_subtimestep_cm;
  double volon_subtimestep_cm;
  double volrunoff_subtimestep_cm;
  double volrech_subtimestep_cm;
  double lower_boundary_flux_for_cache_subtimestep_cm = 0.0;
  double lower_boundary_flux_for_mobile_groundwater_subtimestep_cm = 0.0;
  double lateral_flow_subtimestep_cm;
  double precip_previous_subtimestep_cm;
  double volCRstart_subtimestep_cm;
  double volCRend_subtimestep_cm = volCRend_timestep_cm;
  
  double subtimestep_h = state->lgar_bmi_params.timestep_h;
  int nint = state->lgar_bmi_params.nint;
  double wilting_point_psi_cm = state->lgar_bmi_params.wilting_point_psi_cm;
  double field_capacity_psi_cm = state->lgar_bmi_params.field_capacity_psi_cm;
  double a = state->lgar_bmi_params.a;
  double b = state->lgar_bmi_params.b;
  double frac_to_CR = state->lgar_bmi_params.frac_to_CR;
  double a_slow = state->lgar_bmi_params.a_slow;
  double b_slow = state->lgar_bmi_params.b_slow;
  double frac_slow = state->lgar_bmi_params.frac_slow;
  double CR_fast_discharge_threshold_cm =
    state->lgar_bmi_params.CR_fast_discharge_threshold_cm;
  double CR_slow_discharge_threshold_cm =
    state->lgar_bmi_params.CR_slow_discharge_threshold_cm;
  double spf_factor = state->lgar_bmi_params.spf_factor;
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
    if (state->lgar_bmi_input_params->precipitation_mm_per_h > 10.0 || volon_timestep_cm > 0.0 ) {
      subtimestep_h = state->lgar_bmi_params.minimum_timestep_h;  //case where precip > 1 cm/h, or there is ponded head from the last time step
    }
    else if (state->lgar_bmi_input_params->precipitation_mm_per_h > 0.0) {
      subtimestep_h = state->lgar_bmi_params.minimum_timestep_h * 2.0;  //case where precip is less than 1 cm/h but greater than 0, and there is no ponded head 
    }
    subtimestep_h = fmin(subtimestep_h, state->lgar_bmi_params.forcing_resolution_h);  //just in case the user has specified a minimum time step that would make the subtimestep_h greater than the forcing resolution 
    state->lgar_bmi_params.timestep_h = subtimestep_h;
  }

  bool caching_at_start = state->lgar_mass_balance.cache_fluxes;
  bool switch_caching = FALSE;

  const double max_abs_mobile_dzdt_for_cache =
    lgar_max_abs_mobile_dzdt_for_flux_cache(state->head);
  const double cached_lower_boundary_flux_budget_cm =
    fmax(SMALL_EPS, CACHE_LOWER_BOUNDARY_MBAL_FRACTION * mbal_tol);
  const bool cached_lower_boundary_flux_budget_exceeded =
    (fabs(state->lgar_mass_balance.accumulated_lower_boundary_flux_cm) +
     fabs(state->lgar_mass_balance.previous_lower_boundary_flux_cm)) >=
    cached_lower_boundary_flux_budget_cm;

  if (state->lgar_bmi_params.allow_flux_caching){
    //The idea here is that, during dry periods, AET will become a small fraction of PET and wetting fronts will be very slow moving. In these cases, it is not necessary to compute fluxes for every time step.
    //To save on runtime, and if allow_flux_caching is set to true in the config file, we simply cache computed fluxes to be used for subsequent time steps rather than recomputing them.
    //The current implementation is to not move the wetting fronts under these conditions, but then move them more rapidly once it is time to calculate fluxes again. Also, during these periods, PET will be 0 but made to be larger to conserve mass when it is time to recalculate fluxes.
    //The signed lower-boundary flux is copied from the last timestep, then applied as a cached correction when dynamics are recomputed.
    //There is very little change to the simulation when this is enabled.
    //Note that for NextGen models, it is ultimately desirable that there is technically output for every hour, so simply relaxing the adaptive time step to be coarser than 1 hour isn't the best solution
    if (subtimestep_h == state->lgar_bmi_params.forcing_resolution_h){
      if ( (state->lgar_bmi_input_params->precipitation_mm_per_h < PRECIP_THRESHOLD_MM_PER_H) && ( (state->lgar_mass_balance.previous_AET / (state->lgar_mass_balance.previous_PET + PET_EPSILON )) < AET_PET_RATIO_THRESHOLD) && (state->lgar_bmi_params.cache_count!=NUM_TIMESTEPS_BEFORE_RESET_CACHE) && 
           (volon_timestep_cm<VOLON_TIMESTEP_THRESHOLD_CM) && (fabs(state->lgar_mass_balance.previous_lower_boundary_flux_cm)<BOTTOM_BDY_FLUX_THRESHOLD_CM) ){
        state->lgar_mass_balance.cache_fluxes = TRUE;
      }
      if (max_abs_mobile_dzdt_for_cache > THRESHOLD_DZDT_CM_PER_H ){
        state->lgar_mass_balance.cache_fluxes = FALSE;
      }
      if (cached_lower_boundary_flux_budget_exceeded){
        state->lgar_mass_balance.cache_fluxes = FALSE;
      }
      if (state->lgar_bmi_params.lateral_flow_enabled){
        state->lgar_mass_balance.cache_fluxes = FALSE;
      }
      if (volon_timestep_cm>=VOLON_TIMESTEP_THRESHOLD_CM){
        state->lgar_mass_balance.cache_fluxes = FALSE;
      }
      if ( (state->lgar_bmi_params.cache_count==NUM_TIMESTEPS_BEFORE_RESET_CACHE ) || (state->lgar_bmi_input_params->precipitation_mm_per_h >= PRECIP_THRESHOLD_MM_PER_H) ){
        state->lgar_mass_balance.cache_fluxes = FALSE;
      }
    }
    else {
      state->lgar_mass_balance.cache_fluxes = FALSE;
    }
  }

  if (state->lgar_bmi_params.timesteps < 1){ // decision happens before the current timestep counter is incremented
    state->lgar_mass_balance.cache_fluxes = FALSE;
  }
  if (caching_at_start && !state->lgar_mass_balance.cache_fluxes){
    switch_caching = TRUE;//if you switch from cached to not, you need to add the "missing" PET back into the mass balance and AET calculation 
  }

  if (state->lgar_mass_balance.cache_fluxes){
    state->lgar_bmi_params.cache_count ++;
  }

  const double base_subtimestep_h = state->lgar_bmi_params.timestep_h;
  state->lgar_bmi_params.forcing_interval = int(state->lgar_bmi_params.forcing_resolution_h/base_subtimestep_h+1.0e-08); // add 1.0e-08 to prevent truncation error
  subcycles = state->lgar_bmi_params.forcing_interval;

  if (verbosity.compare("high") == 0) {
    printf("time step size in hours: %lf \n", state->lgar_bmi_params.timestep_h);
  }

  // ensure precip and PET are non-negative
  if (state->lgar_bmi_input_params->precipitation_mm_per_h < 0.0) {
      std::cerr << "Warning: Pr [mm/h] (timestep) is negative ("
                << state->lgar_bmi_input_params->precipitation_mm_per_h
                << "), setting to 0.\n";
      state->lgar_bmi_input_params->precipitation_mm_per_h = 0.0;
  }

  if (state->lgar_bmi_input_params->PET_mm_per_h < 0.0) {
      std::cerr << "Warning: PET [mm/h] (timestep) is negative ("
                << state->lgar_bmi_input_params->PET_mm_per_h
                << "), setting to 0.\n";
      state->lgar_bmi_input_params->PET_mm_per_h = 0.0;
  }

  if (PET_affects_precip){ // if the user wants PET subtracted from precip
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
  
  // subcycling loop (loop over model's timestep). LGARTO can add dynamic event
  // splits at surface layer crossings and mobile TO/GW packet contacts.
  double remaining_forcing_h = subcycles * base_subtimestep_h;
  int cycle = 0;
  int lgarto_event_splits_this_forcing = 0;
  int repeated_mobile_TO_packet_front_num = -1;
  int repeated_mobile_TO_packet_next_front_num = -1;
  int repeated_mobile_TO_packet_layer_num = -1;
  int repeated_mobile_TO_packet_split_count = 0;
  while (remaining_forcing_h > SMALL_EPS) {
    cycle++;
    subtimestep_h = fmin(base_subtimestep_h, remaining_forcing_h);
    if (state->lgar_bmi_params.TO_enabled && !state->lgar_mass_balance.cache_fluxes) {
      const double groundwater_depth_for_event_split_cm =
        state->lgar_bmi_params.mobile_groundwater_level
          ? lgar_effective_groundwater_depth_cm(&state->lgar_bmi_params)
          : -1.0;
      double event_limited_subtimestep_h =
        lgarto_limit_subtimestep_for_surface_boundary_crossing(
          subtimestep_h,
          groundwater_depth_for_event_split_cm,
          state->lgar_bmi_params.num_layers,
          state->lgar_bmi_params.cum_layer_thickness_cm,
          state->head);
      event_limited_subtimestep_h =
        lgarto_limit_subtimestep_for_mobile_TO_packet_overtake(
          event_limited_subtimestep_h,
          state,
          &repeated_mobile_TO_packet_front_num,
          &repeated_mobile_TO_packet_next_front_num,
          &repeated_mobile_TO_packet_layer_num,
          &repeated_mobile_TO_packet_split_count);

      if (event_limited_subtimestep_h + SMALL_EPS < subtimestep_h) {
        lgarto_event_splits_this_forcing++;
        if (lgarto_event_splits_this_forcing > LGARTO_EVENT_SPLIT_MAX_PER_FORCING) {
          fprintf(stderr,
                  "Error: adaptive LGARTO event splitting exceeded split cap.\n"
                  "  splits=%d max_splits=%d base_subtimestep_h=%.17lf "
                  "remaining_forcing_h=%.17lf requested_subtimestep_h=%.17lf "
                  "event_limited_subtimestep_h=%.17lf\n"
                  "  Wetting front list follows:\n",
                  lgarto_event_splits_this_forcing,
                  LGARTO_EVENT_SPLIT_MAX_PER_FORCING,
                  base_subtimestep_h,
                  remaining_forcing_h,
                  subtimestep_h,
                  event_limited_subtimestep_h);
          fflush(stderr);
          listPrint(state->head);
          fflush(stdout);
          abort();
        }

        subtimestep_h = event_limited_subtimestep_h;
      }
    }

    if (subtimestep_h <= 0.0) {
      fprintf(stderr,
              "Error: non-positive subtimestep selected in BMI update.\n"
              "  subtimestep_h=%.17lf remaining_forcing_h=%.17lf base_subtimestep_h=%.17lf\n",
              subtimestep_h,
              remaining_forcing_h,
              base_subtimestep_h);
      abort();
    }

    state->lgar_bmi_params.timestep_h = subtimestep_h;

    bool top_near_sat = false;
    this->state->lgar_bmi_params.time_s    += subtimestep_h * state->units.hr_to_sec;
    this->state->lgar_bmi_params.timesteps ++;

    double precip_for_CR_subtimestep_cm_per_h = 0.0;
    precip_subtimestep_cm_per_h = state->lgar_bmi_input_params->precipitation_mm_per_h * mm_to_cm; // rate [cm/hour]
    
    if (verbosity.compare("high") == 0 || verbosity.compare("low") == 0) {
      std::cerr<<"BMI Update |---------------------------------------------------------------|\n";
      std::cerr<<"BMI Update |Timesteps = "<< state->lgar_bmi_params.timesteps<<", Time [h] = "<<this->state->lgar_bmi_params.time_s / 3600.<<", Subcycle = "<< cycle <<" of "<<subcycles;
      if (cycle > subcycles) {
        std::cerr<<" (adaptive LGARTO event split)";
      }
      std::cerr<<std::endl;
    }

    if( state->state_previous != NULL ){
      listDelete(state->state_previous);
      state->state_previous = NULL;
    }
    state->state_previous = listCopy(state->head);

    double ponded_flux_for_CR = 0.0;

    // allocates some water to conceptual reservoir storage via conditional preferential flow
    if (state->lgar_bmi_params.runoff_in_prev_step){
      double precip_subtimestep_cm_per_h_total = precip_subtimestep_cm_per_h;
      precip_for_CR_subtimestep_cm_per_h = frac_to_CR * precip_subtimestep_cm_per_h_total;
      precip_subtimestep_cm_per_h = (1.0 - frac_to_CR) * precip_subtimestep_cm_per_h_total;

      ponded_flux_for_CR = volon_timestep_cm/subtimestep_h*frac_to_CR;
      volon_timestep_cm = volon_timestep_cm*(1 - frac_to_CR);
    }

    /* Note unit conversion:
       Pr and PET are rates (fluxes) in mm/h
       Pr [mm/h] * 1h/3600sec = Pr [mm/3600sec]
       Model timestep (dt) = 300 sec (5 minutes for example)
       convert rate to amount
       Pr [mm/3600sec] * dt [300 sec] = Pr[mm] * 300/3600.
       in the code below, subtimestep_h is this 300/3600 factor (see initialize from config in lgar.cxx)
    */

    AET_subtimestep_cm            = 0.0;
    volstart_subtimestep_cm       = 0.0;
    volin_subtimestep_cm          = 0.0;
    volrunoff_subtimestep_cm      = 0.0;
    volrech_subtimestep_cm        = 0.0;
    lower_boundary_flux_for_cache_subtimestep_cm = 0.0;
    double temp_rch               = 0.0; //handles case when a fraction of a wetting front technically crosses the lower boundary of the vadose zone
    double creation_excess_gw_flux_subtimestep_cm = 0.0;
    double creation_excess_runoff_subtimestep_cm = 0.0;
    double free_drainage_subtimestep_cm = 0.0;
    lateral_flow_subtimestep_cm = 0.0;
    double lower_boundary_flux_for_CR = 0.0;
    double lower_boundary_CR_exchange_cm = 0.0;
    double mobile_groundwater_explicit_mass_change_subtimestep_cm = 0.0;
    volCRstart_subtimestep_cm =
      state->lgar_mass_balance.CR_fast_storage_cm + state->lgar_mass_balance.CR_slow_storage_cm;

    PET_subtimestep_cm_per_h = state->lgar_bmi_input_params->PET_mm_per_h * mm_to_cm;

    ponded_depth_subtimestep_cm = precip_subtimestep_cm_per_h * subtimestep_h; // the amount of water on the surface before any infiltration and runoff

	    ponded_depth_subtimestep_cm += volon_timestep_cm; // add volume of water on the surface (from the last timestep) to ponded depth as well

	    precip_subtimestep_cm = precip_subtimestep_cm_per_h * subtimestep_h; // rate x dt = amount (portion of the water on the suface for model's timestep [cm])
	    PET_subtimestep_cm = PET_subtimestep_cm_per_h * subtimestep_h;      // potential ET for this subtimestep [cm]

    volstart_subtimestep_cm = lgar_calc_mass_bal(state->lgar_bmi_params.cum_layer_thickness_cm, state->head);
	    const double trace_ponded_start_cm = ponded_depth_subtimestep_cm;
	    const int trace_fronts_start = listLength(state->head);
	    const int trace_surface_fronts_start = listLength_surface(state->head);

	    if (!state->lgar_mass_balance.cache_fluxes){

      //this code makes sure that AET or free drainage will not be extracted in a way that would result in an impossible storage
      double min_storage = 0.0;
      double mass_used_to_check_impossible_storages = lgar_calc_mass_bal(state->lgar_bmi_params.cum_layer_thickness_cm, state->head);
      for (int k = 1; k < num_layers+1; k++) {
        int layer_num_min_check = k;
        int soil_num_min_check = state->lgar_bmi_params.layer_soil_type[layer_num_min_check];
        min_storage += state->soil_properties[soil_num_min_check].theta_r * (state->lgar_bmi_params.cum_layer_thickness_cm[k]-state->lgar_bmi_params.cum_layer_thickness_cm[k-1]);
      }

      int wf_free_drainage_demand =
        state->lgar_bmi_params.mobile_groundwater_level
          ? wetting_front_free_drainage_mobile_groundwater(state->head)
          : wetting_front_free_drainage(state->head);

      double min_water_possible_for_FD_WF = calc_min_water_possible_for_free_drainage_wetting_front(wf_free_drainage_demand,  &state->head, state->lgar_bmi_params.layer_soil_type, state->soil_properties);
      double storage_in_FD_WF = calc_storage_in_free_drainage_wetting_front(wf_free_drainage_demand, &state->head);

      double cached_lower_boundary_flux_correction_cm = 0.0;

      if (switch_caching){
        PET_subtimestep_cm_per_h += state->lgar_mass_balance.accumulated_PET_cm / subtimestep_h;
        state->lgar_mass_balance.accumulated_PET_cm = 0.0;
        cached_lower_boundary_flux_correction_cm += state->lgar_mass_balance.accumulated_lower_boundary_flux_cm;
        state->lgar_mass_balance.accumulated_lower_boundary_flux_cm = 0.0;
        if (verbosity.compare("high") == 0) {
          printf("flux caching increments PET \n");
          printf("PET_subtimestep_cm_per_h: %lf \n", PET_subtimestep_cm_per_h);
        }
      }

      if (state->lgar_bmi_params.free_drainage_enabled){
        struct wetting_front *front = listFindFront(listLength(state->head), state->head, NULL);
        free_drainage_subtimestep_cm += subtimestep_h*front->K_cm_per_h;
        int iter_mass_check_FD = 0;
        while ( (mass_used_to_check_impossible_storages - free_drainage_subtimestep_cm < min_storage) || (storage_in_FD_WF - free_drainage_subtimestep_cm < min_water_possible_for_FD_WF) ){
          free_drainage_subtimestep_cm *= 0.5; //give it a chance to merely become smaller before setting to 0
          if (iter_mass_check_FD > 5){
            free_drainage_subtimestep_cm = 0.0;
            break;
          }
          iter_mass_check_FD ++;
        }
        if (free_drainage_subtimestep_cm<1.E-7){
          free_drainage_subtimestep_cm = 0.0;
        }
        if (front->psi_cm>1.E6){
          free_drainage_subtimestep_cm = 0.0;
        }
      }

      if ((state->lgar_bmi_params.free_drainage_enabled) && verbosity.compare("high") == 0){
        printf("free_drainage_subtimestep_cm: %.10lf \n", free_drainage_subtimestep_cm);
      }

      //using cerr instead of cout due to some cout buffering issues when running in the ngen framework, cerr doesn't buffer so it prints immediately to the sreeen.
      if (verbosity.compare("high") == 0 || verbosity.compare("low") == 0) {

        std::cerr<<"Pr [cm/h], Pr [cm] (subtimestep), subtimestep [h] = "<<state->lgar_bmi_input_params->precipitation_mm_per_h * mm_to_cm <<", "<< precip_subtimestep_cm <<", "<< subtimestep_h<<" ("<<subtimestep_h*3600<<" sec)"<<"\n";
        std::cerr<<"PET [cm/h], PET [cm] (subtimestep) = "<<state->lgar_bmi_input_params->PET_mm_per_h * mm_to_cm <<", "<< PET_subtimestep_cm<<"\n";
      }

      precip_previous_subtimestep_cm = state->lgar_bmi_params.precip_previous_timestep_cm; // creation of a new wetting front depends on previous timestep's rainfall

      num_layers = state->lgar_bmi_params.num_layers;
      double delta_theta;   // the width of a front, such that its volume=depth*delta_theta
      double dry_depth;
      double surf_frac_rz = 0.0;
      std::vector<double> surf_AET_vec(listLength(state->head) + 1, 0.0);
      double groundwater_supported_AET_potential_subtimestep_cm = 0.0;
      const char *trace_branch = "dynamic_no_insert";
      bool trace_create_surficial_front = false;
      bool trace_is_top_wf_saturated = false;
      double trace_projected_TO_storage_release_cm = 0.0;
      double trace_insertion_storage_release_cm = 0.0;
      double trace_surface_creation_gw_capacity_cm = 0.0;
      double trace_surface_creation_post_creation_TO_release_cm = 0.0;
      double trace_insert_raw_fp_cm_per_h = NAN;
      double trace_insert_storage_limit_fp_cm_per_h = NAN;
      double trace_insert_capped_fp_cm_per_h = NAN;

      // Calculate AET from PET if PET is non-zero
      if (PET_subtimestep_cm_per_h > 0.0) {
        const double PET_budget_subtimestep_cm =
          PET_subtimestep_cm_per_h * subtimestep_h;
        const double explicit_aet_lower_boundary_depth_cm =
          lgar_explicit_aet_lower_boundary_depth_cm(&state->lgar_bmi_params);
        groundwater_supported_AET_potential_subtimestep_cm =
          lgar_groundwater_supported_aet_potential_cm(
            &state->lgar_bmi_params, PET_budget_subtimestep_cm);
        if (!state->lgar_bmi_params.TO_enabled) {
          AET_subtimestep_cm = calc_aet(PET_subtimestep_cm_per_h, subtimestep_h, wilting_point_psi_cm,
                                        field_capacity_psi_cm, state->lgar_bmi_params.layer_soil_type,
                                        AET_thresh_Theta, AET_expon, state->head, state->soil_properties);
        }
        else {
          calc_aet(state->lgar_bmi_params.TO_enabled, PET_subtimestep_cm_per_h, subtimestep_h,
                   wilting_point_psi_cm, field_capacity_psi_cm,
                   state->lgar_bmi_params.root_zone_depth_cm, &surf_frac_rz,
                   state->lgar_bmi_params.layer_soil_type, AET_thresh_Theta, AET_expon,
                   state->head, state->soil_properties, surf_AET_vec.data(),
                   explicit_aet_lower_boundary_depth_cm);
          AET_subtimestep_cm = 0.0;
        }
      }

      if (!state->lgar_bmi_params.TO_enabled) {
        int iter_mass_check_AET = 0;
        while ( (mass_used_to_check_impossible_storages - AET_subtimestep_cm < min_storage) || (storage_in_FD_WF - AET_subtimestep_cm < min_water_possible_for_FD_WF) ){
          AET_subtimestep_cm *= 0.5; //give it a chance to merely become smaller before setting to 0
          if (iter_mass_check_AET > 5){
            AET_subtimestep_cm = 0.0;
            break;
          }
          iter_mass_check_AET ++;
        }

        int iter_mass_check_AET_and_FD = 0;
        while ( ( (mass_used_to_check_impossible_storages - AET_subtimestep_cm - free_drainage_subtimestep_cm - cached_lower_boundary_flux_correction_cm) < min_storage) || ( (storage_in_FD_WF - AET_subtimestep_cm - free_drainage_subtimestep_cm - cached_lower_boundary_flux_correction_cm) < min_water_possible_for_FD_WF) ){
          // both should also be checked at the same because while individually these might not make an impossible storage, together they might
          AET_subtimestep_cm *= 0.5;
          free_drainage_subtimestep_cm *= 0.5;
          cached_lower_boundary_flux_correction_cm *=0.5;
          if (iter_mass_check_AET_and_FD > 5){
            AET_subtimestep_cm = 0.0;
            free_drainage_subtimestep_cm = 0.0;
            cached_lower_boundary_flux_correction_cm = 0.0;
            break;
          }
          iter_mass_check_AET_and_FD ++;
        }
      }

      // precip_timestep_cm += precip_subtimestep_cm;
      precip_timestep_cm += precip_subtimestep_cm + precip_for_CR_subtimestep_cm_per_h*subtimestep_h;
      PET_timestep_cm += fmax(PET_subtimestep_cm,0.0); // ensures non-negative PET

      //addressed machine precision issues where volon_timestep_error could be for example -1E-17 or 1.E-20 or smaller
      volon_timestep_cm = fmax(volon_timestep_cm,0.0);
      volon_timestep_cm = volon_timestep_cm > SMALL_EPS ? volon_timestep_cm : 0.0;

      /*----------------------------------------------------------------------*/
      // Should a new wetting front be created?
	      int soil_num = state->lgar_bmi_params.layer_soil_type[state->head->layer_num];
	      double theta_e = state->soil_properties[soil_num].theta_e;
	      bool is_top_wf_saturated = false;
	      struct wetting_front *top_wf_for_preferential_flow = state->head;
	      double theta_e_for_preferential_flow = theta_e;
	      if (!state->lgar_bmi_params.TO_enabled) {
	        is_top_wf_saturated = (state->head->theta + SMALL_EPS) >= theta_e;
	      }
	      else {
	        if (listLength_surface(state->head) > 0) {
	          struct wetting_front *top_most_surface_WF = state->head;
	          while (top_most_surface_WF != NULL && top_most_surface_WF->is_WF_GW) {
	            top_most_surface_WF = top_most_surface_WF->next;
	          }

	          if (top_most_surface_WF != NULL) {
	            const int soil_num_highest_surf =
	              state->lgar_bmi_params.layer_soil_type[top_most_surface_WF->layer_num];
	            const double theta_e_highest_surf = state->soil_properties[soil_num_highest_surf].theta_e;
	            is_top_wf_saturated = (top_most_surface_WF->theta + SMALL_EPS) >= theta_e_highest_surf;
	            top_wf_for_preferential_flow = top_most_surface_WF;
	            theta_e_for_preferential_flow = theta_e_highest_surf;
	          }
	        }
	        else {
	          is_top_wf_saturated = (state->head->theta + SMALL_EPS) >= theta_e;
	        }
	      }
	      double theta_above_which_precip_contribs_to_GW =
	        theta_e_for_preferential_flow * spf_factor;
	      top_near_sat =
	        top_wf_for_preferential_flow != NULL &&
	        top_wf_for_preferential_flow->theta > theta_above_which_precip_contribs_to_GW; //is the top WF near saturation, thus triggering simple preferential flow if enabled

      // checks on creatign a new surficial front
      // 1. check current and previous timestep precipitation
      // bool create_surficial_front = (precip_previous_subtimestep_cm == 0.0 && precip_subtimestep_cm > 0.0);
	      bool create_surficial_front =
	        (precip_previous_subtimestep_cm == 0.0 && precip_subtimestep_cm > 0.0 && volon_timestep_cm == 0.0) ||
	        ((precip_subtimestep_cm > 0.0 || volon_timestep_cm > 0.0) &&
	         (listLength(state->head) == num_layers) &&
	         !(state->lgar_bmi_params.TO_enabled));

	      int mobile_active_surface_front_count = listLength_surface(state->head);
	      if (state->lgar_bmi_params.mobile_groundwater_level &&
	          state->lgar_bmi_params.TO_enabled) {
	        mobile_active_surface_front_count = 0;
	        for (struct wetting_front *front = state->head;
	             front != NULL;
	             front = front->next) {
	          if (front->is_WF_GW == FALSE && front->to_bottom == FALSE) {
	            mobile_active_surface_front_count++;
	          }
	        }
	      }
	      
	      // 2. check soil top wetting front condition (saturated/unsaturated), and surface ponded water
	      if (is_top_wf_saturated)
	        create_surficial_front = false;

	      if (state->lgar_bmi_params.TO_enabled &&
	          (precip_subtimestep_cm > 0.0 || volon_timestep_cm > 0.0) &&
	          (mobile_active_surface_front_count == 0) &&
	          (state->head->theta < (theta_e - SMALL_EPS))) {
	        create_surficial_front = true;
	      }
	      trace_create_surficial_front = create_surficial_front;
	      trace_is_top_wf_saturated = is_top_wf_saturated;

	      if (verbosity.compare("high") == 0 || verbosity.compare("low") == 0) {
        std::string flag        = (create_surficial_front && !is_top_wf_saturated) == true ? "Yes" : "No";
        std::string flag_top_wf = is_top_wf_saturated == true ? "Yes" : "No";
        std::cerr<<"Is top wetting front saturated? "<< flag_top_wf  << "\n";
        std::cerr<<"Create superficial wetting front? "<< flag << "\n";
      }

      /*----------------------------------------------------------------------*/
      /* create a new wetting front if the following is true. Meaning there is no
        wetting front in the top layer to accept the water, must create one. */
	      if(create_surficial_front) {
	        trace_branch = "create_surficial_front";

        double temp_pd = 0.0; // necessary to assign zero precip due to the creation of new wetting front; AET will still be taken out of the layers

        // move the wetting fronts without adding any water; this is done to close the mass balance
        // and also to merge / cross if necessary 
        lgar_limit_upward_TO_dzdt_by_CR_supply(state, subtimestep_h);
        temp_rch = lgar_move_wetting_fronts(subtimestep_h, &free_drainage_subtimestep_cm, &temp_pd, wf_free_drainage_demand, volend_subtimestep_cm, cached_lower_boundary_flux_correction_cm,
              num_layers, &AET_subtimestep_cm, state->lgar_bmi_params.cum_layer_thickness_cm,
              state->lgar_bmi_params.layer_soil_type, state->lgar_bmi_params.frozen_factor,
	              &state->head, state->state_previous, state->soil_properties,
		              state->lgar_bmi_params.TO_enabled ? surf_AET_vec.data() : nullptr,
		              PET_subtimestep_cm_per_h, wilting_point_psi_cm, field_capacity_psi_cm,
		              state->lgar_bmi_params.root_zone_depth_cm, surf_frac_rz,
		              state->lgar_bmi_params.mbal_tol,
		              &lateral_flow_subtimestep_cm,
		              state->lgar_bmi_params.lateral_flow_psi_threshold_cm,
		              state->lgar_bmi_params.lateral_flow_factor,
		              state->lgar_bmi_params.TO_enabled
		                ? lgar_effective_groundwater_depth_cm(&state->lgar_bmi_params)
		                : lgar_fixed_soil_depth_cm(&state->lgar_bmi_params),
		              state->lgar_bmi_params.mobile_groundwater_level);

	        // if (temp_pd != 0.0){ //if temp_pd != 0.0, that means that some water left the model through the lower model bdy. For LGARTO preparation, this has been refactored such that temp_rch handles this now.
	        //   // volrech_subtimestep_cm = temp_pd;
	        //   // volrech_timestep_cm += volrech_subtimestep_cm;
	        //   // temp_pd = 0.0;
	        // }

	        if (state->lgar_bmi_params.TO_enabled) {
	          const double already_booked_lower_boundary_cm =
	            temp_rch + free_drainage_subtimestep_cm;
	          trace_surface_creation_gw_capacity_cm =
	            lgarto_remaining_lower_boundary_capacity_cm(
	              subtimestep_h, num_layers, state->lgar_bmi_params.cum_layer_thickness_cm,
	              state->head, already_booked_lower_boundary_cm);
	          trace_surface_creation_gw_capacity_cm =
	            fmax(trace_surface_creation_gw_capacity_cm,
	                 lgarto_saturated_creation_lower_boundary_capacity_cm(
	                   subtimestep_h, num_layers,
	                   state->lgar_bmi_params.cum_layer_thickness_cm,
	                   state->lgar_bmi_params.layer_soil_type,
	                   state->lgar_bmi_params.frozen_factor,
	                   state->head, state->soil_properties,
	                   already_booked_lower_boundary_cm));
	        }

	        // depth of the surficial front to be created
		        dry_depth = lgar_calc_dry_depth(state->lgar_bmi_params.TO_enabled, use_closed_form_G, nint, subtimestep_h, &delta_theta, state->lgar_bmi_params.layer_soil_type,
	                state->lgar_bmi_params.cum_layer_thickness_cm, state->lgar_bmi_params.frozen_factor,
	                state->head, state->soil_properties);

        if (verbosity.compare("high") == 0) {
          printf("State before moving creating new WF...\n");
          listPrint(state->head);
        }
        
	        double theta_for_new_wf = state->head->theta;
	        struct wetting_front *top_most_surface_WF = state->head;
	        if (listLength_surface(state->head) > 0) {
	          while (top_most_surface_WF != NULL && top_most_surface_WF->is_WF_GW) {
	            top_most_surface_WF = top_most_surface_WF->next;
	          }
	        }
	        if (top_most_surface_WF != NULL && top_most_surface_WF->depth_cm == 0.0) {
	          while (top_most_surface_WF != NULL && top_most_surface_WF->depth_cm == 0.0) {
	            top_most_surface_WF = top_most_surface_WF->next;
	          }
	        }
	        if (top_most_surface_WF != NULL) {
	          theta_for_new_wf = top_most_surface_WF->theta;
	        }

          double creation_excess_gw_flux_cm = 0.0;
          double creation_excess_runoff_cm = 0.0;
		        lgar_create_surficial_front(state->lgar_bmi_params.TO_enabled, num_layers, &ponded_depth_subtimestep_cm, &volin_subtimestep_cm, dry_depth, theta_for_new_wf,
		            state->lgar_bmi_params.layer_soil_type, state->lgar_bmi_params.cum_layer_thickness_cm,
		            state->lgar_bmi_params.frozen_factor, &state->head, state->soil_properties,
		              &creation_excess_gw_flux_cm, &creation_excess_runoff_cm,
		              trace_surface_creation_gw_capacity_cm,
		              state->lgar_bmi_params.TO_enabled
		                ? lgar_effective_groundwater_depth_cm(&state->lgar_bmi_params)
		                : lgar_fixed_soil_depth_cm(&state->lgar_bmi_params));
          creation_excess_gw_flux_subtimestep_cm += creation_excess_gw_flux_cm;
          creation_excess_runoff_subtimestep_cm += creation_excess_runoff_cm;

        if (verbosity.compare("high") == 0) {
          printf("State after moving creating new WF...\n");
          listPrint(state->head);
        }

        if(state->state_previous != NULL ){
          listDelete(state->state_previous);
          state->state_previous = NULL;
        }
        state->state_previous = listCopy(state->head);

        // volin_timestep_cm += volin_subtimestep_cm;

        if (verbosity.compare("high") == 0) {
    std::cerr<<"New wetting front created...\n";
    listPrint(state->head);
        }
      }

      /*----------------------------------------------------------------------*/
      /* infiltrate water based on the infiltration capacity given no new wetting front
        is created and that there is water on the surface (or raining). */

	      if (ponded_depth_subtimestep_cm > 0 && !create_surficial_front) {
	        trace_branch = "lgar_insert_water";
	        double insertion_storage_release_subtimestep_cm = free_drainage_subtimestep_cm;
	        if (state->lgar_bmi_params.TO_enabled) {
          // LGAR includes same-substep free drainage when deciding whether storage is
          // available for infiltration. LGARTO's TO drainage happens later in the
          // substep, so add a read-only positive projection here to avoid routing
          // water to runoff just before TO motion creates room.
          trace_projected_TO_storage_release_cm =
            lgarto_project_TO_motion_lower_boundary_flux_cm(subtimestep_h, num_layers,
                                                            state->lgar_bmi_params.cum_layer_thickness_cm,
                                                            state->lgar_bmi_params.layer_soil_type,
                                                            state->head, state->soil_properties,
                                                            lgar_effective_groundwater_depth_cm(&state->lgar_bmi_params));
	          insertion_storage_release_subtimestep_cm += trace_projected_TO_storage_release_cm;
	        }
	        trace_insertion_storage_release_cm = insertion_storage_release_subtimestep_cm;

	        volrunoff_subtimestep_cm = lgar_insert_water(use_closed_form_G, nint, subtimestep_h, AET_subtimestep_cm, insertion_storage_release_subtimestep_cm, &ponded_depth_subtimestep_cm,
	                &volin_subtimestep_cm, precip_subtimestep_cm_per_h,
	                wf_free_drainage_demand, num_layers,
	                ponded_depth_max_cm, state->lgar_bmi_params.layer_soil_type,
	                state->lgar_bmi_params.cum_layer_thickness_cm,
	                state->lgar_bmi_params.frozen_factor, state->head,
		                state->soil_properties,
		                &trace_insert_raw_fp_cm_per_h,
		                &trace_insert_storage_limit_fp_cm_per_h,
			                &trace_insert_capped_fp_cm_per_h,
		                state->lgar_bmi_params.mobile_groundwater_level);
        // volin_timestep_cm += volin_subtimestep_cm;
        // volrunoff_timestep_cm += volrunoff_subtimestep_cm;
        
        volrech_subtimestep_cm = volin_subtimestep_cm;
        volon_subtimestep_cm = ponded_depth_subtimestep_cm;
        if (volrunoff_subtimestep_cm < 0) {
          std::cerr<<"Runoff is less than 0, which should not happen.\n";
          abort();
        }
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
          volon_subtimestep_cm = ponded_depth_max_cm;
          ponded_depth_subtimestep_cm = ponded_depth_max_cm;
        }
      }
      /*----------------------------------------------------------------------*/

      /* move wetting fronts if no new wetting front is created. Otherwise, movement
        of wetting fronts has already happened at the time of creating surficial front,
        so no need to move them here. */
      if (!create_surficial_front) {
        double volin_subtimestep_cm_temp = volin_subtimestep_cm;  /* passing this for mass balance only, the method modifies it
                    and returns percolated value, so we need to keep its original
                    value stored to copy it back*/
        lgar_limit_upward_TO_dzdt_by_CR_supply(state, subtimestep_h);
	        temp_rch = lgar_move_wetting_fronts(subtimestep_h, &free_drainage_subtimestep_cm, &volin_subtimestep_cm, wf_free_drainage_demand, volend_subtimestep_cm, cached_lower_boundary_flux_correction_cm,
	              num_layers, &AET_subtimestep_cm, state->lgar_bmi_params.cum_layer_thickness_cm,
	              state->lgar_bmi_params.layer_soil_type, state->lgar_bmi_params.frozen_factor,
		              &state->head, state->state_previous, state->soil_properties,
			              state->lgar_bmi_params.TO_enabled ? surf_AET_vec.data() : nullptr,
		              PET_subtimestep_cm_per_h, wilting_point_psi_cm, field_capacity_psi_cm,
		              state->lgar_bmi_params.root_zone_depth_cm, surf_frac_rz,
		              state->lgar_bmi_params.mbal_tol,
		              &lateral_flow_subtimestep_cm,
		              state->lgar_bmi_params.lateral_flow_psi_threshold_cm,
		              state->lgar_bmi_params.lateral_flow_factor,
			              state->lgar_bmi_params.TO_enabled
			                ? lgar_effective_groundwater_depth_cm(&state->lgar_bmi_params)
			                : lgar_fixed_soil_depth_cm(&state->lgar_bmi_params),
		              state->lgar_bmi_params.mobile_groundwater_level);

        // this is the volume of water leaving through the bottom
        volrech_subtimestep_cm = volin_subtimestep_cm;
        // volrech_timestep_cm += volrech_subtimestep_cm;

        volin_subtimestep_cm = volin_subtimestep_cm_temp;
      }

	      const double column_depth_cm =
	        state->lgar_bmi_params.cum_layer_thickness_cm[state->lgar_bmi_params.num_layers];
	      const double surface_cleanup_lower_boundary_cm =
	        state->lgar_bmi_params.TO_enabled
	          ? lgar_effective_groundwater_depth_cm(&state->lgar_bmi_params)
	          : column_depth_cm;
		      lgar_clean_redundant_fronts(&state->head, state->lgar_bmi_params.layer_soil_type,
		                                  state->soil_properties,
		                                  PET_subtimestep_cm_per_h > 0.0,
		                                  state->lgar_bmi_params.cum_layer_thickness_cm,
		                                  column_depth_cm); // deletes redundant WFs; leading zero-depth TO/GW capping is limited to PET-active substeps

	      int correction_type_surf_after_cleanup =
	        lgarto_correction_type_surf(num_layers, state->lgar_bmi_params.cum_layer_thickness_cm,
	                                    &state->head, surface_cleanup_lower_boundary_cm);
	      while (correction_type_surf_after_cleanup == 4) {
	        double cleanup_mass_change_cm = 0.0;
		        lgar_fix_dry_over_wet_wetting_fronts(&cleanup_mass_change_cm,
		                                             state->lgar_bmi_params.cum_layer_thickness_cm,
		                                             state->lgar_bmi_params.layer_soil_type, &state->head,
		                                             state->soil_properties);
	        if (cleanup_mass_change_cm > SMALL_EPS) {
	          temp_rch -= cleanup_mass_change_cm;
	        }
	        else {
	          AET_subtimestep_cm -= cleanup_mass_change_cm;
	        }
	        correction_type_surf_after_cleanup =
	          lgarto_correction_type_surf(num_layers, state->lgar_bmi_params.cum_layer_thickness_cm,
	                                      &state->head, surface_cleanup_lower_boundary_cm);
	      }

      /*----------------------------------------------------------------------*/
      // calculate derivative (dz/dt) for all wetting fronts
      int new_front = -1;
      if (switch_caching){
        if (create_surficial_front){
          new_front = state->head->front_num;
        }
      }
      lgar_dzdt_calc(use_closed_form_G, nint, num_layers, ponded_depth_subtimestep_cm, subtimestep_h, state->lgar_bmi_params.layer_soil_type,
        state->lgar_bmi_params.cum_layer_thickness_cm, state->lgar_bmi_params.frozen_factor,
        state->head, state->soil_properties, switch_caching, state->lgar_bmi_params.cache_count, new_front,
        lgar_effective_groundwater_depth_cm(&state->lgar_bmi_params));

      if (switch_caching){
        state->lgar_bmi_params.cache_count = 1;
      }

      if (state->lgar_bmi_params.TO_enabled &&
          create_surficial_front &&
          creation_excess_runoff_subtimestep_cm > SMALL_EPS) {
        /* Experimental bookkeeping: after creating a saturated surface front,
           dzdt has been recomputed for the updated TO/GW scaffold but that
           scaffold will not physically move again in this creation substep.
           Route only the read-only projected same-substep TO drainage from the
           creation residual to lower-boundary flux, instead of treating all of
           that residual as immediate surface runoff. */
        trace_surface_creation_post_creation_TO_release_cm =
          lgarto_project_TO_motion_lower_boundary_flux_cm(
            subtimestep_h, num_layers,
            state->lgar_bmi_params.cum_layer_thickness_cm,
            state->lgar_bmi_params.layer_soil_type,
            state->head, state->soil_properties,
            lgar_effective_groundwater_depth_cm(&state->lgar_bmi_params));
        const double creation_residual_to_lower_boundary_cm =
          fmin(creation_excess_runoff_subtimestep_cm,
               trace_surface_creation_post_creation_TO_release_cm);
        if (creation_residual_to_lower_boundary_cm > SMALL_EPS) {
          creation_excess_runoff_subtimestep_cm -= creation_residual_to_lower_boundary_cm;
          creation_excess_gw_flux_subtimestep_cm += creation_residual_to_lower_boundary_cm;
        }

        const double already_booked_lower_boundary_cm =
          temp_rch + free_drainage_subtimestep_cm + creation_excess_gw_flux_subtimestep_cm;
        const double saturated_column_capacity_cm =
          lgarto_saturated_creation_lower_boundary_capacity_cm(
            subtimestep_h, num_layers,
            state->lgar_bmi_params.cum_layer_thickness_cm,
            state->lgar_bmi_params.layer_soil_type,
            state->lgar_bmi_params.frozen_factor,
            state->head, state->soil_properties,
            already_booked_lower_boundary_cm);
        const double saturated_column_flux_cm =
          fmin(creation_excess_runoff_subtimestep_cm,
               saturated_column_capacity_cm);
        if (saturated_column_flux_cm > SMALL_EPS) {
          creation_excess_runoff_subtimestep_cm -= saturated_column_flux_cm;
          creation_excess_gw_flux_subtimestep_cm += saturated_column_flux_cm;
          if (verbosity.compare("high") == 0) {
            printf("surface creation saturated-column regime routed %.12e cm "
                   "from runoff to lower-boundary flux "
                   "(remaining capacity %.12e cm).\n",
                   saturated_column_flux_cm,
                   saturated_column_capacity_cm);
          }
        }
      }

	      const double lower_boundary_flux_subtimestep_cm =
	        temp_rch + free_drainage_subtimestep_cm + creation_excess_gw_flux_subtimestep_cm;
	      lower_boundary_flux_for_cache_subtimestep_cm = lower_boundary_flux_subtimestep_cm;
      lower_boundary_flux_for_mobile_groundwater_subtimestep_cm =
        lower_boundary_flux_subtimestep_cm;
	      volrech_subtimestep_cm = 0.0;
      lgar_partition_lower_boundary_flux_for_CR(
        state->lgar_bmi_params.lower_bdy_flux_to_CR,
        lower_boundary_flux_subtimestep_cm,
        &volrech_subtimestep_cm,
	        &lower_boundary_flux_for_CR,
	        &state->lgar_mass_balance.CR_fast_storage_cm,
	        &state->lgar_mass_balance.CR_slow_storage_cm,
        &lower_boundary_CR_exchange_cm);
	      if (state->lgar_bmi_params.lower_bdy_flux_to_CR) {
	        lower_boundary_flux_for_mobile_groundwater_subtimestep_cm =
	          lower_boundary_CR_exchange_cm;
	      }

      if (groundwater_supported_AET_potential_subtimestep_cm > SMALL_EPS) {
        const double remaining_PET_budget_cm =
          fmax(0.0, PET_subtimestep_cm_per_h * subtimestep_h - AET_subtimestep_cm);
        const double groundwater_supported_AET_demand_cm =
          fmin(groundwater_supported_AET_potential_subtimestep_cm,
               remaining_PET_budget_cm);
        const double groundwater_supported_AET_subtimestep_cm =
          lgar_extract_from_CR_storage_fast_then_slow(
            groundwater_supported_AET_demand_cm, state);
        AET_subtimestep_cm += groundwater_supported_AET_subtimestep_cm;
        if (verbosity.compare("high") == 0 &&
            groundwater_supported_AET_subtimestep_cm > SMALL_EPS) {
          printf("Groundwater-supported AET extracted %.17lf cm from CR storage "
                 "(potential=%.17lf demand=%.17lf remaining_PET_budget=%.17lf).\n",
                 groundwater_supported_AET_subtimestep_cm,
                 groundwater_supported_AET_potential_subtimestep_cm,
                 groundwater_supported_AET_demand_cm,
                 remaining_PET_budget_cm);
        }
      }

	      const double runoff_before_creation_excess_cm = volrunoff_subtimestep_cm;
	      volrunoff_subtimestep_cm += creation_excess_runoff_subtimestep_cm;
	      if (creation_excess_runoff_subtimestep_cm > 0.0) {
	        const double accepted_creation_infiltration_cm =
	          volin_subtimestep_cm - creation_excess_runoff_subtimestep_cm;
        if (verbosity.compare("high") == 0) {
          printf("creation-time runoff reclassification updated infiltration bookkeeping: "
                 "volin %.12e -> %.12e cm after subtracting runoff %.12e cm\n",
                 volin_subtimestep_cm, fmax(0.0, accepted_creation_infiltration_cm),
                 creation_excess_runoff_subtimestep_cm);
	        }
	        volin_subtimestep_cm = fmax(0.0, accepted_creation_infiltration_cm);
	      }

	      const bool insert_storage_cap_active =
	        std::isfinite(trace_insert_raw_fp_cm_per_h) &&
	        std::isfinite(trace_insert_capped_fp_cm_per_h) &&
	        trace_insert_capped_fp_cm_per_h < trace_insert_raw_fp_cm_per_h - 1.0e-12;
	      const lgarto_infiltration_limit_trace_row trace_row = {
	        state->lgar_bmi_params.timesteps,
	        this->state->lgar_bmi_params.time_s / 3600.0,
	        cycle,
	        subcycles,
	        subtimestep_h,
	        trace_branch,
	        0,
	        trace_create_surficial_front ? 1 : 0,
	        trace_is_top_wf_saturated ? 1 : 0,
	        trace_fronts_start,
	        trace_surface_fronts_start,
	        listLength(state->head),
	        listLength_surface(state->head),
	        precip_subtimestep_cm,
	        PET_subtimestep_cm,
		        trace_ponded_start_cm,
		        ponded_depth_subtimestep_cm,
		        volin_subtimestep_cm,
	        runoff_before_creation_excess_cm,
	        creation_excess_runoff_subtimestep_cm,
		        volrunoff_subtimestep_cm,
		        creation_excess_gw_flux_subtimestep_cm,
		        trace_surface_creation_gw_capacity_cm,
		        trace_surface_creation_post_creation_TO_release_cm,
		        temp_rch,
	        free_drainage_subtimestep_cm,
	        trace_projected_TO_storage_release_cm,
	        trace_insertion_storage_release_cm,
	        trace_insert_raw_fp_cm_per_h,
	        trace_insert_storage_limit_fp_cm_per_h,
	        trace_insert_capped_fp_cm_per_h,
	        insert_storage_cap_active ? 1 : 0,
	        lower_boundary_flux_subtimestep_cm
	      };
	      lgarto_write_infiltration_limit_trace(trace_row);

	    }

    else {//in this case, we just use cached fluxes in order to save time. No direct computation of fluxes are necessary, and the ones from the last time step are used.
      //also in this case, there will be no runoff due to precipitation partitioning because there is no precipitation 
      volon_subtimestep_cm = volon_timestep_cm;
      PET_timestep_cm += fmax(PET_subtimestep_cm,0.0); // ensures non-negative PET
      state->lgar_mass_balance.accumulated_PET_cm += PET_subtimestep_cm;

      lower_boundary_flux_for_cache_subtimestep_cm = state->lgar_mass_balance.previous_lower_boundary_flux_cm;
      lower_boundary_flux_for_cache_subtimestep_cm =
        lgar_limit_cached_negative_lower_boundary_flux_by_CR_supply(
          state, lower_boundary_flux_for_cache_subtimestep_cm);
      lower_boundary_flux_for_mobile_groundwater_subtimestep_cm =
        lower_boundary_flux_for_cache_subtimestep_cm;
      lgar_partition_lower_boundary_flux_for_CR(
        state->lgar_bmi_params.lower_bdy_flux_to_CR,
        lower_boundary_flux_for_cache_subtimestep_cm,
        &volrech_subtimestep_cm,
        &lower_boundary_flux_for_CR,
        &state->lgar_mass_balance.CR_fast_storage_cm,
        &state->lgar_mass_balance.CR_slow_storage_cm,
        &lower_boundary_CR_exchange_cm);
      if (state->lgar_bmi_params.lower_bdy_flux_to_CR) {
        lower_boundary_flux_for_mobile_groundwater_subtimestep_cm =
          lower_boundary_CR_exchange_cm;
      }
	      state->lgar_mass_balance.accumulated_lower_boundary_flux_cm += lower_boundary_flux_for_cache_subtimestep_cm;

	      const lgarto_infiltration_limit_trace_row trace_row = {
	        state->lgar_bmi_params.timesteps,
	        this->state->lgar_bmi_params.time_s / 3600.0,
	        cycle,
	        subcycles,
	        subtimestep_h,
	        "cached_fluxes",
	        1,
	        0,
	        0,
	        trace_fronts_start,
	        trace_surface_fronts_start,
	        listLength(state->head),
	        listLength_surface(state->head),
	        precip_subtimestep_cm,
	        PET_subtimestep_cm,
		        trace_ponded_start_cm,
		        ponded_depth_subtimestep_cm,
		        volin_subtimestep_cm,
	        volrunoff_subtimestep_cm,
	        0.0,
		        volrunoff_subtimestep_cm,
		        0.0,
		        0.0,
		        0.0,
		        0.0,
		        free_drainage_subtimestep_cm,
	        0.0,
	        0.0,
	        NAN,
	        NAN,
	        NAN,
	        0,
	        lower_boundary_flux_for_cache_subtimestep_cm
	      };
	      lgarto_write_infiltration_limit_trace(trace_row);
	    }

    volend_subtimestep_cm = lgar_calc_mass_bal(state->lgar_bmi_params.cum_layer_thickness_cm, state->head);
    volend_timestep_cm = volend_subtimestep_cm;
    state->lgar_bmi_params.precip_previous_timestep_cm = precip_subtimestep_cm;

    if (verbosity.compare("high") == 0) {
      printf("volCRstart_subtimestep_cm before: %lf \n", volCRstart_subtimestep_cm);
    }

    const double mobile_groundwater_transaction_start_depth_cm =
      lgar_effective_groundwater_depth_cm(&state->lgar_bmi_params);
    double mobile_groundwater_non_CR_exchange_cm = 0.0;

    double volpref_flow_to_CR_subtimestep_cm =
      (precip_for_CR_subtimestep_cm_per_h + ponded_flux_for_CR) * subtimestep_h;
    double vollgarto_domain_to_CR_subtimestep_cm = lower_boundary_flux_for_CR;
    double vollgarto_domain_CR_net_exchange_subtimestep_cm = 0.0;

    const double raw_volin_CR_subtimestep_cm =
      volpref_flow_to_CR_subtimestep_cm + vollgarto_domain_to_CR_subtimestep_cm;
    const double raw_volQ_CR_subtimestep_cm =
      lgar_predict_CR_discharge_cm(
        state,
        subtimestep_h,
        a, a_slow,
        b, b_slow,
        CR_fast_discharge_threshold_cm,
        CR_slow_discharge_threshold_cm,
        frac_slow,
        subtimestep_h > 0.0 ? raw_volin_CR_subtimestep_cm / subtimestep_h : 0.0);
    const double rejected_CR_input_subtimestep_cm =
      lgar_cap_positive_CR_inputs_by_mobile_groundwater_surface(
        state,
        raw_volQ_CR_subtimestep_cm,
        &volpref_flow_to_CR_subtimestep_cm,
        &vollgarto_domain_to_CR_subtimestep_cm,
        "regular CR input");
    if (rejected_CR_input_subtimestep_cm > SMALL_EPS) {
      volrunoff_subtimestep_cm += rejected_CR_input_subtimestep_cm;
    }
    lower_boundary_flux_for_CR = vollgarto_domain_to_CR_subtimestep_cm;
    vollgarto_domain_CR_net_exchange_subtimestep_cm =
      (lower_boundary_CR_exchange_cm < 0.0) ? lower_boundary_CR_exchange_cm : lower_boundary_flux_for_CR;
    if (state->lgar_bmi_params.lower_bdy_flux_to_CR &&
        lower_boundary_flux_for_mobile_groundwater_subtimestep_cm > 0.0) {
      lower_boundary_flux_for_mobile_groundwater_subtimestep_cm =
        lower_boundary_flux_for_CR;
    }

    double volin_CR_subtimestep_cm =
      volpref_flow_to_CR_subtimestep_cm + vollgarto_domain_to_CR_subtimestep_cm;
    double volQ_CR_subtimestep_cm = calc_CR_Q(subtimestep_h, a, a_slow, b, b_slow,
                                              CR_fast_discharge_threshold_cm,
                                              CR_slow_discharge_threshold_cm,
                                              frac_slow,
                                              subtimestep_h > 0.0 ? volin_CR_subtimestep_cm / subtimestep_h : 0.0,
                                              &state->lgar_mass_balance.CR_fast_storage_cm,
                                              &state->lgar_mass_balance.CR_slow_storage_cm);
    if (state->lgar_bmi_params.mobile_groundwater_level) {
      if (state->lgar_bmi_params.lower_bdy_flux_to_CR) {
        mobile_groundwater_explicit_mass_change_subtimestep_cm +=
          lgar_sync_mobile_groundwater_chain_from_CR_storage(state);
      }
      else {
        mobile_groundwater_non_CR_exchange_cm =
          lower_boundary_flux_for_mobile_groundwater_subtimestep_cm +
          (lgar_total_CR_storage_cm(state) - volCRstart_subtimestep_cm);
        lgar_set_mobile_groundwater_depth_from_storage_exchange(
          state,
          mobile_groundwater_transaction_start_depth_cm,
          mobile_groundwater_non_CR_exchange_cm,
          "substep-ledger");
      }
    }
    if (state->lgar_bmi_params.mobile_groundwater_level &&
        state->lgar_bmi_params.TO_enabled) {
      const double mobile_groundwater_post_transaction_depth_cm =
        lgar_effective_groundwater_depth_cm(&state->lgar_bmi_params);
      const double available_recession_rewet_storage_cm =
        fmax(0.0, lgar_total_CR_storage_cm(state));
      if (mobile_groundwater_post_transaction_depth_cm >
            mobile_groundwater_transaction_start_depth_cm + SMALL_EPS &&
          mobile_groundwater_transaction_start_depth_cm <=
            LGARTO_MOBILE_GW_RECESSION_REWET_START_MAX_DEPTH_CM &&
          available_recession_rewet_storage_cm > SMALL_EPS) {
        const double mobile_groundwater_recession_rewet_cm =
          lgarto_rewet_receded_groundwater_zone(
            mobile_groundwater_transaction_start_depth_cm,
            mobile_groundwater_post_transaction_depth_cm,
            num_layers,
            state->lgar_bmi_params.cum_layer_thickness_cm,
            state->lgar_bmi_params.layer_soil_type,
            state->lgar_bmi_params.frozen_factor,
            &state->head,
            state->soil_properties,
            available_recession_rewet_storage_cm);
        if (mobile_groundwater_recession_rewet_cm > SMALL_EPS) {
          const double debited_CR_storage_cm =
            lgar_extract_from_CR_storage_fast_then_slow(
              mobile_groundwater_recession_rewet_cm,
              state);
          vollgarto_domain_CR_net_exchange_subtimestep_cm -=
            debited_CR_storage_cm;

          if (state->lgar_bmi_params.lower_bdy_flux_to_CR) {
            mobile_groundwater_explicit_mass_change_subtimestep_cm +=
              lgar_sync_mobile_groundwater_chain_from_CR_storage(state);
          }
          else {
            mobile_groundwater_non_CR_exchange_cm -= debited_CR_storage_cm;
            lgar_set_mobile_groundwater_depth_from_storage_exchange(
              state,
              mobile_groundwater_transaction_start_depth_cm,
              mobile_groundwater_non_CR_exchange_cm,
              "substep-ledger");
          }

          volend_subtimestep_cm =
            lgar_calc_mass_bal(state->lgar_bmi_params.cum_layer_thickness_cm,
                               state->head);
          volend_timestep_cm = volend_subtimestep_cm;
        }
      }
    }
    if (state->lgar_bmi_params.mobile_groundwater_level &&
        state->lgar_bmi_params.TO_enabled) {
      for (int repair_iter = 0; repair_iter < MAX_NUM_WETTING_FRONTS; repair_iter++) {
        double mobile_groundwater_repair_flux_cm =
          lgarto_submerge_wetting_fronts_below_groundwater(
            lgar_effective_groundwater_depth_cm(&state->lgar_bmi_params),
            num_layers,
            state->lgar_bmi_params.cum_layer_thickness_cm,
            state->lgar_bmi_params.layer_soil_type,
            state->lgar_bmi_params.frozen_factor,
            &state->head,
            state->soil_properties);

        if (fabs(mobile_groundwater_repair_flux_cm) <= SMALL_EPS) {
          break;
        }

        lower_boundary_flux_for_cache_subtimestep_cm +=
          mobile_groundwater_repair_flux_cm;

        double mobile_groundwater_storage_exchange_cm =
          mobile_groundwater_repair_flux_cm;

        if (state->lgar_bmi_params.lower_bdy_flux_to_CR &&
            mobile_groundwater_repair_flux_cm > 0.0) {
          // This repair is discovered after CR discharge was computed for the
          // substep, so add it to storage for subsequent discharge timing.
          double repair_preferential_CR_input_cm = 0.0;
          double repair_domain_CR_input_cm =
            mobile_groundwater_repair_flux_cm;
          const double rejected_repair_CR_input_cm =
            lgar_cap_positive_CR_inputs_by_mobile_groundwater_surface(
              state,
              0.0,
              &repair_preferential_CR_input_cm,
              &repair_domain_CR_input_cm,
              "mobile-GW submergence repair");
          if (rejected_repair_CR_input_cm > SMALL_EPS) {
            volrunoff_subtimestep_cm += rejected_repair_CR_input_cm;
          }
          mobile_groundwater_storage_exchange_cm = repair_domain_CR_input_cm;

          const double slow_input_cm =
            repair_domain_CR_input_cm * frac_slow;
          const double fast_input_cm =
            repair_domain_CR_input_cm - slow_input_cm;
          state->lgar_mass_balance.CR_fast_storage_cm += fast_input_cm;
          state->lgar_mass_balance.CR_slow_storage_cm += slow_input_cm;
          lower_boundary_flux_for_CR += repair_domain_CR_input_cm;
          vollgarto_domain_to_CR_subtimestep_cm +=
            repair_domain_CR_input_cm;
          vollgarto_domain_CR_net_exchange_subtimestep_cm +=
            repair_domain_CR_input_cm;
          volin_CR_subtimestep_cm += repair_domain_CR_input_cm;
        }
        else if (state->lgar_bmi_params.lower_bdy_flux_to_CR &&
                 mobile_groundwater_repair_flux_cm < 0.0) {
          double repair_percolation_cm = 0.0;
          double repair_CR_input_cm = 0.0;
          double repair_CR_exchange_cm = 0.0;
          lgar_partition_lower_boundary_flux_for_CR(
            true,
            mobile_groundwater_repair_flux_cm,
            &repair_percolation_cm,
            &repair_CR_input_cm,
            &state->lgar_mass_balance.CR_fast_storage_cm,
            &state->lgar_mass_balance.CR_slow_storage_cm,
            &repair_CR_exchange_cm);
          mobile_groundwater_storage_exchange_cm = repair_CR_exchange_cm;
          volrech_subtimestep_cm += repair_percolation_cm;
          lower_boundary_flux_for_CR += repair_CR_input_cm;
          vollgarto_domain_to_CR_subtimestep_cm += repair_CR_input_cm;
          vollgarto_domain_CR_net_exchange_subtimestep_cm +=
            repair_CR_exchange_cm;
          volin_CR_subtimestep_cm += repair_CR_input_cm;
        }
        else {
          volrech_subtimestep_cm += mobile_groundwater_repair_flux_cm;
        }

        if (state->lgar_bmi_params.lower_bdy_flux_to_CR) {
          mobile_groundwater_explicit_mass_change_subtimestep_cm +=
            lgar_sync_mobile_groundwater_chain_from_CR_storage(state);
        }
        else {
          mobile_groundwater_non_CR_exchange_cm +=
            mobile_groundwater_storage_exchange_cm;
          lgar_set_mobile_groundwater_depth_from_storage_exchange(
            state,
            mobile_groundwater_transaction_start_depth_cm,
            mobile_groundwater_non_CR_exchange_cm,
            "substep-ledger");
        }
        volend_subtimestep_cm =
          lgar_calc_mass_bal(state->lgar_bmi_params.cum_layer_thickness_cm,
                             state->head);
        volend_timestep_cm = volend_subtimestep_cm;
      }
    }
    if (state->lgar_bmi_params.mobile_groundwater_level &&
        state->lgar_bmi_params.TO_enabled &&
        state->lgar_bmi_params.lower_bdy_flux_to_CR) {
      mobile_groundwater_explicit_mass_change_subtimestep_cm +=
        lgar_sync_mobile_groundwater_chain_from_CR_storage(state);
    }

    state->lgar_mass_balance.volrunoff_CR_cm += volQ_CR_subtimestep_cm;
    volQ_CR_timestep_cm += volQ_CR_subtimestep_cm;
    volpref_flow_to_CR_timestep_cm += volpref_flow_to_CR_subtimestep_cm;
    vollgarto_domain_to_CR_timestep_cm += vollgarto_domain_CR_net_exchange_subtimestep_cm;
    volCRend_subtimestep_cm = state->lgar_mass_balance.CR_fast_storage_cm + state->lgar_mass_balance.CR_slow_storage_cm;
    volCRend_timestep_cm = volCRend_subtimestep_cm;
    mobile_groundwater_explicit_mass_change_timestep_cm +=
      mobile_groundwater_explicit_mass_change_subtimestep_cm;
    if (fabs(mobile_groundwater_explicit_mass_change_subtimestep_cm) > SMALL_EPS) {
      volend_subtimestep_cm =
        lgar_calc_mass_bal(state->lgar_bmi_params.cum_layer_thickness_cm,
                           state->head);
      volend_timestep_cm = volend_subtimestep_cm;
    }

    if (verbosity.compare("high") == 0) {
      printf("volCRstart_subtimestep_cm after: %lf \n", volCRend_timestep_cm);
    }

    // set runoff_in_prev_step for next step
    if ((volrunoff_subtimestep_cm > SMALL_EPS) || (top_near_sat)){
      state->lgar_bmi_params.runoff_in_prev_step = true;
    }
    else {
      state->lgar_bmi_params.runoff_in_prev_step = false;
    }

    //add precip_for_CR_subtimestep_cm_per_h back into precip for mass balance, note that this means volin_CR_timestep is not necessary for mass balance 
    precip_subtimestep_cm += precip_for_CR_subtimestep_cm_per_h * subtimestep_h;

    /*----------------------------------------------------------------------*/
    // mass balance at the subtimestep (local mass balance)

    double local_mb = volstart_subtimestep_cm + precip_subtimestep_cm + volon_timestep_cm + mobile_groundwater_explicit_mass_change_subtimestep_cm - volrunoff_subtimestep_cm - volQ_CR_subtimestep_cm - volCRend_subtimestep_cm + volCRstart_subtimestep_cm
                      - AET_subtimestep_cm - volon_subtimestep_cm - volrech_subtimestep_cm - lateral_flow_subtimestep_cm - volend_subtimestep_cm;

    /*----------------------------------------------------------------------*/

    ///////
    //separating code such that most non substep vars (so xxx_timestep and not xxx_subtimestep) are updated in just one place. not all, because some must be set before substepping.
    volin_timestep_cm += volin_subtimestep_cm;
    volrech_timestep_cm += volrech_subtimestep_cm;

    volrunoff_timestep_cm += volrunoff_subtimestep_cm;
    vollateral_flow_timestep_cm += lateral_flow_subtimestep_cm;

    AET_timestep_cm += AET_subtimestep_cm;
    volon_timestep_cm = volon_subtimestep_cm; // surface ponded water at the end of the timestep 
    ///////
    
    if (verbosity.compare("high") == 0 || verbosity.compare("low") == 0) {
      printf("Printing wetting fronts at this subtimestep... \n");
      listPrint(state->head);
    }

    bool unexpected_local_error = fabs(local_mb) > mbal_tol ? true : false; //1.0E-4 was the default for initial stability testing 
    if (isinf(local_mb)){
      unexpected_local_error = true;
    }
    
    if (verbosity.compare("high") == 0 || verbosity.compare("low") == 0 || unexpected_local_error) {
      if (!state->lgar_bmi_params.frac_to_CR && !state->lgar_bmi_params.lower_bdy_flux_to_CR){
        printf("\nLocal mass balance at this timestep... \n\
        Error         = %14.10f \n\
        Initial water = %14.10f \n\
        Water added   = %14.10f \n\
        Ponded water  = %14.10f \n\
        Infiltration  = %14.10f \n\
        Runoff        = %14.10f \n\
        AET           = %14.10f \n\
        Lateral flow  = %14.10f \n\
        Percolation   = %14.10f \n\
        Final water   = %14.10f \n", local_mb, volstart_subtimestep_cm, precip_subtimestep_cm, volon_subtimestep_cm,
        volin_subtimestep_cm, volrunoff_subtimestep_cm, AET_subtimestep_cm, lateral_flow_subtimestep_cm, volrech_subtimestep_cm,
        volend_subtimestep_cm);
      }
      else {
        printf("\nLocal mass balance at this timestep... \n\
        Error                   = %14.10f \n\
        Initial water (LGAR)    = %14.10f \n\
        Water added (total)     = %14.10f \n\
        Ponded water            = %14.10f \n\
        Infiltration (LGAR)     = %14.10f \n\
        Runoff (total)          = %14.10f \n\
        AET                     = %14.10f \n\
        Lateral flow            = %14.10f \n\
        Percolation             = %14.10f \n\
        Final water (LGAR)      = %14.10f \n\
        Water added (con res)   = %14.10f \n\
        Initial water (con res) = %14.10f \n\
        Final water (con res)   = %14.10f \n\
        Runoff (con res)        = %14.10f \n", local_mb, volstart_subtimestep_cm, precip_subtimestep_cm, volon_subtimestep_cm,
        volin_subtimestep_cm, volrunoff_subtimestep_cm, AET_subtimestep_cm, lateral_flow_subtimestep_cm, volrech_subtimestep_cm,
        volend_subtimestep_cm, volin_CR_subtimestep_cm, volCRstart_subtimestep_cm,
        volCRend_subtimestep_cm, volQ_CR_subtimestep_cm);
      }

      if (unexpected_local_error) {
	printf("Local mass balance (in this timestep) is %14.10f, larger than expected, needs some debugging...\n ",local_mb);
	abort();
      }

    }

    // store local mass balance error to the struct
    state->lgar_mass_balance.local_mass_balance = local_mb;

    const double column_depth_cm =
      state->lgar_bmi_params.cum_layer_thickness_cm[state->lgar_bmi_params.num_layers];
    double vadose_assertion_depth_cm = column_depth_cm;
    if (state->lgar_bmi_params.mobile_groundwater_level &&
        state->lgar_bmi_params.TO_enabled) {
      vadose_assertion_depth_cm =
        fmax(column_depth_cm, lgar_effective_groundwater_depth_cm(&state->lgar_bmi_params));
    }
    // lgarto_abort_if_deferred_gw_flux_mass_balance_correction_exceeded(state->head);
    lgar_assert_wetting_fronts_nonnegative_depth(state->head);
    lgar_assert_wetting_front_depth_order(state->head);
    lgar_assert_wetting_fronts_within_vadose_zone(vadose_assertion_depth_cm, state->head);
    lgar_assert_to_psi_monotonic_with_depth(state->head);
    lgar_assert_zero_depth_TO_supports_drier_than_surface_TO_chain(state->head);
    lgar_assert_boundary_psi_continuity(state->head);
    lgar_assert_surface_fronts_not_partial_to_bottom_scaffold(
      state->lgar_bmi_params.TO_enabled,
      state->lgar_bmi_params.num_layers,
      state->head);
    if (state->lgar_bmi_params.mobile_groundwater_level &&
        state->lgar_bmi_params.TO_enabled &&
        state->lgar_bmi_params.lower_bdy_flux_to_CR) {
      lgarto_assert_mobile_groundwater_CR_chain_consistency(
        lgar_total_CR_storage_cm(state),
        lgar_effective_groundwater_depth_cm(&state->lgar_bmi_params),
        state->lgar_bmi_params.num_layers,
        state->lgar_bmi_params.cum_layer_thickness_cm,
        state->lgar_bmi_params.layer_soil_type,
        state->head,
        state->soil_properties);
    }
    lgar_assert_to_bottom_scaffold(state->lgar_bmi_params.num_layers,
                                   state->lgar_bmi_params.cum_layer_thickness_cm,
                                   state->head);

    remaining_forcing_h -= subtimestep_h;
    if (remaining_forcing_h < SMALL_EPS) {
      remaining_forcing_h = 0.0;
    }

    bool lasam_standalone = true;
#ifdef NGEN
    lasam_standalone = false;
#endif
    // simuation time can't exceed the endtime when running standalone
    if ( (this->state->lgar_bmi_params.time_s >= this->state->lgar_bmi_params.endtime_s) && lasam_standalone)
      break;

  } // end of subcycling

  state->lgar_bmi_params.timestep_h = base_subtimestep_h;

  lgar_assert_wetting_fronts_nonnegative_depth(state->head);
  lgar_assert_wetting_front_depth_order(state->head);
  lgar_assert_zero_depth_TO_supports_drier_than_surface_TO_chain(state->head);
  lgar_assert_to_bottom_scaffold(state->lgar_bmi_params.num_layers,
                                 state->lgar_bmi_params.cum_layer_thickness_cm,
                                 state->head);

  //update giuh at the time step level (was previously updated at the sub time step level)
  volrunoff_giuh_timestep_cm = giuh_convolution_integral(volrunoff_timestep_cm + volQ_CR_timestep_cm + vollateral_flow_timestep_cm, num_giuh_ordinates, giuh_ordinates, giuh_runoff_queue);

  // total mass of water leaving the system, at this time it is the giuh-only, but later will add groundwater component as well.
  // when groundwater component is added, it should probably happen inside of the subcycling loop.
  volQ_timestep_cm = volrunoff_giuh_timestep_cm; //note that volQ_CR_timestep_cm was added to volQ_timestep_cm in the input for giuh_convolution_integral

  /*----------------------------------------------------------------------*/
  // Everything related to lgar state is done at this point, now time to update some dynamic variables

  // update number of wetting fronts
  state->lgar_bmi_params.num_wetting_fronts = listLength(state->head);

  // allocate new memory based on updated wetting fronts; we could make it conditional i.e. create only if no. of wf are changed
  realloc_soil();

  // update thickness/depth and soil moisture of wetting fronts (used for state coupling)
  struct wetting_front *current = state->head;
  int to_bottom_count = 0;
  for (int i=0; i<state->lgar_bmi_params.num_wetting_fronts; i++) {
    if (current->to_bottom){
      to_bottom_count ++;
    }
    if (to_bottom_count>num_layers){
      std::cerr << "Error: too many to_bottom WFs! This should be equal to the number of layers.\n";
      listPrint(state->head);
	    abort();
    }
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
  state->lgar_mass_balance.volprecip_timestep_cm  = precip_timestep_cm;
  state->lgar_mass_balance.volin_timestep_cm      = volin_timestep_cm;
  state->lgar_mass_balance.volon_timestep_cm      = volon_timestep_cm;
  state->lgar_mass_balance.volend_timestep_cm     = volend_timestep_cm;
  state->lgar_mass_balance.volCRend_timestep_cm   = volCRend_timestep_cm;
  state->lgar_mass_balance.volAET_timestep_cm     = AET_timestep_cm;
  state->lgar_mass_balance.volrech_timestep_cm    = volrech_timestep_cm;
  state->lgar_mass_balance.volrunoff_timestep_cm  = volrunoff_timestep_cm;
  state->lgar_mass_balance.volQ_timestep_cm       = volQ_timestep_cm;
  state->lgar_mass_balance.volQ_CR_timestep_cm    = volQ_CR_timestep_cm;
  state->lgar_mass_balance.volpref_flow_to_CR_timestep_cm = volpref_flow_to_CR_timestep_cm;
  state->lgar_mass_balance.vollgarto_domain_to_CR_timestep_cm = vollgarto_domain_to_CR_timestep_cm;
  state->lgar_mass_balance.vollateral_flow_timestep_cm = vollateral_flow_timestep_cm;
  state->lgar_mass_balance.volPET_timestep_cm     = PET_timestep_cm;
  state->lgar_mass_balance.volrunoff_giuh_timestep_cm = volrunoff_giuh_timestep_cm;

  //for caching
  state->lgar_mass_balance.previous_AET         = AET_subtimestep_cm;
  state->lgar_mass_balance.previous_PET         = PET_subtimestep_cm;
  if (!state->lgar_mass_balance.cache_fluxes){
    state->lgar_mass_balance.previous_lower_boundary_flux_cm = lower_boundary_flux_for_cache_subtimestep_cm;
  }

  // add to mass balance accumulated variables
  state->lgar_mass_balance.volprecip_cm  += precip_timestep_cm;
  state->lgar_mass_balance.volin_cm      += volin_timestep_cm;
  state->lgar_mass_balance.volon_cm       = volon_timestep_cm;
  state->lgar_mass_balance.volend_cm      = volend_timestep_cm;
  state->lgar_mass_balance.volCRend_cm    = volCRend_timestep_cm;
  state->lgar_mass_balance.volAET_cm     += AET_timestep_cm;
  state->lgar_mass_balance.volrech_cm    += volrech_timestep_cm;
  state->lgar_mass_balance.volrunoff_cm  += volrunoff_timestep_cm;
  state->lgar_mass_balance.volQ_cm       += volQ_timestep_cm;
  state->lgar_mass_balance.volQ_CR_cm    += volQ_CR_timestep_cm;
  state->lgar_mass_balance.volpref_flow_to_CR_cm += volpref_flow_to_CR_timestep_cm;
  state->lgar_mass_balance.vollgarto_domain_to_CR_cm += vollgarto_domain_to_CR_timestep_cm;
  state->lgar_mass_balance.vollateral_flow_cm += vollateral_flow_timestep_cm;
  state->lgar_mass_balance.volPET_cm     += PET_timestep_cm;
  state->lgar_mass_balance.volrunoff_giuh_cm  += volrunoff_giuh_timestep_cm;
  state->lgar_mass_balance.volchange_calib_cm +=
    volchange_calib_cm + mobile_groundwater_explicit_mass_change_timestep_cm;
 
  // converted values, a struct local to the BMI and has bmi output variables
  bmi_unit_conv.mass_balance_m        = state->lgar_mass_balance.local_mass_balance * state->units.cm_to_m;
  bmi_unit_conv.volprecip_timestep_m  = precip_timestep_cm * state->units.cm_to_m;
  bmi_unit_conv.volin_timestep_m      = volin_timestep_cm * state->units.cm_to_m;
  bmi_unit_conv.volend_timestep_m     = volend_timestep_cm * state->units.cm_to_m;
  bmi_unit_conv.volCRend_timestep_m   = volCRend_timestep_cm * state->units.cm_to_m;
  bmi_unit_conv.volAET_timestep_m     = AET_timestep_cm * state->units.cm_to_m;
  bmi_unit_conv.volrech_timestep_m    = volrech_timestep_cm * state->units.cm_to_m;
  bmi_unit_conv.volrunoff_timestep_m  = volrunoff_timestep_cm * state->units.cm_to_m;
  bmi_unit_conv.volQ_timestep_m       = volQ_timestep_cm * state->units.cm_to_m;
  bmi_unit_conv.volQ_CR_timestep_m    = volQ_CR_timestep_cm * state->units.cm_to_m;
  bmi_unit_conv.volpref_flow_to_CR_timestep_m = volpref_flow_to_CR_timestep_cm * state->units.cm_to_m;
  bmi_unit_conv.vollgarto_domain_to_CR_timestep_m = vollgarto_domain_to_CR_timestep_cm * state->units.cm_to_m;
  bmi_unit_conv.vollateral_flow_timestep_m = vollateral_flow_timestep_cm * state->units.cm_to_m;
  bmi_unit_conv.volPET_timestep_m     = PET_timestep_cm * state->units.cm_to_m;
  bmi_unit_conv.volrunoff_giuh_timestep_m = volrunoff_giuh_timestep_cm * state->units.cm_to_m;
  
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
update_calibratable_parameters()
{
  int soil, layer_num;
  struct wetting_front *current = state->head;

  if (verbosity.compare("high") == 0)
    listPrint(state->head);
  
  double volstart_before = lgar_calc_mass_bal(state->lgar_bmi_params.cum_layer_thickness_cm, state->head);

  // first we update the parameters that depend on soil layer, for each layer. 
  // This no longer relies on arrays (for compatability with ngen-cal) and instead supports calibration of up to 3 soil layers with scalar parameters.
  for (int i=0; i<state->lgar_bmi_params.num_wetting_fronts; i++) {
    layer_num  = current->layer_num;
    soil = state->lgar_bmi_params.layer_soil_type[layer_num];
    
    assert (current != NULL);

    if (verbosity.compare("high") == 0 || verbosity.compare("low") == 0) {
      std::cerr<<"----------- Calibratable parameters depending on soil layer (initial values) ----------- \n";
      std::cerr<<"| soil_type = "<< soil <<", layer = "<<layer_num
	       <<", smcmax = "   << state->soil_properties[soil].theta_e
	       <<", smcmin = "   << state->soil_properties[soil].theta_r
	       <<", vg_n = "     << state->soil_properties[soil].vg_n
	       <<", vg_alpha = " << state->soil_properties[soil].vg_alpha_per_cm
	       <<", Ksat = "     << state->soil_properties[soil].Ksat_cm_per_h
	       <<", theta = "    << current->theta <<"\n";
    }
    
    // state->soil_properties[soil].theta_e = state->lgar_calib_params.theta_e[layer_num-1];
    // state->soil_properties[soil].theta_r = state->lgar_calib_params.theta_r[layer_num-1];
    // state->soil_properties[soil].vg_n    = state->lgar_calib_params.vg_n[layer_num-1];
    // state->soil_properties[soil].vg_m    = 1.0 - 1.0/state->soil_properties[soil].vg_n;
    // state->soil_properties[soil].vg_alpha_per_cm = state->lgar_calib_params.vg_alpha[layer_num-1];
    // state->soil_properties[soil].Ksat_cm_per_h   = state->lgar_calib_params.Ksat[layer_num-1];

    if (layer_num==1){
      state->soil_properties[soil].theta_e = state->lgar_calib_params.theta_e_1;
      state->soil_properties[soil].theta_r = state->lgar_calib_params.theta_r_1;
      state->soil_properties[soil].vg_n    = state->lgar_calib_params.vg_n_1;
      state->soil_properties[soil].vg_m    = 1.0 - 1.0/state->soil_properties[soil].vg_n;
      state->soil_properties[soil].vg_alpha_per_cm = state->lgar_calib_params.vg_alpha_1;
      state->soil_properties[soil].Ksat_cm_per_h   = state->lgar_calib_params.Ksat_1;
      if (state->lgar_bmi_params.log_mode){
        state->soil_properties[soil].vg_alpha_per_cm = pow(10.0, state->lgar_calib_params.vg_alpha_1);
        state->soil_properties[soil].Ksat_cm_per_h   = pow(10.0, state->lgar_calib_params.Ksat_1);
      }
    }

    if (layer_num==2){
      state->soil_properties[soil].theta_e = state->lgar_calib_params.theta_e_2;
      state->soil_properties[soil].theta_r = state->lgar_calib_params.theta_r_2;
      state->soil_properties[soil].vg_n    = state->lgar_calib_params.vg_n_2;
      state->soil_properties[soil].vg_m    = 1.0 - 1.0/state->soil_properties[soil].vg_n;
      state->soil_properties[soil].vg_alpha_per_cm = state->lgar_calib_params.vg_alpha_2;
      state->soil_properties[soil].Ksat_cm_per_h   = state->lgar_calib_params.Ksat_2;
      if (state->lgar_bmi_params.log_mode){
        state->soil_properties[soil].vg_alpha_per_cm = pow(10.0, state->lgar_calib_params.vg_alpha_2);
        state->soil_properties[soil].Ksat_cm_per_h   = pow(10.0, state->lgar_calib_params.Ksat_2);
      }
    }

    if (layer_num==3){
      state->soil_properties[soil].theta_e = state->lgar_calib_params.theta_e_3;
      state->soil_properties[soil].theta_r = state->lgar_calib_params.theta_r_3;
      state->soil_properties[soil].vg_n    = state->lgar_calib_params.vg_n_3;
      state->soil_properties[soil].vg_m    = 1.0 - 1.0/state->soil_properties[soil].vg_n;
      state->soil_properties[soil].vg_alpha_per_cm = state->lgar_calib_params.vg_alpha_3;
      state->soil_properties[soil].Ksat_cm_per_h   = state->lgar_calib_params.Ksat_3;
      if (state->lgar_bmi_params.log_mode){
        state->soil_properties[soil].vg_alpha_per_cm = pow(10.0, state->lgar_calib_params.vg_alpha_3);
        state->soil_properties[soil].Ksat_cm_per_h   = pow(10.0, state->lgar_calib_params.Ksat_3);
      }
    }
    
    current->theta = calc_theta_from_h(current->psi_cm, state->soil_properties[soil].vg_alpha_per_cm,
				       state->soil_properties[soil].vg_m, state->soil_properties[soil].vg_n,
				       state->soil_properties[soil].theta_e, state->soil_properties[soil].theta_r);

    if (verbosity.compare("high") == 0 || verbosity.compare("low") == 0) {
      std::cerr<<"----------- Calibratable parameters depending on soil layer (updated values) ----------- \n";
      std::cerr<<"| soil_type = "<< soil <<", layer = "<<layer_num
	       <<", smcmax = "   << state->soil_properties[soil].theta_e
	       <<", smcmin = "   << state->soil_properties[soil].theta_r
	       <<", vg_n = "     << state->soil_properties[soil].vg_n
	       <<", vg_alpha = " << state->soil_properties[soil].vg_alpha_per_cm
	       <<", Ksat = "     << state->soil_properties[soil].Ksat_cm_per_h
	       <<", theta = "    << current->theta <<"\n";
    }
    
    current = current->next;
  }

  //next we update the parameters that apply to the whole model domain and do not depend on soil layer
  if (verbosity.compare("high") == 0 || verbosity.compare("low") == 0) {
    std::cerr<<"----------- Calibratable parameters independent of soil layer (initial values) ----------- \n";
    std::cerr<<"field_capacity_psi = "   << state->lgar_bmi_params.field_capacity_psi_cm
      <<", ponded_depth_max = "     << state->lgar_bmi_params.ponded_depth_max_cm
      <<", a = "     << state->lgar_bmi_params.a
      <<", b = "     << state->lgar_bmi_params.b
      <<", frac_to_CR = "     << state->lgar_bmi_params.frac_to_CR
      <<", spf_factor = "     << state->lgar_bmi_params.spf_factor <<
      "\n";
    if (state->lgar_bmi_params.frac_slow){
      std::cerr
      <<", a_slow = "     << state->lgar_bmi_params.a_slow
      <<", b_slow = "     << state->lgar_bmi_params.b_slow
      <<", frac__slow = "     << state->lgar_bmi_params.frac_slow <<
      "\n";
    }
  }



  state->lgar_bmi_params.field_capacity_psi_cm = state->lgar_calib_params.field_capacity_psi;
  state->lgar_bmi_params.ponded_depth_max_cm   = state->lgar_calib_params.ponded_depth_max;
  state->lgar_bmi_params.a                     = state->lgar_calib_params.a;
  state->lgar_bmi_params.b                     = state->lgar_calib_params.b;
  state->lgar_bmi_params.frac_to_CR            = state->lgar_calib_params.frac_to_CR;
  state->lgar_bmi_params.spf_factor            = state->lgar_calib_params.spf_factor;

  if (state->lgar_bmi_params.frac_slow){
    state->lgar_bmi_params.a_slow              = state->lgar_calib_params.a_slow;
    state->lgar_bmi_params.b_slow              = state->lgar_calib_params.b_slow;
    state->lgar_bmi_params.frac_slow           = state->lgar_calib_params.frac_slow;
  }

  if (state->lgar_bmi_params.log_mode){
    state->lgar_bmi_params.a = pow(10.0, state->lgar_calib_params.a);
    if (state->lgar_bmi_params.frac_slow){
      state->lgar_bmi_params.a_slow = pow(10.0, state->lgar_calib_params.a_slow);
    }
  }

  if (verbosity.compare("high") == 0 || verbosity.compare("low") == 0) {
    std::cerr<<"----------- Calibratable parameters independent of soil layer (updated values) ----------- \n";
    std::cerr<<"field_capacity_psi = "   << state->lgar_bmi_params.field_capacity_psi_cm
      <<", ponded_depth_max = "     << state->lgar_bmi_params.ponded_depth_max_cm
      <<", a = "     << state->lgar_bmi_params.a
      <<", b = "     << state->lgar_bmi_params.b
      <<", frac_to_CR = "     << state->lgar_bmi_params.frac_to_CR
      <<", spf_factor = "     << state->lgar_bmi_params.spf_factor <<
      "\n";
    if (state->lgar_bmi_params.frac_slow){
      std::cerr
      <<", a_slow = "     << state->lgar_bmi_params.a_slow
      <<", b_slow = "     << state->lgar_bmi_params.b_slow
      <<", frac__slow = "     << state->lgar_bmi_params.frac_slow << 
      "\n";
    }
  }
  
  if (verbosity.compare("high") == 0)
    listPrint(state->head);
  
  double volstart_after = lgar_calc_mass_bal(state->lgar_bmi_params.cum_layer_thickness_cm, state->head);

  if (verbosity.compare("high") == 0 || verbosity.compare("low") == 0)
    std::cerr<<"Mass of water (before and after) = "<< volstart_before<<", "<< volstart_after <<"\n";
  
  return volstart_after - volstart_before;
}

void BmiLGAR::
Finalize()
{
  global_mass_balance();
  if (state->lgar_bmi_params.TO_enabled && state->lgar_bmi_params.allow_flux_caching) {
    printf("Warning: TO mode ran with allow_flux_caching=true. "
           "LGARTO flux caching is experimental; check global mass balance and hydrograph outputs.\n");
    fflush(stdout);
  }
  listDelete(state->head);
  listDelete(state->state_previous);

  delete [] state->soil_properties;

  delete [] state->lgar_bmi_params.soil_depth_wetting_fronts;
  delete [] state->lgar_bmi_params.soil_moisture_wetting_fronts;

  delete [] state->lgar_bmi_params.soil_temperature;
  delete [] state->lgar_bmi_params.soil_temperature_z;
  delete [] state->lgar_bmi_params.layer_soil_type;

  // no longer needed as per layer calibratable parameters are no longer handled with arrays
  // delete [] state->lgar_calib_params.theta_e;
  // delete [] state->lgar_calib_params.theta_r;
  // delete [] state->lgar_calib_params.vg_n;
  // delete [] state->lgar_calib_params.vg_alpha;
  // delete [] state->lgar_calib_params.Ksat;

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
  else if (name.compare("surface_runoff") == 0 || name.compare("giuh_runoff") == 0 || name.compare("a") == 0 || name.compare("b") == 0 || name.compare("frac_to_CR") == 0 || name.compare("spf_factor") == 0
	   || name.compare("a_slow") == 0 || name.compare("b_slow") == 0 || name.compare("frac_slow") == 0 || name.compare("soil_storage") == 0 || name.compare("field_capacity") == 0 || name.compare("ponded_depth_max") == 0)// double
    return 1;
  else if (name.compare("total_discharge") == 0 || name.compare("infiltration") == 0
	   || name.compare("percolation") == 0  || name.compare("conceptual_reservoir_to_stream_discharge") == 0
	   || name.compare("preferential_flow_to_conceptual_reservoir") == 0
	   || name.compare("lgarto_domain_to_conceptual_reservoir") == 0
	   || name.compare("lateral_flow") == 0) // double
    return 1;
  else if (name.compare("mass_balance") == 0)
    return 1;
  else if (name.compare("smcmax_1") == 0 || name.compare("smcmin_1") == 0 // per layer parameters are now handled with scalars and not arrays
	   || name.compare("van_genuchten_m_1") == 0 || name.compare("van_genuchten_alpha_1") == 0 || name.compare("van_genuchten_n_1") == 0 
	   || name.compare("hydraulic_conductivity_1") == 0)
    return 1;
  else if (name.compare("smcmax_2") == 0 || name.compare("smcmin_2") == 0 // per layer parameters are now handled with scalars and not arrays
	   || name.compare("van_genuchten_m_2") == 0 || name.compare("van_genuchten_alpha_2") == 0 || name.compare("van_genuchten_n_2") == 0 
	   || name.compare("hydraulic_conductivity_2") == 0)
    return 1;
  else if (name.compare("smcmax_2") == 0 || name.compare("smcmin_2") == 0 // per layer parameters are now handled with scalars and not arrays
	   || name.compare("van_genuchten_m_2") == 0 || name.compare("van_genuchten_alpha_2") == 0 || name.compare("van_genuchten_n_2") == 0 
	   || name.compare("hydraulic_conductivity_2") == 0)
    return 1;
  // else if (name.compare("soil_depth_layers") == 0  || name.compare("smcmax") == 0 || name.compare("smcmin") == 0
	//    || name.compare("van_genuchten_m") == 0 || name.compare("van_genuchten_alpha") == 0 || name.compare("van_genuchten_n") == 0 
	//    || name.compare("hydraulic_conductivity") == 0) // array of doubles (fixed length)
  else if (name.compare("soil_depth_layers") == 0 )// array of doubles (fixed length)
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
  else if (name.compare("surface_runoff") == 0 || name.compare("giuh_runoff") == 0 
	   || name.compare("soil_storage") == 0) // double
    return "m";
  else if (name.compare("total_discharge") == 0 || name.compare("infiltration") == 0
	   || name.compare("percolation") == 0) // double
    return "m";
  else if (name.compare("mass_balance") == 0 || name.compare("conceptual_reservoir_to_stream_discharge") == 0
	   || name.compare("preferential_flow_to_conceptual_reservoir") == 0
	   || name.compare("lgarto_domain_to_conceptual_reservoir") == 0
	   || name.compare("lateral_flow") == 0)
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
  if (name.compare("precipitation_rate") == 0 || name.compare("precipitation") == 0 ||
      name.compare("potential_evapotranspiration") == 0 || name.compare("potential_evapotranspiration_rate") == 0
      || name.compare("actual_evapotranspiration") == 0) // double
    return "node";
  else if (name.compare("surface_runoff") == 0 || name.compare("giuh_runoff") == 0 || name.compare("a") == 0 || name.compare("b") == 0 || name.compare("frac_to_CR") == 0 || name.compare("spf_factor") == 0
	   || name.compare("a_slow") == 0 || name.compare("b_slow") == 0 || name.compare("frac_slow") == 0 || name.compare("soil_storage") == 0) // double
    return "node";
   else if (name.compare("total_discharge") == 0 || name.compare("infiltration") == 0
	    || name.compare("percolation") == 0 || name.compare("conceptual_reservoir_to_stream_discharge") == 0
	    || name.compare("preferential_flow_to_conceptual_reservoir") == 0
	    || name.compare("lgarto_domain_to_conceptual_reservoir") == 0
	    || name.compare("lateral_flow") == 0) // double
    return "node";
  else if (name.compare("soil_moisture_wetting_fronts") == 0) // array of doubles
    return "node";
  else if (name.compare("mass_balance") == 0)
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
  else if (name.compare("giuh_runoff") == 0)
    return (void*)(&bmi_unit_conv.volrunoff_giuh_timestep_m);
  else if (name.compare("soil_storage") == 0)
    return (void*)(&bmi_unit_conv.volend_timestep_m);
  else if (name.compare("conceptual_reservoir_storage") == 0)
    return (void*)(&bmi_unit_conv.volCRend_timestep_m);
  else if (name.compare("total_discharge") == 0)
    return (void*)(&bmi_unit_conv.volQ_timestep_m);
  else if (name.compare("infiltration") == 0)
    return (void*)(&bmi_unit_conv.volin_timestep_m);
  else if (name.compare("percolation") == 0)
    return (void*)(&bmi_unit_conv.volrech_timestep_m);
  else if (name.compare("conceptual_reservoir_to_stream_discharge") == 0)
    return (void*)(&bmi_unit_conv.volQ_CR_timestep_m);
  else if (name.compare("preferential_flow_to_conceptual_reservoir") == 0)
    return (void*)(&bmi_unit_conv.volpref_flow_to_CR_timestep_m);
  else if (name.compare("lgarto_domain_to_conceptual_reservoir") == 0)
    return (void*)(&bmi_unit_conv.vollgarto_domain_to_CR_timestep_m);
  else if (name.compare("lateral_flow") == 0)
    return (void*)(&bmi_unit_conv.vollateral_flow_timestep_m);
  else if (name.compare("mass_balance") == 0)
    return (void*)(&bmi_unit_conv.mass_balance_m);
  else if (name.compare("soil_depth_layers") == 0)
    return (void*)this->state->lgar_bmi_params.cum_layer_thickness_cm;  // this too and, if needed, change soil_moisture_layers to soil_thickness_layers
  else if (name.compare("soil_moisture_wetting_fronts") == 0)
    return (void*)this->state->lgar_bmi_params.soil_moisture_wetting_fronts;
  else if (name.compare("soil_depth_wetting_fronts") == 0)
    return (void*)this->state->lgar_bmi_params.soil_depth_wetting_fronts;
  else if (name.compare("soil_num_wetting_fronts") == 0)
    return (void*)(&state->lgar_bmi_params.num_wetting_fronts);
  else if (name.compare("soil_temperature_profile") == 0)
    return (void*)this->state->lgar_bmi_params.soil_temperature;
  // else if (name.compare("smcmax") == 0)
  //   return (void*)this->state->lgar_calib_params.theta_e;
  // else if (name.compare("smcmin") == 0)
  //   return (void*)this->state->lgar_calib_params.theta_r;
  // else if (name.compare("van_genuchten_n") == 0)
  //   return (void*)this->state->lgar_calib_params.vg_n;
  // else if (name.compare("van_genuchten_alpha") == 0)
  //   return (void*)this->state->lgar_calib_params.vg_alpha;
  // else if (name.compare("hydraulic_conductivity") == 0)
  //   return (void*)this->state->lgar_calib_params.Ksat;

  // per layer calibratable params are now scalars and not arrays
  else if (name.compare("smcmax_1") == 0)
    return (void*)&this->state->lgar_calib_params.theta_e_1;
  else if (name.compare("smcmin_1") == 0)
    return (void*)&this->state->lgar_calib_params.theta_r_1;
  else if (name.compare("van_genuchten_n_1") == 0)
    return (void*)&this->state->lgar_calib_params.vg_n_1;
  else if (name.compare("van_genuchten_alpha_1") == 0)
    return (void*)&this->state->lgar_calib_params.vg_alpha_1;
  else if (name.compare("hydraulic_conductivity_1") == 0)
    return (void*)&this->state->lgar_calib_params.Ksat_1;
  else if (name.compare("smcmax_2") == 0)
    return (void*)&this->state->lgar_calib_params.theta_e_2;
  else if (name.compare("smcmin_2") == 0)
    return (void*)&this->state->lgar_calib_params.theta_r_2;
  else if (name.compare("van_genuchten_n_2") == 0)
    return (void*)&this->state->lgar_calib_params.vg_n_2;
  else if (name.compare("van_genuchten_alpha_2") == 0)
    return (void*)&this->state->lgar_calib_params.vg_alpha_2;
  else if (name.compare("hydraulic_conductivity_2") == 0)
    return (void*)&this->state->lgar_calib_params.Ksat_2;
  else if (name.compare("smcmax_3") == 0)
    return (void*)&this->state->lgar_calib_params.theta_e_3;
  else if (name.compare("smcmin_3") == 0)
    return (void*)&this->state->lgar_calib_params.theta_r_3;
  else if (name.compare("van_genuchten_n_3") == 0)
    return (void*)&this->state->lgar_calib_params.vg_n_3;
  else if (name.compare("van_genuchten_alpha_3") == 0)
    return (void*)&this->state->lgar_calib_params.vg_alpha_3;
  else if (name.compare("hydraulic_conductivity_3") == 0)
    return (void*)&this->state->lgar_calib_params.Ksat_3;

  else if (name.compare("ponded_depth_max") == 0)
    return (void*)&this->state->lgar_calib_params.ponded_depth_max;
  else if (name.compare("field_capacity") == 0)
    return (void*)&this->state->lgar_calib_params.field_capacity_psi;
  else if (name.compare("a") == 0)
    return (void*)&this->state->lgar_calib_params.a;
  else if (name.compare("b") == 0)
    return (void*)&this->state->lgar_calib_params.b;
  else if (name.compare("frac_to_CR") == 0)
    return (void*)&this->state->lgar_calib_params.frac_to_CR;
  else if (name.compare("a_slow") == 0)
    return (void*)&this->state->lgar_calib_params.a_slow;
  else if (name.compare("b_slow") == 0)
    return (void*)&this->state->lgar_calib_params.b_slow;
  else if (name.compare("frac_slow") == 0)
    return (void*)&this->state->lgar_calib_params.frac_slow;
  else if (name.compare("spf_factor") == 0)
    return (void*)&this->state->lgar_calib_params.spf_factor;
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

// helper functions that enable GIUH queue saving
double* BmiLGAR::get_giuh_runoff_queue()
{
    return giuh_runoff_queue;
}

int BmiLGAR::get_num_giuh_ordinates()
{
    return num_giuh_ordinates;
}

#endif
