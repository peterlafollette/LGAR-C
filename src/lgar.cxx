#ifndef LGAR_CXX_INCLUDED
#define LGAR_CXX_INCLUDED

#include "../include/all.hxx"
#include <iostream>
#include <fstream>
#include <string.h>
#include <sstream>

using namespace std;


//#####################################################################################
/* authors : Ahmad Jan, Fred Ogden, and Peter La Follette
   year    : 2022
   email   : ahmad.jan@noaa.gov, peter.lafollette@noaa.gov
   - The file constains lgar subroutines */
//#####################################################################################


/*##################################################*/
//--------------------------------------------------
//
//     LL          GGG
//     LL        GGG GGG       AAA  A    RR   RRRR
//     LL       GG    GGG    AAA AAAA     RR RR  RR
//     LL       GG     GG   AA     AA     RRR
//     LL      GGG    GGG  AAA     AA     RR
//     LL       GG  GG GG   AA     AAA    RR
//     LL        GGGG  GG    AAA  AA A    RR
//     LL              GG      AAAA   AA  RR
//     LL              GG
//     LLLLLLLL  GG   GG
//                 GGGG
//
// The lgar specific functions based on soil physics
//-----------------------------------------------------
//
//            SKETCH SHOWING 3 SOIL LAYERS AND 4 WETTING FRONTS
//
//    theta_r                                  theta1             theta_e
//      --------------------------------------------------------------   --depth_cm = 0    -------> theta
//      r                                        f                   e
//      r                                        f                   e
//      r                                --------1                   e   --depth_cm(f1)
//      r                               f         \                  e
//      r     1st soil layer            f          \                 e
//      r                               f wetting front number       e
//      r                               f /                          e
//      r                               f/                           e
//      --------------------------------2-----------------------------  -- depth_cm(f2)
//         r                                      f              e
//         r     2nd soil layer                   f              e
//         r                                      f              e
//  |      ---------------------------------------3---------------      -- depth_cm(f3)
//  |        r                                          f     e
//  |        r                                          f     e
//  |        r       3rd  soil layer                    f     e
//  |        r                                          f     e
//  V        -------------------------------------------4------         -- depth_cm(f4)
// depth

#define THRESHOLD_NO_MOISTURE_DIFF 1.E-15      // threshold that will be used to check if adjacent WFs are redundant
#define CREATION_COLOCATED_TOLERANCE_CM 0.0//1.E-8
#define MBAL_ITERATIVE_TOLERANCE 1.E-10        // in the loops that close mass balance across multiple layers, the before and after masses (considering fluxes as well) must match by this number or less
#define MAX_ITER_MBAL_LOOP 1.E5                // the loop that adjusts theta after WFs move (shich that psi will be equal across soil layer boundaries) will iterate this many times before accepting a mass balance error.
#define MAX_ITER_SATURATION_MBAL_LOOP 1.E4     // the loop that adjusts the depth of a saturated WF after the WF becomes saturated to conserve mass has this maximum number of iterations before it just accepts that there will be a mass balance error
#define TRUNCATION_DEPTH 1.E-9                 // when a WF exceeds the lower boundary, we want it to only slightly do this in order to keep the lower boundary condition effectively no flow but also correctly set psi for WFs that it passed
#define LOWER_BOUNDARY_FINAL_TOL_CM 1.E-6      // end-of-step lower-boundary assertion tolerance; ignores intentional epsilon overshoots from internal boundary handling
#define FACTOR_LIMITS_LAYER_CROSSING_SPEED 2.0 // when a WF crosses a layer boundary, it shouldn't go too far into the next layer -- for example in the case of sand over clay, a WF in sand might have a large dzdt value that leads to crossing to an unrealistic depth in the clay below
#define DEPTH_AVOIDS_SAME_WF_DEPTH 1.E-6       // in the event that multiple WFs all would cross a layer boundary and would each have their depth in the new layer limited by FACTOR_LIMITS_LAYER_CROSSING_SPEED, this just prevents these WFs from being exactly at the same depth.
#define PSI_UPPER_LIM 1.E7                     // in loops that close the mass balance by iterating theta and psi, we impose an upper limit on capillary head because some values are just not physically realistic
#define ZERO_DEPTH_TO_PSI_CAP_CM 1.E6          // zero-depth TO fronts with larger finite capillary heads are capped before leaving wetting-front motion
#define ZERO_DEPTH_TO_DELETE_DEPTH_TOL_CM 1.E-10
#define ZERO_DEPTH_TO_DELETE_MASS_TOL_CM 1.E-10
#define FINITE_SAME_LAYER_TO_DELETE_MASS_TOL_CM 1.E-4 // finite duplicate deletes can re-snap nearby to_bottom fronts by roundoff-scale amounts
#define ZERO_DEPTH_TO_DRY_SUPPORT_MAX_SE 0.1   // only delete zero-depth TO supports when they are close to residual saturation
#define TO_PSI_GAP_REFINEMENT_MIN_GAP_CM 5.0
#define TO_PSI_GAP_REFINEMENT_RELATIVE_GAP_FACTOR 0.25
#define TO_PSI_GAP_REFINEMENT_MASS_TOLERANCE_CM 1.E-8
#define TO_PSI_GAP_REFINEMENT_ZERO_DEPTH_DRY_FACTOR 5.0
#define TO_PSI_GAP_REFINEMENT_ZERO_DEPTH_MAX_CM 1.E-6
#define NO_SURFACE_TO_PSI_SPAN_INCREMENT_FACTOR 0.25
#define NO_SURFACE_TO_PSI_SPAN_MAX_INSERTIONS 4
#define NO_SURFACE_TO_PSI_SPAN_TOLERANCE_CM 1.E-8
#define NO_SURFACE_TO_PSI_GAP_REFINEMENT_MAX_INSERTIONS 3
#define SURFACE_PRESENT_TO_PSI_GAP_REFINEMENT_MAX_INSERTIONS 3
#define TO_HYDROSTATIC_DAMPING_MIN_DEPTH_FRACTION 0.8
#define TO_HYDROSTATIC_DAMPING_FRACTION 0.5
#define GW_FLUX_MBAL_CORRECTION_DEBUG_THRESHOLD_CM 1.0
#define BOUNDARY_NEAR_SATURATION_PSI_SNAP_MAX_CM 1.E-2
#define BOUNDARY_STORAGE_NEUTRAL_THETA_TOL 1.E-12
#define MOBILE_GROUNDWATER_SUBMERGENCE_TOL_CM 1.E-6

static bool deferred_gw_flux_mbal_correction_violation = false;
static double deferred_gw_flux_mbal_correction_cm = 0.0;
static const char *deferred_gw_flux_mbal_correction_routine = NULL;
static const char *deferred_gw_flux_mbal_correction_reason = NULL;

static void lgarto_clear_deferred_gw_flux_mass_balance_correction(void)
{
  deferred_gw_flux_mbal_correction_violation = false;
  deferred_gw_flux_mbal_correction_cm = 0.0;
  deferred_gw_flux_mbal_correction_routine = NULL;
  deferred_gw_flux_mbal_correction_reason = NULL;
}

static void lgarto_cap_zero_depth_TO_psi(int num_layers,
                                         int *soil_type,
                                         double *frozen_factor,
                                         struct soil_properties_ *soil_properties,
                                         struct wetting_front *head);
static double lgar_surface_creation_downstream_TO_psi_bound_cm(struct wetting_front *repair_target);
static double lgar_theta_mass_balance_correction_with_min_psi(
  bool use_dry_over_wet, int front_num, double prior_mass, struct wetting_front** head,
  double *cum_layer_thickness_cm, int *soil_type,
  struct soil_properties_ *soil_properties, double min_psi_cm);
static double lgarto_repair_negative_depth_fronts_to_lower_boundary_flux(
  int num_layers, double *cum_layer_thickness_cm, int *soil_type,
  double *frozen_factor, struct soil_properties_ *soil_properties,
  struct wetting_front **head, const char *context);
static double lgarto_surface_fronts_cross_layer_boundary_upward(
  int num_layers, double *cum_layer_thickness_cm, int *soil_type,
  double *frozen_factor, struct soil_properties_ *soil_properties,
  struct wetting_front **head, const char *context);

static double lgar_boundary_roundtrip_psi_tolerance_cm(double psi_a_cm, double psi_b_cm)
{
  const double psi_scale_cm = fmax(1.0, fmax(std::fabs(psi_a_cm), std::fabs(psi_b_cm)));

  return fmax(1.E-3, 1.E-5 * psi_scale_cm);
}

static bool lgar_boundary_snap_preserves_upper_storage(const struct wetting_front *upper,
                                                       const struct wetting_front *lower,
                                                       int *soil_type,
                                                       struct soil_properties_ *soil_properties)
{
  if (upper == NULL || lower == NULL || !upper->to_bottom ||
      upper->layer_num == lower->layer_num || lower->psi_cm < 0.0 ||
      lower->psi_cm > BOUNDARY_NEAR_SATURATION_PSI_SNAP_MAX_CM ||
      !std::isfinite(lower->psi_cm)) {
    return false;
  }

  const int soil_num = soil_type[upper->layer_num];
  const double theta_e = soil_properties[soil_num].theta_e;
  const double theta_r = soil_properties[soil_num].theta_r;
  const double vg_a = soil_properties[soil_num].vg_alpha_per_cm;
  const double vg_m = soil_properties[soil_num].vg_m;
  const double vg_n = soil_properties[soil_num].vg_n;
  const double theta_tol =
    fmax(BOUNDARY_STORAGE_NEUTRAL_THETA_TOL * fmax(1.0, fabs(theta_e - theta_r)),
         1.E-14);

  /* Some soil water retention curves are numerically flat right at saturation:
     a tiny but nonzero psi can evaluate to theta_e exactly in the upper soil,
     while the lower soil still resolves it as barely unsaturated. In that case
     preserving the lower front's psi on the upper to_bottom front restores the
     layer-boundary continuity invariant without changing upper-layer storage. */
  const double upper_theta_at_lower_psi =
    calc_theta_from_h(lower->psi_cm, vg_a, vg_m, vg_n, theta_e, theta_r);

  return fabs(upper->theta - theta_e) <= theta_tol &&
         fabs(upper_theta_at_lower_psi - theta_e) <= theta_tol &&
         fabs(upper_theta_at_lower_psi - upper->theta) <= theta_tol;
}

static bool lgar_should_snap_boundary_psi_to_lower(const struct wetting_front *upper,
                                                   const struct wetting_front *lower,
                                                   int *soil_type,
                                                   struct soil_properties_ *soil_properties)
{
  const double psi_mismatch_cm = std::fabs(upper->psi_cm - lower->psi_cm);
  const double roundtrip_psi_tol_cm =
    lgar_boundary_roundtrip_psi_tolerance_cm(upper->psi_cm, lower->psi_cm);

  return psi_mismatch_cm <= roundtrip_psi_tol_cm ||
         lgar_boundary_snap_preserves_upper_storage(upper, lower, soil_type, soil_properties);
}

static void lgarto_assert_gw_flux_mass_balance_correction_within_debug_threshold(
  double correction_cm,
  const char *routine,
  const char *reason,
  struct wetting_front *head)
{
  (void)head;

  if (!std::isfinite(correction_cm) ||
      std::fabs(correction_cm) > GW_FLUX_MBAL_CORRECTION_DEBUG_THRESHOLD_CM) {
    if (!deferred_gw_flux_mbal_correction_violation) {
      deferred_gw_flux_mbal_correction_violation = true;
      deferred_gw_flux_mbal_correction_cm = correction_cm;
      deferred_gw_flux_mbal_correction_routine = routine;
      deferred_gw_flux_mbal_correction_reason = reason;
    }
  }
}

extern void lgarto_abort_if_deferred_gw_flux_mass_balance_correction_exceeded(
  struct wetting_front *head)
{
  if (!deferred_gw_flux_mbal_correction_violation) {
    return;
  }

  fprintf(stderr,
          "Error: wetting-front editing routine routed a large mass-balance "
          "correction/residual to lower-boundary flux.\n"
          "  Deferred until end of substep after first threshold exceedance.\n"
          "  routine=%s reason=%s\n"
          "  |flux_correction_cm|=%.16g threshold_cm=%.16g flux_correction_cm=%.16g\n"
          "  Wetting front list follows:\n",
          deferred_gw_flux_mbal_correction_routine != NULL ?
            deferred_gw_flux_mbal_correction_routine : "(unknown)",
          deferred_gw_flux_mbal_correction_reason != NULL ?
            deferred_gw_flux_mbal_correction_reason : "(unspecified)",
          std::fabs(deferred_gw_flux_mbal_correction_cm),
          GW_FLUX_MBAL_CORRECTION_DEBUG_THRESHOLD_CM,
          deferred_gw_flux_mbal_correction_cm);
  fflush(stderr);
  listPrint(head);
  fflush(stdout);
  abort();
}

struct lgarto_surface_TO_merge_front_state
{
  struct wetting_front *front;
  double psi_cm;
  double theta;
  double K_cm_per_h;
};

static void lgarto_restore_front_hydraulic_states(
  const std::vector<lgarto_surface_TO_merge_front_state> &states)
{
  for (size_t i = 0; i < states.size(); i++) {
    states[i].front->psi_cm = states[i].psi_cm;
    states[i].front->theta = states[i].theta;
    states[i].front->K_cm_per_h = states[i].K_cm_per_h;
  }
}

static bool lgarto_update_front_from_shared_psi(struct wetting_front *front,
                                                double psi_cm,
                                                int num_layers,
                                                int *soil_type,
                                                double *frozen_factor,
                                                struct soil_properties_ *soil_properties)
{
  if (front == NULL || front->layer_num < 1 || front->layer_num > num_layers ||
      !std::isfinite(psi_cm) || psi_cm < 0.0) {
    return false;
  }

  const int soil_num = soil_type[front->layer_num];
  const double theta_e = soil_properties[soil_num].theta_e;
  const double theta_r = soil_properties[soil_num].theta_r;
  const double vg_a = soil_properties[soil_num].vg_alpha_per_cm;
  const double vg_m = soil_properties[soil_num].vg_m;
  const double vg_n = soil_properties[soil_num].vg_n;
  const double frozen = frozen_factor != NULL ? frozen_factor[front->layer_num] : 1.0;
  const double Ksat_cm_per_h = soil_properties[soil_num].Ksat_cm_per_h * frozen;

  front->psi_cm = psi_cm;
  front->theta = calc_theta_from_h(psi_cm, vg_a, vg_m, vg_n, theta_e, theta_r);
  const double Se = calc_Se_from_theta(front->theta, theta_e, theta_r);
  front->K_cm_per_h = calc_K_from_Se(Se, Ksat_cm_per_h, vg_m);

  return std::isfinite(front->theta) && std::isfinite(front->K_cm_per_h) &&
         front->theta >= theta_r - 1.0e-12 && front->theta <= theta_e + 1.0e-12;
}

static bool lgarto_apply_shared_psi_to_fronts(
  const std::vector<struct wetting_front *> &fronts,
  double psi_cm,
  int num_layers,
  int *soil_type,
  double *frozen_factor,
  struct soil_properties_ *soil_properties)
{
  for (size_t i = 0; i < fronts.size(); i++) {
    if (!lgarto_update_front_from_shared_psi(fronts[i], psi_cm, num_layers, soil_type,
                                             frozen_factor, soil_properties)) {
      return false;
    }
  }
  return true;
}

static bool lgarto_psis_match_for_local_chain(double psi_a_cm, double psi_b_cm)
{
  const double psi_scale_cm = fmax(1.0, fmax(fabs(psi_a_cm), fabs(psi_b_cm)));
  return fabs(psi_a_cm - psi_b_cm) <= fmax(1.0e-8, 1.0e-8 * psi_scale_cm);
}

static std::vector<struct wetting_front *>
lgarto_collect_surface_TO_merge_psi_repair_chain(struct wetting_front *head,
                                                 struct wetting_front *surface_front)
{
  std::vector<struct wetting_front *> fronts;
  if (surface_front == NULL) {
    return fronts;
  }

  struct wetting_front *previous = NULL;
  for (struct wetting_front *probe = head; probe != NULL && probe != surface_front;
       probe = probe->next) {
    previous = probe;
  }

  /* Include the adjacent upper to_bottom scaffold when it already shares psi
     with the boundary-pinned surface front. Otherwise a later theta snap sees
     that connection and changes storage outside this mass solve. */
  if (previous != NULL && previous->to_bottom && previous->layer_num != surface_front->layer_num &&
      lgarto_psis_match_for_local_chain(previous->psi_cm, surface_front->psi_cm)) {
    fronts.push_back(previous);
  }

  fronts.push_back(surface_front);

  /* If this surface front already shares a capillary head with an adjacent
     to_bottom/lower-front boundary chain, keep that local chain tied together.
     In the common failure mode this list contains only the boundary-pinned
     surface front; the separate to_bottom/lower-front pair already has its own
     continuity and should not be broadened unless it is actually connected by psi. */
  struct wetting_front *current = surface_front;
  const double shared_psi_cm = surface_front->psi_cm;
  while (current->next != NULL && current->next->to_bottom &&
         lgarto_psis_match_for_local_chain(shared_psi_cm, current->next->psi_cm)) {
    current = current->next;
    fronts.push_back(current);
  }

  if (current != surface_front && current->next != NULL && !current->next->to_bottom &&
      lgarto_psis_match_for_local_chain(shared_psi_cm, current->next->psi_cm)) {
    fronts.push_back(current->next);
  }

  return fronts;
}

static bool lgarto_surface_TO_merge_psi_repair_preserves_local_ordering(
  struct wetting_front *head)
{
  for (struct wetting_front *current = head; current != NULL && current->next != NULL;
       current = current->next) {
    struct wetting_front *next = current->next;
    if (current->layer_num != next->layer_num) {
      continue;
    }

    if (!current->is_WF_GW && !next->is_WF_GW &&
        current->theta + 1.0e-10 < next->theta) {
      return false;
    }

    if (current->is_WF_GW && next->is_WF_GW &&
        current->psi_cm + 1.0e-10 < next->psi_cm) {
      return false;
    }
  }

  return true;
}

extern double lgarto_restore_surface_TO_merge_mass_via_boundary_psi(
  double target_mass_cm,
  int boundary_pinned_surface_front_num,
  int num_layers,
  double *cum_layer_thickness_cm,
  int *soil_type,
  double *frozen_factor,
  struct soil_properties_ *soil_properties,
  struct wetting_front *head)
{
  double residual_cm = target_mass_cm - lgar_calc_mass_bal(cum_layer_thickness_cm, head);
  if (residual_cm <= MBAL_ITERATIVE_TOLERANCE || head == NULL ||
      boundary_pinned_surface_front_num < 1) {
    return residual_cm;
  }

  struct wetting_front *surface_front =
    listFindFront(boundary_pinned_surface_front_num, head, NULL);
  if (surface_front == NULL || surface_front->is_WF_GW || surface_front->to_bottom ||
      surface_front->layer_num < 1 || surface_front->layer_num > num_layers ||
      !std::isfinite(surface_front->psi_cm) || surface_front->psi_cm <= 0.0) {
    return residual_cm;
  }

  const double layer_bottom_cm = cum_layer_thickness_cm[surface_front->layer_num];
  const double boundary_tol_cm = fmax(1.0e-8, 10.0 * TRUNCATION_DEPTH);
  if (fabs(surface_front->depth_cm - layer_bottom_cm) > boundary_tol_cm) {
    return residual_cm;
  }

  std::vector<struct wetting_front *> repair_fronts =
    lgarto_collect_surface_TO_merge_psi_repair_chain(head, surface_front);
  if (repair_fronts.empty()) {
    return residual_cm;
  }

  std::vector<lgarto_surface_TO_merge_front_state> saved_states;
  for (size_t i = 0; i < repair_fronts.size(); i++) {
    lgarto_surface_TO_merge_front_state saved = {
      repair_fronts[i],
      repair_fronts[i]->psi_cm,
      repair_fronts[i]->theta,
      repair_fronts[i]->K_cm_per_h
    };
    saved_states.push_back(saved);
  }

  const double original_psi_cm = surface_front->psi_cm;

  /* The depth repair has exhausted geometry: the surface front is already at
     the layer bottom. Lowering psi wets that local profile volume without
     moving deeper TO fronts or inventing a lower-boundary flux correction. */
  if (!lgarto_apply_shared_psi_to_fronts(repair_fronts, 0.0, num_layers, soil_type,
                                         frozen_factor, soil_properties) ||
      !lgarto_surface_TO_merge_psi_repair_preserves_local_ordering(head)) {
    lgarto_restore_front_hydraulic_states(saved_states);
    return residual_cm;
  }

  const double saturated_residual_cm =
    target_mass_cm - lgar_calc_mass_bal(cum_layer_thickness_cm, head);
  if (saturated_residual_cm > MBAL_ITERATIVE_TOLERANCE) {
    if (verbosity.compare("high") == 0) {
      printf("surface/TO merge boundary-psi repair could not close residual %.12e cm "
             "by saturating front %d; remaining residual at saturation %.12e cm\n",
             residual_cm,
             surface_front->front_num,
             saturated_residual_cm);
    }
    lgarto_restore_front_hydraulic_states(saved_states);
    return residual_cm;
  }

  double psi_wet_cm = 0.0;
  double psi_dry_cm = original_psi_cm;
  double best_psi_cm = 0.0;
  double best_abs_residual_cm = fabs(saturated_residual_cm);

  for (int iter = 0; iter < 80; iter++) {
    const double probe_psi_cm = 0.5 * (psi_wet_cm + psi_dry_cm);
    if (!lgarto_apply_shared_psi_to_fronts(repair_fronts, probe_psi_cm, num_layers,
                                           soil_type, frozen_factor, soil_properties) ||
        !lgarto_surface_TO_merge_psi_repair_preserves_local_ordering(head)) {
      psi_wet_cm = probe_psi_cm;
      continue;
    }

    const double probe_residual_cm =
      target_mass_cm - lgar_calc_mass_bal(cum_layer_thickness_cm, head);
    if (fabs(probe_residual_cm) < best_abs_residual_cm) {
      best_abs_residual_cm = fabs(probe_residual_cm);
      best_psi_cm = probe_psi_cm;
    }

    if (fabs(probe_residual_cm) <= MBAL_ITERATIVE_TOLERANCE) {
      best_psi_cm = probe_psi_cm;
      break;
    }

    if (probe_residual_cm > 0.0) {
      psi_dry_cm = probe_psi_cm;
    }
    else {
      psi_wet_cm = probe_psi_cm;
    }
  }

  if (!lgarto_apply_shared_psi_to_fronts(repair_fronts, best_psi_cm, num_layers,
                                         soil_type, frozen_factor, soil_properties) ||
      !lgarto_surface_TO_merge_psi_repair_preserves_local_ordering(head)) {
    lgarto_restore_front_hydraulic_states(saved_states);
    return residual_cm;
  }

  const double repaired_residual_cm =
    target_mass_cm - lgar_calc_mass_bal(cum_layer_thickness_cm, head);
  if (fabs(repaired_residual_cm) > fabs(residual_cm)) {
    lgarto_restore_front_hydraulic_states(saved_states);
    return residual_cm;
  }

  if (verbosity.compare("high") == 0) {
    printf("surface/TO merge boundary-psi repair updated front %d in layer %d: "
           "psi %.17lf cm -> %.17lf cm, theta %.17lf, residual %.12e cm "
           "(connected fronts adjusted: %zu)\n",
           surface_front->front_num,
           surface_front->layer_num,
           original_psi_cm,
           surface_front->psi_cm,
           surface_front->theta,
           repaired_residual_cm,
           repair_fronts.size());
  }

  return repaired_residual_cm;
}

static double lgarto_restore_surface_TO_merge_mass_via_surface_depths(double target_mass_cm,
                                                                      int num_layers,
                                                                      double *cum_layer_thickness_cm,
                                                                      struct wetting_front *head,
                                                                      int *boundary_pinned_surface_front_num)
{
  if (boundary_pinned_surface_front_num != NULL) {
    *boundary_pinned_surface_front_num = -1;
  }

  double residual_cm = target_mass_cm - lgar_calc_mass_bal(cum_layer_thickness_cm, head);
  if (std::fabs(residual_cm) <= MBAL_ITERATIVE_TOLERANCE || head == NULL) {
    return residual_cm;
  }

  for (int pass = 0; pass < listLength(head); pass++) {
    bool repaired_this_pass = false;

    for (struct wetting_front *candidate = head; candidate != NULL; candidate = candidate->next) {
      if (candidate->is_WF_GW || candidate->to_bottom ||
          candidate->layer_num < 1 || candidate->layer_num > num_layers) {
        continue;
      }

      double lower_depth_cm = cum_layer_thickness_cm[candidate->layer_num - 1];
      double upper_depth_cm = cum_layer_thickness_cm[candidate->layer_num];

      struct wetting_front *previous = NULL;
      for (struct wetting_front *probe = head; probe != NULL && probe != candidate; probe = probe->next) {
        previous = probe;
      }

      if (previous != NULL && previous->layer_num == candidate->layer_num) {
        lower_depth_cm = fmax(lower_depth_cm, previous->depth_cm);
      }
      if (candidate->next != NULL && candidate->next->layer_num == candidate->layer_num) {
        upper_depth_cm = fmin(upper_depth_cm, candidate->next->depth_cm);
      }

      if (upper_depth_cm <= lower_depth_cm) {
        continue;
      }

      const double original_depth_cm = candidate->depth_cm;
      if (original_depth_cm < lower_depth_cm || original_depth_cm > upper_depth_cm) {
        continue;
      }

      const double endpoint_depth_cm = residual_cm > 0.0 ? upper_depth_cm : lower_depth_cm;
      if (endpoint_depth_cm == original_depth_cm) {
        continue;
      }

      candidate->depth_cm = endpoint_depth_cm;
      const double endpoint_residual_cm =
        target_mass_cm - lgar_calc_mass_bal(cum_layer_thickness_cm, head);
      candidate->depth_cm = original_depth_cm;

      if (std::fabs(endpoint_residual_cm) >= std::fabs(residual_cm)) {
        continue;
      }

      double bracket_lo_cm = fmin(original_depth_cm, endpoint_depth_cm);
      double bracket_hi_cm = fmax(original_depth_cm, endpoint_depth_cm);
      double best_depth_cm = endpoint_depth_cm;
      double best_abs_residual_cm = std::fabs(endpoint_residual_cm);

      for (int iter = 0; iter < 80; iter++) {
        const double probe_depth_cm = 0.5 * (bracket_lo_cm + bracket_hi_cm);
        candidate->depth_cm = probe_depth_cm;
        const double probe_residual_cm =
          target_mass_cm - lgar_calc_mass_bal(cum_layer_thickness_cm, head);

        if (std::fabs(probe_residual_cm) < best_abs_residual_cm) {
          best_abs_residual_cm = std::fabs(probe_residual_cm);
          best_depth_cm = probe_depth_cm;
        }

        if (std::fabs(probe_residual_cm) <= MBAL_ITERATIVE_TOLERANCE) {
          best_depth_cm = probe_depth_cm;
          break;
        }

        if ((residual_cm > 0.0 && probe_residual_cm > 0.0) ||
            (residual_cm < 0.0 && probe_residual_cm < 0.0)) {
          bracket_lo_cm = residual_cm > 0.0 ? probe_depth_cm : bracket_lo_cm;
          bracket_hi_cm = residual_cm < 0.0 ? probe_depth_cm : bracket_hi_cm;
        }
        else {
          bracket_lo_cm = residual_cm < 0.0 ? probe_depth_cm : bracket_lo_cm;
          bracket_hi_cm = residual_cm > 0.0 ? probe_depth_cm : bracket_hi_cm;
        }
      }

      candidate->depth_cm = best_depth_cm;
      residual_cm = target_mass_cm - lgar_calc_mass_bal(cum_layer_thickness_cm, head);
      repaired_this_pass = true;

      const double layer_bottom_cm = cum_layer_thickness_cm[candidate->layer_num];
      const double boundary_tol_cm = fmax(1.0e-8, 10.0 * TRUNCATION_DEPTH);
      if (residual_cm > MBAL_ITERATIVE_TOLERANCE &&
          fabs(candidate->depth_cm - layer_bottom_cm) <= boundary_tol_cm &&
          boundary_pinned_surface_front_num != NULL) {
        *boundary_pinned_surface_front_num = candidate->front_num;
      }

      if (verbosity.compare("high") == 0) {
        printf("surface/TO merge depth-repair updated front %d in layer %d: "
               "depth %.17lf cm -> %.17lf cm, residual %.12e cm\n",
               candidate->front_num,
               candidate->layer_num,
               original_depth_cm,
               candidate->depth_cm,
               residual_cm);
      }

      if (std::fabs(residual_cm) <= MBAL_ITERATIVE_TOLERANCE) {
        return residual_cm;
      }
    }

    if (!repaired_this_pass) {
      break;
    }
  }

  return residual_cm;
}

static double lgar_limit_TO_dzdt_near_hydrostatic_target(double dzdt_cm_per_h,
                                                         double current_depth_cm,
                                                         double target_depth_cm,
                                                         double column_depth_cm,
                                                         double subtimestep_h)
{
  if (subtimestep_h <= 0.0) {
    return dzdt_cm_per_h;
  }

  if (current_depth_cm < TO_HYDROSTATIC_DAMPING_MIN_DEPTH_FRACTION * column_depth_cm) {
    return dzdt_cm_per_h;
  }

  if (target_depth_cm < 0.0 || target_depth_cm > column_depth_cm) {
    return dzdt_cm_per_h;
  }

  const double projected_depth_cm = current_depth_cm + dzdt_cm_per_h * subtimestep_h;
  const double current_offset_cm = current_depth_cm - target_depth_cm;
  const double projected_offset_cm = projected_depth_cm - target_depth_cm;

  if (std::fabs(current_offset_cm) <= TRUNCATION_DEPTH) {
    return dzdt_cm_per_h;
  }

  if (current_offset_cm * projected_offset_cm < 0.0) {
    // Approach the hydrostatic target smoothly without allowing a one-step sign flip.
    return TO_HYDROSTATIC_DAMPING_FRACTION * (target_depth_cm - current_depth_cm) / subtimestep_h;
  }

  return dzdt_cm_per_h;
}

static bool lgarto_is_last_layer_GW_overshoot_pair(struct wetting_front *current,
                                                   int num_layers)
{
  if (current == NULL || current->next == NULL) {
    return false;
  }

  struct wetting_front *next = current->next;
  return current->is_WF_GW && !current->to_bottom &&
         next->is_WF_GW && next->to_bottom &&
         current->layer_num == num_layers &&
         next->layer_num == num_layers &&
         current->layer_num == next->layer_num &&
         current->depth_cm > next->depth_cm;
}

extern double lgarto_truncate_last_layer_GW_overshoot(double *cum_layer_thickness_cm,
                                                      int num_layers,
                                                      struct wetting_front **head,
                                                      int *soil_type,
                                                      struct soil_properties_ *soil_properties)
{
  if (cum_layer_thickness_cm == NULL || head == NULL || *head == NULL) {
    return 0.0;
  }

  struct wetting_front *overshooting_front = NULL;
  for (struct wetting_front *current = *head; current != NULL; current = current->next) {
    if (lgarto_is_last_layer_GW_overshoot_pair(current, num_layers)) {
      overshooting_front = current;
    }
  }

  if (overshooting_front == NULL) {
    return 0.0;
  }

  struct wetting_front *bottom_front = overshooting_front->next;
  const double lower_boundary_depth_cm = cum_layer_thickness_cm[num_layers];
  const double geometric_bottom_flux_cm =
    (overshooting_front->depth_cm - bottom_front->depth_cm) *
    (overshooting_front->theta - bottom_front->theta);

  struct lgarto_saved_to_bottom_state
  {
    struct wetting_front *front;
    bool is_WF_GW;
    double theta;
    double psi_cm;
    double K_cm_per_h;
  };

  std::vector<lgarto_saved_to_bottom_state> upstream_boundary_states;
  for (struct wetting_front *current = *head;
       current != NULL && current != overshooting_front;
       current = current->next) {
    if (current->to_bottom) {
      upstream_boundary_states.push_back({
        current,
        current->is_WF_GW,
        current->theta,
        current->psi_cm,
        current->K_cm_per_h
      });
    }
  }

  const double mass_before_cm = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);

  if (verbosity.compare("high") == 0) {
    printf("TO/general correction type 9 truncating last-layer GW front %d into bottom front %d.\n",
           overshooting_front->front_num, bottom_front->front_num);
    printf("  overshoot_depth_cm=%.15f bottom_depth_cm=%.15f lower_boundary_depth_cm=%.15f\n",
           overshooting_front->depth_cm, bottom_front->depth_cm, lower_boundary_depth_cm);
    printf("  overshoot_theta=%.15f bottom_theta=%.15f geometric_bottom_flux_cm=%.15f\n",
           overshooting_front->theta, bottom_front->theta, geometric_bottom_flux_cm);
  }

  bottom_front->depth_cm = lower_boundary_depth_cm;
  bottom_front->theta = overshooting_front->theta;
  bottom_front->psi_cm = overshooting_front->psi_cm;
  bottom_front->K_cm_per_h = overshooting_front->K_cm_per_h;
  bottom_front->is_WF_GW = TRUE;
  bottom_front->to_bottom = TRUE;

  listDeleteFront(overshooting_front->front_num, head, soil_type, soil_properties);

  /* listDeleteFront globally re-snaps every upstream to_bottom front to the
     next front's psi. For this correction, the overshooting last-layer front is
     a local layer-crossing artifact; deleting it should not dry or wet unrelated
     upper boundary fronts. */
  for (size_t i = 0; i < upstream_boundary_states.size(); i++) {
    upstream_boundary_states[i].front->is_WF_GW = upstream_boundary_states[i].is_WF_GW;
    upstream_boundary_states[i].front->theta = upstream_boundary_states[i].theta;
    upstream_boundary_states[i].front->psi_cm = upstream_boundary_states[i].psi_cm;
    upstream_boundary_states[i].front->K_cm_per_h = upstream_boundary_states[i].K_cm_per_h;
  }

  const double bottom_flux_cm =
    mass_before_cm - lgar_calc_mass_bal(cum_layer_thickness_cm, *head);

  if (verbosity.compare("high") == 0) {
    printf("  routed_bottom_flux_cm=%.15f\n", bottom_flux_cm);
    printf("State after truncating last-layer GW overshoot into the bottom boundary...\n");
    listPrint(*head);
  }

  return bottom_flux_cm;
}

extern double lgarto_clip_final_layer_GW_overshoot_to_vadose_boundary(double *cum_layer_thickness_cm,
                                                                      int num_layers,
                                                                      struct wetting_front **head)
{
  if (cum_layer_thickness_cm == NULL || head == NULL || *head == NULL || num_layers < 1) {
    return 0.0;
  }

  const double lower_boundary_depth_cm = cum_layer_thickness_cm[num_layers];
  if (!std::isfinite(lower_boundary_depth_cm) || lower_boundary_depth_cm <= 0.0) {
    return 0.0;
  }

  const double final_tol_cm = LOWER_BOUNDARY_FINAL_TOL_CM;
  const double max_accepted_overshoot_cm =
    fmax(final_tol_cm, 0.001 * lower_boundary_depth_cm);
  double bottom_flux_cm = 0.0;

  for (struct wetting_front *current = *head; current != NULL; current = current->next) {
    if (!current->is_WF_GW || current->to_bottom || current->layer_num != num_layers ||
        !std::isfinite(current->depth_cm) ||
        current->depth_cm <= lower_boundary_depth_cm + final_tol_cm) {
      continue;
    }

    const double overshoot_cm = current->depth_cm - lower_boundary_depth_cm;
    if (overshoot_cm > max_accepted_overshoot_cm + final_tol_cm) {
      continue;
    }

    const double original_depth_cm = current->depth_cm;
    const double mass_before_cm = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
    current->depth_cm = lower_boundary_depth_cm;
    double clipped_flux_cm = mass_before_cm - lgar_calc_mass_bal(cum_layer_thickness_cm, *head);

    if (clipped_flux_cm < -MBAL_ITERATIVE_TOLERANCE) {
      current->depth_cm = original_depth_cm;
      if (verbosity.compare("high") == 0) {
        printf("Final lower-boundary TO/GW clip rejected for front %d because clipping "
               "would imply negative recharge %.12e cm.\n",
               current->front_num, clipped_flux_cm);
      }
      continue;
    }

    if (clipped_flux_cm < 0.0) {
      clipped_flux_cm = 0.0;
    }
    bottom_flux_cm += clipped_flux_cm;

    if (verbosity.compare("high") == 0) {
      printf("Final lower-boundary TO/GW clip moved front %d from %.17lf cm to "
             "%.17lf cm and routed %.12e cm to bottom flux.\n",
             current->front_num,
             original_depth_cm,
             current->depth_cm,
             clipped_flux_cm);
    }
  }

  return bottom_flux_cm;
}

// ############################################################################################
/*
  @param soil_moisture_wetting_fronts  : 1D array of thetas (soil moisture content) per wetting front;
                                         output to other models (e.g. soil freeze-thaw)
  @param soil_depth_wetting_fronts : 1D array of absolute depths of the wetting fronts [meters];
					 output to other models (e.g. soil freeze-thaw)
*/
// ############################################################################################
extern void lgar_initialize(string config_file, struct model_state *state)
{
  int soil;
  
  InitFromConfigFile(config_file, state);
  state->lgar_bmi_params.shape[0] = state->lgar_bmi_params.num_layers;
  state->lgar_bmi_params.shape[1] = state->lgar_bmi_params.num_wetting_fronts;

  state->lgar_bmi_params.soil_depth_wetting_fronts    = new double[state->lgar_bmi_params.num_wetting_fronts];
  state->lgar_bmi_params.soil_moisture_wetting_fronts = new double[state->lgar_bmi_params.num_wetting_fronts];

  // initialize array for holding calibratable parameters
  // calibratabale parameters are scalars now not arrays
  // state->lgar_calib_params.theta_e  = new double[state->lgar_bmi_params.num_layers];
  // state->lgar_calib_params.theta_r  = new double[state->lgar_bmi_params.num_layers];
  // state->lgar_calib_params.vg_n     = new double[state->lgar_bmi_params.num_layers];
  // state->lgar_calib_params.vg_alpha = new double[state->lgar_bmi_params.num_layers];
  // state->lgar_calib_params.Ksat     = new double[state->lgar_bmi_params.num_layers];
  
  // initialize thickness/depth and soil moisture of wetting fronts (used for model coupling)
  // also initialize calibratable parameters
  state->lgar_calib_params.field_capacity_psi = state->lgar_bmi_params.field_capacity_psi_cm;
  state->lgar_calib_params.ponded_depth_max = state->lgar_bmi_params.ponded_depth_max_cm;
  state->lgar_calib_params.a = state->lgar_bmi_params.a;
  state->lgar_calib_params.b = state->lgar_bmi_params.b;
  state->lgar_calib_params.frac_to_CR = state->lgar_bmi_params.frac_to_CR;
  state->lgar_calib_params.a_slow = state->lgar_bmi_params.a_slow;
  state->lgar_calib_params.b_slow = state->lgar_bmi_params.b_slow;
  state->lgar_calib_params.frac_slow = state->lgar_bmi_params.frac_slow;
  state->lgar_calib_params.spf_factor = state->lgar_bmi_params.spf_factor;

  struct wetting_front *current = state->head;
  bool layer_1_params_set = false;
  bool layer_2_params_set = false;
  bool layer_3_params_set = false;
  for (int i=0; i<state->lgar_bmi_params.num_wetting_fronts; i++) {
    assert (current != NULL);
    
    soil = state->lgar_bmi_params.layer_soil_type[current->layer_num];

    state->lgar_bmi_params.soil_moisture_wetting_fronts[i] = current->theta;
    state->lgar_bmi_params.soil_depth_wetting_fronts[i]    = current->depth_cm * state->units.cm_to_m;

    // // we now handle calibration of layered parameters with scalars
    // state->lgar_calib_params.theta_e[i]  = state->soil_properties[soil].theta_e;
    // state->lgar_calib_params.theta_r[i]  = state->soil_properties[soil].theta_r;
    // state->lgar_calib_params.vg_n[i]     = state->soil_properties[soil].vg_n;
    // state->lgar_calib_params.vg_alpha[i] = state->soil_properties[soil].vg_alpha_per_cm;
    // state->lgar_calib_params.Ksat[i]     = state->soil_properties[soil].Ksat_cm_per_h;

    if (current->layer_num == 1 && !layer_1_params_set){
      state->lgar_calib_params.theta_e_1  = state->soil_properties[soil].theta_e;
      state->lgar_calib_params.theta_r_1  = state->soil_properties[soil].theta_r;
      state->lgar_calib_params.vg_n_1     = state->soil_properties[soil].vg_n;
      state->lgar_calib_params.vg_alpha_1 = state->soil_properties[soil].vg_alpha_per_cm;
      state->lgar_calib_params.Ksat_1     = state->soil_properties[soil].Ksat_cm_per_h;
      layer_1_params_set = true;
    }

    if (current->layer_num == 2 && !layer_2_params_set){
      state->lgar_calib_params.theta_e_2  = state->soil_properties[soil].theta_e;
      state->lgar_calib_params.theta_r_2  = state->soil_properties[soil].theta_r;
      state->lgar_calib_params.vg_n_2     = state->soil_properties[soil].vg_n;
      state->lgar_calib_params.vg_alpha_2 = state->soil_properties[soil].vg_alpha_per_cm;
      state->lgar_calib_params.Ksat_2     = state->soil_properties[soil].Ksat_cm_per_h;
      layer_2_params_set = true;
    }

    if (current->layer_num == 3 && !layer_3_params_set){
      state->lgar_calib_params.theta_e_3  = state->soil_properties[soil].theta_e;
      state->lgar_calib_params.theta_r_3  = state->soil_properties[soil].theta_r;
      state->lgar_calib_params.vg_n_3     = state->soil_properties[soil].vg_n;
      state->lgar_calib_params.vg_alpha_3 = state->soil_properties[soil].vg_alpha_per_cm;
      state->lgar_calib_params.Ksat_3     = state->soil_properties[soil].Ksat_cm_per_h;
      layer_3_params_set = true;
    }
    
    current = current->next;
  }


  /* initialize bmi input variables to -1.0 (on purpose), this should be assigned (non-negative) and if not,
     the code will throw an error in the Update method */
  state->lgar_bmi_input_params->precipitation_mm_per_h = -1.0;
  state->lgar_bmi_input_params->PET_mm_per_h           = -1.0;

  // // initialize all global mass balance variables to zero
  // state->lgar_mass_balance.volprecip_cm              = 0.0;
  // state->lgar_mass_balance.volin_cm                  = 0.0;
  // state->lgar_mass_balance.volend_cm                 = 0.0;
  // state->lgar_mass_balance.volCRend_cm               = 0.0;
  // state->lgar_mass_balance.volAET_cm                 = 0.0;
  // state->lgar_mass_balance.volrech_cm                = 0.0;
  // state->lgar_mass_balance.volrunoff_cm              = 0.0;
  // state->lgar_mass_balance.volrunoff_giuh_cm         = 0.0;
  // state->lgar_mass_balance.volQ_cm                   = 0.0;
  // state->lgar_mass_balance.volQ_CR_cm                = 0.0;
  // state->lgar_mass_balance.volPET_cm                 = 0.0;
  // state->lgar_mass_balance.volon_cm                  = 0.0;
  // state->lgar_mass_balance.volon_timestep_cm         = 0.0; /* setting volon and precip at the initial time to 0.0
	// 						       as they determine the creation of surficail wetting front */
  // state->lgar_bmi_params.precip_previous_timestep_cm = 0.0;
  // state->lgar_mass_balance.volchange_calib_cm        = 0.0;

  bool non_vadose_restart_loaded =
    !state->lgar_bmi_params.init_state_path.empty() &&
    !state->lgar_bmi_params.init_non_vadose_state_path.empty();

  state->lgar_mass_balance.volprecip_cm      = 0.0;
  state->lgar_mass_balance.volin_cm          = 0.0;
  state->lgar_mass_balance.volend_cm         = 0.0;
  state->lgar_mass_balance.volAET_cm         = 0.0;
  state->lgar_mass_balance.volrech_cm        = 0.0;
  state->lgar_mass_balance.volrunoff_cm      = 0.0;
  state->lgar_mass_balance.volrunoff_giuh_cm = 0.0;
  state->lgar_mass_balance.volQ_cm           = 0.0;
  state->lgar_mass_balance.volQ_CR_cm        = 0.0;
  state->lgar_mass_balance.volpref_flow_to_CR_cm = 0.0;
  state->lgar_mass_balance.vollgarto_domain_to_CR_cm = 0.0;
  state->lgar_mass_balance.volPET_cm         = 0.0;
  state->lgar_mass_balance.volon_cm          = 0.0;
  state->lgar_mass_balance.volchange_calib_cm = 0.0;

  if (!non_vadose_restart_loaded) {
    state->lgar_mass_balance.CR_fast_storage_cm =
      state->lgar_bmi_params.initial_CR_fast_storage_cm;
    state->lgar_mass_balance.CR_slow_storage_cm =
      state->lgar_bmi_params.initial_CR_slow_storage_cm;
    const double initial_CR_storage_cm =
      state->lgar_mass_balance.CR_fast_storage_cm +
      state->lgar_mass_balance.CR_slow_storage_cm;

    state->lgar_mass_balance.volCRend_cm = initial_CR_storage_cm;
    state->lgar_mass_balance.volCRend_timestep_cm = initial_CR_storage_cm;
    state->lgar_mass_balance.volpref_flow_to_CR_timestep_cm = 0.0;
    state->lgar_mass_balance.vollgarto_domain_to_CR_timestep_cm = 0.0;
    state->lgar_mass_balance.volon_timestep_cm = 0.0;
    state->lgar_bmi_params.precip_previous_timestep_cm = 0.0;
    state->lgar_bmi_params.runoff_in_prev_step = false;
    state->lgar_mass_balance.volCRstart_cm = initial_CR_storage_cm;
  }
  else {
    state->lgar_mass_balance.volCRend_cm =
        state->lgar_mass_balance.CR_fast_storage_cm +
        state->lgar_mass_balance.CR_slow_storage_cm;

    state->lgar_mass_balance.volCRend_timestep_cm =
        state->lgar_mass_balance.CR_fast_storage_cm +
        state->lgar_mass_balance.CR_slow_storage_cm;
    state->lgar_mass_balance.volCRstart_cm =
      state->lgar_mass_balance.CR_fast_storage_cm +
      state->lgar_mass_balance.CR_slow_storage_cm;
  }

}

struct TOPsiGapRefinementCandidate {
  struct wetting_front *upper_front;
  struct wetting_front *lower_front;
  double psi_gap_cm;
};

static int lgar_layer_num_from_depth_cm(int num_layers,
                                        double *cum_layer_thickness_cm,
                                        double depth_cm)
{
  for (int layer_num = 1; layer_num <= num_layers; layer_num++) {
    if (depth_cm <= cum_layer_thickness_cm[layer_num] &&
        depth_cm > cum_layer_thickness_cm[layer_num - 1]) {
      return layer_num;
    }
  }

  return -1;
}

static struct wetting_front *lgar_insert_groundwater_front_before(double depth_cm,
                                                                  double psi_new_cm,
                                                                  struct wetting_front *lower_front,
                                                                  int num_layers,
                                                                  int *soil_type,
                                                                  double *cum_layer_thickness_cm,
                                                                  double *frozen_factor,
                                                                  struct soil_properties_ *soil_properties,
                                                                  struct wetting_front **head)
{
  if (lower_front == NULL || head == NULL || *head == NULL) {
    return NULL;
  }

  const int layer_num = lgar_layer_num_from_depth_cm(num_layers, cum_layer_thickness_cm, depth_cm);
  if (layer_num < 1) {
    return NULL;
  }

  const int soil_num = soil_type[layer_num];
  const double theta_e = soil_properties[soil_num].theta_e;
  const double theta_r = soil_properties[soil_num].theta_r;
  const double vg_a = soil_properties[soil_num].vg_alpha_per_cm;
  const double vg_m = soil_properties[soil_num].vg_m;
  const double vg_n = soil_properties[soil_num].vg_n;
  const double theta_new = calc_theta_from_h(psi_new_cm, vg_a, vg_m, vg_n, theta_e, theta_r);

  struct wetting_front *inserted = NULL;
  if (lower_front == *head) {
    listInsertFirst(depth_cm, theta_new, 1, layer_num, FALSE, head);
    inserted = *head;
  }
  else {
    struct wetting_front *previous = *head;
    while (previous != NULL && previous->next != lower_front) {
      previous = previous->next;
    }

    if (previous == NULL) {
      return NULL;
    }

    inserted = (struct wetting_front *) malloc(sizeof(struct wetting_front));
    inserted->depth_cm = depth_cm;
    inserted->theta = theta_new;
    inserted->front_num = lower_front->front_num;
    inserted->layer_num = layer_num;
    inserted->to_bottom = FALSE;
    inserted->dzdt_cm_per_h = 0.0;
    inserted->is_WF_GW = TRUE;
    inserted->next = lower_front;
    previous->next = inserted;

    for (struct wetting_front *current = lower_front; current != NULL; current = current->next) {
      current->front_num++;
    }
  }

  if (inserted == NULL) {
    return NULL;
  }

  inserted->is_WF_GW = TRUE;
  inserted->to_bottom = FALSE;
  inserted->psi_cm = psi_new_cm;
  inserted->dzdt_cm_per_h = 0.0;

  const double Se = calc_Se_from_theta(inserted->theta, theta_e, theta_r);
  const double Ksat_cm_per_h = soil_properties[soil_num].Ksat_cm_per_h * frozen_factor[layer_num];
  inserted->K_cm_per_h = calc_K_from_Se(Se, Ksat_cm_per_h, vg_m);

  return inserted;
}

static struct wetting_front *lgar_find_previous_front(struct wetting_front *head,
                                                      struct wetting_front *target)
{
  if (head == NULL || target == NULL || head == target) {
    return NULL;
  }

  struct wetting_front *previous = head;
  while (previous != NULL && previous->next != target) {
    previous = previous->next;
  }

  return previous;
}

static bool lgarto_insert_zero_depth_groundwater_front(double psi_new_cm,
                                                       int *soil_type,
                                                       double *frozen_factor,
                                                       struct soil_properties_ *soil_properties,
                                                       struct wetting_front **head)
{
  if (head == NULL || *head == NULL || listLength_surface(*head) > 0) {
    return false;
  }

  const int layer_num = 1;
  const int soil_num = soil_type[layer_num];
  const double theta_e = soil_properties[soil_num].theta_e;
  const double theta_r = soil_properties[soil_num].theta_r;
  const double vg_a = soil_properties[soil_num].vg_alpha_per_cm;
  const double vg_m = soil_properties[soil_num].vg_m;
  const double vg_n = soil_properties[soil_num].vg_n;
  const double theta_new = calc_theta_from_h(psi_new_cm, vg_a, vg_m, vg_n, theta_e, theta_r);

  listInsertFirst(0.0, theta_new, 1, layer_num, FALSE, head);
  struct wetting_front *inserted = *head;
  inserted->is_WF_GW = TRUE;
  inserted->to_bottom = FALSE;
  inserted->psi_cm = psi_new_cm;
  inserted->dzdt_cm_per_h = 0.0;

  const double Se = calc_Se_from_theta(theta_new, theta_e, theta_r);
  const double Ksat_cm_per_h = soil_properties[soil_num].Ksat_cm_per_h * frozen_factor[layer_num];
  inserted->K_cm_per_h = calc_K_from_Se(Se, Ksat_cm_per_h, vg_m);
  return true;
}

static int lgarto_ensure_no_surface_to_psi_span(int num_layers,
                                                double *cum_layer_thickness_cm,
                                                int *soil_type,
                                                double *frozen_factor,
                                                struct soil_properties_ *soil_properties,
                                                struct wetting_front **head)
{
  if (head == NULL || *head == NULL || listLength_surface(*head) > 0) {
    return 0;
  }

  if (!(*head)->is_WF_GW) {
    return 0;
  }

  const double column_depth_cm = cum_layer_thickness_cm[num_layers];
  if (column_depth_cm <= 0.0) {
    return 0;
  }

  const double psi_increment_cm = NO_SURFACE_TO_PSI_SPAN_INCREMENT_FACTOR * column_depth_cm;
  // Avoid creating a zero-depth support just to close machine-roundoff gaps
  // between the current dry-end psi and the column-depth psi target.
  const double psi_span_tolerance_cm =
    fmax(NO_SURFACE_TO_PSI_SPAN_TOLERANCE_CM,
         1.0e-12 * column_depth_cm);
  int insertions = 0;

  while (*head != NULL &&
         (*head)->is_WF_GW &&
         (*head)->psi_cm + psi_span_tolerance_cm < column_depth_cm &&
         insertions < NO_SURFACE_TO_PSI_SPAN_MAX_INSERTIONS) {
    if ((*head)->psi_cm >= PSI_UPPER_LIM) {
      break;
    }

    const double psi_target_cm = fmin(column_depth_cm, (*head)->psi_cm + psi_increment_cm);
    if (psi_target_cm <= (*head)->psi_cm + psi_span_tolerance_cm) {
      break;
    }

    if (!lgarto_insert_zero_depth_groundwater_front(psi_target_cm, soil_type, frozen_factor,
                                                    soil_properties, head)) {
      break;
    }

    insertions++;
  }

  if (insertions > 0 && verbosity.compare("high") == 0) {
    printf("Inserted %d zero-depth TO support front(s) to extend the no-surface TO psi span toward %.6f cm.\n",
           insertions, column_depth_cm);
    listPrint(*head);
  }

  return insertions;
}

static double lgarto_restore_to_psi_gap_mass_via_depth(struct wetting_front *upper_front,
                                                       struct wetting_front *inserted,
                                                       struct wetting_front *lower_front,
                                                       double target_mass,
                                                       int num_layers,
                                                       double *cum_layer_thickness_cm,
                                                       struct wetting_front *head)
{
  if (upper_front == NULL || inserted == NULL || lower_front == NULL || head == NULL) {
    return target_mass;
  }

  if (upper_front->layer_num != inserted->layer_num ||
      inserted->layer_num != lower_front->layer_num) {
    return target_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, head);
  }

  const double original_upper_depth_cm = upper_front->depth_cm;
  const double original_lower_depth_cm = lower_front->depth_cm;
  if (original_lower_depth_cm <= original_upper_depth_cm + 2.0 * CREATION_COLOCATED_TOLERANCE_CM) {
    return target_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, head);
  }

  const int layer_num = upper_front->layer_num;
  double minimum_upper_depth_cm = cum_layer_thickness_cm[layer_num - 1] + CREATION_COLOCATED_TOLERANCE_CM;
  struct wetting_front *previous = lgar_find_previous_front(head, upper_front);
  if (previous != NULL) {
    minimum_upper_depth_cm =
      fmax(minimum_upper_depth_cm, previous->depth_cm + CREATION_COLOCATED_TOLERANCE_CM);
  }

  if (minimum_upper_depth_cm >= original_upper_depth_cm - CREATION_COLOCATED_TOLERANCE_CM) {
    return target_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, head);
  }

  auto apply_geometry = [&](double upper_depth_cm) {
    upper_front->depth_cm = upper_depth_cm;
    inserted->depth_cm = 0.5 * (upper_depth_cm + original_lower_depth_cm);
  };

  apply_geometry(original_upper_depth_cm);
  const double mass_at_original_geometry = lgar_calc_mass_bal(cum_layer_thickness_cm, head);

  apply_geometry(minimum_upper_depth_cm);
  const double mass_at_shallow_limit = lgar_calc_mass_bal(cum_layer_thickness_cm, head);

  if (mass_at_shallow_limit + TO_PSI_GAP_REFINEMENT_MASS_TOLERANCE_CM < target_mass) {
    apply_geometry(original_upper_depth_cm);
    return target_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, head);
  }

  double bracket_shallow_depth_cm = minimum_upper_depth_cm;
  double bracket_deep_depth_cm = original_upper_depth_cm;
  double best_upper_depth_cm = original_upper_depth_cm;
  double best_mass_error_cm = fabs(target_mass - mass_at_original_geometry);

  for (int iter = 0; iter < 80; iter++) {
    const double probe_upper_depth_cm =
      0.5 * (bracket_shallow_depth_cm + bracket_deep_depth_cm);
    apply_geometry(probe_upper_depth_cm);
    const double current_mass_cm = lgar_calc_mass_bal(cum_layer_thickness_cm, head);
    const double current_mass_error_cm = target_mass - current_mass_cm;

    if (fabs(current_mass_error_cm) < best_mass_error_cm) {
      best_mass_error_cm = fabs(current_mass_error_cm);
      best_upper_depth_cm = probe_upper_depth_cm;
    }

    if (fabs(current_mass_error_cm) <= TO_PSI_GAP_REFINEMENT_MASS_TOLERANCE_CM) {
      best_upper_depth_cm = probe_upper_depth_cm;
      break;
    }

    if (current_mass_error_cm > 0.0) {
      bracket_deep_depth_cm = probe_upper_depth_cm;
    }
    else {
      bracket_shallow_depth_cm = probe_upper_depth_cm;
    }
  }

  apply_geometry(best_upper_depth_cm);
  return target_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, head);
}

static bool lgarto_find_large_to_psi_gap_candidate(int num_layers,
                                                   double *cum_layer_thickness_cm,
                                                   struct wetting_front *head,
                                                   double minimum_candidate_depth_cm,
                                                   TOPsiGapRefinementCandidate *candidate)
{
  if (head == NULL || candidate == NULL) {
    return false;
  }

  const double dry_zero_depth_psi_cm =
    TO_PSI_GAP_REFINEMENT_ZERO_DEPTH_DRY_FACTOR * cum_layer_thickness_cm[num_layers];
  bool found = false;
  double best_excess_gap_cm = 0.0;

  for (struct wetting_front *current = head; current != NULL && current->next != NULL;
       current = current->next) {
    struct wetting_front *next = current->next;

    if (current->is_WF_GW == FALSE || next->is_WF_GW == FALSE) {
      continue;
    }

    if (current->to_bottom == TRUE) {
      continue;
    }

    // Near-zero TO/GW supports are metadata for the dry top of the TO chain,
    // not spatial resolution. Do not subdivide a dry near-surface support pair.
    if (current->depth_cm <= TO_PSI_GAP_REFINEMENT_ZERO_DEPTH_MAX_CM &&
        next->depth_cm <= TO_PSI_GAP_REFINEMENT_ZERO_DEPTH_MAX_CM &&
        current->psi_cm > dry_zero_depth_psi_cm &&
        next->psi_cm > dry_zero_depth_psi_cm) {
      continue;
    }

    if (current->depth_cm <= CREATION_COLOCATED_TOLERANCE_CM) {
      continue;
    }

    if (current->depth_cm <= minimum_candidate_depth_cm) {
      continue;
    }

    if (current->layer_num != next->layer_num) {
      continue;
    }

    // Leave the boundary-support interval alone. Refining the gap immediately above a
    // to_bottom scaffold front can create deep, high-psi TO structure that changes
    // recharge behavior much more than it helps missing-gradient coverage.
    if (next->to_bottom == TRUE) {
      continue;
    }

    if (next->depth_cm <= current->depth_cm + CREATION_COLOCATED_TOLERANCE_CM) {
      continue;
    }

    const double psi_gap_cm = current->psi_cm - next->psi_cm;
    const double local_psi_midpoint_cm =
      0.5 * (fmax(0.0, current->psi_cm) + fmax(0.0, next->psi_cm));
    const double threshold_cm =
      fmax(TO_PSI_GAP_REFINEMENT_MIN_GAP_CM,
           TO_PSI_GAP_REFINEMENT_RELATIVE_GAP_FACTOR * local_psi_midpoint_cm);
    const double excess_gap_cm = psi_gap_cm - threshold_cm;
    if (excess_gap_cm <= best_excess_gap_cm) {
      continue;
    }

    found = true;
    best_excess_gap_cm = excess_gap_cm;
    candidate->upper_front = current;
    candidate->lower_front = next;
    candidate->psi_gap_cm = psi_gap_cm;
  }

  return found;
}

static bool lgarto_try_refine_large_to_psi_gap(int num_layers,
                                               double *cum_layer_thickness_cm,
                                               int *soil_type,
                                               double *frozen_factor,
                                               struct soil_properties_ *soil_properties,
                                               struct wetting_front **head,
                                               double minimum_candidate_depth_cm)
{
  TOPsiGapRefinementCandidate candidate = {};
  if (!lgarto_find_large_to_psi_gap_candidate(num_layers, cum_layer_thickness_cm, *head,
                                              minimum_candidate_depth_cm, &candidate)) {
    return false;
  }

  const double prior_mass = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
  struct wetting_front *snapshot = listCopy(*head, NULL);

  const double insert_depth_cm =
    0.5 * (candidate.upper_front->depth_cm + candidate.lower_front->depth_cm);
  const double inserted_psi_cm =
    0.5 * (candidate.upper_front->psi_cm + candidate.lower_front->psi_cm);
  const double dry_zero_depth_psi_cm =
    TO_PSI_GAP_REFINEMENT_ZERO_DEPTH_DRY_FACTOR * cum_layer_thickness_cm[num_layers];

  if (insert_depth_cm <= TO_PSI_GAP_REFINEMENT_ZERO_DEPTH_MAX_CM &&
      candidate.upper_front->psi_cm > dry_zero_depth_psi_cm &&
      candidate.lower_front->psi_cm > dry_zero_depth_psi_cm) {
    if (snapshot != NULL) {
      listDelete(snapshot);
    }
    return false;
  }

  struct wetting_front *inserted =
    lgar_insert_groundwater_front_before(insert_depth_cm, inserted_psi_cm, candidate.lower_front,
                                         num_layers, soil_type, cum_layer_thickness_cm,
                                         frozen_factor, soil_properties, head);
  if (inserted == NULL) {
    if (snapshot != NULL) {
      listDelete(snapshot);
    }
    return false;
  }

  const double residual_mass_error_cm =
    fabs(lgarto_restore_to_psi_gap_mass_via_depth(candidate.upper_front, inserted, candidate.lower_front,
                                                  prior_mass, num_layers, cum_layer_thickness_cm, *head));
  if (residual_mass_error_cm > TO_PSI_GAP_REFINEMENT_MASS_TOLERANCE_CM) {
    if (*head != NULL) {
      listDelete(*head);
    }
    *head = snapshot;
    if (verbosity.compare("high") == 0) {
      printf("Rolling back TO psi-gap refinement due to residual mass error %.17lf cm.\n",
             residual_mass_error_cm);
    }
    return false;
  }

  if (snapshot != NULL) {
    listDelete(snapshot);
  }

  if (verbosity.compare("high") == 0) {
    printf("Inserted corrective TO front at depth %.6f cm with psi %.6f cm to refine a psi gap of %.6f cm.\n",
           inserted->depth_cm, inserted->psi_cm, candidate.psi_gap_cm);
    listPrint(*head);
  }

  return true;
}

static int lgarto_refine_large_to_psi_gaps(int num_layers,
                                           double *cum_layer_thickness_cm,
                                           int *soil_type,
                                           double *frozen_factor,
                                           struct soil_properties_ *soil_properties,
                                           struct wetting_front **head,
                                           int max_insertions,
                                           double minimum_candidate_depth_cm)
{
  int insertions = 0;
  while (insertions < max_insertions &&
         lgarto_try_refine_large_to_psi_gap(num_layers, cum_layer_thickness_cm, soil_type,
                                            frozen_factor, soil_properties, head,
                                            minimum_candidate_depth_cm)) {
    insertions++;
  }
  return insertions;
}

// #############################################################################################################################
/*
  Read and initialize values from a configuration file
  @param verbosity              : supress all outputs, file writing if 'none', screen output and write file if 'high'
  @param layer_thickness_cm     : 1D (double) array of layer thicknesses in cm, read from config file
  @param layer_soil_type        : 1D (int) array of layers soil type, read from config file, each integer represent a soil type
  @param num_layers             : number of actual soil layers
  @param num_wetting_fronts     : number of wetting fronts
  @param num_cells_temp         : number of cells of the discretized soil temperature profile
  @param cum_layer_thickness_cm : 1D (double) array of cumulative thickness of layers, allocate memory at run time
  @param soil_depth_cm          : depth of the computational domain (i.e., depth of the last/deepest soil layer from the surface)
  @param initial_psi_cm         : model initial (psi) condition
  @param timestep_h             : model timestep in hours
  @param forcing_resolution_h   : forcing resolution in hours
  @param forcing_interval       : factor equals to forcing_resolution_h/timestep_h (used to determine model subtimestep's forcings)
  @param num_soil_types         : number of soil types; must be less than or equal to MAX_NUM_SOIL_TYPES
  @param AET_cm                 : actual evapotranspiration in cm
  @param soil_temperature       : 1D (double) array of soil temperature [K]; bmi input for coupling lasam to soil freeze thaw model
  @param soil_temperature_z     : 1D (double) array of soil discretization associated with temperature profile [m];
                                  depth from the surface in meters
  @param frozen_factor          : frozen factor causing the hydraulic conductivity to decrease due to frozen soil
                                  (when coupled to soil freeze thaw model)
  @param wilting_point_psi_cm   : wilting point (the amount of water not available for plants or not accessible by plants)
  @param field_capacity_psi_cm  : field capacity, represented with a capillary head (head above which drainage is much faster)
  @param ponded_depth_cm        : amount of water on the surface not available for surface drainage (initialized to zero)
  @param ponded_depth_max cm    : maximum amount of water on the surface not available for surface drainage (default is zero)
  @param nint                   : number of trapezoids used in integrating the Geff function (set to 120)
  @param time_s                 : current time [s] (initially set to zero)
  @param sft_coupled            : model coupling flag. if true, lasam is coupled to soil freeze thaw model; default is uncoupled version
  @param giuh_ordinates         : geomorphological instantaneous unit hydrograph
  @param num_giuh_ordinates     : number of giuh ordinates
*/

// #############################################################################################################################
extern void InitFromConfigFile(string config_file, struct model_state *state)
{

  ifstream fp; //FILE *fp = fopen(config_file.c_str(),"r");
  fp.open(config_file);
  //struct wetting_front* head = state->head;
  
  // loop over the variables in the file to see if verbosity is provided, if not default is "none" (prints nothing)
  while (fp) {
    string line;
    string param_key, param_value, param_unit;

    getline(fp, line);

    int loc_eq = line.find("=") + 1;
    int loc_u = line.find("[");
    param_key = line.substr(0,line.find("="));

    param_value = line.substr(loc_eq,loc_u - loc_eq);

    if (param_key == "verbosity") {
      verbosity = param_value;
      if (verbosity.compare("none") != 0) {
	std::cerr<<"Verbosity is set to \' "<<verbosity<<"\' \n";
	std::cerr<<"          *****         \n";
      }

      fp.clear();
      break;
    }
  }

  // seek to beginning of input after searching for 'verbosity'
  fp.clear();
  fp.seekg(0, fp.beg);

  if (verbosity.compare("none") != 0) {
    std::cerr<<"------------- Initialization from config file ---------------------- \n";
  }

  // setting these options to false (defualt) 
  state->lgar_bmi_params.sft_coupled           = false;
  state->lgar_bmi_params.use_closed_form_G     = false;
  state->lgar_bmi_params.adaptive_timestep     = false;
  state->lgar_bmi_params.runoff_in_prev_step   = false;
  state->lgar_bmi_params.PET_affects_precip    = false;
  state->lgar_bmi_params.allow_flux_caching    = false;
  state->lgar_bmi_params.log_mode              = false;
  state->lgar_bmi_params.TO_enabled            = false;
  state->lgar_bmi_params.free_drainage_enabled = false;
  state->lgar_bmi_params.lower_bdy_flux_to_CR  = false;
  state->lgar_bmi_params.mobile_groundwater_level = false;
  state->lgar_bmi_params.groundwater_depth_cm = 0.0;
  state->lgar_bmi_params.CR_fast_discharge_threshold_cm = 0.0;
  state->lgar_bmi_params.CR_slow_discharge_threshold_cm = 0.0;
  state->lgar_bmi_params.initial_CR_fast_storage_cm = 0.0;
  state->lgar_bmi_params.initial_CR_slow_storage_cm = 0.0;
  state->lgar_bmi_params.CR_capillary_supply_threshold_cm = 0.1;
  // setting mass balance tolerance to be large by default; this can be specified in the config file
  state->lgar_bmi_params.mbal_tol = 1.E1;
  
  bool is_layer_thickness_set       = false;
  bool is_initial_psi_set           = false;
  bool is_timestep_set              = false;
  bool is_endtime_set               = false;
  bool is_forcing_resolution_set    = false;
  bool is_layer_soil_type_set       = false;
  bool is_wilting_point_psi_cm_set  = false;
  bool is_field_capacity_psi_cm_set = false;
  bool is_root_zone_depth_cm_set    = false;
  bool is_a_set                     = false;
  bool is_b_set                     = false;
  bool is_frac_to_CR_set            = false;
  bool is_a_slow_set                = false;
  bool is_b_slow_set                = false;
  bool is_frac_slow_set             = false;
  bool is_soil_params_file_set      = false;
  bool is_max_valid_soil_types_set  = false;
  bool is_giuh_ordinates_set        = false;
  bool is_soil_z_set                = false;
  bool is_ponded_depth_max_cm_set   = false;
  bool is_state_path_set            = false;
  bool is_non_vadose_state_path_set = false;
  bool is_giuh_state_path_set       = false;

  string soil_params_file;

  // a temporary array to store the original (hourly based) giuh values
  std::vector<double> giuh_ordinates_temp;
 
  while (fp) {

    string line;
    string param_key, param_value, param_unit;

    getline(fp, line);

    int loc_eq = line.find("=") + 1;
    int loc_u = line.find("[");
    param_key = line.substr(0,line.find("="));

    bool is_unit = line.find("[") != string::npos;

    if (is_unit)
      param_unit = line.substr(loc_u,line.find("]")+1);
    else
      param_unit = "";

    param_value = line.substr(loc_eq,loc_u - loc_eq);
    
    if (param_key == "layer_thickness") {
      vector<double> vec = ReadVectorData(param_value);

      state->lgar_bmi_params.layer_thickness_cm = new double[vec.size()+1];
      state->lgar_bmi_params.cum_layer_thickness_cm = new double[vec.size()+1];

      state->lgar_bmi_params.layer_thickness_cm[0] = 0.0; // the value at index 0 is never used
      // calculate the cumulative (absolute) depth from land surface to bottom of each soil layer
      state->lgar_bmi_params.cum_layer_thickness_cm[0] = 0.0;

      for (unsigned int layer=1; layer <= vec.size(); layer++) {
      	state->lgar_bmi_params.layer_thickness_cm[layer] = vec[layer-1];
	state->lgar_bmi_params.cum_layer_thickness_cm[layer] = state->lgar_bmi_params.cum_layer_thickness_cm[layer-1] + vec[layer-1];
      }

	      state->lgar_bmi_params.num_layers = vec.size();

	      state->lgar_bmi_params.soil_depth_cm = state->lgar_bmi_params.cum_layer_thickness_cm[state->lgar_bmi_params.num_layers];
	      state->lgar_bmi_params.groundwater_depth_cm = state->lgar_bmi_params.soil_depth_cm;
	      is_layer_thickness_set = true;

      if (verbosity.compare("high") == 0) {
	std::cerr<<"Number of layers : "<<state->lgar_bmi_params.num_layers<<"\n";
	for (int i=1; i<=state->lgar_bmi_params.num_layers; i++)
	  std::cerr<<"Thickness, cum. depth : "<<state->lgar_bmi_params.layer_thickness_cm[i]<<" , "
		   <<state->lgar_bmi_params.cum_layer_thickness_cm[i]<<"\n";
	std::cerr<<"          *****         \n";
      }

      continue;
    }
    else if (param_key == "layer_soil_type") {
      vector<double> vec = ReadVectorData(param_value);

      state->lgar_bmi_params.layer_soil_type = new int[vec.size()+1];

      for (unsigned int layer=1; layer <= vec.size(); layer++)
      	state->lgar_bmi_params.layer_soil_type[layer] = vec[layer-1];

      is_layer_soil_type_set = true;

      continue;
    }
    else if (param_key == "giuh_ordinates") {
      vector<double> vec = ReadVectorData(param_value);

      giuh_ordinates_temp.resize(vec.size()+1);
 
      for (unsigned int i=1; i <= vec.size(); i++)
	giuh_ordinates_temp[i] = vec[i-1];

      is_giuh_ordinates_set = true;

      if (verbosity.compare("high") == 0) {
	for (unsigned int i=1; i <= vec.size(); i++)
	  std::cerr<<"GIUH ordinates (hourly) : "<<giuh_ordinates_temp[i]<<"\n";

	std::cerr<<"          *****         \n";
      }

      continue;
    }
    else if (param_key == "soil_z") {
      vector<double> vec = ReadVectorData(param_value);

      state->lgar_bmi_params.soil_temperature_z = new double[vec.size()];

      for (unsigned int i=0; i < vec.size(); i++)
      	state->lgar_bmi_params.soil_temperature_z[i] = vec[i];

      state->lgar_bmi_params.num_cells_temp = vec.size();

      is_soil_z_set = true;

      if (verbosity.compare("high") == 0) {
	for (int i=0; i<state->lgar_bmi_params.num_cells_temp; i++)
	  std::cerr<<"Soil z (temperature resolution) : "<<state->lgar_bmi_params.soil_temperature_z[i]<<"\n";

	std::cerr<<"          *****         \n";
      }

      continue;
    }
    else if (param_key == "initial_psi") {
      state->lgar_bmi_params.initial_psi_cm = stod(param_value);
      is_initial_psi_set = true;

      if (verbosity.compare("high") == 0) {
	std::cerr<<"Initial Psi : "<<state->lgar_bmi_params.initial_psi_cm<<"\n";
	std::cerr<<"          *****         \n";
      }

      continue;
    }
    else if (param_key == "max_valid_soil_types") {
      state->lgar_bmi_params.num_soil_types = std::min(stoi(param_value), MAX_NUM_SOIL_TYPES);
      is_max_valid_soil_types_set = true;
      continue;
    }
    else if (param_key == "soil_params_file") {
      soil_params_file = param_value;
      is_soil_params_file_set = true;

      if (verbosity.compare("high") == 0) {
	std::cerr<<"Soil paramaters file : "<<soil_params_file<<"\n";
	std::cerr<<"          *****         \n";
      }

      continue;
    }
    else if (param_key == "wilting_point_psi") {
      state->lgar_bmi_params.wilting_point_psi_cm = stod(param_value);
      is_wilting_point_psi_cm_set = true;

      if (verbosity.compare("high") == 0) {
	std::cerr<<"Wilting point Psi [cm] : "<<state->lgar_bmi_params.wilting_point_psi_cm<<"\n";
	std::cerr<<"          *****         \n";
      }

      continue;
    }
    else if (param_key == "field_capacity_psi") {
      state->lgar_bmi_params.field_capacity_psi_cm = stod(param_value);
      is_field_capacity_psi_cm_set = true;

      if (verbosity.compare("high") == 0) {
	std::cerr<<"Field capacity Psi [cm] : "<<state->lgar_bmi_params.field_capacity_psi_cm<<"\n";
	std::cerr<<"          *****         \n";
      }

      continue;
    }
    else if (param_key == "root_zone_depth") {
      state->lgar_bmi_params.root_zone_depth_cm = stod(param_value);

      if (state->lgar_bmi_params.root_zone_depth_cm < 0.0) {
        printf("root_zone_depth is less than 0 \n");
        abort();
      }

      is_root_zone_depth_cm_set = true;

      if (verbosity.compare("high") == 0) {
        std::cerr<<"root zone depth [cm] : "<<state->lgar_bmi_params.root_zone_depth_cm<<"\n";
        std::cerr<<"          *****         \n";
      }

      continue;
    }
    else if (param_key == "a") {
      state->lgar_bmi_params.a = stod(param_value);
      is_a_set = true;

      if (verbosity.compare("high") == 0) {
	std::cerr<<"a : "<<state->lgar_bmi_params.a<<"\n";
	std::cerr<<"          *****         \n";
      }
      continue;
    }
    else if (param_key == "b") {
      state->lgar_bmi_params.b = stod(param_value);
      is_b_set = true;

      if (verbosity.compare("high") == 0) {
	std::cerr<<"b : "<<state->lgar_bmi_params.b<<"\n";
	std::cerr<<"          *****         \n";
      }
      continue;
    }
    else if (param_key == "frac_to_CR") {
      state->lgar_bmi_params.frac_to_CR = stod(param_value);
      is_frac_to_CR_set = true;

      if (verbosity.compare("high") == 0) {
	std::cerr<<"frac_to_CR : "<<state->lgar_bmi_params.frac_to_CR<<"\n";
	std::cerr<<"          *****         \n";
      }
      continue;
    }
    else if (param_key == "a_slow") {
      state->lgar_bmi_params.a_slow = stod(param_value);
      is_a_slow_set = true;

      if (verbosity.compare("high") == 0) {
	std::cerr<<"a_slow : "<<state->lgar_bmi_params.a_slow<<"\n";
	std::cerr<<"          *****         \n";
      }
      continue;
    }
    else if (param_key == "b_slow") {
      state->lgar_bmi_params.b_slow = stod(param_value);
      is_b_slow_set = true;

      if (verbosity.compare("high") == 0) {
	std::cerr<<"b_slow : "<<state->lgar_bmi_params.b_slow<<"\n";
	std::cerr<<"          *****         \n";
      }
      continue;
    }
    else if (param_key == "frac_slow") {
      state->lgar_bmi_params.frac_slow = stod(param_value);
      is_frac_slow_set = true;

      if (verbosity.compare("high") == 0) {
	std::cerr<<"frac_slow : "<<state->lgar_bmi_params.frac_slow<<"\n";
	std::cerr<<"          *****         \n";
      }
      continue;
    }
    else if (param_key == "CR_fast_discharge_threshold" ||
             param_key == "CR_fast_discharge_threshold_cm") {
      state->lgar_bmi_params.CR_fast_discharge_threshold_cm = stod(param_value);

      if (state->lgar_bmi_params.CR_fast_discharge_threshold_cm < 0.0) {
        std::cerr<<"Invalid option: CR_fast_discharge_threshold must be >= 0. \n";
        abort();
      }

      if (verbosity.compare("high") == 0) {
	std::cerr<<"CR_fast_discharge_threshold : "
                 <<state->lgar_bmi_params.CR_fast_discharge_threshold_cm<<"\n";
	std::cerr<<"          *****         \n";
      }
      continue;
    }
    else if (param_key == "CR_slow_discharge_threshold" ||
             param_key == "CR_slow_discharge_threshold_cm") {
      state->lgar_bmi_params.CR_slow_discharge_threshold_cm = stod(param_value);

      if (state->lgar_bmi_params.CR_slow_discharge_threshold_cm < 0.0) {
        std::cerr<<"Invalid option: CR_slow_discharge_threshold must be >= 0. \n";
        abort();
      }

      if (verbosity.compare("high") == 0) {
	std::cerr<<"CR_slow_discharge_threshold : "
                 <<state->lgar_bmi_params.CR_slow_discharge_threshold_cm<<"\n";
	std::cerr<<"          *****         \n";
      }
      continue;
    }
    else if (param_key == "initial_CR_storage" ||
             param_key == "initial_CR_storage_cm" ||
             param_key == "initial_CR_fast_storage" ||
             param_key == "initial_CR_fast_storage_cm") {
      state->lgar_bmi_params.initial_CR_fast_storage_cm = stod(param_value);

      if (state->lgar_bmi_params.initial_CR_fast_storage_cm < 0.0) {
        std::cerr<<"Invalid option: initial_CR_fast_storage must be >= 0. \n";
        abort();
      }

      if (verbosity.compare("high") == 0) {
	std::cerr<<"initial_CR_fast_storage : "
                 <<state->lgar_bmi_params.initial_CR_fast_storage_cm<<"\n";
	std::cerr<<"          *****         \n";
      }
      continue;
    }
    else if (param_key == "initial_CR_slow_storage" ||
             param_key == "initial_CR_slow_storage_cm") {
      state->lgar_bmi_params.initial_CR_slow_storage_cm = stod(param_value);

      if (state->lgar_bmi_params.initial_CR_slow_storage_cm < 0.0) {
        std::cerr<<"Invalid option: initial_CR_slow_storage must be >= 0. \n";
        abort();
      }

      if (verbosity.compare("high") == 0) {
	std::cerr<<"initial_CR_slow_storage : "
                 <<state->lgar_bmi_params.initial_CR_slow_storage_cm<<"\n";
	std::cerr<<"          *****         \n";
      }
      continue;
    }
    else if (param_key == "CR_capillary_supply_threshold" ||
             param_key == "CR_capillary_supply_threshold_cm") {
      state->lgar_bmi_params.CR_capillary_supply_threshold_cm = stod(param_value);

      if (state->lgar_bmi_params.CR_capillary_supply_threshold_cm < 0.0) {
        std::cerr<<"Invalid option: CR_capillary_supply_threshold must be >= 0. \n";
        abort();
      }

      if (verbosity.compare("high") == 0) {
	std::cerr<<"CR_capillary_supply_threshold : "
                 <<state->lgar_bmi_params.CR_capillary_supply_threshold_cm<<"\n";
	std::cerr<<"          *****         \n";
      }
      continue;
    }
    else if (param_key == "frac_to_GW") {
      state->lgar_bmi_params.frac_to_CR = stod(param_value);
      is_frac_to_CR_set = true;

      if (verbosity.compare("high") == 0) {
	std::cerr<<"frac_to_CR (using old name in config): "<<state->lgar_bmi_params.frac_to_CR<<"\n";
	std::cerr<<"          *****         \n";
      }
      continue;
    }
    else if (param_key == "spf_factor") {
      state->lgar_bmi_params.spf_factor = stod(param_value);

      if (verbosity.compare("high") == 0) {
	std::cerr<<"spf_factor : "<<state->lgar_bmi_params.spf_factor<<"\n";
	std::cerr<<"          *****         \n";
      }
      continue;
    }
    else if (param_key == "use_closed_form_G") { 
      if (param_value == "false") {
        state->lgar_bmi_params.use_closed_form_G = false;
      }
      else if (param_value == "true") {
        state->lgar_bmi_params.use_closed_form_G = true;
      }
      else {
	std::cerr<<"Invalid option: use_closed_form_G must be true or false. \n";
        abort();
      }

      continue;
    }
  else if (param_key == "init_state_path") {
      // remove potential whitespace
      param_value.erase(0, param_value.find_first_not_of(" \t"));
      param_value.erase(param_value.find_last_not_of(" \t") + 1);

      state->lgar_bmi_params.init_state_path = param_value;
      is_state_path_set = true;

      if (verbosity.compare("high") == 0) {
          std::cerr << "init_state_path set to: "
                    << state->lgar_bmi_params.init_state_path << "\n";
          std::cerr << "          *****         \n";
      }
      continue;
  }
  else if (param_key == "init_non_vadose_state_path") {
      // remove potential whitespace
      param_value.erase(0, param_value.find_first_not_of(" \t"));
      param_value.erase(param_value.find_last_not_of(" \t") + 1);

      state->lgar_bmi_params.init_non_vadose_state_path = param_value;
      is_non_vadose_state_path_set = true;

      if (verbosity.compare("high") == 0) {
          std::cerr << "init_non_vaodse_state_path set to: "
                    << state->lgar_bmi_params.init_non_vadose_state_path << "\n";
          std::cerr << "          *****         \n";
      }
      continue;
  }
  else if (param_key == "init_giuh_state_path") {
      // remove potential whitespace
      param_value.erase(0, param_value.find_first_not_of(" \t"));
      param_value.erase(param_value.find_last_not_of(" \t") + 1);

      state->lgar_bmi_params.init_giuh_state_path = param_value;
      is_giuh_state_path_set = true;

      if (verbosity.compare("high") == 0) {
          std::cerr << "init_giuh_state_path set to: "
                    << state->lgar_bmi_params.init_giuh_state_path << "\n";
          std::cerr << "          *****         \n";
      }
      continue;
  }
    else if (param_key == "free_drainage_enabled") { 
      if (param_value == "false") {
        state->lgar_bmi_params.free_drainage_enabled = false;
      }
      else if (param_value == "true") {
        state->lgar_bmi_params.free_drainage_enabled = true;
      }
      else {
	std::cerr<<"Invalid option: free_drainage_enabled must be true or false, or left unspecified (defaulting to false). \n";
        abort();
      }

      continue;
    }
    else if (param_key == "lower_bdy_flux_to_CR" || param_key == "free_drainage_to_CR") {
      if (param_value == "false") {
        state->lgar_bmi_params.lower_bdy_flux_to_CR = false;
      }
      else if (param_value == "true") {
        state->lgar_bmi_params.lower_bdy_flux_to_CR = true;
      }
      else {
	std::cerr<<"Invalid option: lower_bdy_flux_to_CR must be true or false, or left unspecified (defaulting to false). \n";
        abort();
      }

      continue;
    }
    else if (param_key == "mobile_groundwater_level") {
      if (param_value == "false") {
        state->lgar_bmi_params.mobile_groundwater_level = false;
      }
      else if (param_value == "true") {
        state->lgar_bmi_params.mobile_groundwater_level = true;
      }
      else {
		std::cerr<<"Invalid option: mobile_groundwater_level must be true or false, or left unspecified (defaulting to false). \n";
        abort();
      }

      continue;
    }
    else if (param_key == "PET_affects_precip") { 
      if (param_value == "false") {
        state->lgar_bmi_params.PET_affects_precip = false;
      }
      else if (param_value == "true") {
        state->lgar_bmi_params.PET_affects_precip = true;
      }
      else {
	std::cerr<<"Invalid option: PET_affects_precip must be true or false, or left unspecified (defaulting to false). \n";
        abort();
      }

      continue;
    }
    else if (param_key == "allow_flux_caching") { 
      if (param_value == "false") {
        state->lgar_bmi_params.allow_flux_caching = false;
      }
      else if (param_value == "true") {
        state->lgar_bmi_params.allow_flux_caching = true;
      }
      else {
	std::cerr<<"Invalid option: allow_flux_caching must be true or false, or left unspecified (defaulting to false). \n";
        abort();
      }

      continue;
    }
    else if (param_key == "log_mode") { 
      if ((param_value == "false") || (param_value == "0")) {
        state->lgar_bmi_params.log_mode = false;
      }
      else if ( (param_value == "true") || (param_value == "1")) {
        state->lgar_bmi_params.log_mode = true;
        if (verbosity.compare("high") == 0) {
          printf("log_mode enabled. So K_s for each layer, alpha for each layer, and a for the nonlinear reservoir(s) will use the log of their input values. \n");
        }
      }
      else {
	std::cerr<<"Invalid option: log_mode must be true or false. \n";
        abort();
      }

      continue;
    }
    else if (param_key == "mbal_tol") {
      state->lgar_bmi_params.mbal_tol = stod(param_value);

      if (verbosity.compare("high") == 0) {
	std::cerr<<"Mass balance tolerance [cm] : "<<state->lgar_bmi_params.mbal_tol<<"\n";
	std::cerr<<"          *****         \n";
      }

      continue;
    }
    else if (param_key == "adaptive_timestep") { 
      if ((param_value == "false") || (param_value == "0")) {
        state->lgar_bmi_params.adaptive_timestep = false;
      }
      else if ( (param_value == "true") || (param_value == "1")) {
        state->lgar_bmi_params.adaptive_timestep = true;
      }
      else {
	std::cerr<<"Invalid option: adaptive_timestep must be true or false. \n";
        abort();
      }

      continue;
    }
    else if (param_key == "TO_enabled") {
      if ((param_value == "false") || (param_value == "0")) {
        state->lgar_bmi_params.TO_enabled = false;
      }
      else if ((param_value == "true") || (param_value == "1")) {
        state->lgar_bmi_params.TO_enabled = true;
      }
      else {
	std::cerr<<"Invalid option: TO_enabled must be true or false. \n";
        abort();
      }

      continue;
    }
    else if (param_key == "timestep") {
      state->lgar_bmi_params.timestep_h = stod(param_value);

      if (param_unit == "[s]" || param_unit == "[sec]" || param_unit == "") // defalut time unit is seconds
	state->lgar_bmi_params.timestep_h /= 3600; // convert to hours
      else if (param_unit == "[min]" || param_unit == "[minute]")
	state->lgar_bmi_params.timestep_h /= 60; // convert to hours
      else if (param_unit == "[h]" || param_unit == "[hr]")
	state->lgar_bmi_params.timestep_h /= 1.0; // convert to hours

      assert (state->lgar_bmi_params.timestep_h > 0);
      is_timestep_set = true;

      state->lgar_bmi_params.minimum_timestep_h = state->lgar_bmi_params.timestep_h;

      if (verbosity.compare("high") == 0) {
	std::cerr<<"Model timestep [hours,seconds]: "<<state->lgar_bmi_params.timestep_h<<" , "
		 <<state->lgar_bmi_params.timestep_h*3600<<"\n";
	std::cerr<<"          *****         \n";
      }

      continue;
    }
    else if (param_key == "endtime") {

      if (param_unit == "[s]" || param_unit == "[sec]" || param_unit == "") // defalut time unit is seconds
	state->lgar_bmi_params.endtime_s = stod(param_value);
      else if (param_unit == "[min]" || param_unit == "[minute]")
	state->lgar_bmi_params.endtime_s = stod(param_value) * 60.0;
      else if (param_unit == "[h]" || param_unit == "[hr]")
	state->lgar_bmi_params.endtime_s = stod(param_value) * 3600.0;
      else if (param_unit == "[d]" || param_unit == "[day]")
	state->lgar_bmi_params.endtime_s = stod(param_value) * 86400.0;

      assert (state->lgar_bmi_params.endtime_s > 0);
      is_endtime_set = true;

      if (verbosity.compare("high") == 0) {
	std::cerr<<"Endtime [days, hours]: "<< state->lgar_bmi_params.endtime_s/86400.0 <<" , "
		 << state->lgar_bmi_params.endtime_s/3600.0<<"\n";
	std::cerr<<"          *****         \n";
      }

      continue;
    }
    else if (param_key == "forcing_resolution") {
      state->lgar_bmi_params.forcing_resolution_h = stod(param_value);

      if (param_unit == "[s]" || param_unit == "[sec]" || param_unit == "") // defalut time unit is seconds
	state->lgar_bmi_params.forcing_resolution_h /= 3600;                // convert to hours
      else if (param_unit == "[min]" || param_unit == "[minute]")
	state->lgar_bmi_params.forcing_resolution_h /= 60;                 // convert to hours
      else if (param_unit == "[h]" || param_unit == "[hr]")
	state->lgar_bmi_params.forcing_resolution_h /= 1.0;               // convert to hours

      assert (state->lgar_bmi_params.forcing_resolution_h > 0);
      is_forcing_resolution_set = true;

      if (verbosity.compare("high") == 0) {
	std::cerr<<"Forcing resolution [hours]: "<<state->lgar_bmi_params.forcing_resolution_h<<"\n";
	std::cerr<<"          *****         \n";
      }

      continue;
    }
    else if (param_key == "sft_coupled") {
      if (param_value == "true") {
	state->lgar_bmi_params.sft_coupled = 1;
      }
      else if (param_value == "false") {
	state->lgar_bmi_params.sft_coupled = 0; // false
      }
      else {
	std::cerr<<"Invalid option: sft_coupled must be true or false. \n";
        abort();
      }
      
      continue;
    }
    else if (param_key == "ponded_depth_max") {
      state->lgar_bmi_params.ponded_depth_max_cm = fmax(stod(param_value), 0.0);
      is_ponded_depth_max_cm_set = true;

      if (verbosity.compare("high") == 0) {
	std::cerr<<"Maximum ponded depth [cm] : "<<state->lgar_bmi_params.ponded_depth_max_cm<<"\n";
	std::cerr<<"          *****         \n";
      }

      continue;
    }
    else if (param_key == "calib_params") {
      if (param_value == "true") {
	state->lgar_bmi_params.calib_params_flag = 1;
      }
      else if (param_value == "false") {
	state->lgar_bmi_params.calib_params_flag = 0; // false
      }
      else {
	std::cerr<<"Invalid option: calib_params must be true or false. \n";
        abort();
      }
      
      continue;
    }
  }

  fp.close();

  if (verbosity.compare("high") == 0) {
    std::string flag = state->lgar_bmi_params.use_closed_form_G == true ? "Yes" : "No";
    std::cerr<<"Using closed_form_G? "<< flag <<"\n";
    std::cerr<<"          *****         \n";
  }

  if (verbosity.compare("high") == 0) {
    std::string flag = state->lgar_bmi_params.PET_affects_precip == true ? "Yes" : "No";
    std::cerr<<"Does AET reduce precip? "<< flag <<"\n";
    std::cerr<<"          *****         \n";
  }

  if (verbosity.compare("high") == 0) {
    std::string flag = state->lgar_bmi_params.allow_flux_caching == true ? "Yes" : "No";
    std::cerr<<"Will fluxes be cached and used for subsequent time steps rather than computed during dry conditions? "<< flag <<"\n";
    std::cerr<<"          *****         \n";
  }

  if (verbosity.compare("high") == 0) {
    std::string flag = state->lgar_bmi_params.sft_coupled == true ? "Yes" : "No";
    std::cerr<<"Coupled to SoilFreezeThaw? "<< flag <<"\n";
    std::cerr<<"          *****         \n";
  }
  
  if(!is_max_valid_soil_types_set)
     state->lgar_bmi_params.num_soil_types = MAX_NUM_SOIL_TYPES;     // maximum number of valid soil types defaults to 15

  if (verbosity.compare("high") == 0) {
    std::cerr<<"Maximum number of soil types: "<<state->lgar_bmi_params.num_soil_types<<"\n";
    std::cerr<<"          *****         \n";
  }

  if (!is_layer_soil_type_set) {
    stringstream errMsg;
    errMsg << "The configuration file \'" << config_file <<"\' does not set layer_soil_type. \n";
    throw runtime_error(errMsg.str());
  }

  if (state->lgar_bmi_params.log_mode){
    state->lgar_bmi_params.a = pow(10.0, state->lgar_bmi_params.a);
    if (is_a_slow_set){
      state->lgar_bmi_params.a_slow = pow(10.0, state->lgar_bmi_params.a_slow);
    }
  }
    
  if(is_soil_params_file_set) {
    //allocate memory to create an array of structures to hold the soils properties data.
    //state->soil_properties = (struct soil_properties_*) malloc((state->lgar_bmi_params.num_layers+1)*sizeof(struct soil_properties_));


    state->soil_properties = new soil_properties_[state->lgar_bmi_params.num_soil_types+1];
    int num_soil_types = state->lgar_bmi_params.num_soil_types;
    double wilting_point_psi_cm = state->lgar_bmi_params.wilting_point_psi_cm;
    lgar_read_vG_param_file(soil_params_file.c_str(), num_soil_types,
						    wilting_point_psi_cm, state->soil_properties, state->lgar_bmi_params.log_mode);

    // check if soil layers provided are within the range
    state->lgar_bmi_params.is_invalid_soil_type = false; // model not valid for soil types = waterbody, glacier, lava, etc.
    for (int layer=1; layer <= state->lgar_bmi_params.num_layers; layer++) {
      //assert (state->lgar_bmi_params.layer_soil_type[layer] <= state->lgar_bmi_params.num_soil_types);
      //assert (state->lgar_bmi_params.layer_soil_type[layer] <= max_num_soil_in_file);
      if (state->lgar_bmi_params.layer_soil_type[layer] > state->lgar_bmi_params.num_soil_types) {
	state->lgar_bmi_params.is_invalid_soil_type = true;
	if (verbosity.compare("high") == 0) {
	  std::cerr << "Invalid soil type: "
		    << state->lgar_bmi_params.layer_soil_type[layer]
		    <<". Model returns input_precip = ouput_Qout. \n";
	}
	break;
      }
    }

    if (verbosity.compare("high") == 0) {
      for (int layer=1; layer<=state->lgar_bmi_params.num_layers; layer++) {
	int soil = state->lgar_bmi_params.layer_soil_type[layer];
	std::cerr<<"Soil type/name : "<<state->lgar_bmi_params.layer_soil_type[layer]
		 <<" "<<state->soil_properties[soil].soil_name<<"\n";
      }
      std::cerr<<"          *****         \n";
    }
  }
  else {
    stringstream errMsg;
    errMsg << "The configuration file \'" << config_file <<"\' does not set soil_params_file. \n";
    throw runtime_error(errMsg.str());
  }
  
  if (!is_layer_thickness_set) {
    stringstream errMsg;
    errMsg << "The configuration file \'" << config_file <<"\' does not set layer_thickness. \n";
    throw runtime_error(errMsg.str());
  }

  if (!is_initial_psi_set) {
    stringstream errMsg;
    errMsg << "The configuration file \'" << config_file <<"\' does not set initial_psi. \n";
    throw runtime_error(errMsg.str());
  }

  if (!is_timestep_set) {
    stringstream errMsg;
    errMsg << "The configuration file \'" << config_file <<"\' does not set timestep. \n";
    throw runtime_error(errMsg.str());
  }

  if (!is_endtime_set) {
    stringstream errMsg;
    errMsg << "The configuration file \'" << config_file <<"\' does not set endtime. \n";
    throw runtime_error(errMsg.str());
  }

  if(!is_wilting_point_psi_cm_set) {
    stringstream errMsg;
    errMsg << "The configuration file \'" << config_file <<"\' does not set wilting_point_psi. \n Recommended value of 15495.0[cm], corresponding to 15 atm. \n";
    throw runtime_error(errMsg.str());
  }

  if(!is_field_capacity_psi_cm_set) {
    stringstream errMsg;
    errMsg << "The configuration file \'" << config_file <<"\' does not set field_capacity_psi. \n Recommended value of 340.9[cm] for most soils, corresponding to 1/3 atm, or 103.3[cm] for sands, corresponding to 1/10 atm. \n";
    throw runtime_error(errMsg.str());
  }

  if(!is_root_zone_depth_cm_set && state->lgar_bmi_params.TO_enabled) {
    stringstream errMsg;
    errMsg << "root zone depth not set in the config file while TO mode is enabled "<< config_file << "\n";
    throw runtime_error(errMsg.str());
  }

  if (! ( (is_a_set == is_b_set) && (is_frac_to_CR_set == is_b_set)) ){
    //in this case, it must be either the case that all of these have been set (the user wants a nonlinear reservoir), or that none of these are set (the user does not want this).
    //it can not be the case that only one or two of these three have been set.
    stringstream errMsg;
    errMsg << "The configuration file \'" << config_file <<"\' does not correctly set a, b, and frac_to_CR. Either all or none must be set. a and b must be 0 or greater and frac_to_CR must be between 0 and 1. \n";
    throw runtime_error(errMsg.str());
  }

  if (! ( (is_a_slow_set == is_b_slow_set) && (is_frac_slow_set == is_b_slow_set)) ){
    //in this case, it must be either the case that all of these have been set (the user wants a second nonlinear reservoir), or that none of these are set (the user does not want this).
    //technically you can set the "slow" reservoir and not the other one -- in either case it amounts to 1 nonlinear reservoir.
    //it can not be the case that only one or two of these three have been set.
    stringstream errMsg;
    errMsg << "The configuration file \'" << config_file <<"\' does not correctly set a_slow, b_slow, and frac_slow. Either all or none must be set. a_slow and b_slow must be 0 or greater and frac_slow must be between 0 and 1 (but greater than 0). \n";
    throw runtime_error(errMsg.str());
  }

  if ( ((is_non_vadose_state_path_set && !is_state_path_set) || (!is_non_vadose_state_path_set && is_state_path_set)) && (is_a_set) ){
    stringstream errMsg;
    errMsg << "The configuration file \'" << config_file <<"\' sets one of init_non_vadose_state_path or init_state_path but not both, while a nonlinear reservoir is desired. Either both or none of init_non_vadose_state_path and init_state_path must be set in this case. \n";
    throw runtime_error(errMsg.str());
  }

  if (state->lgar_bmi_params.lower_bdy_flux_to_CR &&
      !state->lgar_bmi_params.free_drainage_enabled &&
      !state->lgar_bmi_params.TO_enabled){
    stringstream errMsg;
    errMsg << "The configuration file \'" << config_file <<"\' sets lower_bdy_flux_to_CR as true but sets both free_drainage_enabled and TO_enabled as false. lower_bdy_flux_to_CR requires either LGAR free drainage or LGARTO lower-boundary exchange. \n";
    throw runtime_error(errMsg.str());
  }

  if (state->lgar_bmi_params.TO_enabled && state->lgar_bmi_params.allow_flux_caching) {
    std::cerr << "Warning: allow_flux_caching=true with TO_enabled=true. "
              << "LGARTO flux caching is experimental; check global mass balance and hydrograph outputs.\n";
  }

  if (is_giuh_ordinates_set && !is_giuh_state_path_set && is_state_path_set){
    stringstream errMsg;
    errMsg << "The configuration file \'" << config_file <<"\' sets GIUH ordinates and is loading wetting fronts but is not loading a GIUH queue. Either do not set GIUH ordinates, provide a GIUH file to load from, or do not load a soil moisture profile. \n";
    throw runtime_error(errMsg.str());
  }

  if (!is_forcing_resolution_set) {
    stringstream errMsg;
    errMsg << "The configuration file \'" << config_file <<"\' does not set forcing_resolution. \n";
    throw runtime_error(errMsg.str());
  }

  if (is_giuh_ordinates_set) {
    int factor = int(1.0/state->lgar_bmi_params.forcing_resolution_h);

    state->lgar_bmi_params.num_giuh_ordinates = factor * (giuh_ordinates_temp.size() - 1);
    state->lgar_bmi_params.giuh_ordinates = new double[state->lgar_bmi_params.num_giuh_ordinates+1];
    
    for (unsigned int i=0; i<giuh_ordinates_temp.size()-1; i++) {
      for (int j=0; j<factor; j++) {
        int index = j + i * factor + 1;
        state->lgar_bmi_params.giuh_ordinates[index] = giuh_ordinates_temp[i+1]/double(factor);
      }
    }
    
    if (verbosity.compare("high") == 0) {
      for (int i=1; i<=state->lgar_bmi_params.num_giuh_ordinates; i++)
	      std::cerr<<"GIUH ordinates (scaled) : "<<state->lgar_bmi_params.giuh_ordinates[i]<<"\n";
      
      std::cerr<<"          *****         \n";
    }
    giuh_ordinates_temp.clear();
  }

  else if (!is_giuh_ordinates_set) {
    stringstream errMsg;
    errMsg << "The configuration file \'" << config_file <<"\' does not set giuh_ordinates. \n";
    throw runtime_error(errMsg.str());
  }

  if (state->lgar_bmi_params.sft_coupled) {
    state->lgar_bmi_params.soil_temperature = new double[state->lgar_bmi_params.num_cells_temp]();
    if (!is_soil_z_set) {
      stringstream errMsg;
      errMsg << "The configuration file \'" << config_file <<"\' does not set soil_z. \n";
      throw runtime_error(errMsg.str());
    }
  }
  else {
    //NJF FIXME these arrays should be allocated based on num_cells_temp...
    state->lgar_bmi_params.soil_temperature   = new double[1]();
    state->lgar_bmi_params.soil_temperature_z = new double[1]();
    state->lgar_bmi_params.num_cells_temp     = 1;
  }

  if (!is_ponded_depth_max_cm_set)
    state->lgar_bmi_params.ponded_depth_max_cm = 0.0; // default maximum ponded depth is set to zero (i.e. no surface ponding)


  state->lgar_bmi_params.forcing_interval = int(state->lgar_bmi_params.forcing_resolution_h/state->lgar_bmi_params.timestep_h+1.0e-08); // add 1.0e-08 to prevent truncation error

  // initialize frozen factor array to 1.
  state->lgar_bmi_params.frozen_factor = new double[state->lgar_bmi_params.num_layers+1];
  for (int i=0; i <= state->lgar_bmi_params.num_layers; i++)
    state->lgar_bmi_params.frozen_factor[i] = 1.0;

  if (!is_state_path_set){
    InitializeWettingFronts(state->lgar_bmi_params.TO_enabled, state->lgar_bmi_params.num_layers,
          state->lgar_bmi_params.initial_psi_cm, state->lgar_bmi_params.layer_soil_type,
          state->lgar_bmi_params.cum_layer_thickness_cm, state->lgar_bmi_params.layer_thickness_cm,
          state->lgar_bmi_params.frozen_factor, &state->head, state->soil_properties);
  }
  else {
    InitializeWettingFrontsFromCSV( //note that loading can yield a small mass balance error if theta and psi values were not recorded to high precision
        //will load the first line of the .csv 
        state->lgar_bmi_params.num_layers,
        state->lgar_bmi_params.init_state_path.c_str(),
        state->lgar_bmi_params.layer_soil_type,
        state->lgar_bmi_params.cum_layer_thickness_cm,
        state->lgar_bmi_params.frozen_factor,
        &state->head,
        state->soil_properties
    );
    if (!state->lgar_bmi_params.init_non_vadose_state_path.empty()) {
      InitializenonvadoseStateFromCSV(
          state->lgar_bmi_params.init_non_vadose_state_path.c_str(),
          state);
    }
  }
  
  if (verbosity.compare("none") != 0) {
    std::cerr<<"--- Initial state/conditions --- \n";
    listPrint(state->head);
    std::cerr<<"          *****         \n";
  }

  // initial mass in the system
  state->lgar_mass_balance.volstart_cm      = lgar_calc_mass_bal(state->lgar_bmi_params.cum_layer_thickness_cm, state->head);

  state->lgar_bmi_params.ponded_depth_cm    = 0.0; // initially we start with a dry surface (no surface ponding)
  state->lgar_bmi_params.nint               = 120; // hacked, not needed to be an input option

  // state->lgar_bmi_params.num_wetting_fronts = state->lgar_bmi_params.num_layers; // only if using InitializeWettingFronts and not InitializeWettingFrontsFromCSV
  state->lgar_bmi_params.num_wetting_fronts = listLength(state->head);
  // assert (state->lgar_bmi_params.num_layers == listLength(state->head)); // only if InitializeWettingFronts used and not InitializeWettingFrontsFromCSV

  if (verbosity.compare("high") == 0) {
    std::cerr<<"Initial ponded depth is set to zero. \n";
    std::cerr<<"No. of spatial intervals used in trapezoidal integration to compute G : "<<state->lgar_bmi_params.nint<<"\n";
  }

  state->lgar_bmi_input_params     = new lgar_bmi_input_parameters;
  state->lgar_bmi_params.time_s    = 0.0;
  state->lgar_bmi_params.timesteps = 0.0;

  if (verbosity.compare("none") != 0) {
    std::cerr<<"------------- Initialization done! ---------------------- \n";
    std::cerr<<"--------------------------------------------------------- \n";
  }

}

//##############################################################################
/*
  calculates initial theta (soil moisture content) and hydraulic conductivity
  from the prescribed psi value for each of the soil layers
*/
// #############################################################################
extern void InitializeWettingFronts(bool TO_enabled, int num_layers, double initial_psi_cm, int *layer_soil_type,
				    double *cum_layer_thickness_cm, double *layer_thickness_cm, double *frozen_factor,
				    struct wetting_front** head, struct soil_properties_ *soil_properties)
{
  if (TO_enabled) {
    int soil;
    int layer = 1;
    double Se, theta_init;
    bool bottom_flag;
    double Ksat_cm_per_h;
    struct wetting_front *current;
    const int number_of_WFs_per_layer = 16;
    bool switch_to_next_layer_flag = false;
    double prior_psi_cm = cum_layer_thickness_cm[num_layers];
    double new_wf_depth;
    int wf_in_layer = 1;
    double extra_moisture_factor_cm = 0.0;
    double extra_height_factor = 0.0;

    for (int front = 1; front <= (num_layers * number_of_WFs_per_layer); front++) {
      soil = layer_soil_type[layer];
      double total_depth = cum_layer_thickness_cm[num_layers];

      if ((front % number_of_WFs_per_layer) == 0) {
        wf_in_layer = 1;
        initial_psi_cm = (layer_thickness_cm[layer] / number_of_WFs_per_layer) - extra_moisture_factor_cm;
        if (layer < num_layers) {
          initial_psi_cm =
            (total_depth -
             (cum_layer_thickness_cm[layer - 1] +
              (number_of_WFs_per_layer - 1) * layer_thickness_cm[layer] / number_of_WFs_per_layer)) -
            extra_moisture_factor_cm;
        }
        new_wf_depth = cum_layer_thickness_cm[layer - 1] + layer_thickness_cm[layer];
        prior_psi_cm = initial_psi_cm;
      }
      else {
        if (wf_in_layer == 1) {
          initial_psi_cm = prior_psi_cm;
        }
        else {
          initial_psi_cm =
            (total_depth -
             (cum_layer_thickness_cm[layer - 1] +
              (wf_in_layer - 1) * layer_thickness_cm[layer] / number_of_WFs_per_layer)) -
            extra_moisture_factor_cm;
        }
        new_wf_depth =
          (cum_layer_thickness_cm[layer - 1] +
           wf_in_layer * layer_thickness_cm[layer] / number_of_WFs_per_layer) -
          extra_height_factor;
      }

      if (initial_psi_cm < 0.0) {
        initial_psi_cm = 0.0;
      }

      theta_init = calc_theta_from_h(initial_psi_cm, soil_properties[soil].vg_alpha_per_cm,
				     soil_properties[soil].vg_m, soil_properties[soil].vg_n,
				     soil_properties[soil].theta_e, soil_properties[soil].theta_r);

      if (verbosity.compare("high") == 0) {
        printf("layer, theta, psi, alpha, m, n, theta_e, theta_r = %d, %6.6f, %6.6f, %6.6f, %6.6f, %6.6f, %6.6f, %6.6f \n",
	       layer, theta_init, initial_psi_cm, soil_properties[soil].vg_alpha_per_cm,
	       soil_properties[soil].vg_m, soil_properties[soil].vg_n,
	       soil_properties[soil].theta_e, soil_properties[soil].theta_r);
      }

      bottom_flag = ((front % number_of_WFs_per_layer) == 0);

      if (new_wf_depth < cum_layer_thickness_cm[layer - 1]) {
        new_wf_depth = cum_layer_thickness_cm[layer - 1];
      }

      current = listInsertFront(new_wf_depth, theta_init, front, layer, bottom_flag, head);
      current->psi_cm = initial_psi_cm;
      current->is_WF_GW = true;

      Se = calc_Se_from_theta(current->theta, soil_properties[soil].theta_e, soil_properties[soil].theta_r);
      Ksat_cm_per_h = frozen_factor[layer] * soil_properties[soil].Ksat_cm_per_h;
      current->K_cm_per_h = calc_K_from_Se(Se, Ksat_cm_per_h, soil_properties[soil].vg_m);

      if (switch_to_next_layer_flag) {
        layer++;
        switch_to_next_layer_flag = false;
      }
      else {
        wf_in_layer++;
        if (((front + 1) % number_of_WFs_per_layer) == 0) {
          switch_to_next_layer_flag = true;
        }
      }
    }
  }
  else {
    int soil;
    int front = 0;
    double Se, theta_init;
    bool bottom_flag;
    double Ksat_cm_per_h;
    struct wetting_front *current;

    for(int layer=1;layer<=num_layers;layer++) {
      front++;

      soil = layer_soil_type[layer];
      theta_init = calc_theta_from_h(initial_psi_cm,soil_properties[soil].vg_alpha_per_cm,
				     soil_properties[soil].vg_m,soil_properties[soil].vg_n,
				     soil_properties[soil].theta_e,soil_properties[soil].theta_r);

      if (verbosity.compare("high") == 0) {
        printf("layer, theta, psi, alpha, m, n, theta_e, theta_r = %d, %6.6f, %6.6f, %6.6f, %6.6f, %6.6f, %6.6f, %6.6f \n",
	       layer, theta_init, initial_psi_cm, soil_properties[soil].vg_alpha_per_cm, soil_properties[soil].vg_m,
	       soil_properties[soil].vg_n,soil_properties[soil].theta_e,soil_properties[soil].theta_r);
      }

      bottom_flag = true;
      current = listInsertFront(cum_layer_thickness_cm[layer],theta_init,front,layer,bottom_flag, head);

      current->psi_cm = initial_psi_cm;
      current->is_WF_GW = false;
      Se = calc_Se_from_theta(current->theta,soil_properties[soil].theta_e,soil_properties[soil].theta_r);

      Ksat_cm_per_h = frozen_factor[layer] * soil_properties[soil].Ksat_cm_per_h;
      current->K_cm_per_h = calc_K_from_Se(Se, Ksat_cm_per_h , soil_properties[soil].vg_m);
    }
  }
}

// ##################################################################################
/*
  Reads 1D data from the config file
  - used for reading soil discretization (1D)
  - used for reading layers depth from the surface if model `layered` is chosen
*/
// ##################################################################################
extern vector<double>
ReadVectorData(string key)
{
  int pos = 0;
  string delimiter = ",";
  vector<double> value(0.0);
  string z1 = key;

  if (z1.find(delimiter) == string::npos) {
    double v = stod(z1);
    value.push_back(v);

  }
  else {
    while (z1.find(delimiter) != string::npos) {
      pos = z1.find(delimiter);
      string z_v = z1.substr(0, pos);

      value.push_back(stod(z_v.c_str()));

      z1.erase(0, pos + delimiter.length());
      if (z1.find(delimiter) == string::npos)
	value.push_back(stod(z1));
    }
  }

  return value;
}

// ############################################################################################
/*
  calculates frozen factor based on L. Wang et al. (www.hydrol-earth-syst-sci.net/14/557/2010/)
  uses layered-average soil temperatures and an exponential function to compute frozen fraction
  for each layer
 */
// ############################################################################################
extern void frozen_factor_hydraulic_conductivity(struct lgar_bmi_parameters lgar_bmi_params)
{

  int c = 0, count;
  double layer_temp;
  double factor;

  for (int layer=1; layer<=lgar_bmi_params.num_layers; layer++) {
    layer_temp = 0.0;
    count = 0;

    while (lgar_bmi_params.soil_temperature_z[c] <= lgar_bmi_params.cum_layer_thickness_cm[layer]) {
      layer_temp +=lgar_bmi_params.soil_temperature[c];
      c++;
      count++;

      if (c == lgar_bmi_params.num_cells_temp) // to ensure we don't access out of bound values
	break;
    }


    // assert (layer_temp > 100.0); // just a check to ensure the while loop executes at least once
    assert (count > 0); // just a check to ensure the while loop executes at least once

    layer_temp /= count;  // layer-averaged temperature

    factor = exp(-10 * (273.15 - layer_temp)); /* Eq. 6 (L. Wang et al.,Frozen soil parameterization in a distributed
                                                biosphere hydrological model, www.hydrol-earth-syst-sci.net/14/557/2010/)
						and Eq. 22 (Bao et al., An enthalpy-based frozen model) */

    factor = fmax(fmin(factor,1.0), 0.05); // 0.05 <= factor <= 1.0
    lgar_bmi_params.frozen_factor[layer] = factor;

  }

  if (verbosity.compare("high") == 0) {
    for (int i=1; i <= lgar_bmi_params.num_layers; i++)
      std::cerr<<"frozen factor = "<< lgar_bmi_params.frozen_factor[i]<<"\n";
  }

}

/*
extern void lgar_update(struct model_state *state)
{ if we ever decided to run this version without the bmi then we simply need to copy `update method` from the bmi here.}
*/


// #########################################################################################
/*
  calculates global mass balance at the end of simulation
*/
// #########################################################################################
extern void lgar_global_mass_balance(struct model_state *state, double *giuh_runoff_queue_cm)
{
  double volstart           = state->lgar_mass_balance.volstart_cm;
  double volprecip          = state->lgar_mass_balance.volprecip_cm;
  double volrunoff          = state->lgar_mass_balance.volrunoff_cm;
  double volrunoff_CR       = state->lgar_mass_balance.volrunoff_CR_cm;
  double volAET             = state->lgar_mass_balance.volAET_cm;
  double volPET             = state->lgar_mass_balance.volPET_cm;
  double volon              = state->lgar_mass_balance.volon_cm;
  double volin              = state->lgar_mass_balance.volin_cm;
  double volrech            = state->lgar_mass_balance.volrech_cm;
  double volend             = state->lgar_mass_balance.volend_cm;
  double volCRend           = state->lgar_mass_balance.volCRend_cm;
  double volCRstart         = state->lgar_mass_balance.volCRstart_cm;
  double volrunoff_giuh     = state->lgar_mass_balance.volrunoff_giuh_cm;
  double volend_giuh_cm     = 0.0;
  double total_Q_cm         = state->lgar_mass_balance.volQ_cm;
  double volchange_calib_cm = state->lgar_mass_balance.volchange_calib_cm;
  
  //check if the giuh queue have some water left at the end of simulaiton; needs to be included in the global mass balance
  // hold on; this is probably not needed as we have volrunoff in the balance; revist AJK
  for(int i=0; i <= state->lgar_bmi_params.num_giuh_ordinates; i++)
    volend_giuh_cm += giuh_runoff_queue_cm[i];

  double global_error_cm = volstart + volprecip - volrunoff - volAET - volon - volrech - volend + volchange_calib_cm - volrunoff_CR - volCRend + volCRstart;
  
  printf("\n********************************************************* \n");
  printf("-------------------- Simulation Summary ----------------- \n");
  //printf("Time (sec)                 = %6.10f \n", elapsed);
  printf("------------------------ Mass balance ------------------- \n");
  printf("Initial water in soil       = %14.10f cm\n", volstart);
  printf("Total precipitation         = %14.10f cm\n", volprecip);
  printf("Total infiltration (matrix) = %14.10f cm\n", volin);
  printf("Final water in soil         = %14.10f cm\n", volend);
  printf("Surface ponded water        = %14.10f cm\n", volon);
  printf("Surface runoff              = %14.10f cm\n", volrunoff);
  printf("GIUH runoff                 = %14.10f cm\n", volrunoff_giuh);
  printf("GIUH water (in array)       = %14.10f cm\n", volend_giuh_cm);
  printf("Total percolation           = %14.10f cm\n", volrech);
  printf("Total AET                   = %14.10f cm\n", volAET);
  printf("Total PET                   = %14.10f cm\n", volPET);
  if (state->lgar_bmi_params.frac_to_CR || state->lgar_bmi_params.lower_bdy_flux_to_CR){
    printf("storage in reservoir        = %14.10f cm\n", state->lgar_mass_balance.CR_fast_storage_cm);
    if (state->lgar_bmi_params.frac_slow){
      printf("storage in slow reservoir   = %14.10f cm\n", state->lgar_mass_balance.CR_slow_storage_cm);
      printf("runoff from CRs             = %14.10f cm\n", volrunoff_CR);
    }
    else {
    printf("runoff from CR              = %14.10f cm\n", volrunoff_CR);
    }
  }
  printf("Total discharge (Q)         = %14.10f cm\n", total_Q_cm);
  printf("Vol change (calibration)    = %14.10f cm\n", volchange_calib_cm);
  printf("Global balance              =   %.6e cm\n", global_error_cm);

}

// ############################################################################################
/*
 finds the wetting front that corresponds to psi (head) value closest to zero
 (i.e., saturation in terms of psi). This is the wetting front that experiences infiltration
 and actual ET based on precipitatona and PET, respectively. For example, the actual ET
 is extracted from this wetting front plus the wetting fronts above it.
 Note: the free_drainage name came from its python version, which is probably not the correct name.
 */
// ############################################################################################
extern int wetting_front_free_drainage(struct wetting_front* head) {

  int wf_that_supplies_free_drainage_demand = 1;
  struct wetting_front *current;

  //current = head;
  int number_of_wetting_fronts = listLength(head);

  for(current = head; current != NULL; current = current->next)
  {
    if (current->next != NULL) {
      if ((current->layer_num == current->next->layer_num) && (current->is_WF_GW==0))
	break;
      else
	wf_that_supplies_free_drainage_demand++;

    }
  }

  if (wf_that_supplies_free_drainage_demand > number_of_wetting_fronts)
    wf_that_supplies_free_drainage_demand--;

  if (verbosity.compare("high") == 0) {
    printf("wetting_front_free_drainage = %d \n", wf_that_supplies_free_drainage_demand);
  }

  return  wf_that_supplies_free_drainage_demand;
}

static void lgar_apply_surface_depth_update_with_event_limit(struct wetting_front *current,
                                                             struct wetting_front *next,
                                                             double timestep_h,
                                                             double column_depth_cm,
                                                             int num_layers,
                                                             bool lgarto_active)
{
  if (current == NULL) {
    return;
  }

  const double original_depth_cm = current->depth_cm;
  const double original_dzdt_cm_per_h = current->dzdt_cm_per_h;
  const double projected_depth_cm = original_depth_cm + original_dzdt_cm_per_h * timestep_h;

  (void)next;
  (void)num_layers;

  double limited_depth_cm = projected_depth_cm;

  bool limit_at_domain_overshoot = false;
  if (lgarto_active && timestep_h > 0.0 && !current->is_WF_GW && !current->to_bottom &&
      original_dzdt_cm_per_h > 0.0 && column_depth_cm > 0.0) {
    const double max_surface_depth_cm = 1.10 * column_depth_cm;
    if (limited_depth_cm > max_surface_depth_cm) {
      limited_depth_cm = max_surface_depth_cm;
      limit_at_domain_overshoot = true;
    }
  }

  if (!limit_at_domain_overshoot) {
    current->depth_cm = projected_depth_cm;
    return;
  }

  current->depth_cm = limited_depth_cm;
  current->dzdt_cm_per_h = (limited_depth_cm - original_depth_cm) / timestep_h;

  if (limit_at_domain_overshoot && verbosity.compare("high") == 0) {
    printf("Domain-limited LGARTO surface wetting front %d from projected depth %.17lf cm "
           "to %.17lf cm (110%% of column depth %.17lf cm; old_dzdt=%.17lf cm/h, "
           "new_dzdt=%.17lf cm/h).\n",
           current->front_num,
           projected_depth_cm,
           current->depth_cm,
           column_depth_cm,
           original_dzdt_cm_per_h,
           current->dzdt_cm_per_h);
  }
}

static void lgarto_convert_surface_fronts_drier_than_TO_below(struct wetting_front **head,
                                                              const char *correction_context)
{
  if (head == NULL || *head == NULL) {
    return;
  }

  bool converted_any = false;
  const int max_passes = listLength(*head) + 1;

  for (int pass = 0; pass < max_passes; pass++) {
    bool converted_this_pass = false;
    struct wetting_front *temp_WF = *head;

    while (temp_WF != NULL && temp_WF->next != NULL) {
      if (temp_WF->is_WF_GW == FALSE && temp_WF->next->is_WF_GW == TRUE &&
          temp_WF->theta < temp_WF->next->theta &&
          temp_WF->layer_num == temp_WF->next->layer_num) {
        if (verbosity.compare("high") == 0) {
          printf("%s converting front %d to GW because surface theta %.15f "
                 "is drier than GW theta %.15f below it in layer %d.\n",
                 correction_context, temp_WF->front_num, temp_WF->theta,
                 temp_WF->next->theta, temp_WF->layer_num);
        }
        temp_WF->is_WF_GW = TRUE;
        converted_any = true;
        converted_this_pass = true;
      }
      else if (temp_WF->to_bottom == TRUE && temp_WF->is_WF_GW == FALSE &&
               temp_WF->next->is_WF_GW == TRUE) {
        if (verbosity.compare("high") == 0) {
          printf("%s converting front %d to GW because a to_bottom surface front "
                 "is immediately above a GW front.\n",
                 correction_context, temp_WF->front_num);
        }
        temp_WF->is_WF_GW = TRUE;
        converted_any = true;
        converted_this_pass = true;
      }

      temp_WF = temp_WF->next;
    }

    if (!converted_this_pass) {
      break;
    }
  }

  if (converted_any && verbosity.compare("high") == 0) {
    printf("%s converted surf wf to TO \n", correction_context);
    listPrint(*head);
  }
}

static bool lgarto_surface_front_overtook_surface_front_above_TO_chain(const struct wetting_front *current)
{
  if (current == NULL || current->next == NULL || current->next->next == NULL) {
    return false;
  }

  const struct wetting_front *next = current->next;
  const struct wetting_front *next_to_next = next->next;

  return current->is_WF_GW == FALSE &&
         next->is_WF_GW == FALSE &&
         next_to_next->is_WF_GW == TRUE &&
         current->depth_cm > next->depth_cm &&
         current->theta > next->theta &&
         current->theta <= next_to_next->theta &&
         current->layer_num == next->layer_num &&
         !next->to_bottom;
}

static bool lgarto_convert_overtaken_surface_front_above_TO_chain(struct wetting_front **head,
                                                                  const char *correction_context)
{
  if (head == NULL || *head == NULL) {
    return false;
  }

  for (struct wetting_front *current = *head; current != NULL; current = current->next) {
    if (!lgarto_surface_front_overtook_surface_front_above_TO_chain(current)) {
      continue;
    }

    struct wetting_front *next = current->next;
    struct wetting_front *next_to_next = next->next;
    if (verbosity.compare("high") == 0) {
      printf("%s converting overtaking surface front %d and overtaken surface front %d "
             "to TO above an existing TO chain (overtaking depth %.15f theta %.15f; "
             "overtaken depth %.15f theta %.15f; TO below front %d theta %.15f).\n",
             correction_context,
             current->front_num,
             next->front_num,
             current->depth_cm,
             current->theta,
             next->depth_cm,
             next->theta,
             next_to_next->front_num,
             next_to_next->theta);
    }

    current->is_WF_GW = TRUE;
    next->is_WF_GW = TRUE;
    /*
     * These fronts already moved as surface fronts in this substep.  In the
     * movement loop, dzdt == 0 for a non-to_bottom front is the existing
     * one-substep "do not move this newly created/converted front" latch.
     * The next lgar_dzdt_calc call will recompute mobile TO/GW dzdt values.
     */
    current->dzdt_cm_per_h = 0.0;
    next->dzdt_cm_per_h = 0.0;
    return true;
  }

  return false;
}

static double lgarto_resolve_mixed_surface_surface_TO_overtake(int num_layers,
                                                               double *cum_layer_thickness_cm,
                                                               int *soil_type,
                                                               double *frozen_factor,
                                                               struct wetting_front **head,
                                                               struct soil_properties_ *soil_properties)
{
  const double mass_before_mixed_correction =
    lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
  const bool converted =
    lgarto_convert_overtaken_surface_front_above_TO_chain(head, "mixed surface/surface/TO correction");
  if (!converted) {
    return 0.0;
  }

  (void) num_layers;
  (void) soil_type;
  (void) frozen_factor;
  (void) soil_properties;
  double mass_balance_flux_correction_cm =
    mass_before_mixed_correction - lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
  if (verbosity.compare("high") == 0) {
    printf("mixed surface/surface/TO correction residual mass correction: %.17lf cm "
           "(flag-only conversion)\n",
           mass_balance_flux_correction_cm);
  }

  return mass_balance_flux_correction_cm;
}

// #######################################################################################################
/*
  the function moves wetting fronts, merge wetting fronts and does the mass balance correction when needed
  @param current : wetting front pointing to the current node of the current state
  @param next    : wetting front pointing to the next node of the current state
  @param previous    : wetting front pointing to the previous node of the current state
  @param current_old : wetting front pointing to the current node of the previous state
  @param next_old    : wetting front pointing to the next node of the previous state
  @param head : pointer to the first wetting front in the list of the current state

  Note: '_old' denotes the wetting_front or variables at the previous timestep (or state)
*/
// #######################################################################################################
extern double lgar_move_wetting_fronts(double timestep_h, double *free_drainage_subtimestep_cm, double *volin_cm, int wf_free_drainage_demand,
						     double old_mass, double cached_lower_boundary_flux_correction_cm, int num_layers, double *AET_demand_cm, double *cum_layer_thickness_cm,
						     int *soil_type, double *frozen_factor, struct wetting_front** head,
						     struct wetting_front* state_previous, struct soil_properties_ *soil_properties,
						     const double *surf_AET_vec, double PET_timestep_cm, double wilting_point_psi_cm,
						     double field_capacity_psi_cm, double root_zone_depth_cm, double surf_frac_rz,
						     double lgar_global_theta_snap_mass_tolerance_cm,
						     double groundwater_depth_cm)
{

  lgarto_clear_deferred_gw_flux_mass_balance_correction();

  if (verbosity.compare("high") == 0) {
    printf("State before moving wetting fronts...\n");
    listPrint(*head);
  }

  struct wetting_front *current;
  struct wetting_front *next;
  struct wetting_front *previous;

  struct wetting_front *current_old;
  struct wetting_front *next_old;

  double column_depth = cum_layer_thickness_cm[num_layers];
  double surface_lower_boundary_depth_cm = column_depth;
  if (std::isfinite(groundwater_depth_cm) && groundwater_depth_cm > 0.0 &&
      groundwater_depth_cm < column_depth) {
    surface_lower_boundary_depth_cm = groundwater_depth_cm;
  }

  previous = *head;
  double theta_e,theta_r;
  double vg_a, vg_m, vg_n;
  int layer_num, soil_num;

  int number_of_wetting_fronts = listLength(*head);
  int number_of_surface_WFs = 0;
  int number_of_TO_WFs_above_surface_WFs = 0;
  bool encountered_surface_WF = false;
  bool lgarto_active = false;

  for (struct wetting_front *count_front = *head; count_front != NULL; count_front = count_front->next) {
    if (count_front->is_WF_GW) {
      lgarto_active = true;
      if (!encountered_surface_WF) {
        number_of_TO_WFs_above_surface_WFs++;
      }
    }
    else {
      encountered_surface_WF = true;
      number_of_surface_WFs++;
    }
  }

  if (number_of_surface_WFs == 0) {
    number_of_TO_WFs_above_surface_WFs = 0;
  }

  current = *head;

  int last_wetting_front_index = number_of_wetting_fronts;
  int deepest_surface_front_index = number_of_surface_WFs + number_of_TO_WFs_above_surface_WFs;
  int top_most_surface_front_index = number_of_surface_WFs > 0 ? (number_of_TO_WFs_above_surface_WFs + 1) : 1;
  int layer_num_above, layer_num_below;
  double deepest_surf_depth_at_start = 0.0;

  for (struct wetting_front *depth_front = *head; depth_front != NULL; depth_front = depth_front->next) {
    if (!depth_front->is_WF_GW) {
      deepest_surf_depth_at_start = depth_front->depth_cm;
    }
  }

  double precip_mass_to_add = (*volin_cm); // water to be added to the soil

	  double bottom_boundary_flux_cm = 0.0; // water leaving the system through the bottom boundary
	  const bool use_TO_surface_AET = (surf_AET_vec != nullptr);
	  const double AET_pet_budget_cm = fmax(0.0, PET_timestep_cm * timestep_h);
	  double accepted_positive_AET_bookkeeping_cm = 0.0;

  auto add_AET_with_pet_budget = [&](double requested_increment_cm,
                                     const char *context,
                                     bool route_overflow_to_bottom_boundary) -> double {
    if (requested_increment_cm <= 0.0) {
      *AET_demand_cm += requested_increment_cm;
      return requested_increment_cm;
    }

    const double remaining_AET_budget_cm = fmax(0.0, AET_pet_budget_cm - *AET_demand_cm);
	    const double accepted_increment_cm = fmin(requested_increment_cm, remaining_AET_budget_cm);
	    const double overflow_cm = requested_increment_cm - accepted_increment_cm;
	    *AET_demand_cm += accepted_increment_cm;
	    if (use_TO_surface_AET) {
	      accepted_positive_AET_bookkeeping_cm += accepted_increment_cm;
	    }

    if (overflow_cm > MBAL_ITERATIVE_TOLERANCE) {
      if (verbosity.compare("high") == 0) {
        printf("Capping AET increment for %s from %.17lf cm to %.17lf cm "
               "because PET budget is %.17lf cm and current AET is %.17lf cm.\n",
               context, requested_increment_cm, accepted_increment_cm,
               AET_pet_budget_cm, *AET_demand_cm);
      }

      if (route_overflow_to_bottom_boundary) {
        lgarto_assert_gw_flux_mass_balance_correction_within_debug_threshold(
          overflow_cm,
          "lgar_move_wetting_fronts",
          context,
          *head);
        bottom_boundary_flux_cm += overflow_cm;
      }
    }

    return accepted_increment_cm;
  };

  *volin_cm = 0.0; // assuming that all the water can fit in, if not then re-assign the left over water at the end. now handled from the returned value from this function
  double free_drainage_demand = *free_drainage_subtimestep_cm;

  /* ************************************************************ */
  // main loop advancing only surface wetting fronts and doing the mass balance
  // groundwater wetting fronts are updated in a separate TO-specific block below
  // wf denotes wetting front

  for (int wf = deepest_surface_front_index; wf != 0; wf--) {

    current = listFindFront(wf, *head, NULL);
    if (current == NULL) {
      break;
    }

    if (current->is_WF_GW) {
      break;
    }

    if (verbosity.compare("high") == 0) {
      printf("Moving |******** Wetting Front = %d *********| \n", wf);
    }

    if (wf == 1 && number_of_wetting_fronts >0) {
      current = listFindFront(wf, *head, NULL);
      next = listFindFront(wf+1, *head, NULL);
      previous = NULL;

      current_old = listFindFront(wf, *head, state_previous);
      next_old = listFindFront(wf+1, *head, state_previous);
    }
    else if (wf < number_of_wetting_fronts) {
      current = listFindFront(wf, *head, NULL);
      next = listFindFront(wf+1, *head, NULL);
      previous = listFindFront(wf-1, *head, NULL);

      current_old = listFindFront(wf, *head, state_previous);
      next_old = listFindFront(wf+1, *head, state_previous);
    }
    else if (wf == number_of_wetting_fronts) {
      current = listFindFront(wf, *head, NULL);
      next = NULL;
      previous = listFindFront(wf-1, *head, NULL);

      current_old = listFindFront(wf, *head, state_previous);
      next_old = NULL;
    }

    layer_num   = current->layer_num;
    soil_num    = soil_type[layer_num];
    theta_e     = soil_properties[soil_num].theta_e;
    theta_r     = soil_properties[soil_num].theta_r;
    vg_a        = soil_properties[soil_num].vg_alpha_per_cm;
    vg_m        = soil_properties[soil_num].vg_m;
    vg_n        = soil_properties[soil_num].vg_n;

    // find indices of above and below layers
    layer_num_above = (wf == 1) ? layer_num : previous->layer_num;
    layer_num_below = (wf == last_wetting_front_index) ? layer_num + 1 : next->layer_num;

    if (verbosity.compare("high") == 0) {
       printf ("Layers (current, above, below) == %d %d %d \n", layer_num, layer_num_above, layer_num_below);
       listPrint(*head);
    }

    double actual_ET_demand = use_TO_surface_AET ? surf_AET_vec[wf] : *AET_demand_cm;
    if (use_TO_surface_AET) {
      actual_ET_demand =
        add_AET_with_pet_budget(actual_ET_demand, "surface AET demand", false);
    }

    // case to check if the wetting front is at the interface, i.e. deepest wetting front within a layer
    // psi of the layer below is already known/updated, so we that psi to compute the theta of the deepest current layer
    // todo. this condition can be replace by current->to_depth = FALSE && l<last_wetting_front_index
    /*             _____________
       layer_above             |
                            ___|
			   |
                   ________|____    <----- wetting fronts at the interface have same psi value
		           |
       layer current    ___|
                       |
                   ____|________
       layer_below     |
                  _____|________
    */
    /*************************************************************************************/
    if ( (wf < last_wetting_front_index) && (layer_num_below != layer_num) ) {
      
      if (verbosity.compare("high") == 0) {
	printf("case (deepest wetting front within layer) : layer_num (%d) != layer_num_below (%d) \n", layer_num, layer_num_below);
      }

      current->theta = calc_theta_from_h(next->psi_cm, vg_a,vg_m, vg_n, theta_e, theta_r);
      current->psi_cm = next->psi_cm;
    }

    // case to check if the number of wetting fronts are equal to the number of layers, i.e., one wetting front per layer
    /*************************************************************************************/
    /* For example, 3 layers and 3 wetting fronts in a state. psi profile is constant, and theta profile is non-uniform due
       to different van Genuchten parameters
                theta profile       psi profile  (constant head)
               _____________       ______________
                         |                   |
               __________|__       __________|___
	                     |                     |
               ________|____       __________|___
                   |                         |
               ____|________       __________|___
    */

    if (wf == number_of_wetting_fronts && layer_num_below != layer_num && number_of_wetting_fronts == num_layers) {

      if (verbosity.compare("high") == 0) {
	printf("case (number_of_wetting_fronts equal to num_layers) : l (%d) == num_layers (%d) == num_wetting_fronts(%d) \n", wf, num_layers,number_of_wetting_fronts);
      }

      // local variables
      double vg_a_k, vg_m_k, vg_n_k;
      double theta_e_k, theta_r_k;

      lgar_apply_surface_depth_update_with_event_limit(current, next, timestep_h, column_depth, num_layers,
                                                       lgarto_active); // this is probably not needed, as dz/dt = 0 for the deepest wetting front

      double *delta_thetas = (double *) malloc(sizeof(double)*(layer_num+1));
      double *delta_thickness = (double *) malloc(sizeof(double)*(layer_num+1));

      double psi_cm_old = current_old->psi_cm;
      //double psi_cm_below_old = 0.0;

      double psi_cm = current->psi_cm;
      //double psi_cm_below = 0.0;

      // mass = delta(depth) * delta(theta)
      double prior_mass = (current_old->depth_cm - cum_layer_thickness_cm[layer_num-1]) * (current_old->theta - 0.0); // 0.0 = next_old->theta

      double new_mass = (current->depth_cm - cum_layer_thickness_cm[layer_num-1]) * (current->theta - 0.0); // 0.0 = next->theta;

      for (int k=1; k<layer_num; k++) {
	int soil_num_k  = soil_type[k];
	theta_e_k = soil_properties[soil_num_k].theta_e;
	theta_r_k = soil_properties[soil_num_k].theta_r;
	vg_a_k    = soil_properties[soil_num_k].vg_alpha_per_cm;
	vg_m_k    = soil_properties[soil_num_k].vg_m;
	vg_n_k    = soil_properties[soil_num_k].vg_n;

	// using psi_cm_old for all layers because the psi is constant across layers in this particular case
	double theta_old             = calc_theta_from_h(psi_cm_old, vg_a_k, vg_m_k, vg_n_k, theta_e_k,theta_r_k);
	double theta_below_old       = 0.0;
	double local_delta_theta_old = theta_old - theta_below_old;
	double layer_thickness       = cum_layer_thickness_cm[k] - cum_layer_thickness_cm[k-1];

	prior_mass += (layer_thickness * local_delta_theta_old);

	double theta       = calc_theta_from_h(psi_cm, vg_a_k, vg_m_k, vg_n_k, theta_e_k, theta_r_k);
	double theta_below = 0.0;

	new_mass += layer_thickness * (theta - theta_below);
  //NJF theta_below is always 0, so all delta_thetas are always 0...
  //does this really need a dynamic array in this case???
	delta_thetas[k] = theta_below;
	delta_thickness[k] = layer_thickness;
      }

      delta_thetas[layer_num] = 0.0;
      delta_thickness[layer_num] = current->depth_cm - cum_layer_thickness_cm[layer_num-1];

      // double free_drainage_demand = 0;

      if (use_TO_surface_AET) {
        if (wf_free_drainage_demand == wf) {
		  prior_mass += precip_mass_to_add - (free_drainage_demand + cached_lower_boundary_flux_correction_cm);
        }
        prior_mass -= actual_ET_demand;
      }
      else if (wf_free_drainage_demand == wf) {
		prior_mass += precip_mass_to_add - (free_drainage_demand + cached_lower_boundary_flux_correction_cm + actual_ET_demand);
      }

      // theta mass balance computes new theta that conserves the mass; new theta is assigned to the current wetting front

      const double AET_before_mass_balance_cm = *AET_demand_cm;
      double theta_new = lgar_theta_mass_balance(layer_num, soil_num, psi_cm, new_mass, prior_mass, precip_mass_to_add, AET_demand_cm,
						 delta_thetas, delta_thickness, soil_type, soil_properties,
                         !use_TO_surface_AET);
      if (use_TO_surface_AET) {
        actual_ET_demand += (*AET_demand_cm - AET_before_mass_balance_cm);
      }
      else {
        actual_ET_demand = *AET_demand_cm;
      }
      //done with delta_thetas and delta_thickness, cleanup memory
      free(delta_thetas);
      free(delta_thickness);
      current->theta = fmax(theta_r, fmin(theta_new, theta_e));

      double Se = calc_Se_from_theta(current->theta,theta_e,theta_r);
      current->psi_cm = calc_h_from_Se(Se, vg_a, vg_m, vg_n);

      /* note: theta and psi of the current wetting front are updated here based on the wetting front's mass balance,
	 upper wetting fronts will be updated later in the lgar_merge_ module (the place where all state
	 variables are updated before proceeding to the next timestep. */

    }


    // case to check if the 'current' wetting front is within the layer and not at the layer's interface
    // layer_num == layer_num_below means there is another wetting front below the current wetting front
    // and they both belong to the same layer (in simple words, wetting fronts not at the interface)
    // l < last_wetting_front_index means that the current wetting front is not the deepest wetting front in the domain
    /*************************************************************************************/
  
    if ( (wf < last_wetting_front_index) && (layer_num == layer_num_below) ) {


      if (verbosity.compare("high") == 0) {
	printf("case (wetting front within a layer) : layer_num (%d) == layer_num_below (%d) \n", layer_num,layer_num_below);
      }

      // if wetting front is in the most surficial layer
      if (layer_num == 1) {

	// double free_drainage_demand = 0;
	// prior mass = mass contained in the current old wetting front
	double prior_mass = current_old->depth_cm * (current_old->theta -  next_old->theta);

	if (use_TO_surface_AET) {
	  if (wf_free_drainage_demand == wf) {
		    prior_mass += precip_mass_to_add - (free_drainage_demand + cached_lower_boundary_flux_correction_cm);
	  }
	  prior_mass -= actual_ET_demand;
	}
	else if (wf_free_drainage_demand == wf) {
		  prior_mass += precip_mass_to_add - (free_drainage_demand + cached_lower_boundary_flux_correction_cm + actual_ET_demand);
	}

			lgar_apply_surface_depth_update_with_event_limit(current, next, timestep_h, column_depth, num_layers,
	                                                     lgarto_active);

			/* condition to bound the wetting front depth, if depth of a wf, at this timestep,
			   gets greater than the domain depth, it will be merge anyway as it is passing
		   the layer depth */

  if (current->depth_cm > surface_lower_boundary_depth_cm && listLength_surface(*head)==listLength(*head)) {
    if (surface_lower_boundary_depth_cm >= column_depth - LOWER_BOUNDARY_FINAL_TOL_CM) {
	    current->depth_cm = column_depth + TRUNCATION_DEPTH; //we want WFs to exceed the lower boundary in the event that they must be partially truncated and then WFs above this one will correctly have their moisture corrected, but also want WFs to not exceed the lower boundary much
    }
	  }

	if (current->dzdt_cm_per_h == 0.0 && current->to_bottom == FALSE) // a new front was just created, so don't update it.
	  current->theta = current->theta;
	else {
      if ((prior_mass/current->depth_cm + next->theta)<theta_r){
        if (verbosity.compare("high") == 0) {
          printf("Deleting WF (%d) that will go below theta_r (before)...\n", current->front_num);
          listPrint(*head);
        }

        //the idea here is that in some cases, the reduction in theta via WF movement or AET will be intense enough such that theta goes below theta_r.
        //it requires a fairly unusual soil, which I encountered during random parameter sampling.
        double mass_before_theta_went_below_theta_r = lgar_calc_mass_bal(cum_layer_thickness_cm, *head) - current->depth_cm*(current->theta - (prior_mass/current->depth_cm + next->theta));
        current = listDeleteFront(current->front_num, head, soil_type, soil_properties);
        current = next;
        double mass_after_theta_went_below_theta_r = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
        const double reduced_aet_cm =
          fabs(mass_before_theta_went_below_theta_r - mass_after_theta_went_below_theta_r);
        *AET_demand_cm = *AET_demand_cm - reduced_aet_cm;
        if (use_TO_surface_AET) {
          actual_ET_demand -= reduced_aet_cm;
        }
        else {
          actual_ET_demand = *AET_demand_cm;
        }
        if (verbosity.compare("high") == 0) {
          printf("Deleting WF that will go below theta_r (after)...\n");
          listPrint(*head);
        }
      }
      else {//This is the case where theta>theta_r, which will be almost all of the time 
	      current->theta = fmax(theta_r, fmin(theta_e, prior_mass/current->depth_cm + next->theta));
      }
    }

      }
      else {

	/*
	  this note is copied from Python version:
	  "However, calculation of theta via mass balance is a bit trickier. This is because each wetting front
	  in deeper layers can be thought of as extending all the way to the surface, in terms of psi values.
	  For example, a wetting front in layer 2 with a theta value of 0.4 will in reality extend to layer
	  1 with a theta value that is different (usually smaller) due to different soil hydraulic properties.
	  But, the theta value of this extended wetting front is not recorded in current or previous states.
          So, simply from states, the mass balance of a wetting front that, in terms of psi, extends between
	  multiple layers cannot be calculated. Therefore, the theta values that the current wetting front *would*
	  have in above layers is calculated from the psi value of the current wetting front, with the assumption
	  that the hydraulic head of this wetting front is the same all the way up to the surface.

	  - LGAR paper (https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2022WR033742) has a better description, using diagrams, of the mass balance of wetting fronts
	*/

	double vg_a_k, vg_m_k, vg_n_k;
	double theta_e_k, theta_r_k;

			lgar_apply_surface_depth_update_with_event_limit(current, next, timestep_h, column_depth, num_layers,
	                                                     lgarto_active);

	  if (current->depth_cm > surface_lower_boundary_depth_cm && listLength_surface(*head)==listLength(*head)) {
      if (surface_lower_boundary_depth_cm >= column_depth - LOWER_BOUNDARY_FINAL_TOL_CM) {
		    current->depth_cm = column_depth + TRUNCATION_DEPTH; //we want WFs to exceed the lower boundary in the event that they must be partially truncated and then WFs above this one will correctly have their moisture corrected, but also want WFs to not exceed the lower boundary much
      }
	  }

	double *delta_thetas    = (double *)malloc(sizeof(double)*(layer_num+1));
	double *delta_thickness = (double *)malloc(sizeof(double)*(layer_num+1));


	double psi_cm_old = current_old->psi_cm;
	double psi_cm_below_old = current_old->next->psi_cm;

	double psi_cm = current->psi_cm;
	double psi_cm_below = next->psi_cm;
	double theta_mass_balance_psi_upper_limit_cm = PSI_UPPER_LIM;
	if (std::isfinite(psi_cm_below) && psi_cm_below > PSI_UPPER_LIM) {
	  /*
	   * A surface front in a lower layer can sit above an extremely dry TO/GW
	   * front.  In that local column the lower front, not the generic dry cap,
	   * is the physically relevant dry bound for mass closure.
	   */
	  theta_mass_balance_psi_upper_limit_cm = psi_cm_below;
	}

	// mass = delta(depth) * delta(theta)
	//      = difference in current and next wetting front thetas times depth of the current wetting front
	double prior_mass = (current_old->depth_cm - cum_layer_thickness_cm[layer_num-1]) * (current_old->theta - next_old->theta);
	double new_mass = (current->depth_cm - cum_layer_thickness_cm[layer_num-1]) * (current->theta - next->theta);

	// compute mass in the layers above the current wetting front
	// use the psi of the current wetting front and van Genuchten parameters of
	// the respective layers to get the total mass above the current wetting front
	for (int k=1; k<layer_num; k++) {
	  int soil_num_k  = soil_type[k];
	  theta_e_k = soil_properties[soil_num_k].theta_e;
	  theta_r_k = soil_properties[soil_num_k].theta_r;
	  vg_a_k    = soil_properties[soil_num_k].vg_alpha_per_cm;
	  vg_m_k    = soil_properties[soil_num_k].vg_m;
	  vg_n_k    = soil_properties[soil_num_k].vg_n;

	  double theta_old = calc_theta_from_h(psi_cm_old, vg_a_k, vg_m_k, vg_n_k, theta_e_k,theta_r_k);
	  double theta_below_old = calc_theta_from_h(psi_cm_below_old, vg_a_k, vg_m_k, vg_n_k, theta_e_k, theta_r_k);
	  double local_delta_theta_old = theta_old - theta_below_old;
	  double layer_thickness = (cum_layer_thickness_cm[k] - cum_layer_thickness_cm[k-1]);

	  prior_mass += (layer_thickness * local_delta_theta_old);

	  //-------------------------------------------
	  // do the same for the current state
	  double theta = calc_theta_from_h(psi_cm, vg_a_k, vg_m_k, vg_n_k, theta_e_k, theta_r_k);

	  double theta_below = calc_theta_from_h(psi_cm_below, vg_a_k, vg_m_k, vg_n_k, theta_e_k, theta_r_k);

	  new_mass += layer_thickness * (theta - theta_below);

	  delta_thetas[k] = theta_below;
	  delta_thickness[k] = layer_thickness;
	}

	delta_thetas[layer_num] = next->theta;
	delta_thickness[layer_num] = current->depth_cm - cum_layer_thickness_cm[layer_num-1];

	// double free_drainage_demand = 0;
  

	if (use_TO_surface_AET) {
	  if (wf_free_drainage_demand == wf) {
		    prior_mass += precip_mass_to_add - (free_drainage_demand + cached_lower_boundary_flux_correction_cm);
	  }
	  prior_mass -= actual_ET_demand;
	}
	else if (wf_free_drainage_demand == wf) {
		  prior_mass += precip_mass_to_add - (free_drainage_demand + cached_lower_boundary_flux_correction_cm + actual_ET_demand);
	}
  // theta mass balance computes new theta that conserves the mass; new theta is assigned to the current wetting front
	const double AET_before_mass_balance_cm = *AET_demand_cm;
	double theta_new = lgar_theta_mass_balance(layer_num, soil_num, psi_cm, new_mass, prior_mass, precip_mass_to_add, AET_demand_cm,
						   delta_thetas, delta_thickness, soil_type, soil_properties,
                           !use_TO_surface_AET,
                           theta_mass_balance_psi_upper_limit_cm);
  if (use_TO_surface_AET) {
    actual_ET_demand += (*AET_demand_cm - AET_before_mass_balance_cm);
  }
  else {
    actual_ET_demand = *AET_demand_cm;
  }
  //done with delta_thetas and delta_thickness, cleanup memory
  free(delta_thetas);
  free(delta_thickness);
	current->theta = fmax(theta_r, fmin(theta_new, theta_e));

      }

      double Se = calc_Se_from_theta(current->theta,theta_e,theta_r);
      current->psi_cm = calc_h_from_Se(Se, vg_a, vg_m, vg_n);

    }
  

    // if f_p (predicted infiltration) causes theta > theta_e, mass correction is needed.
    // depth of the wetting front is increased to close the mass balance when theta > theta_e.
    // l == 1 is the last iteration (top most wetting front), so do a check on the mass balance)
    // this part should be moved out of here to a subroutine; add a call to that subroutine
    if (wf == top_most_surface_front_index) { 

      wf_free_drainage_demand = wetting_front_free_drainage(*head);
      struct wetting_front *wf_free_drainage = listFindFront(wf_free_drainage_demand, *head, NULL);
    // if ((wf == wf_free_drainage_demand) && (current->theta>=theta_e) ) {
      int soil_num_k1  = soil_type[wf_free_drainage->layer_num];
      double theta_e_k1 = soil_properties[soil_num_k1].theta_e;

      const double total_ET_demand_cm = use_TO_surface_AET ? *AET_demand_cm : actual_ET_demand;
	      double mass_timestep = (old_mass + precip_mass_to_add) - (total_ET_demand_cm + free_drainage_demand + cached_lower_boundary_flux_correction_cm);

      assert (old_mass > 0.0);
      
      if (fabs(wf_free_drainage->theta - theta_e_k1) < 1E-15) {
	
      double current_mass = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);

      double mass_balance_error = fabs(current_mass - mass_timestep); // mass error

      double factor = 1.0;
      if (wf_free_drainage->next!=NULL){
        if (fabs(wf_free_drainage->theta - wf_free_drainage->next->theta)<0.01){// if two adjacent theta values are quite close, an initial factor of 1.0 might not make the mass balance close within 10000 iterations
          factor = factor / fabs(wf_free_drainage->theta - wf_free_drainage->next->theta);
        }
      }

      // double factor = fmax(1,current->psi_cm/100); speed optimization should look at optimal factor values 
      bool switched = false;

      double depth_new = wf_free_drainage->depth_cm;
      bool use_lgarto_saturated_depth_min_bound = false;
      bool lgarto_saturated_depth_bound_hit = false;
      double lgarto_saturated_depth_min_cm = TRUNCATION_DEPTH;
      if (lgarto_active && !wf_free_drainage->is_WF_GW && !wf_free_drainage->to_bottom &&
          wf_free_drainage->layer_num >= 1 && wf_free_drainage->layer_num <= num_layers) {
        const int bounded_layer = wf_free_drainage->layer_num;
        lgarto_saturated_depth_min_cm = bounded_layer == 1
          ? TRUNCATION_DEPTH
          : cum_layer_thickness_cm[bounded_layer - 1] + TRUNCATION_DEPTH;

        struct wetting_front *previous_free_drainage = NULL;
        if (wf_free_drainage->front_num > 1) {
          previous_free_drainage = listFindFront(wf_free_drainage->front_num - 1, *head, NULL);
        }
        if (previous_free_drainage != NULL &&
            previous_free_drainage->layer_num == bounded_layer) {
          lgarto_saturated_depth_min_cm =
            fmax(lgarto_saturated_depth_min_cm,
                 previous_free_drainage->depth_cm + DEPTH_AVOIDS_SAME_WF_DEPTH);
        }

        use_lgarto_saturated_depth_min_bound = true;
      }

      // loop to adjust the depth for mass balance
      int iter = 0;
      bool iter_aug_flag = FALSE;
      bool break_flag = FALSE;
      int speedup_thresh = MAX_ITER_SATURATION_MBAL_LOOP/10;

      if (fabs(mass_balance_error) > MBAL_ITERATIVE_TOLERANCE){
        if (verbosity.compare("high") == 0) {
          printf("start WF depth adjustment due to saturation");
          listPrint(*head);
        }
      }

      while (fabs(mass_balance_error) > MBAL_ITERATIVE_TOLERANCE) {
        iter++;
	        if (iter>MAX_ITER_SATURATION_MBAL_LOOP) {
	          break_flag = TRUE;
	          const double accepted_AET_increment_cm =
	            add_AET_with_pet_budget(mass_balance_error,
	                                    "saturated free-drainage depth adjustment residual",
	                                    false);
	          if (use_TO_surface_AET) {
	            actual_ET_demand += accepted_AET_increment_cm;
	          }
	          else {
	            actual_ET_demand = *AET_demand_cm;
          }
          break;
        }

        if ((iter>speedup_thresh) && (!iter_aug_flag)){
          factor = factor * 100;
          iter_aug_flag = TRUE;
        }

        if (current_mass < mass_timestep) {
          depth_new += 0.01 * factor;
          switched = false;
        }
        else {
          if (!switched) {
            switched = true;
            factor = factor * 0.1;
          }
          depth_new -= 0.01 * factor;

        }

        if ( (wf_free_drainage->to_bottom==TRUE) && (wf_free_drainage->layer_num==num_layers) ){
          depth_new = cum_layer_thickness_cm[num_layers];
        }

        // LGARTO event limiting can pin a saturated surface WF to a TO/to_bottom
        // boundary after infiltration was accepted. If the old saturated-depth
        // mass solve then wants to pull that front above its assigned layer,
        // keep the geometry valid and expose the remaining signed residual.
        if (use_lgarto_saturated_depth_min_bound &&
            depth_new < lgarto_saturated_depth_min_cm) {
          depth_new = lgarto_saturated_depth_min_cm;
          wf_free_drainage->depth_cm = depth_new;
          current_mass = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
          mass_balance_error = fabs(current_mass - mass_timestep);
          break_flag = TRUE;
          lgarto_saturated_depth_bound_hit = true;
          if (verbosity.compare("high") == 0) {
            printf("LGARTO saturated free-drainage depth adjustment hit lower "
                   "layer-bound %.17lf cm for front %d; "
                   "remaining mass residual %.17e cm will be routed as signed "
                   "lower-boundary correction.\n",
                   depth_new,
                   wf_free_drainage->front_num,
                   mass_timestep - current_mass);
          }
          break;
        }

        wf_free_drainage->depth_cm = depth_new;

        current_mass = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
        mass_balance_error = fabs(current_mass - mass_timestep);

      }

      if (depth_new<TRUNCATION_DEPTH){ // extremely rare error where, if the WFs below this one are extremely dry (psi values > 1.E7) and WF below have layer n values close to 1 (say 1.02 or so),
                                       // theta below this WF will sometimes change very slightly. If the next WF is thick enough, the current WF is thin enough, and the added infiltraiton is small enough, then 
                                       // technically mass balance will need a negative WF depth. Instead, we just accept the very small mass balance error and move on.
        // double mass_before = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
        depth_new=TRUNCATION_DEPTH;
        wf_free_drainage->depth_cm = depth_new;
        // double mass_after = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
        // mass_balance_error = fabs(mass_after - mass_before); //optionally could send mass balance error to AET to force closure, but this could yield an invalid AET value
        // *AET_demand_cm = *AET_demand_cm + mass_balance_error;
        // actual_ET_demand = *AET_demand_cm;
      }

      if (isinf(wf_free_drainage->depth_cm)){ // there is a rare case where the psi-theta relationship, for psi very close to 0, is not 1:1, so psi can technically change and theta=theta_e for either psi. 
                                              // Then, it can be that the WF below WF that accepts infiltration has a theta value equal to theta_e for its layer but a psi value that is very slightly above 0.
                                              // In this case, no adjustment of depth will close the mass balance so the depth will become infinite. 
        wf_free_drainage->depth_cm = cum_layer_thickness_cm[num_layers] + 1.E-6;
        break_flag = TRUE;
      }

      //there is a general class of problem where a very small psi value that is greater than 0 (say 1e-3 or so) will for some but not all soils mathematically yield theta = theta_e, even though theta should be slightly less than theta_e.
      //in layered soils, this can cause a mass balance error. It is fairly rare and only seems to impact cases where the model domain is entirely saturated, which shouldn't happen when LGAR is applied in the correct environment / with sufficient layer thicknesses.
      if (break_flag) {
        current_mass = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
        const double total_ET_demand_cm = use_TO_surface_AET ? *AET_demand_cm : actual_ET_demand;
	        mass_timestep = (old_mass + precip_mass_to_add) - (total_ET_demand_cm + free_drainage_demand + cached_lower_boundary_flux_correction_cm);
        mass_balance_error = mass_timestep - current_mass;
        lgarto_assert_gw_flux_mass_balance_correction_within_debug_threshold(
          mass_balance_error,
          "lgar_move_wetting_fronts",
          lgarto_saturated_depth_bound_hit
            ? "bounded saturated free-drainage depth residual"
            : "saturated free-drainage depth adjustment residual",
          *head);
        bottom_boundary_flux_cm += mass_balance_error;
      }

      }
    }
    
  }
  
  /*******************************************************************/
  // end of the for loop
  /*******************************************************************/


  if (verbosity.compare("high") == 0) {
    printf("State after moving but before merging wetting fronts...\n");
    listPrint(*head);
  }

  // ********************** MERGING AND CROSSING WETTING FRONT ****************************
  
  /* In the python version of LGAR, wetting front merging, layer boundary crossing, and lower boundary crossing
     all occur in a loop that happens after wetting fronts have been moved. This prevents the model from crashing,
     as there are rare but possible cases where multiple merging / layer boundary crossing events will happen in
     the same time step. For example, if two wetting fronts cross a layer boundary in the same time step, it will
     be necessary for merging to occur before layer boundary crossing. LGAR-C iteratively checks for these cases
     and addresses all until none are left via lgarto_correction_type_surf. */

  double mass_change = 0.0;

  int correction_type_surf =
    lgarto_correction_type_surf(num_layers, cum_layer_thickness_cm, head,
                                surface_lower_boundary_depth_cm);
  int surface_correction_iterations = 0;
  const int max_surface_correction_iterations = 1000;

  while (correction_type_surf!=0){
    surface_correction_iterations++;

    if (correction_type_surf==1){
      lgar_merge_wetting_fronts(soil_type, frozen_factor, head, soil_properties);
    }

    if (correction_type_surf==2){
      double mass_before_bdy_crossing = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
      lgar_wetting_fronts_cross_layer_boundary(num_layers, cum_layer_thickness_cm, soil_type, frozen_factor, head, soil_properties);
	      double mass_after_bdy_crossing  = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
	      if (fabs(mass_before_bdy_crossing - mass_after_bdy_crossing)>100.*MBAL_ITERATIVE_TOLERANCE){//the inclusion of 100.*MBAL_ITERATIVE_TOLERANCE is due to the fact that mass_before_bdy_crossing > mass_after_bdy_crossing might be true, but only within the mass balance tolerance, in which case we should not run this
	        (void) add_AET_with_pet_budget(mass_before_bdy_crossing - mass_after_bdy_crossing,
	                                       "surface layer-boundary crossing mass residual",
	                                       true);
	      }
	    }

    if (correction_type_surf==3){
      // bottom_boundary_flux_cm += lgar_wetting_front_cross_domain_boundary(TO_enabled, cum_layer_thickness_cm, soil_type, frozen_factor, head, soil_properties);
      bottom_boundary_flux_cm +=
        lgar_wetting_front_cross_domain_boundary(surface_lower_boundary_depth_cm,
                                                 soil_type, frozen_factor, head,
                                                 soil_properties);
          if (isnan(bottom_boundary_flux_cm)){
            bottom_boundary_flux_cm = 0.0;
          }
    }

	    if (correction_type_surf==4){
	      mass_change = 0.0;
	      lgar_fix_dry_over_wet_wetting_fronts(&mass_change, cum_layer_thickness_cm, soil_type, head, soil_properties);
	      (void) add_AET_with_pet_budget(-mass_change,
	                                     "surface dry-over-wet mass residual",
	                                     true);
	    }

    if (correction_type_surf==5){
      lgarto_convert_surface_fronts_drier_than_TO_below(head, "surface correction type 5");
    }

    if (correction_type_surf==6){
      const double mass_balance_flux_correction_cm =
        lgarto_resolve_mixed_surface_surface_TO_overtake(num_layers, cum_layer_thickness_cm,
                                                         soil_type, frozen_factor, head,
                                                         soil_properties);
      if (fabs(mass_balance_flux_correction_cm) > MBAL_ITERATIVE_TOLERANCE) {
        lgarto_assert_gw_flux_mass_balance_correction_within_debug_threshold(
          mass_balance_flux_correction_cm,
          "lgar_move_wetting_fronts",
          "mixed surface/surface/TO correction residual",
          *head);
        bottom_boundary_flux_cm += mass_balance_flux_correction_cm;
      }
    }

    if (correction_type_surf==7){
      bool merged_in_non_top_layer = false;
      const double mass_before_surface_TO_merge = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
      merged_in_non_top_layer =
        lgar_merge_surface_and_TO_wetting_fronts(merged_in_non_top_layer, num_layers,
                                                 cum_layer_thickness_cm, soil_type,
                                                 frozen_factor, soil_properties, head,
                                                 true);
      const bool did_a_WF_have_negative_depth = lgarto_correct_negative_depths(head);
      lgarto_cleanup_after_surface_TO_merging_in_layer_below_top(merged_in_non_top_layer, soil_type,
                                                                 soil_properties, head);
      const bool close_psis = correct_close_psis(soil_type, soil_properties, head);

      double mass_balance_flux_correction_cm =
        mass_before_surface_TO_merge - lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
      int boundary_pinned_surface_front_num = -1;
      if (fabs(mass_balance_flux_correction_cm) > MBAL_ITERATIVE_TOLERANCE) {
        mass_balance_flux_correction_cm =
          lgarto_restore_surface_TO_merge_mass_via_surface_depths(mass_before_surface_TO_merge,
                                                                  num_layers,
                                                                  cum_layer_thickness_cm,
                                                                  *head,
                                                                  &boundary_pinned_surface_front_num);
      }
      if (mass_balance_flux_correction_cm > MBAL_ITERATIVE_TOLERANCE &&
          boundary_pinned_surface_front_num > -1) {
        mass_balance_flux_correction_cm =
          lgarto_restore_surface_TO_merge_mass_via_boundary_psi(mass_before_surface_TO_merge,
                                                                boundary_pinned_surface_front_num,
                                                                num_layers,
                                                                cum_layer_thickness_cm,
                                                                soil_type,
                                                                frozen_factor,
                                                                soil_properties,
                                                                *head);
      }
      if (did_a_WF_have_negative_depth || close_psis ||
          fabs(mass_balance_flux_correction_cm) > MBAL_ITERATIVE_TOLERANCE) {
        lgarto_assert_gw_flux_mass_balance_correction_within_debug_threshold(
          mass_balance_flux_correction_cm,
          "lgar_move_wetting_fronts",
          "pre-TO-motion surface/TO merge cleanup mass residual",
          *head);
        bottom_boundary_flux_cm += mass_balance_flux_correction_cm;
      }
    }

    if (correction_type_surf==9){
      const double mass_balance_flux_correction_cm =
        lgarto_surface_fronts_cross_layer_boundary_upward(
          num_layers, cum_layer_thickness_cm, soil_type, frozen_factor,
          soil_properties, head, "surface correction loop");
      if (fabs(mass_balance_flux_correction_cm) > MBAL_ITERATIVE_TOLERANCE) {
        lgarto_assert_gw_flux_mass_balance_correction_within_debug_threshold(
          mass_balance_flux_correction_cm,
          "lgar_move_wetting_fronts",
          "surface upward layer-crossing residual",
          *head);
        bottom_boundary_flux_cm += mass_balance_flux_correction_cm;
      }
    }

    correction_type_surf =
      lgarto_correction_type_surf(num_layers, cum_layer_thickness_cm, head,
                                  surface_lower_boundary_depth_cm);
    if (verbosity.compare("high") == 0) {
      printf("correction_type_surf at end of iteration in while loop: %d \n", correction_type_surf);
    }
    if (correction_type_surf != 0 &&
        surface_correction_iterations >= max_surface_correction_iterations) {
      fprintf(stderr,
              "Error: surface wetting-front correction loop exceeded iteration cap.\n"
              "  While not technically physically impossible for a very long, complex wetting front list, it much more likely means that a correciton type was detected but not resolved, resulting in an infinite loop.\n"
              "  iterations=%d max_iterations=%d remaining_correction_type_surf=%d\n"
              "  Wetting front list follows:\n",
              surface_correction_iterations,
              max_surface_correction_iterations,
              correction_type_surf);
      fflush(stderr);
      listPrint(*head);
      fflush(stdout);
      abort();
    }
	  }

	  const double negative_depth_repair_flux_cm =
	    lgarto_repair_negative_depth_fronts_to_lower_boundary_flux(
	      num_layers, cum_layer_thickness_cm, soil_type, frozen_factor,
	      soil_properties, head, "post-surface correction loop");
	  if (fabs(negative_depth_repair_flux_cm) > MBAL_ITERATIVE_TOLERANCE) {
	    lgarto_assert_gw_flux_mass_balance_correction_within_debug_threshold(
	      negative_depth_repair_flux_cm,
	      "lgar_move_wetting_fronts",
	      "post-surface negative-depth repair residual",
	      *head);
	    bottom_boundary_flux_cm += negative_depth_repair_flux_cm;
	  }

		  const double bottom_boundary_flux_cm_from_surf_WFs = bottom_boundary_flux_cm;
		  double bottom_boundary_flux_above_surface_WFs_cm = 0.0;

  bool has_TO_fronts = false;
  for (struct wetting_front *count_front = *head; count_front != NULL; count_front = count_front->next) {
    if (count_front->is_WF_GW) {
      has_TO_fronts = true;
      break;
    }
  }

  if (has_TO_fronts) {
    if (use_TO_surface_AET) {
      if (verbosity.compare("high") == 0) {
        printf("before calc AET from TO WFs: \n");
        listPrint(*head);
      }

      const double cumulative_ET_from_TO_WFs_cm =
        lgarto_calc_aet_from_TO_WFs(num_layers, deepest_surf_depth_at_start, root_zone_depth_cm,
                                    PET_timestep_cm, timestep_h, surf_frac_rz,
	                                    precip_mass_to_add <= 1.0e-12,
	                                    wilting_point_psi_cm, field_capacity_psi_cm, soil_type,
	                                    cum_layer_thickness_cm, frozen_factor, soil_properties, head);
	      (void) add_AET_with_pet_budget(cumulative_ET_from_TO_WFs_cm,
	                                     "TO AET extraction",
	                                     true);

      if (verbosity.compare("high") == 0) {
        printf("cumulative_ET_from_TO_WFs_cm: %10.16lf \n", cumulative_ET_from_TO_WFs_cm);
        printf("after calc AET from TO WFs: \n");
        listPrint(*head);
      }
    }

    if (verbosity.compare("high") == 0) {
      printf("State before TO WF depth update via dzdt...\n");
      listPrint(*head);
    }

	    const int updated_wetting_front_count = listLength(*head);
	    bool TO_WFs_above_surface_WFs_flag = false;
	    for (int wf = updated_wetting_front_count - 1; wf >= 1; wf--) {
	      current = listFindFront(wf, *head, NULL);
	      if (current == NULL) {
	        continue;
	      }

	      if (!current->is_WF_GW) {
	        TO_WFs_above_surface_WFs_flag = true;
	        continue;
	      }

	      if (current->to_bottom) {
	        continue;
	      }

      struct wetting_front *next_to_use = current->next;
      while (next_to_use != NULL && !next_to_use->is_WF_GW) {
        next_to_use = next_to_use->next;
      }

      if (next_to_use == NULL) {
        continue;
      }

      double delta_depth = current->dzdt_cm_per_h * timestep_h;
      double dzdt_cap_boundary_cm = column_depth;
      double allowed_boundary_overshoot_cm = 0.001 * column_depth;
      if (current->layer_num >= 1 && current->layer_num <= num_layers) {
        dzdt_cap_boundary_cm = fmin(column_depth, cum_layer_thickness_cm[current->layer_num]);
      }
      if (current->layer_num >= 1 && current->layer_num < num_layers) {
        const double next_layer_thickness_cm =
          cum_layer_thickness_cm[current->layer_num + 1] - cum_layer_thickness_cm[current->layer_num];
        allowed_boundary_overshoot_cm = 0.1 * next_layer_thickness_cm;
      }
      const double max_TO_dzdt_depth_cm = dzdt_cap_boundary_cm + allowed_boundary_overshoot_cm;
      if (timestep_h > 0.0 && current->dzdt_cm_per_h > 0.0 &&
          current->depth_cm + delta_depth > max_TO_dzdt_depth_cm) {
        const double old_dzdt_cm_per_h = current->dzdt_cm_per_h;
        delta_depth = fmax(0.0, max_TO_dzdt_depth_cm - current->depth_cm);
        current->dzdt_cm_per_h = delta_depth / timestep_h;

        if (verbosity.compare("high") == 0) {
          printf("Capping TO dzdt depth update for front %d from %.17lf cm/h to %.17lf cm/h "
                 "to limit assigned-boundary overshoot "
                 "(old_projected_depth_cm=%.17lf capped_projected_depth_cm=%.17lf "
                 "lower_boundary_cm=%.17lf allowed_overshoot_cm=%.17lf).\n",
                 current->front_num, old_dzdt_cm_per_h, current->dzdt_cm_per_h,
                 current->depth_cm + old_dzdt_cm_per_h * timestep_h,
                 current->depth_cm + delta_depth, dzdt_cap_boundary_cm,
                 allowed_boundary_overshoot_cm);
        }
      }

      struct wetting_front *next = current->next;
      const int layer_num = current->layer_num;
      const double upper_layer_boundary_cm =
        layer_num > 1 ? cum_layer_thickness_cm[layer_num - 1] : 0.0;
      double projected_depth_cm = current->depth_cm + delta_depth;

      if (timestep_h > 0.0 && current->dzdt_cm_per_h < 0.0) {
        double limiting_surface_depth_cm = -1.0;
        for (struct wetting_front *candidate = *head; candidate != NULL; candidate = candidate->next) {
          if (candidate->is_WF_GW || candidate->to_bottom ||
              candidate->layer_num != current->layer_num) {
            continue;
          }

          if (candidate->depth_cm <= current->depth_cm &&
              candidate->depth_cm > current->depth_cm + delta_depth) {
            limiting_surface_depth_cm = fmax(limiting_surface_depth_cm, candidate->depth_cm);
          }
        }

        if (limiting_surface_depth_cm >= 0.0) {
          const double old_dzdt_cm_per_h = current->dzdt_cm_per_h;
          const double old_projected_depth_cm = current->depth_cm + delta_depth;
          const double capped_projected_depth_cm =
            fmax(upper_layer_boundary_cm,
                 limiting_surface_depth_cm - DEPTH_AVOIDS_SAME_WF_DEPTH);
          /* Land just above the active surface front instead of exactly on it.
             Exact contact can be invisible to the strict > merge/correction
             checks, while a tiny event overshoot lets the existing correction
             loop resolve the surface/TO contact before the substep ends. */
          delta_depth = capped_projected_depth_cm - current->depth_cm;
          current->dzdt_cm_per_h = delta_depth / timestep_h;
          projected_depth_cm = current->depth_cm + delta_depth;

          if (verbosity.compare("high") == 0) {
            printf("Capping upward TO dzdt depth update for front %d from %.17lf cm/h to %.17lf cm/h "
                   "to land just above an active surface front "
                   "(old_projected_depth_cm=%.17lf capped_projected_depth_cm=%.17lf "
                   "surface_depth_cm=%.17lf layer=%d).\n",
                   current->front_num, old_dzdt_cm_per_h, current->dzdt_cm_per_h,
                   old_projected_depth_cm, current->depth_cm + delta_depth,
                   limiting_surface_depth_cm, current->layer_num);
          }
        }
      }

      // Use the bounded top-layer TO motion for both the depth update and flux bookkeeping.
      if (timestep_h > 0.0 && current->dzdt_cm_per_h < 0.0 && layer_num == 1 &&
          projected_depth_cm < upper_layer_boundary_cm) {
        const double old_dzdt_cm_per_h = current->dzdt_cm_per_h;
        const double old_projected_depth_cm = projected_depth_cm;
        delta_depth = upper_layer_boundary_cm - current->depth_cm;
        current->dzdt_cm_per_h = delta_depth / timestep_h;
        projected_depth_cm = current->depth_cm + delta_depth;

        if (verbosity.compare("high") == 0) {
          printf("Capping upward TO dzdt depth update for front %d from %.17lf cm/h to %.17lf cm/h "
                 "at the soil surface "
                 "(old_projected_depth_cm=%.17lf capped_projected_depth_cm=%.17lf).\n",
                 current->front_num, old_dzdt_cm_per_h, current->dzdt_cm_per_h,
                 old_projected_depth_cm, projected_depth_cm);
        }
      }

      if (timestep_h > 0.0 && current->dzdt_cm_per_h < 0.0 && layer_num > 1 &&
          current->depth_cm >= upper_layer_boundary_cm &&
          projected_depth_cm < upper_layer_boundary_cm && next != NULL) {
        const int soil_num = soil_type[layer_num];
        const double theta_e = soil_properties[soil_num].theta_e;
        const double theta_r = soil_properties[soil_num].theta_r;
        const double vg_a = soil_properties[soil_num].vg_alpha_per_cm;
        const double vg_m = soil_properties[soil_num].vg_m;
        const double vg_n = soil_properties[soil_num].vg_n;

        const int soil_num_above = soil_type[layer_num - 1];
        const double above_theta_e = soil_properties[soil_num_above].theta_e;
        const double above_theta_r = soil_properties[soil_num_above].theta_r;
        const double above_vg_a = soil_properties[soil_num_above].vg_alpha_per_cm;
        const double above_vg_m = soil_properties[soil_num_above].vg_m;
        const double above_vg_n = soil_properties[soil_num_above].vg_n;

        const double current_theta = fmin(theta_e, current->theta);
        const double Se = calc_Se_from_theta(current->theta, theta_e, theta_r);
        const double current_psi_cm = calc_h_from_Se(Se, vg_a, vg_m, vg_n);
        const double theta_new =
          calc_theta_from_h(current_psi_cm, above_vg_a, above_vg_m, above_vg_n, above_theta_e,
                            above_theta_r);
        const double new_boundary_theta =
          calc_theta_from_h(next->psi_cm, above_vg_a, above_vg_m, above_vg_n, above_theta_e,
                            above_theta_r);
        const double theta_denominator = theta_new - new_boundary_theta;
        const double remap_ratio =
          (next->theta - current_theta) / theta_denominator;
        const double projected_overshoot_cm = upper_layer_boundary_cm - projected_depth_cm;
        const double projected_remap_depth_cm =
          cum_layer_thickness_cm[layer_num - 1] + projected_overshoot_cm * remap_ratio;
        const double max_remap_depth_cm = cum_layer_thickness_cm[layer_num - 1];
        const int remap_layer_num = layer_num - 1;
        const double layer_min_remap_depth_cm = cum_layer_thickness_cm[layer_num - 2];
        double min_remap_depth_cm = layer_min_remap_depth_cm;

        for (struct wetting_front *candidate = *head; candidate != NULL; candidate = candidate->next) {
          if (candidate->is_WF_GW || candidate->to_bottom ||
              candidate->layer_num != remap_layer_num ||
              candidate->depth_cm < layer_min_remap_depth_cm ||
              candidate->depth_cm > max_remap_depth_cm) {
            continue;
          }

          min_remap_depth_cm = fmax(min_remap_depth_cm, candidate->depth_cm);
        }

        if (!std::isfinite(projected_remap_depth_cm) ||
            projected_remap_depth_cm < min_remap_depth_cm ||
            projected_remap_depth_cm > max_remap_depth_cm) {
          double capped_overshoot_cm = 0.0;
          if (std::isfinite(remap_ratio) && remap_ratio < 0.0) {
            capped_overshoot_cm = (min_remap_depth_cm - max_remap_depth_cm) / remap_ratio;
          }

          capped_overshoot_cm = fmax(0.0, fmin(capped_overshoot_cm, projected_overshoot_cm));
          const double capped_projected_depth_cm = upper_layer_boundary_cm - capped_overshoot_cm;
          const double old_dzdt_cm_per_h = current->dzdt_cm_per_h;
          delta_depth = capped_projected_depth_cm - current->depth_cm;
          current->dzdt_cm_per_h = delta_depth / timestep_h;

          if (verbosity.compare("high") == 0) {
            char projected_remap_depth_text[64];
            char remap_ratio_text[64];
            if (std::isfinite(projected_remap_depth_cm)) {
              snprintf(projected_remap_depth_text, sizeof(projected_remap_depth_text),
                       "%.17lf", projected_remap_depth_cm);
            }
            else {
              snprintf(projected_remap_depth_text, sizeof(projected_remap_depth_text),
                       "undefined");
            }
            if (std::isfinite(remap_ratio)) {
              snprintf(remap_ratio_text, sizeof(remap_ratio_text), "%.17lf", remap_ratio);
            }
            else {
              snprintf(remap_ratio_text, sizeof(remap_ratio_text), "undefined");
            }

            printf("Capping upward TO dzdt depth update for front %d from %.17lf cm/h to %.17lf cm/h "
                   "to keep moving-up boundary remap inside the layer above "
                   "(old_projected_depth_cm=%.17lf capped_projected_depth_cm=%.17lf "
                   "projected_remap_depth_cm=%s min_remap_depth_cm=%.17lf "
                   "max_remap_depth_cm=%.17lf remap_ratio=%s remap_layer=%d).\n",
                   current->front_num, old_dzdt_cm_per_h, current->dzdt_cm_per_h,
                   projected_depth_cm, capped_projected_depth_cm, projected_remap_depth_text,
                   min_remap_depth_cm, max_remap_depth_cm, remap_ratio_text, remap_layer_num);
          }
        }
      }

      double delta_theta = 0.0;

      if (current->layer_num == next_to_use->layer_num) {
        delta_theta = next_to_use->theta - current->theta;
      }
      else {
        const int soil_num_current = soil_type[current->layer_num];
        const double theta_e_current = soil_properties[soil_num_current].theta_e;
        const double theta_r_current = soil_properties[soil_num_current].theta_r;
        const double vg_a_current = soil_properties[soil_num_current].vg_alpha_per_cm;
        const double vg_m_current = soil_properties[soil_num_current].vg_m;
        const double vg_n_current = soil_properties[soil_num_current].vg_n;
        const double equiv_next_theta =
          calc_theta_from_h(next_to_use->psi_cm, vg_a_current, vg_m_current, vg_n_current,
                            theta_e_current, theta_r_current);
        delta_theta = equiv_next_theta - current->theta;
      }

	      if (!TO_WFs_above_surface_WFs_flag) {
	        current->depth_cm += delta_depth;
	      }
	      else {
	        bottom_boundary_flux_above_surface_WFs_cm += delta_depth * delta_theta;
	      }
	      bottom_boundary_flux_cm += delta_depth * delta_theta;
	    }

	    if (verbosity.compare("high") == 0) {
	      printf("State after TO WF depth update via dzdt...\n");
	      listPrint(*head);
	      printf("Bottom boundary flux after TO WF movement = %lf \n", bottom_boundary_flux_cm);
	    }

	    if (listLength_surface(*head) > 0) {
	      bottom_boundary_flux_cm = lgarto_extract_TO_GW_flux_from_surface_WFs(
	        &bottom_boundary_flux_above_surface_WFs_cm, bottom_boundary_flux_cm, AET_demand_cm,
	        cum_layer_thickness_cm, soil_type, soil_properties, head);
	      const double target_mass_after_fluxes_cm =
	        (old_mass + precip_mass_to_add) -
	        (*AET_demand_cm + free_drainage_demand +
	         cached_lower_boundary_flux_correction_cm + bottom_boundary_flux_cm);
	      lgar_global_theta_update(bottom_boundary_flux_above_surface_WFs_cm,
	                               target_mass_after_fluxes_cm,
	                               cum_layer_thickness_cm, soil_type, soil_properties,
	                               head, lgar_global_theta_snap_mass_tolerance_cm);
	      lgar_global_psi_update(soil_type, soil_properties, head);
	    }

	    mass_change = 0.0;
	    if (lgar_check_dry_over_wet_wetting_fronts(*head)) {
	      lgar_fix_dry_over_wet_wetting_fronts(&mass_change, cum_layer_thickness_cm, soil_type, head, soil_properties);
	      (void) add_AET_with_pet_budget(-mass_change,
	                                     "TO dry-over-wet mass residual",
	                                     true);
	    }

    if (listLength_surface(*head) == 0 && PET_timestep_cm <= 1.0e-12) {
      lgarto_ensure_no_surface_to_psi_span(num_layers, cum_layer_thickness_cm, soil_type,
                                           frozen_factor, soil_properties, head);
    }

    int surface_front_count = listLength_surface(*head);
    double deepest_surface_depth_cm = 0.0;
    if (surface_front_count > 0) {
      for (struct wetting_front *front = *head; front != NULL; front = front->next) {
        if (front->is_WF_GW == FALSE) {
          deepest_surface_depth_cm = fmax(deepest_surface_depth_cm, front->depth_cm);
        }
      }
    }

    // With active surface fronts, refine only the TO/GW profile beneath the deepest
    // surface front. The surface-supported bridge and zero-depth TO metadata stay out
    // of this resolution pass.
    const int psi_gap_refinement_count =
      lgarto_refine_large_to_psi_gaps(num_layers, cum_layer_thickness_cm, soil_type,
                                      frozen_factor, soil_properties, head,
                                      surface_front_count == 0 ?
                                        NO_SURFACE_TO_PSI_GAP_REFINEMENT_MAX_INSERTIONS :
                                        SURFACE_PRESENT_TO_PSI_GAP_REFINEMENT_MAX_INSERTIONS,
                                      surface_front_count == 0 ?
                                        0.0 : deepest_surface_depth_cm + DEPTH_AVOIDS_SAME_WF_DEPTH);
    if (psi_gap_refinement_count > 0 && verbosity.compare("high") == 0) {
      printf("State after %d TO psi-gap refinement iteration(s) and before TO/general correction loop...\n",
             psi_gap_refinement_count);
      listPrint(*head);
    }

    bool merged_in_non_top_layer = false;
    int correction_type =
      lgarto_correction_type(num_layers, cum_layer_thickness_cm, head,
                             surface_lower_boundary_depth_cm);

    if (verbosity.compare("high") == 0) {
      printf("correction_type before TO/general correction loop: %d \n", correction_type);
      printf("mass before TO/general correction loop: %lf \n", lgar_calc_mass_bal(cum_layer_thickness_cm, *head));
    }

    double net_last_layer_GW_overshoot_truncation_cm = 0.0;
    while (correction_type != 0) {
      if (verbosity.compare("high") == 0) {
        printf("correction_type at start of TO/general correction iteration: %d \n", correction_type);
      }

      if (correction_type == 1) {
        const double mass_before_surface_TO_merge = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
        merged_in_non_top_layer =
          lgar_merge_surface_and_TO_wetting_fronts(merged_in_non_top_layer, num_layers,
                                                   cum_layer_thickness_cm, soil_type,
                                                   frozen_factor, soil_properties, head);
        const bool did_a_WF_have_negative_depth = lgarto_correct_negative_depths(head);
        lgarto_cleanup_after_surface_TO_merging_in_layer_below_top(merged_in_non_top_layer, soil_type,
                                                                   soil_properties, head);
        const bool close_psis = correct_close_psis(soil_type, soil_properties, head);

        double mass_balance_flux_correction_cm =
          mass_before_surface_TO_merge - lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
        int boundary_pinned_surface_front_num = -1;
        if (fabs(mass_balance_flux_correction_cm) > MBAL_ITERATIVE_TOLERANCE) {
          mass_balance_flux_correction_cm =
            lgarto_restore_surface_TO_merge_mass_via_surface_depths(mass_before_surface_TO_merge,
                                                                    num_layers,
                                                                    cum_layer_thickness_cm,
                                                                    *head,
                                                                    &boundary_pinned_surface_front_num);
        }
        if (mass_balance_flux_correction_cm > MBAL_ITERATIVE_TOLERANCE &&
            boundary_pinned_surface_front_num > -1) {
          mass_balance_flux_correction_cm =
            lgarto_restore_surface_TO_merge_mass_via_boundary_psi(mass_before_surface_TO_merge,
                                                                  boundary_pinned_surface_front_num,
                                                                  num_layers,
                                                                  cum_layer_thickness_cm,
                                                                  soil_type,
                                                                  frozen_factor,
                                                                  soil_properties,
                                                                  *head);
        }
        if (did_a_WF_have_negative_depth || close_psis ||
            fabs(mass_balance_flux_correction_cm) > MBAL_ITERATIVE_TOLERANCE) {
          lgarto_assert_gw_flux_mass_balance_correction_within_debug_threshold(
            mass_balance_flux_correction_cm,
            "lgar_move_wetting_fronts",
            "surface/TO merge cleanup mass residual",
            *head);
          bottom_boundary_flux_cm += mass_balance_flux_correction_cm;
        }
      }

      if (correction_type == 2) {
        const double target_mass_after_fluxes =
          (old_mass + precip_mass_to_add) -
	          (*AET_demand_cm + free_drainage_demand + cached_lower_boundary_flux_correction_cm +
           bottom_boundary_flux_cm);
        const double mass_balance_flux_correction_cm =
          lgarto_TO_WFs_merge_via_depth(target_mass_after_fluxes, column_depth, cum_layer_thickness_cm, head,
                                        soil_type, soil_properties);
        lgarto_assert_gw_flux_mass_balance_correction_within_debug_threshold(
          mass_balance_flux_correction_cm,
          "lgarto_TO_WFs_merge_via_depth",
          "TO depth merge residual",
          *head);
        bottom_boundary_flux_cm += mass_balance_flux_correction_cm;
        lgar_global_psi_update(soil_type, soil_properties, head);
      }

      if (correction_type == 3) {
        const double target_mass_after_fluxes =
          (old_mass + precip_mass_to_add) -
	          (*AET_demand_cm + free_drainage_demand + cached_lower_boundary_flux_correction_cm +
           bottom_boundary_flux_cm);
        const double mass_balance_flux_correction_cm =
          lgarto_TO_WFs_merge_via_theta(target_mass_after_fluxes, column_depth, cum_layer_thickness_cm, head,
                                        soil_type, soil_properties);
        lgarto_assert_gw_flux_mass_balance_correction_within_debug_threshold(
          mass_balance_flux_correction_cm,
          "lgarto_TO_WFs_merge_via_theta",
          "TO theta merge residual",
          *head);
        bottom_boundary_flux_cm += mass_balance_flux_correction_cm;
        lgar_global_psi_update(soil_type, soil_properties, head);
      }

      if (correction_type == 4) {
        int front_num_with_negative_depth = -1;
        double mass_balance_flux_correction_cm = 0.0;
        lgar_TO_wetting_fronts_cross_layer_boundary(&front_num_with_negative_depth, num_layers,
                                                    cum_layer_thickness_cm, soil_type, frozen_factor, head,
                                                    soil_properties,
                                                    &mass_balance_flux_correction_cm);
        lgar_global_psi_update(soil_type, soil_properties, head);
        if (std::fabs(mass_balance_flux_correction_cm) > MBAL_ITERATIVE_TOLERANCE) {
          lgarto_assert_gw_flux_mass_balance_correction_within_debug_threshold(
            mass_balance_flux_correction_cm,
            "lgar_TO_wetting_fronts_cross_layer_boundary",
            "TO layer-boundary remap residual",
            *head);
          bottom_boundary_flux_cm += mass_balance_flux_correction_cm;
        }

        // if (front_num_with_negative_depth > -1) {
        //   const double mass_before_negative_depth_delete = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
        //   listDeleteFront(front_num_with_negative_depth, head, soil_type, soil_properties);
        //   const double mass_after_negative_depth_delete = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
        //   bottom_boundary_flux_cm +=
        //     (mass_before_negative_depth_delete - mass_after_negative_depth_delete);
        // }
      }

      if (correction_type == 5) {
        lgar_merge_wetting_fronts(soil_type, frozen_factor, head, soil_properties);
      }

      if (correction_type == 6) {
        const double mass_before_surface_bdy_crossing = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
        const bool close_psis = correct_close_psis(soil_type, soil_properties, head);
        if (close_psis) {
          const double mass_balance_flux_correction_cm =
            mass_before_surface_bdy_crossing - lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
          lgarto_assert_gw_flux_mass_balance_correction_within_debug_threshold(
            mass_balance_flux_correction_cm,
            "lgar_move_wetting_fronts",
            "surface boundary-crossing close-psi residual",
            *head);
          bottom_boundary_flux_cm += mass_balance_flux_correction_cm;
        }
        const double mass_before_surface_bdy_crossing_remap = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
        lgar_wetting_fronts_cross_layer_boundary(num_layers, cum_layer_thickness_cm, soil_type, frozen_factor, head,
                                                 soil_properties);
        const double mass_balance_flux_correction_cm =
          mass_before_surface_bdy_crossing_remap - lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
        if (fabs(mass_balance_flux_correction_cm) > MBAL_ITERATIVE_TOLERANCE) {
          lgarto_assert_gw_flux_mass_balance_correction_within_debug_threshold(
            mass_balance_flux_correction_cm,
            "lgar_move_wetting_fronts",
            "surface boundary-crossing remap residual",
            *head);
          bottom_boundary_flux_cm += mass_balance_flux_correction_cm;
        }
      }

      if (correction_type == 7) {
        bottom_boundary_flux_cm +=
          lgar_wetting_front_cross_domain_boundary(surface_lower_boundary_depth_cm, soil_type, frozen_factor,
                                                   head, soil_properties);
        if (isnan(bottom_boundary_flux_cm)) {
          bottom_boundary_flux_cm = 0.0;
        }
      }

      if (correction_type == 8) {
        lgarto_convert_surface_fronts_drier_than_TO_below(head, "correction type 8");
      }

      if (correction_type == 9) {
        const double mass_balance_flux_correction_cm =
          lgarto_truncate_last_layer_GW_overshoot(cum_layer_thickness_cm, num_layers,
                                                  head, soil_type, soil_properties);
        net_last_layer_GW_overshoot_truncation_cm += mass_balance_flux_correction_cm;
        bottom_boundary_flux_cm += mass_balance_flux_correction_cm;
      }

      if (correction_type == 10) {
        const double mass_balance_flux_correction_cm =
          lgarto_surface_fronts_cross_layer_boundary_upward(
            num_layers, cum_layer_thickness_cm, soil_type, frozen_factor,
            soil_properties, head, "TO/general correction loop");
        if (fabs(mass_balance_flux_correction_cm) > MBAL_ITERATIVE_TOLERANCE) {
          lgarto_assert_gw_flux_mass_balance_correction_within_debug_threshold(
            mass_balance_flux_correction_cm,
            "lgar_move_wetting_fronts",
            "general surface upward layer-crossing residual",
            *head);
          bottom_boundary_flux_cm += mass_balance_flux_correction_cm;
        }
      }

      correction_type =
        lgarto_correction_type(num_layers, cum_layer_thickness_cm, head,
                               surface_lower_boundary_depth_cm);
      if (verbosity.compare("high") == 0) {
        printf("correction_type at end of TO/general correction iteration: %d \n", correction_type);
      }
    }

    lgarto_assert_gw_flux_mass_balance_correction_within_debug_threshold(
      net_last_layer_GW_overshoot_truncation_cm,
      "lgarto_truncate_last_layer_GW_overshoot",
      "net last-layer TO/GW overshoot truncation across correction loop",
      *head);

    if (verbosity.compare("high") == 0) {
      printf("mass after TO/general correction loop: %lf \n", lgar_calc_mass_bal(cum_layer_thickness_cm, *head));
      printf("bottom_boundary_flux_cm after TO/general correction loop = %lf \n", bottom_boundary_flux_cm);
      printf("bottom_boundary_flux_cm_from_surf_WFs retained as %lf \n", bottom_boundary_flux_cm_from_surf_WFs);
    }
  }

  /***********************************************/
  // make sure all psi values are updated
  // there is a rare error where, after a wetting front crosses a layer boundary to a soil layer with a much more sensitive soil water retention curve, and the wetting front is near but not at saturation,
  // the barely unsatruated wetting front in the new layer will update its psi value as 0. This causes unequal psi values among adjacent wetting fronts in different layers and then a mass balance error.
  // For example, consider the following: (1.0/(pow(1.0+pow(0.001039*0.00296620706633,2.938410),0.659680))*(0.618625-0.052623)+0.052623). This expression, using specific values in calc_theta_from_h, will yield a theta of theta_e (0.618625), while having a psi value that is slightly greater than 0, which is 0.00296620706633 in this case.
  // While that is ok for this soil layer in particular, adjacent wetting fronts above this one with a less sensitive soil water retention curve will yield a non-theta_e value for the psi value that is slightly above 0.
  // A solution is to either not run the following code, or to not run it when the wetting front is very close to saturation with a very sensitive soil water retention curve. Adding code that runs the following only for psi>1. 

  current = *head;

  for (int wf=1; wf != listLength(*head) + 1; wf++) { // shifted loop bound to listLength(*head) + 1 so that the very bottom WF can have its K_cm_per_h updated, important for free drainage. 

    int soil_num_k    = soil_type[current->layer_num];

    double theta_e_k   = soil_properties[soil_num_k].theta_e;
    double theta_r_k   = soil_properties[soil_num_k].theta_r;
    double vg_a_k      = soil_properties[soil_num_k].vg_alpha_per_cm;
    double vg_m_k      = soil_properties[soil_num_k].vg_m;
    double vg_n_k      = soil_properties[soil_num_k].vg_n;

    double Ksat_cm_per_h_k  = frozen_factor[current->layer_num] * soil_properties[soil_num_k].Ksat_cm_per_h;

    double Se = calc_Se_from_theta(current->theta,theta_e_k,theta_r_k);
    // if (current->psi_cm>1.0){
      current->psi_cm = calc_h_from_Se(Se, vg_a_k, vg_m_k, vg_n_k); 
    // }

    current->K_cm_per_h = calc_K_from_Se(Se, Ksat_cm_per_h_k, vg_m_k);

    current = current->next;

  }

  for (current = *head; current != NULL && current->next != NULL; current = current->next) {
    struct wetting_front *next = current->next;
    if (!current->to_bottom || current->layer_num == next->layer_num) {
      continue;
    }

    if (!lgar_should_snap_boundary_psi_to_lower(current, next, soil_type, soil_properties)) {
      continue;
    }

    current->psi_cm = next->psi_cm;
    const int soil_num_k = soil_type[current->layer_num];
    const double theta_e_k = soil_properties[soil_num_k].theta_e;
    const double theta_r_k = soil_properties[soil_num_k].theta_r;
    const double vg_a_k = soil_properties[soil_num_k].vg_alpha_per_cm;
    const double vg_m_k = soil_properties[soil_num_k].vg_m;
    const double vg_n_k = soil_properties[soil_num_k].vg_n;
    const double Ksat_cm_per_h_k =
      frozen_factor[current->layer_num] * soil_properties[soil_num_k].Ksat_cm_per_h;

    current->theta = calc_theta_from_h(current->psi_cm, vg_a_k, vg_m_k, vg_n_k, theta_e_k, theta_r_k);
    const double Se = calc_Se_from_theta(current->theta, theta_e_k, theta_r_k);
    current->K_cm_per_h = calc_K_from_Se(Se, Ksat_cm_per_h_k, vg_m_k);
  }

  lgarto_cap_zero_depth_TO_psi(num_layers, soil_type, frozen_factor, soil_properties, *head);

  bottom_boundary_flux_cm +=
    lgarto_clip_final_layer_GW_overshoot_to_vadose_boundary(cum_layer_thickness_cm,
                                                            num_layers,
                                                            head);


  if (verbosity.compare("high") == 0){
    printf("Moving/merging wetting fronts done... \n");
    listPrint(*head);
  }

  //Just a check to make sure that, when there is only 1 layer, than the existing wetting front is at the correct depth.
  //This might have been fixed with other debugging related to scenarios with just 1 layer where the wetting front is completely satruated. Not sure this is necessary.
	  if (listLength(*head)==1) {
	    current = *head;
	    if (current->depth_cm != cum_layer_thickness_cm[1]) {
	      current->depth_cm = cum_layer_thickness_cm[1];
	    }
	  }

	  const double final_mass_cm = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
	  const double storage_implied_AET_cm =
	    fmax(0.0, old_mass + precip_mass_to_add -
	                (free_drainage_demand + cached_lower_boundary_flux_correction_cm +
	                 bottom_boundary_flux_cm + final_mass_cm));
	  const double applied_AET_bookkeeping_cm =
	    fmin(accepted_positive_AET_bookkeeping_cm, storage_implied_AET_cm);
		  if (use_TO_surface_AET &&
		      fabs(*AET_demand_cm - applied_AET_bookkeeping_cm) > MBAL_ITERATIVE_TOLERANCE) {
	    if (verbosity.compare("high") == 0) {
	      printf("Reconciling LGARTO AET bookkeeping from %.17lf cm to %.17lf cm "
	             "after profile storage and fluxes implied a different accepted AET amount "
	             "(accepted_positive_AET_bookkeeping_cm=%.17lf, storage_implied_AET_cm=%.17lf).\n",
	             *AET_demand_cm, applied_AET_bookkeeping_cm,
	             accepted_positive_AET_bookkeeping_cm, storage_implied_AET_cm);
	    }
		    *AET_demand_cm = applied_AET_bookkeeping_cm;
		  }

			  return(bottom_boundary_flux_cm);

	}


// ############################################################################################
/*
  the function merges wetting fronts; called from lgar_move_wetting_fronts.
*/
// ############################################################################################

extern void lgar_merge_wetting_fronts(int *soil_type, double *frozen_factor, struct wetting_front** head,
				      struct soil_properties_ *soil_properties)
{
  
  struct wetting_front *current;
  struct wetting_front *next;
  struct wetting_front *next_to_next;
  current = *head;

  if (verbosity.compare("high") == 0) {
    printf("State before merging wetting fronts...\n");
    listPrint(*head);
    printf("Merging wetting fronts... \n");
  }

  // local variables
  double theta_e,theta_r;
  double vg_a, vg_m, vg_n;
  double Se, Ksat_cm_per_h;
  int layer_num, soil_num;
    
  for (int wf=1; wf != listLength(*head); wf++) {
    
    if (verbosity.compare("high") == 0) {
      printf("Merge | ********* Wetting Front = %d *********\n", wf);
    }


    next = current->next;
    next_to_next = current->next->next;

    // case : wetting front passing another wetting front within a layer
    /**********************************************************/
	    // 'current->depth_cm > next->depth_cm' ensures that merging is needed
	    // 'current->layer_num == next->layer_num' ensures wetting fronts are in the same layer
	    // '!next->to_bottom' ensures that the next wetting front is not the deepest wetting front in the layer
	    if ( (current->depth_cm > next->depth_cm) && (current->layer_num == next->layer_num) && !next->to_bottom) {
        if (lgarto_surface_front_overtook_surface_front_above_TO_chain(current)) {
          (void) lgarto_convert_overtaken_surface_front_above_TO_chain(
            head, "mixed surface/surface/TO merge guard");
          break;
        }

	      double current_mass_this_layer = current->depth_cm * (current->theta - next->theta) + next->depth_cm*(next->theta - next_to_next->theta);
	      current->depth_cm = current_mass_this_layer / (current->theta - next_to_next->theta);

	      // assert (current->depth_cm > 0.0); //in theory only required for LGAR only mode

      layer_num = current->layer_num;
      soil_num  = soil_type[layer_num];
      theta_e   = soil_properties[soil_num].theta_e;
      theta_r   = soil_properties[soil_num].theta_r;
      vg_a      = soil_properties[soil_num].vg_alpha_per_cm;
      vg_m      = soil_properties[soil_num].vg_m;
      vg_n      = soil_properties[soil_num].vg_n;
      Se        = calc_Se_from_theta(current->theta,theta_e,theta_r);

      Ksat_cm_per_h  = soil_properties[soil_num].Ksat_cm_per_h * frozen_factor[current->layer_num];

      current->psi_cm     = calc_h_from_Se(Se, vg_a, vg_m, vg_n);
      current->K_cm_per_h = calc_K_from_Se(Se, Ksat_cm_per_h, vg_m);
      
      if (verbosity.compare("high") == 0) {
        printf ("Deleting wetting front (before)... \n");
        listPrint(*head);
      }
      
      next = listDeleteFront(next->front_num, head, soil_type, soil_properties);;
      
      if (verbosity.compare("high") == 0) {
        printf ("Deleting wetting front (after) ... \n");
        listPrint(*head);
      }
    }
    
    current = current->next;
  }

  if (verbosity.compare("high") == 0) {
    printf("State after merging wetting fronts...\n");
    listPrint(*head);
  }

}



// ############################################################################################
/*
  the function lets wetting fronts of a sufficient depth cross layer boundaries; called from lgar_move_wetting_fronts.
*/
// ############################################################################################

extern void lgar_wetting_fronts_cross_layer_boundary(int num_layers,
						     double* cum_layer_thickness_cm, int *soil_type,
						     double *frozen_factor, struct wetting_front** head,
						     struct soil_properties_ *soil_properties)
{
  struct wetting_front *current;
  struct wetting_front *next;
  struct wetting_front *next_to_next;
  current = *head; 
  bool cross_necessary = false;
  bool theta_correction_necessary = false;
  int front_for_cross = -1;
  bool lgarto_active = false;

  for (struct wetting_front *front = *head; front != NULL; front = front->next) {
    if (front->is_WF_GW) {
      lgarto_active = true;
      break;
    }
  }

  if (verbosity.compare("high") == 0) {
    printf("Layer boundary crossing... \n");
    printf("States before wetting fronts cross layer boundary...\n");
    listPrint(*head);
  }

  for (int wf=1; wf != listLength(*head); wf++) {
    
    // local variables
    double theta_e,theta_r;
    double vg_a, vg_m, vg_n;
    int layer_num, soil_num;


    layer_num   = current->layer_num;
    soil_num    = soil_type[layer_num];
    theta_e     = soil_properties[soil_num].theta_e;
    theta_r     = soil_properties[soil_num].theta_r;
    vg_a        = soil_properties[soil_num].vg_alpha_per_cm;
    vg_m        = soil_properties[soil_num].vg_m;
    vg_n        = soil_properties[soil_num].vg_n;
    double Ksat_cm_per_h  = soil_properties[soil_num].Ksat_cm_per_h * frozen_factor[current->layer_num]; //PTL addition to make K_cm_per_h for this conditon to be correct
    
    next = current->next;
    next_to_next = current->next->next;

    if (current->depth_cm > cum_layer_thickness_cm[layer_num] && (next->depth_cm == cum_layer_thickness_cm[layer_num]) && (next->to_bottom)
	&& (layer_num!=num_layers) ) {

      if (verbosity.compare("high") == 0) {
        printf("Boundary Crossing | ******* Wetting Front = %d ****** \n", wf);
      }

      double current_theta = fmin(theta_e, current->theta);
      double overshot_depth = current->depth_cm - next->depth_cm;
      int soil_num_next = soil_type[layer_num+1];

      double next_theta_e   = soil_properties[soil_num_next].theta_e;
      double next_theta_r   = soil_properties[soil_num_next].theta_r;
      double next_vg_a      = soil_properties[soil_num_next].vg_alpha_per_cm;
      double next_vg_m      = soil_properties[soil_num_next].vg_m;
      double next_vg_n      = soil_properties[soil_num_next].vg_n;
      //double next_Ksat_cm_per_h  = soil_properties[soil_num_next].Ksat_cm_per_h * frozen_factor[current->layer_num]; 

      double Se = calc_Se_from_theta(current->theta,theta_e, theta_r);
      current->psi_cm = calc_h_from_Se(Se, vg_a, vg_m, vg_n);

      current->K_cm_per_h = calc_K_from_Se(Se, Ksat_cm_per_h, vg_m);
      
      // current psi with van Gunechten properties of the next layer to get new theta
      double theta_new = calc_theta_from_h(current->psi_cm, next_vg_a, next_vg_m, next_vg_n, next_theta_e, next_theta_r);

      double mbal_correction = overshot_depth * (current_theta - next->theta);
      double mbal_Z_correction = mbal_correction / (theta_new - next_to_next->theta); // this is the new wetting front depth
      if (isinf(mbal_Z_correction)){//in some rare cases, due to (very) different shapes of soil water rentention curves in adjacent soil layers near saturation, the WF that crossed a layer bdy can yield a theta value identical to the one below it, resulting in division by 0 and an infinite WF depth. 
        mbal_Z_correction = cum_layer_thickness_cm[layer_num+1]/(1e3);
      }

      double depth_new = cum_layer_thickness_cm[layer_num] + mbal_Z_correction; // this is the new wetting front absolute depth
      if (lgarto_active) {
        const double max_surface_remap_depth_cm = 1.10 * cum_layer_thickness_cm[num_layers];
        if (!std::isfinite(depth_new) || depth_new > max_surface_remap_depth_cm) {
          if (verbosity.compare("high") == 0) {
            printf("Domain-limited LGARTO surface layer-boundary remap for front %d "
                   "from depth_new %.17lf cm to %.17lf cm "
                   "(110%% of column depth %.17lf cm).\n",
                   current->front_num,
                   depth_new,
                   max_surface_remap_depth_cm,
                   cum_layer_thickness_cm[num_layers]);
          }
          depth_new = max_surface_remap_depth_cm;
        }
      }

      // A previous attempt capped ordinary surface layer-boundary remaps at the
      // first TO front in the target layer. That was too broad: synth 1 needs
      // some of these remaps to proceed so that recharge stays in family.
      // If the resulting surface/TO drier merge would create an invalid local
      // remap, lgar_merge_surface_and_TO_wetting_fronts now limits that merge
      // directly.
      //
      // double limiting_TO_depth_cm = HUGE_VAL;
      // for (struct wetting_front *candidate = *head; candidate != NULL; candidate = candidate->next) {
      //   if (candidate == current || candidate == next || !candidate->is_WF_GW ||
      //       candidate->to_bottom || candidate->layer_num != layer_num + 1) {
      //     continue;
      //   }
      //
      //   limiting_TO_depth_cm = fmin(limiting_TO_depth_cm, candidate->depth_cm);
      // }
      //
      // if (std::isfinite(limiting_TO_depth_cm) && depth_new > limiting_TO_depth_cm) {
      //   // Old mass-conservative layer-crossing remap can jump a surface front
      //   // far past an existing TO front when theta_new is very close to the
      //   // lower front's theta:
      //   // double depth_new = cum_layer_thickness_cm[layer_num] + mbal_Z_correction;
      //   depth_new = limiting_TO_depth_cm + TRUNCATION_DEPTH;
      // }

      current->depth_cm = cum_layer_thickness_cm[layer_num];
      
      next->theta = theta_new;
      next->psi_cm = current->psi_cm;
      next->depth_cm = depth_new;
      next->layer_num = layer_num + 1;
      next->dzdt_cm_per_h = current->dzdt_cm_per_h;
      current->dzdt_cm_per_h = 0;
      current->to_bottom = TRUE;
      next->to_bottom = FALSE;
      current->is_WF_GW = 0;
      next->is_WF_GW = 0;
      cross_necessary = true;
      front_for_cross = next->front_num;
      if (depth_new>cum_layer_thickness_cm[num_layers] && next->layer_num==num_layers){
        theta_correction_necessary = true;
      }
      
    }
    
    current = current->next;
  }

  if (cross_necessary){
    for (int wf = listLength(*head)-1; wf != 0; wf--) {
      struct wetting_front *current_temp = listFindFront(wf, *head, NULL);
      struct wetting_front *next_temp = current_temp->next;
      if ( (current_temp->to_bottom) ){
        current_temp->psi_cm = next_temp->psi_cm;

        int soil_num_k1 = soil_type[current_temp->layer_num]; 
        double theta_e_k   = soil_properties[soil_num_k1].theta_e;
        double theta_r_k   = soil_properties[soil_num_k1].theta_r;
        double vg_a_k      = soil_properties[soil_num_k1].vg_alpha_per_cm;
        double vg_m_k      = soil_properties[soil_num_k1].vg_m;
        double vg_n_k      = soil_properties[soil_num_k1].vg_n;
        current_temp->theta = calc_theta_from_h(current_temp->psi_cm, vg_a_k, vg_m_k, vg_n_k,theta_e_k,theta_r_k);
        double Ksat_cm_per_h_k  = soil_properties[soil_num_k1].Ksat_cm_per_h;
        double Se = calc_Se_from_theta(current_temp->theta,theta_e_k, theta_r_k);
        current_temp->K_cm_per_h = calc_K_from_Se(Se, Ksat_cm_per_h_k, vg_m_k);
      }
      else{
        int soil_num_k1 = soil_type[current_temp->layer_num]; 
        double theta_e_k   = soil_properties[soil_num_k1].theta_e;
        double theta_r_k   = soil_properties[soil_num_k1].theta_r;
        double vg_m_k      = soil_properties[soil_num_k1].vg_m;
        double Ksat_cm_per_h_k  = soil_properties[soil_num_k1].Ksat_cm_per_h;
        double Se = calc_Se_from_theta(current_temp->theta,theta_e_k, theta_r_k);
        current_temp->K_cm_per_h = calc_K_from_Se(Se, Ksat_cm_per_h_k, vg_m_k);
      }
    }
  if (verbosity.compare("high") == 0) {
    printf("States after wetting fronts cross layer boundary before theta correction...\n");
    listPrint(*head);
  }

    if (theta_correction_necessary){
      for (int wf = listLength(*head)-1; wf != 0; wf--) {
      struct wetting_front *current = listFindFront(wf, *head, NULL);
        double prior_mass = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
        if (current->depth_cm>cum_layer_thickness_cm[num_layers] && current->layer_num==num_layers && current->front_num==front_for_cross){
          current->depth_cm = cum_layer_thickness_cm[num_layers] + 1.E-6;
          int front_num_correction = current->front_num;
          // first, lgar_theta_mass_balance_correction will attempt to close the mass balance by adjusting the theta value of the WF that crossed the layer boundary and other WFs sharing a psi value with it.
          lgar_theta_mass_balance_correction(false, front_num_correction, prior_mass, head, cum_layer_thickness_cm, soil_type, soil_properties); 
            if (verbosity.compare("high") == 0) {
              printf("States after wetting fronts cross layer boundary and after theta correction...\n");
              listPrint(*head);
            }
          double current_mass = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
          if (prior_mass > (current_mass + 100.*MBAL_ITERATIVE_TOLERANCE)) { //the inclusion of 100.*MBAL_ITERATIVE_TOLERANCE is due to the fact that prior_mass > current_mass might be true, but only within the mass balance tolerance, in which case we should not run this
            // if lgar_theta_mass_balance_correction reached theta_e or got very close, then mass balance closure was not possible (theta could not be increased enough so the current mass is too low), which is uncommon but can happen if multiple correction types are necessary in the same time step
            // in this case, the depth that closes the mass balance is searched for by finding two depths for current->depth_cm, one that makes the storage too low, and one that makes it too high, so the answer is somewhere in between 
            double depth_increment = 0.001;                //initial guess
            double max_increment = 1./TRUNCATION_DEPTH;    // Prevent inf depth
            double max_depth = cum_layer_thickness_cm[num_layers] + 1./TRUNCATION_DEPTH; 

            double increment = depth_increment;
            bool found_upper_bound = false;
            int iter_one_direction = 0;

            // Exponential increase to find the depth that leads to a mass that is too large
            while (prior_mass > current_mass && increment < max_increment && current->depth_cm < max_depth) {
              iter_one_direction++;
              if (iter_one_direction > MAX_ITER_MBAL_LOOP){ 
                break;
              }
              current->depth_cm += increment;
              current_mass = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);

              if (current_mass >= prior_mass) {
                  found_upper_bound = true;
                  break;
              }
              increment *= 2.0;  
            }

            // Now that the depth that closes the mass balance will be between these two depths
            if (found_upper_bound) {
              double low = current->depth_cm - increment;
              double high = current->depth_cm;
              double mid;

              double current_mass = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);

              while (fabs(current_mass - prior_mass) > MBAL_ITERATIVE_TOLERANCE) {
                iter_one_direction++;
                mid = 0.5 * (low + high);
                current->depth_cm = mid;
                current_mass = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);

                if (current_mass >= prior_mass) {
                  high = mid;
                } else {
                  low = mid;
                }

                if (iter_one_direction > MAX_ITER_MBAL_LOOP){
                  break;
                }
              }
            }
          }
        }
      }
    } 
  }

  if (verbosity.compare("high") == 0) {
    printf("States after wetting fronts cross layer boundary...\n");
    listPrint(*head);
  }

}

// ############################################################################################
/*
  If TO fronts above surface fronts demand groundwater flux, subtract that mass from the
  right-most surface wetting-front region before the TO correction loop runs.
*/
// ############################################################################################

extern double lgarto_extract_TO_GW_flux_from_surface_WFs(double *bottom_boundary_flux_above_surface_WFs_cm,
							 double bottom_boundary_flux_cm,
							 double *AET_demand_cm,
							 double *cum_layer_thickness_cm,
							 int *soil_type,
							 struct soil_properties_ *soil_properties,
							 struct wetting_front **head)
{
  if (*bottom_boundary_flux_above_surface_WFs_cm == 0.0) {
    return bottom_boundary_flux_cm;
  }

  const double mass_at_start_of_extraction = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
  const double bottom_bdy_flux_at_start_of_extraction = bottom_boundary_flux_cm;
  const int wf_free_drainage_demand = wetting_front_free_drainage(*head);

  if (verbosity.compare("high") == 0) {
    printf("before lgarto_extract_TO_GW_flux_from_surface_WFs \n");
    listPrint(*head);
    printf("associated mass: %lf \n", mass_at_start_of_extraction);
    printf("bottom_boundary_flux_above_surface_WFs_cm: %lf \n", *bottom_boundary_flux_above_surface_WFs_cm);
  }

  if (listLength_surface(*head) > 0 &&
      *bottom_boundary_flux_above_surface_WFs_cm < -MBAL_ITERATIVE_TOLERANCE) {
    // Negative residuals are signed lower-boundary flux corrections, not water
    // to extract from surface fronts.  Keep profile storage unchanged and carry
    // the signed amount in the timestep flux bookkeeping.
    bottom_boundary_flux_cm -= *bottom_boundary_flux_above_surface_WFs_cm;
    if (verbosity.compare("high") == 0) {
      printf("Routed negative TO/GW above-surface residual %.17lf cm to bottom-boundary flux; "
             "updated bottom_boundary_flux_cm = %.17lf.\n",
             *bottom_boundary_flux_above_surface_WFs_cm, bottom_boundary_flux_cm);
    }
    *bottom_boundary_flux_above_surface_WFs_cm = 0.0;
    return bottom_boundary_flux_cm;
  }

  struct wetting_front *current = listFindFront(wf_free_drainage_demand, *head, NULL);
  if (current == NULL) {
    return bottom_boundary_flux_cm;
  }

  const int wf_from_which_to_extract = current->front_num;
  const double original_theta = current->theta;
  const double original_psi_cm = current->psi_cm;
  struct wetting_front *next = current->next;

  if (listLength_surface(*head) > 0) {
    if (current->layer_num == 1) {
      if (current->depth_cm > 0.0) {
        const double theta_reduction = *bottom_boundary_flux_above_surface_WFs_cm / current->depth_cm;
        const double theta_new = current->theta - theta_reduction;
        current->theta = fmax(soil_properties[soil_type[current->layer_num]].theta_r + 1.E-12,
                              fmin(theta_new, soil_properties[soil_type[current->layer_num]].theta_e));

        const int soil_num_k = soil_type[current->layer_num];
        const double theta_e_k = soil_properties[soil_num_k].theta_e;
        const double theta_r_k = soil_properties[soil_num_k].theta_r;
        const double vg_a_k = soil_properties[soil_num_k].vg_alpha_per_cm;
        const double vg_m_k = soil_properties[soil_num_k].vg_m;
        const double vg_n_k = soil_properties[soil_num_k].vg_n;
        const double Se = calc_Se_from_theta(current->theta, theta_e_k, theta_r_k);
        current->psi_cm = calc_h_from_Se(Se, vg_a_k, vg_m_k, vg_n_k);

        if (current->theta == theta_r_k + 1.E-12) {
          const double requested_flux_cm = *bottom_boundary_flux_above_surface_WFs_cm;
          const double actual_extracted_flux_cm =
            mass_at_start_of_extraction - lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
          lgarto_assert_gw_flux_mass_balance_correction_within_debug_threshold(
            actual_extracted_flux_cm - requested_flux_cm,
            "lgarto_extract_TO_GW_flux_from_surface_WFs",
            "top-layer extraction theta_r residual adjustment",
            *head);
          bottom_boundary_flux_cm -= requested_flux_cm;
          bottom_boundary_flux_cm += actual_extracted_flux_cm;
          *bottom_boundary_flux_above_surface_WFs_cm = actual_extracted_flux_cm;
        }
      }
    }
    else if (next != NULL) {
      int layer_num = current->layer_num;
      int soil_num = soil_type[layer_num];

      double *delta_thetas = (double *)malloc(sizeof(double) * (layer_num + 1));
      double *delta_thickness = (double *)malloc(sizeof(double) * (layer_num + 1));

      const double psi_cm = current->psi_cm;
      const double psi_cm_below = next->psi_cm;

      double new_mass =
        (current->depth_cm - cum_layer_thickness_cm[layer_num - 1]) * (current->theta - next->theta);
      double prior_mass = new_mass;

      for (int k = 1; k < layer_num; k++) {
        const int soil_num_k = soil_type[k];
        const double theta_e_k = soil_properties[soil_num_k].theta_e;
        const double theta_r_k = soil_properties[soil_num_k].theta_r;
        const double vg_a_k = soil_properties[soil_num_k].vg_alpha_per_cm;
        const double vg_m_k = soil_properties[soil_num_k].vg_m;
        const double vg_n_k = soil_properties[soil_num_k].vg_n;
        const double layer_thickness = cum_layer_thickness_cm[k] - cum_layer_thickness_cm[k - 1];

        const double theta = calc_theta_from_h(psi_cm, vg_a_k, vg_m_k, vg_n_k, theta_e_k, theta_r_k);
        const double theta_below =
          calc_theta_from_h(psi_cm_below, vg_a_k, vg_m_k, vg_n_k, theta_e_k, theta_r_k);

        new_mass += layer_thickness * (theta - theta_below);
        prior_mass += layer_thickness * (theta - theta_below);

        delta_thetas[k] = theta_below;
        delta_thickness[k] = layer_thickness;
      }

      delta_thetas[layer_num] = next->theta;
      delta_thickness[layer_num] = current->depth_cm - cum_layer_thickness_cm[layer_num - 1];

      prior_mass -= *bottom_boundary_flux_above_surface_WFs_cm;

      double theta_new = current->theta;
      if (*bottom_boundary_flux_above_surface_WFs_cm != 0.0) {
        theta_new = lgar_theta_mass_balance(layer_num, soil_num, current->psi_cm, new_mass, prior_mass, 0.0,
                                            AET_demand_cm, delta_thetas, delta_thickness, soil_type,
                                            soil_properties, true);
      }

      current->theta = fmax(soil_properties[soil_type[current->layer_num]].theta_r,
                            fmin(theta_new, soil_properties[soil_type[current->layer_num]].theta_e));

      layer_num = current->layer_num;
      soil_num = soil_type[layer_num];
      const double theta_e_k = soil_properties[soil_num].theta_e;
      const double theta_r_k = soil_properties[soil_num].theta_r;
      const double vg_a_k = soil_properties[soil_num].vg_alpha_per_cm;
      const double vg_m_k = soil_properties[soil_num].vg_m;
      const double vg_n_k = soil_properties[soil_num].vg_n;

      if (current->theta == theta_r_k) {
        current->theta = current->theta + 1.E-12;
        const double requested_flux_cm = *bottom_boundary_flux_above_surface_WFs_cm;
        const double actual_extracted_flux_cm =
          mass_at_start_of_extraction - lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
        lgarto_assert_gw_flux_mass_balance_correction_within_debug_threshold(
          actual_extracted_flux_cm - requested_flux_cm,
          "lgarto_extract_TO_GW_flux_from_surface_WFs",
          "deeper-layer extraction theta_r residual adjustment",
          *head);
        bottom_boundary_flux_cm -= requested_flux_cm;
        bottom_boundary_flux_cm += actual_extracted_flux_cm;
        *bottom_boundary_flux_above_surface_WFs_cm = actual_extracted_flux_cm;
      }

      const double Se = calc_Se_from_theta(current->theta, theta_e_k, theta_r_k);
      if (*bottom_boundary_flux_above_surface_WFs_cm != 0.0) {
        current->psi_cm = calc_h_from_Se(Se, vg_a_k, vg_m_k, vg_n_k);
        *bottom_boundary_flux_above_surface_WFs_cm -=
          (mass_at_start_of_extraction - lgar_calc_mass_bal(cum_layer_thickness_cm, *head));
      }

      free(delta_thetas);
      free(delta_thickness);
    }
  }
  else {
    bottom_boundary_flux_cm -= *bottom_boundary_flux_above_surface_WFs_cm;
  }

  listSendToTop(*head);

  current = listFindFront(listLength(*head) - 1, *head, NULL);
  if (current == NULL) {
    current = *head;
  }
  next = current != NULL ? current->next : NULL;
  if (next != NULL) {
    for (int wf = listLength(*head) - 1; wf != 1; wf--) {
      if (current->to_bottom == TRUE) {
        current->psi_cm = next->psi_cm;
      }
      current = listFindFront(wf - 1, *head, NULL);
      next = current != NULL ? current->next : NULL;
      if (current == NULL) {
        break;
      }
    }
  }

  current = *head;
  next = current != NULL ? current->next : NULL;
  bool deepest_surface_WF_became_TO_in_layer_deeper_than_top = false;
  for (int wf = 1; wf != listLength(*head); wf++) {
    if (current == NULL || next == NULL) {
      break;
    }
    if ((current->is_WF_GW == 0) && (next->is_WF_GW == 1) &&
        (current->layer_num == next->layer_num) && (current->theta <= next->theta)) {
      current->is_WF_GW = 1;
      if (current->layer_num > 1) {
        deepest_surface_WF_became_TO_in_layer_deeper_than_top = true;
      }
    }
    current = current->next;
    next = current != NULL ? current->next : NULL;
  }

  current = listFindFront(listLength(*head), *head, NULL);
  if (deepest_surface_WF_became_TO_in_layer_deeper_than_top) {
    for (int wf = listLength(*head); wf != 0; wf--) {
      if (current == NULL) {
        break;
      }
      if (current->next != NULL) {
        if (current->to_bottom == TRUE && (current->next->is_WF_GW == TRUE)) {
          current->is_WF_GW = 1;
        }
      }
      current = listFindFront(current->front_num - 1, *head, NULL);
    }
  }

  current = *head;
  // if (current != NULL && current->next != NULL) {
  //   if ((current->depth_cm == 0.0) && (current->next->depth_cm != 0.0)) {
  //     listDeleteFront(1, head, soil_type, soil_properties);
  //   }
  // }

  current = listFindFront(wf_from_which_to_extract, *head, NULL);
  if (current != NULL) {
    const int layer_num = current->layer_num;
    const int soil_num = soil_type[layer_num];
    const double theta_r = soil_properties[soil_num].theta_r;
    bool all_WFs_to_GW = false;
    if (current->theta < theta_r) {
      current->is_WF_GW = TRUE;
      current->theta = original_theta;
      current->psi_cm = original_psi_cm;

      bottom_boundary_flux_cm = bottom_bdy_flux_at_start_of_extraction + *bottom_boundary_flux_above_surface_WFs_cm;
      *bottom_boundary_flux_above_surface_WFs_cm = 0.0;
      all_WFs_to_GW = true;
    }
    if (all_WFs_to_GW) {
      current = *head;
      for (int wf = 1; wf != listLength(*head); wf++) {
        if (current == NULL) {
          break;
        }
        current->is_WF_GW = TRUE;
        current = current->next;
      }
    }
  }

  if (verbosity.compare("high") == 0) {
    printf("After TO demand subtracted from rightmost surface WF ... \n");
    printf("bottom_boundary_flux_above_surface_WFs_cm: %.17lf \n", *bottom_boundary_flux_above_surface_WFs_cm);
    listPrint(*head);
    printf("associated mass: %lf \n", lgar_calc_mass_bal(cum_layer_thickness_cm, *head));
  }

  return bottom_boundary_flux_cm;
}

// ############################################################################################
/*
  Updates theta globally from psi after TO surface-flux extraction modifies psi continuity.
*/
// ############################################################################################

struct lgarto_global_theta_snap_state
{
  struct wetting_front *front;
  double theta;
  double psi_cm;
};

static std::vector<lgarto_global_theta_snap_state>
lgarto_snapshot_theta_snap_state(struct wetting_front *head)
{
  std::vector<lgarto_global_theta_snap_state> states;
  for (struct wetting_front *current = head; current != NULL; current = current->next) {
    lgarto_global_theta_snap_state state;
    state.front = current;
    state.theta = current->theta;
    state.psi_cm = current->psi_cm;
    states.push_back(state);
  }
  return states;
}

static void lgarto_restore_theta_snap_state(
  const std::vector<lgarto_global_theta_snap_state> &states)
{
  for (size_t i = 0; i < states.size(); i++) {
    states[i].front->theta = states[i].theta;
    states[i].front->psi_cm = states[i].psi_cm;
  }
}

static const lgarto_global_theta_snap_state *
lgarto_find_theta_snap_state(const std::vector<lgarto_global_theta_snap_state> &states,
                             const struct wetting_front *front)
{
  for (size_t i = 0; i < states.size(); i++) {
    if (states[i].front == front) {
      return &states[i];
    }
  }
  return NULL;
}

static bool lgarto_front_changed_during_theta_snap(
  const std::vector<lgarto_global_theta_snap_state> &states,
  const struct wetting_front *front)
{
  const lgarto_global_theta_snap_state *state =
    lgarto_find_theta_snap_state(states, front);
  if (state == NULL || front == NULL) {
    return false;
  }

  return fabs(front->theta - state->theta) > THRESHOLD_NO_MOISTURE_DIFF ||
         fabs(front->psi_cm - state->psi_cm) > 1.0e-12;
}

static std::vector<struct wetting_front *>
lgarto_collect_changed_to_bottom_snap_scaffold(
  struct wetting_front *head,
  const std::vector<lgarto_global_theta_snap_state> &pre_snap_states)
{
  std::vector<struct wetting_front *> scaffold;
  bool changed_scaffold_front = false;
  struct wetting_front *previous = NULL;

  for (struct wetting_front *current = head; current != NULL; current = current->next) {
    const bool previous_boundary =
      previous != NULL && previous->to_bottom && current->to_bottom &&
      previous->layer_num != current->layer_num;
    const bool next_boundary =
      current->next != NULL && current->to_bottom && current->next->to_bottom &&
      current->layer_num != current->next->layer_num;
    const bool in_to_bottom_boundary_scaffold =
      current->to_bottom && (previous_boundary || next_boundary);

    if (in_to_bottom_boundary_scaffold) {
      scaffold.push_back(current);
      if (lgarto_front_changed_during_theta_snap(pre_snap_states, current)) {
        changed_scaffold_front = true;
      }
    }

    previous = current;
  }

  if (!changed_scaffold_front) {
    scaffold.clear();
  }
  return scaffold;
}

static bool lgarto_apply_shared_boundary_snap_psi(
  const std::vector<struct wetting_front *> &fronts,
  double psi_cm,
  int *soil_type,
  struct soil_properties_ *soil_properties)
{
  if (!std::isfinite(psi_cm) || psi_cm < 0.0) {
    return false;
  }

  for (size_t i = 0; i < fronts.size(); i++) {
    struct wetting_front *front = fronts[i];
    if (front == NULL) {
      return false;
    }

    const int soil_num = soil_type[front->layer_num];
    const double theta_e = soil_properties[soil_num].theta_e;
    const double theta_r = soil_properties[soil_num].theta_r;
    const double vg_a = soil_properties[soil_num].vg_alpha_per_cm;
    const double vg_m = soil_properties[soil_num].vg_m;
    const double vg_n = soil_properties[soil_num].vg_n;

    front->psi_cm = psi_cm;
    front->theta = calc_theta_from_h(psi_cm, vg_a, vg_m, vg_n, theta_e, theta_r);
    if (!std::isfinite(front->theta) ||
        front->theta < theta_r - 1.0e-12 ||
        front->theta > theta_e + 1.0e-12) {
      return false;
    }
  }

  return true;
}

static void lgarto_apply_to_bottom_boundary_scaffold_theta_snap(
  struct wetting_front *head,
  int *soil_type,
  struct soil_properties_ *soil_properties)
{
  if (head == NULL || soil_type == NULL || soil_properties == NULL) {
    return;
  }

  std::vector<struct wetting_front *> fronts;
  for (struct wetting_front *current = head; current != NULL; current = current->next) {
    fronts.push_back(current);
  }

  if (fronts.size() < 2) {
    return;
  }

  // Surface-flux extraction can copy psi onto a to_bottom boundary scaffold
  // before this routine runs.  Walk bottom-up so every upper boundary front
  // inherits the final lower psi, and update theta immediately for that soil.
  for (int i = static_cast<int>(fronts.size()) - 2; i >= 0; i--) {
    struct wetting_front *current = fronts[i];
    struct wetting_front *next = fronts[i + 1];
    if (current == NULL || next == NULL ||
        current->to_bottom != TRUE || next->to_bottom != TRUE ||
        current->layer_num == next->layer_num ||
        !std::isfinite(next->psi_cm) || next->psi_cm < 0.0) {
      continue;
    }

    const int soil_num_k = soil_type[current->layer_num];
    const double theta_e_k = soil_properties[soil_num_k].theta_e;
    const double theta_r_k = soil_properties[soil_num_k].theta_r;
    const double vg_a_k = soil_properties[soil_num_k].vg_alpha_per_cm;
    const double vg_m_k = soil_properties[soil_num_k].vg_m;
    const double vg_n_k = soil_properties[soil_num_k].vg_n;

    current->psi_cm = next->psi_cm;
    current->theta = calc_theta_from_h(current->psi_cm, vg_a_k, vg_m_k, vg_n_k, theta_e_k, theta_r_k);
  }
}

static bool lgarto_mass_aware_to_bottom_snap(
  double target_mass_cm,
  double *cum_layer_thickness_cm,
  int *soil_type,
  struct soil_properties_ *soil_properties,
  struct wetting_front **head,
  const std::vector<struct wetting_front *> &snap_scaffold)
{
  if (head == NULL || *head == NULL || snap_scaffold.empty()) {
    return false;
  }

  const std::vector<lgarto_global_theta_snap_state> simple_snap_states =
    lgarto_snapshot_theta_snap_state(*head);

  double psi_low_cm = 0.0;
  double psi_high_cm = 0.0;
  for (size_t i = 0; i < snap_scaffold.size(); i++) {
    if (snap_scaffold[i] != NULL && std::isfinite(snap_scaffold[i]->psi_cm)) {
      psi_high_cm = fmax(psi_high_cm, snap_scaffold[i]->psi_cm);
    }
  }
  psi_high_cm = fmax(psi_high_cm, 1.0e-12);

  if (!lgarto_apply_shared_boundary_snap_psi(snap_scaffold, psi_low_cm,
                                             soil_type, soil_properties)) {
    lgarto_restore_theta_snap_state(simple_snap_states);
    return false;
  }
  double mass_low_cm = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);

  if (!lgarto_apply_shared_boundary_snap_psi(snap_scaffold, psi_high_cm,
                                             soil_type, soil_properties)) {
    lgarto_restore_theta_snap_state(simple_snap_states);
    return false;
  }
  double mass_high_cm = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);

  while (mass_high_cm > target_mass_cm && psi_high_cm < PSI_UPPER_LIM) {
    psi_low_cm = psi_high_cm;
    mass_low_cm = mass_high_cm;
    psi_high_cm = fmin(PSI_UPPER_LIM, fmax(psi_high_cm * 2.0, psi_high_cm + 1.0e-12));

    if (!lgarto_apply_shared_boundary_snap_psi(snap_scaffold, psi_high_cm,
                                               soil_type, soil_properties)) {
      lgarto_restore_theta_snap_state(simple_snap_states);
      return false;
    }
    mass_high_cm = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
  }

  const double bracket_tol_cm = fmax(MBAL_ITERATIVE_TOLERANCE, 1.0e-12);
  if (target_mass_cm > mass_low_cm + bracket_tol_cm ||
      target_mass_cm < mass_high_cm - bracket_tol_cm) {
    lgarto_restore_theta_snap_state(simple_snap_states);
    return false;
  }

  double best_psi_cm = psi_low_cm;
  double best_mass_error_cm = fabs(mass_low_cm - target_mass_cm);
  if (fabs(mass_high_cm - target_mass_cm) < best_mass_error_cm) {
    best_psi_cm = psi_high_cm;
    best_mass_error_cm = fabs(mass_high_cm - target_mass_cm);
  }

  for (int iter = 0; iter < 100; iter++) {
    const double psi_mid_cm = 0.5 * (psi_low_cm + psi_high_cm);
    if (!lgarto_apply_shared_boundary_snap_psi(snap_scaffold, psi_mid_cm,
                                               soil_type, soil_properties)) {
      lgarto_restore_theta_snap_state(simple_snap_states);
      return false;
    }

    const double mass_mid_cm = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
    const double mass_error_cm = mass_mid_cm - target_mass_cm;
    if (fabs(mass_error_cm) < best_mass_error_cm) {
      best_psi_cm = psi_mid_cm;
      best_mass_error_cm = fabs(mass_error_cm);
    }
    if (fabs(mass_error_cm) <= MBAL_ITERATIVE_TOLERANCE) {
      best_psi_cm = psi_mid_cm;
      break;
    }

    if (mass_mid_cm > target_mass_cm) {
      psi_low_cm = psi_mid_cm;
    }
    else {
      psi_high_cm = psi_mid_cm;
    }
  }

  if (!lgarto_apply_shared_boundary_snap_psi(snap_scaffold, best_psi_cm,
                                             soil_type, soil_properties)) {
    lgarto_restore_theta_snap_state(simple_snap_states);
    return false;
  }

  return true;
}

extern void lgar_global_theta_update(double bottom_boundary_flux_above_surface_WFs_cm,
				     double target_mass_after_fluxes_cm,
				     double *cum_layer_thickness_cm,
				     int *soil_type,
				     struct soil_properties_ *soil_properties,
				     struct wetting_front **head,
				     double mass_balance_tolerance_cm)
{
  const double mass_before_snap_cm =
    (cum_layer_thickness_cm != NULL && head != NULL && *head != NULL)
      ? lgar_calc_mass_bal(cum_layer_thickness_cm, *head)
      : 0.0;
  const std::vector<lgarto_global_theta_snap_state> pre_snap_states =
    (head != NULL && *head != NULL)
      ? lgarto_snapshot_theta_snap_state(*head)
      : std::vector<lgarto_global_theta_snap_state>();

  if (verbosity.compare("high") == 0) {
    printf("before lgar_global_theta_update \n");
    listPrint(*head);
  }

  const int wf_free_drainage_demand = wetting_front_free_drainage(*head);
  for (int wf = wf_free_drainage_demand; wf > 1; wf--) {
    struct wetting_front *current = listFindFront(wf, *head, NULL);
    struct wetting_front *previous = listFindFront(wf - 1, *head, NULL);

    if ((previous != NULL) && (current != NULL)) {
      if ((bottom_boundary_flux_above_surface_WFs_cm != 0.0) && (current->layer_num != previous->layer_num)) {
        previous->psi_cm = current->psi_cm;
        const int soil_num_k = soil_type[previous->layer_num];
        const double theta_e_k = soil_properties[soil_num_k].theta_e;
        const double theta_r_k = soil_properties[soil_num_k].theta_r;
        const double vg_a_k = soil_properties[soil_num_k].vg_alpha_per_cm;
        const double vg_m_k = soil_properties[soil_num_k].vg_m;
        const double vg_n_k = soil_properties[soil_num_k].vg_n;
        previous->theta = calc_theta_from_h(previous->psi_cm, vg_a_k, vg_m_k, vg_n_k, theta_e_k, theta_r_k);
      }
    }
  }

  if (bottom_boundary_flux_above_surface_WFs_cm != 0.0) {
    lgarto_apply_to_bottom_boundary_scaffold_theta_snap(*head, soil_type, soil_properties);
  }

  for (int wf = wf_free_drainage_demand - 1; wf > 1; wf--) {
    struct wetting_front *current = listFindFront(wf, *head, NULL);
    if ((current != NULL) && (current->next != NULL)) {
      if ((current->to_bottom == TRUE) && (current->next->is_WF_GW == 0)) {
        current->psi_cm = current->next->psi_cm;
        const int soil_num_k = soil_type[current->layer_num];
        const double theta_e_k = soil_properties[soil_num_k].theta_e;
        const double theta_r_k = soil_properties[soil_num_k].theta_r;
        const double vg_a_k = soil_properties[soil_num_k].vg_alpha_per_cm;
        const double vg_m_k = soil_properties[soil_num_k].vg_m;
        const double vg_n_k = soil_properties[soil_num_k].vg_n;
        current->theta = calc_theta_from_h(current->psi_cm, vg_a_k, vg_m_k, vg_n_k, theta_e_k, theta_r_k);
      }
    }
  }

  if (verbosity.compare("high") == 0) {
    printf("after lgar_global_theta_update \n");
    listPrint(*head);
  }

  if (cum_layer_thickness_cm == NULL || head == NULL || *head == NULL ||
      !std::isfinite(mass_balance_tolerance_cm) ||
      mass_balance_tolerance_cm <= 0.0) {
    return;
  }

  const double mass_after_simple_snap_cm =
    lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
  const double simple_snap_mass_change_cm =
    mass_after_simple_snap_cm - mass_before_snap_cm;
  const bool use_target_mass =
    std::isfinite(target_mass_after_fluxes_cm) && target_mass_after_fluxes_cm > 0.0;
  const double target_mass_cm =
    use_target_mass ? target_mass_after_fluxes_cm : mass_before_snap_cm;
  const double simple_snap_target_error_cm =
    mass_after_simple_snap_cm - target_mass_cm;

  // Some simple boundary snaps are the intended storage response to an extracted TO flux.
  // Only replace the snap when its final storage would miss the timestep mass target.
  if (fabs(simple_snap_target_error_cm) <= mass_balance_tolerance_cm) {
    return;
  }

  const std::vector<struct wetting_front *> snap_scaffold =
    lgarto_collect_changed_to_bottom_snap_scaffold(*head, pre_snap_states);
  if (snap_scaffold.empty()) {
    return;
  }

  if (lgarto_mass_aware_to_bottom_snap(target_mass_cm, cum_layer_thickness_cm,
                                       soil_type, soil_properties, head,
                                       snap_scaffold)) {
    if (verbosity.compare("high") == 0) {
      const double mass_after_solve_cm =
        lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
      printf("Mass-aware to_bottom psi snap replaced simple snap "
             "(simple_mass_change_cm=%.17g, simple_target_error_cm=%.17g, "
             "solved_target_error_cm=%.17g, target_mass_cm=%.17g, trigger_tol_cm=%.17g).\n",
             simple_snap_mass_change_cm,
             simple_snap_target_error_cm,
             mass_after_solve_cm - target_mass_cm,
             target_mass_cm,
             mass_balance_tolerance_cm);
      listPrint(*head);
    }
  }
}

// ############################################################################################
/*
  Updates psi globally from theta after TO corrections that modify theta directly.
*/
// ############################################################################################

extern void lgar_global_psi_update(int *soil_type, struct soil_properties_ *soil_properties,
				   struct wetting_front **head)
{
  struct wetting_front *current = *head;
  for (int wf = 1; wf != listLength(*head); wf++) {
    if (current == NULL) {
      break;
    }

    const int soil_num = soil_type[current->layer_num];
    const double theta_e = soil_properties[soil_num].theta_e;
    const double theta_r = soil_properties[soil_num].theta_r;
    const double vg_a = soil_properties[soil_num].vg_alpha_per_cm;
    const double vg_m = soil_properties[soil_num].vg_m;
    const double vg_n = soil_properties[soil_num].vg_n;

    const double Se = calc_Se_from_theta(current->theta, theta_e, theta_r);
    current->psi_cm = calc_h_from_Se(Se, vg_a, vg_m, vg_n);

    current = current->next;
  }

  for (current = *head; current != NULL && current->next != NULL; current = current->next) {
    struct wetting_front *next = current->next;
    if (!current->to_bottom || current->layer_num == next->layer_num) {
      continue;
    }

    if (!lgar_should_snap_boundary_psi_to_lower(current, next, soil_type, soil_properties)) {
      continue;
    }

    current->psi_cm = next->psi_cm;
    const int soil_num = soil_type[current->layer_num];
    const double theta_e = soil_properties[soil_num].theta_e;
    const double theta_r = soil_properties[soil_num].theta_r;
    const double vg_a = soil_properties[soil_num].vg_alpha_per_cm;
    const double vg_m = soil_properties[soil_num].vg_m;
    const double vg_n = soil_properties[soil_num].vg_n;
    current->theta = calc_theta_from_h(current->psi_cm, vg_a, vg_m, vg_n, theta_e, theta_r);
  }

  if (verbosity.compare("high") == 0) {
    printf("states after lgar_global_psi_update, which comes after TO merging and layer bdy crossing ... \n");
    listPrint(*head);
  }
}

// ############################################################################################
/*
  Checks that no wetting front persists below the vadose-zone lower boundary once a
  substep or timestep has completed.
*/
// ############################################################################################

extern void lgar_assert_wetting_fronts_within_vadose_zone(double domain_depth_cm,
                                                          struct wetting_front *head)
{
  const double depth_tol_cm = LOWER_BOUNDARY_FINAL_TOL_CM;

  if (head == NULL) {
    return;
  }

  for (struct wetting_front *current = head; current != NULL; current = current->next) {
    if (current->depth_cm <= domain_depth_cm + depth_tol_cm) {
      continue;
    }

    fprintf(stderr,
            "Error: wetting front persisted below the vadose-zone lower boundary.\n"
            "  offending WF front_num=%d layer=%d depth_cm=%.16g written_depth_mm=%.16g theta=%.16g psi_cm=%.16g dzdt_cm_per_h=%.16g is_WF_GW=%d\n"
            "  lower_boundary_depth_cm=%.16g written_lower_boundary_mm=%.16g\n"
            "  depth_excess_cm=%.16g > %.16g\n"
            "  Wetting front list follows:\n",
            current->front_num,
            current->layer_num,
            current->depth_cm,
            current->depth_cm * 10.0,
            current->theta,
            current->psi_cm,
            current->dzdt_cm_per_h,
            current->is_WF_GW,
            domain_depth_cm,
            domain_depth_cm * 10.0,
            current->depth_cm - domain_depth_cm,
            depth_tol_cm);
    fflush(stderr);
    listPrint(head);
    fflush(stdout);
    abort();
  }
}

// ############################################################################################
/*
  Checks that no wetting front keeps a meaningfully negative depth once a substep or timestep
  has completed.
*/
// ############################################################################################

extern void lgar_assert_wetting_fronts_nonnegative_depth(struct wetting_front *head)
{
  const double min_depth_cm = -1.0E-6;

  if (head == NULL) {
    return;
  }

  for (struct wetting_front *current = head; current != NULL; current = current->next) {
    if (current->depth_cm >= min_depth_cm) {
      continue;
    }

    fprintf(stderr,
            "Error: wetting front persisted with a negative depth.\n"
            "  offending WF front_num=%d layer=%d depth_cm=%.16g written_depth_mm=%.16g theta=%.16g psi_cm=%.16g dzdt_cm_per_h=%.16g is_WF_GW=%d to_bottom=%d\n"
            "  min_allowed_depth_cm=%.16g written_min_allowed_depth_mm=%.16g\n"
            "  Wetting front list follows:\n",
            current->front_num,
            current->layer_num,
            current->depth_cm,
            current->depth_cm * 10.0,
            current->theta,
            current->psi_cm,
            current->dzdt_cm_per_h,
            current->is_WF_GW,
            current->to_bottom,
            min_depth_cm,
            min_depth_cm * 10.0);
    fflush(stderr);
    listPrint(head);
    fflush(stdout);
    abort();
  }
}

// ############################################################################################
/*
  Checks that adjacent TO wetting fronts in the same layer become no wetter in psi as depth
  increases. Equal-depth support fronts are ignored, and fronts on opposite sides of a soil
  layer boundary are intentionally not checked here.
*/
// ############################################################################################

static double lgar_psi_assertion_tolerance_cm(double psi_a_cm, double psi_b_cm)
{
  const double abs_psi_tol_cm = 1.E-3;
  const double rel_psi_tol = 1.E-9;
  const double near_saturation_psi_cm = 0.1;
  const double near_saturation_psi_tol_cm = 0.1;
  const double very_dry_psi_cm = 1.E5;
  const double very_dry_rel_psi_tol = 1.E-5;
  const double psi_scale_cm = fmax(1.0, fmax(std::fabs(psi_a_cm), std::fabs(psi_b_cm)));
  double psi_tol_cm = fmax(abs_psi_tol_cm, rel_psi_tol * psi_scale_cm);

  /* Near saturation, tiny nonzero psi values can round-trip to effectively the
     same theta as psi=0. At very large psi, the water-retention curve is flat
     enough that tiny relative psi reversals are not meaningful storage states. */
  if (fmax(std::fabs(psi_a_cm), std::fabs(psi_b_cm)) <= near_saturation_psi_cm) {
    psi_tol_cm = fmax(psi_tol_cm, near_saturation_psi_tol_cm);
  }

  if (psi_scale_cm >= very_dry_psi_cm) {
    psi_tol_cm = fmax(psi_tol_cm, very_dry_rel_psi_tol * psi_scale_cm);
  }

  return psi_tol_cm;
}

struct lgarto_surface_TO_support_ordering_state
{
  bool active;
  bool has_zero_depth_TO_support;
  double min_zero_depth_TO_psi_cm;
  struct wetting_front *first_TO_below_surface;
};

static lgarto_surface_TO_support_ordering_state
lgarto_find_surface_TO_support_ordering_state(struct wetting_front *head)
{
  lgarto_surface_TO_support_ordering_state state = {};
  state.min_zero_depth_TO_psi_cm = HUGE_VAL;

  bool saw_surface_front = false;
  for (struct wetting_front *current = head; current != NULL; current = current->next) {
    if (!saw_surface_front && current->is_WF_GW == TRUE &&
        current->to_bottom == FALSE &&
        current->layer_num == 1 &&
        std::fabs(current->depth_cm) <= ZERO_DEPTH_TO_DELETE_DEPTH_TOL_CM &&
        std::isfinite(current->psi_cm)) {
      state.has_zero_depth_TO_support = true;
      state.min_zero_depth_TO_psi_cm =
        fmin(state.min_zero_depth_TO_psi_cm, current->psi_cm);
      continue;
    }

    if (current->is_WF_GW == FALSE) {
      saw_surface_front = true;
      continue;
    }

    if (saw_surface_front && current->is_WF_GW == TRUE) {
      state.first_TO_below_surface = current;
      state.active = state.has_zero_depth_TO_support &&
                     std::isfinite(state.min_zero_depth_TO_psi_cm) &&
                     std::isfinite(current->psi_cm);
      return state;
    }
  }

  return state;
}

static bool lgarto_zero_depth_TO_support_ordering_is_valid(
  const lgarto_surface_TO_support_ordering_state &state)
{
  if (!state.active || state.first_TO_below_surface == NULL) {
    return true;
  }

  const double psi_tol_cm =
    lgar_psi_assertion_tolerance_cm(state.min_zero_depth_TO_psi_cm,
                                    state.first_TO_below_surface->psi_cm);
  return state.min_zero_depth_TO_psi_cm + psi_tol_cm >=
         state.first_TO_below_surface->psi_cm;
}

static bool lgarto_same_layer_TO_reversal_is_surface_supported_bridge(
  struct wetting_front *head,
  const struct wetting_front *current)
{
  const lgarto_surface_TO_support_ordering_state state =
    lgarto_find_surface_TO_support_ordering_state(head);
  if (!state.active || state.first_TO_below_surface != current) {
    return false;
  }

  const double psi_tol_cm =
    lgar_psi_assertion_tolerance_cm(state.min_zero_depth_TO_psi_cm,
                                    current->psi_cm);
  return std::fabs(current->psi_cm - state.min_zero_depth_TO_psi_cm) <= psi_tol_cm;
}

extern void lgar_assert_to_psi_monotonic_with_depth(struct wetting_front *head)
{
  const double depth_tol_cm = 1.E-10;

  for (struct wetting_front *current = head; current != NULL && current->next != NULL; current = current->next) {
    struct wetting_front *next = current->next;

    if (!current->is_WF_GW || !next->is_WF_GW) {
      continue;
    }

    if (current->layer_num != next->layer_num) {
      continue;
    }

    if (next->depth_cm <= current->depth_cm + depth_tol_cm) {
      continue;
    }

    const double psi_tol_cm = lgar_psi_assertion_tolerance_cm(current->psi_cm, next->psi_cm);
    if (next->psi_cm <= current->psi_cm + psi_tol_cm) {
      continue;
    }

    if (lgarto_same_layer_TO_reversal_is_surface_supported_bridge(head, current)) {
      continue;
    }

    fprintf(stderr,
            "Error: adjacent TO wetting fronts reverse psi ordering with depth in the same soil layer.\n"
            "  shallower front_num=%d layer=%d depth_cm=%.16g written_depth_mm=%.16g theta=%.16g psi_cm=%.16g dzdt_cm_per_h=%.16g to_bottom=%d\n"
            "  deeper   front_num=%d layer=%d depth_cm=%.16g written_depth_mm=%.16g theta=%.16g psi_cm=%.16g dzdt_cm_per_h=%.16g to_bottom=%d\n"
            "  delta_depth_cm=%.16g delta_psi_cm=%.16g > %.16g\n"
            "  Note: this check intentionally skips adjacent TO fronts that straddle a soil layer boundary.\n"
            "  Wetting front list follows:\n",
            current->front_num,
            current->layer_num,
            current->depth_cm,
            current->depth_cm * 10.0,
            current->theta,
            current->psi_cm,
            current->dzdt_cm_per_h,
            current->to_bottom,
            next->front_num,
            next->layer_num,
            next->depth_cm,
            next->depth_cm * 10.0,
            next->theta,
            next->psi_cm,
            next->dzdt_cm_per_h,
            next->to_bottom,
            next->depth_cm - current->depth_cm,
            next->psi_cm - current->psi_cm,
            psi_tol_cm);
    fflush(stderr);
    listPrint(head);
    fflush(stdout);
    abort();
  }
}

extern void lgar_assert_zero_depth_TO_supports_drier_than_surface_TO_chain(struct wetting_front *head)
{
  const lgarto_surface_TO_support_ordering_state state =
    lgarto_find_surface_TO_support_ordering_state(head);

  if (lgarto_zero_depth_TO_support_ordering_is_valid(state)) {
    return;
  }

  const double psi_tol_cm =
    lgar_psi_assertion_tolerance_cm(state.min_zero_depth_TO_psi_cm,
                                    state.first_TO_below_surface->psi_cm);
  fprintf(stderr,
          "Error: leading zero-depth TO/GW supports are wetter than the first TO front below active surface fronts.\n"
          "  first_below_surface front_num=%d layer=%d depth_cm=%.16g written_depth_mm=%.16g theta=%.16g psi_cm=%.16g dzdt_cm_per_h=%.16g to_bottom=%d\n"
          "  min_zero_depth_TO_psi_cm=%.16g psi_tolerance_cm=%.16g\n"
          "  Rule: with active surface fronts, leading zero-depth TO supports must have psi >= the first TO front below the surface stack.\n"
          "  Wetting front list follows:\n",
          state.first_TO_below_surface->front_num,
          state.first_TO_below_surface->layer_num,
          state.first_TO_below_surface->depth_cm,
          state.first_TO_below_surface->depth_cm * 10.0,
          state.first_TO_below_surface->theta,
          state.first_TO_below_surface->psi_cm,
          state.first_TO_below_surface->dzdt_cm_per_h,
          state.first_TO_below_surface->to_bottom,
          state.min_zero_depth_TO_psi_cm,
          psi_tol_cm);
  fflush(stderr);
  listPrint(head);
  fflush(stdout);
  abort();
}

// ############################################################################################
/*
  Warns when psi is discontinuous across active soil-layer boundaries represented by a to_bottom
  wetting front and the front immediately below it.
*/
// ############################################################################################

extern void lgar_assert_boundary_psi_continuity(struct wetting_front *head)
{
  for (struct wetting_front *current = head; current != NULL && current->next != NULL; current = current->next) {
    struct wetting_front *next = current->next;

    if (!current->to_bottom) {
      continue;
    }

    if (next->layer_num == current->layer_num) {
      continue;
    }

    const double psi_mismatch_cm = std::fabs(current->psi_cm - next->psi_cm);
    const double psi_tol_cm = lgar_boundary_roundtrip_psi_tolerance_cm(current->psi_cm, next->psi_cm);
    if (psi_mismatch_cm > psi_tol_cm && (current->psi_cm > 0.1 && next->psi_cm > 0.1)) { //idea is that for small n values, nonzero but small psi will yield a theta
	                                                                                         //of exactly theta_e.
      std::cerr << "WARNING: psi mismatch across soil layer boundary exceeds tolerance. "
                << "This may indicate a real TO boundary-continuity error, or it may be "
                << "a storage-neutral roundoff case where very small nonzero psi values "
                << "evaluate to the same theta near saturation, or very large dry-tail "
                << "psi values differ in psi but evaluate to nearly identical theta.\n"
                << "  upper front_num=" << current->front_num
                << " layer=" << current->layer_num
                << " depth_cm=" << current->depth_cm
                << " psi_cm=" << current->psi_cm << "\n"
                << "  lower front_num=" << next->front_num
                << " layer=" << next->layer_num
                << " depth_cm=" << next->depth_cm
                << " psi_cm=" << next->psi_cm << "\n"
                << "  |delta_psi_cm|=" << psi_mismatch_cm
                << " > " << psi_tol_cm << "\n";
    }
  }
}

// ############################################################################################
/*
  Checks that the to_bottom scaffold contains exactly one wetting front per soil layer, and
  that each such front sits at that layer's lower boundary.
*/
// ############################################################################################

extern void lgar_assert_to_bottom_scaffold(int num_layers,
                                           double *cum_layer_thickness_cm,
                                           struct wetting_front *head)
{
  const double depth_tol_cm = 1.E-10;

  if (head == NULL) {
    fprintf(stderr,
            "Error: to_bottom scaffold check failed because the wetting-front list is empty.\n");
    fflush(stderr);
    abort();
  }

  std::vector<int> counts(num_layers + 1, 0);
  std::vector<struct wetting_front *> first_to_bottom(num_layers + 1, NULL);

  for (struct wetting_front *current = head; current != NULL; current = current->next) {
    if (!current->to_bottom) {
      continue;
    }

    if (current->layer_num < 1 || current->layer_num > num_layers) {
      fprintf(stderr,
              "Error: to_bottom wetting front has an invalid layer number.\n"
              "  offending front_num=%d layer=%d depth_cm=%.16g written_depth_mm=%.16g theta=%.16g psi_cm=%.16g dzdt_cm_per_h=%.16g is_WF_GW=%d\n"
              "  valid_layer_range=[1,%d]\n"
              "  Wetting front list follows:\n",
              current->front_num,
              current->layer_num,
              current->depth_cm,
              current->depth_cm * 10.0,
              current->theta,
              current->psi_cm,
              current->dzdt_cm_per_h,
              current->is_WF_GW,
              num_layers);
      fflush(stderr);
      listPrint(head);
      fflush(stdout);
      abort();
    }

    const double expected_depth_cm = cum_layer_thickness_cm[current->layer_num];
    const double depth_error_cm = std::fabs(current->depth_cm - expected_depth_cm);
    if (depth_error_cm > depth_tol_cm) {
      fprintf(stderr,
              "Error: to_bottom wetting front is not located at its soil-layer lower boundary.\n"
              "  offending front_num=%d layer=%d depth_cm=%.16g written_depth_mm=%.16g theta=%.16g psi_cm=%.16g dzdt_cm_per_h=%.16g is_WF_GW=%d\n"
              "  expected_boundary_depth_cm=%.16g written_expected_boundary_mm=%.16g\n"
              "  |depth_error_cm|=%.16g > %.16g\n"
              "  Wetting front list follows:\n",
              current->front_num,
              current->layer_num,
              current->depth_cm,
              current->depth_cm * 10.0,
              current->theta,
              current->psi_cm,
              current->dzdt_cm_per_h,
              current->is_WF_GW,
              expected_depth_cm,
              expected_depth_cm * 10.0,
              depth_error_cm,
              depth_tol_cm);
      fflush(stderr);
      listPrint(head);
      fflush(stdout);
      abort();
    }

    counts[current->layer_num] += 1;
    if (first_to_bottom[current->layer_num] == NULL) {
      first_to_bottom[current->layer_num] = current;
      continue;
    }

    struct wetting_front *previous = first_to_bottom[current->layer_num];
    fprintf(stderr,
            "Error: multiple to_bottom wetting fronts were found for the same soil layer.\n"
            "  layer=%d expected_exactly_one=1 found=%d\n"
            "  first offending pair:\n"
            "    front_num=%d depth_cm=%.16g written_depth_mm=%.16g theta=%.16g psi_cm=%.16g dzdt_cm_per_h=%.16g is_WF_GW=%d\n"
            "    front_num=%d depth_cm=%.16g written_depth_mm=%.16g theta=%.16g psi_cm=%.16g dzdt_cm_per_h=%.16g is_WF_GW=%d\n"
            "  Wetting front list follows:\n",
            current->layer_num,
            counts[current->layer_num],
            previous->front_num,
            previous->depth_cm,
            previous->depth_cm * 10.0,
            previous->theta,
            previous->psi_cm,
            previous->dzdt_cm_per_h,
            previous->is_WF_GW,
            current->front_num,
            current->depth_cm,
            current->depth_cm * 10.0,
            current->theta,
            current->psi_cm,
            current->dzdt_cm_per_h,
            current->is_WF_GW);
    fflush(stderr);
    listPrint(head);
    fflush(stdout);
    abort();
  }

  int total_to_bottom_count = 0;
  for (int layer = 1; layer <= num_layers; layer++) {
    total_to_bottom_count += counts[layer];
  }

  for (int layer = 1; layer <= num_layers; layer++) {
    if (counts[layer] == 1) {
      continue;
    }

    fprintf(stderr,
            "Error: the to_bottom scaffold does not contain exactly one boundary front per soil layer.\n"
            "  missing_or_invalid_layer=%d expected_exactly_one=1 found=%d total_to_bottom_found=%d expected_total=%d\n"
            "  expected_boundary_depth_cm=%.16g written_expected_boundary_mm=%.16g\n"
            "  Wetting front list follows:\n",
            layer,
            counts[layer],
            total_to_bottom_count,
            num_layers,
            cum_layer_thickness_cm[layer],
            cum_layer_thickness_cm[layer] * 10.0);
    fflush(stderr);
    listPrint(head);
    fflush(stdout);
    abort();
  }
}

// ############################################################################################
/*
  Clips any negative TO wetting-front depths back to zero after surface/TO correction steps.
*/
// ############################################################################################

extern bool lgarto_correct_negative_depths(struct wetting_front **head)
{
  struct wetting_front *current = *head;
  bool did_a_WF_have_negative_depth = false;

  for (int wf = 1; wf != listLength(*head); wf++) {
    if (current == NULL) {
      break;
    }

    if (current->depth_cm < 0.0) {
      current->depth_cm = 0.0;
      did_a_WF_have_negative_depth = true;
    }

    current = current->next;
  }

  return did_a_WF_have_negative_depth;
}

// ############################################################################################
/*
  the function lets groundwater/TO wetting fronts of a sufficient depth cross layer boundaries.
  This is the old LGARTO up/down crossing logic, adapted to the current branch.
*/
// ############################################################################################

extern bool lgar_TO_wetting_fronts_cross_layer_boundary(int *front_num_with_negative_depth,
							int num_layers, double *cum_layer_thickness_cm,
							int *soil_type, double *frozen_factor,
							struct wetting_front **head,
							struct soil_properties_ *soil_properties,
							double *mass_balance_flux_correction_cm)
{
  bool crossed = false;

  if (front_num_with_negative_depth != NULL) {
    *front_num_with_negative_depth = -1;
  }
  if (mass_balance_flux_correction_cm != NULL) {
    *mass_balance_flux_correction_cm = 0.0;
  }

  for (int wf = listLength(*head); wf != 0; wf--) {
    struct wetting_front *current = listFindFront(wf, *head, NULL);
    if (current == NULL) {
      continue;
    }

    struct wetting_front *previous = listFindFront(wf - 1, *head, NULL);
    struct wetting_front *next = current->next;
    const int layer_num = current->layer_num;

    if (next != NULL && current->depth_cm > cum_layer_thickness_cm[layer_num] && next->to_bottom &&
        current->is_WF_GW && (layer_num != num_layers)) {
      if (verbosity.compare("high") == 0) {
        printf("Inside layer boundary crossing for TO wetting front, moving down case ... \n");
        listPrint(*head);
      }

      const double mass_before_moving_down_crossing =
        lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
      const int soil_num = soil_type[layer_num];
      const double theta_e = soil_properties[soil_num].theta_e;
      const double theta_r = soil_properties[soil_num].theta_r;
      const double vg_a = soil_properties[soil_num].vg_alpha_per_cm;
      const double vg_m = soil_properties[soil_num].vg_m;
      const double vg_n = soil_properties[soil_num].vg_n;
      const double Ksat_cm_per_h =
        soil_properties[soil_num].Ksat_cm_per_h * frozen_factor[current->layer_num];

      struct wetting_front *next_to_next = current->next->next;
      if (next_to_next == NULL) {
        next_to_next = next;
      }

      const double current_theta = fmin(theta_e, current->theta);
      const double overshot_depth = current->depth_cm - next->depth_cm;
      const int soil_num_next = soil_type[layer_num + 1];

      const double next_theta_e = soil_properties[soil_num_next].theta_e;
      const double next_theta_r = soil_properties[soil_num_next].theta_r;
      const double next_vg_a = soil_properties[soil_num_next].vg_alpha_per_cm;
      const double next_vg_m = soil_properties[soil_num_next].vg_m;
      const double next_vg_n = soil_properties[soil_num_next].vg_n;

      const double Se = calc_Se_from_theta(current->theta, theta_e, theta_r);
      current->psi_cm = calc_h_from_Se(Se, vg_a, vg_m, vg_n);
      current->K_cm_per_h = calc_K_from_Se(Se, Ksat_cm_per_h, vg_m);

      double theta_new =
        calc_theta_from_h(current->psi_cm, next_vg_a, next_vg_m, next_vg_n, next_theta_e, next_theta_r);

      double mbal_correction = overshot_depth * (current_theta - next->theta);
      double denominator = theta_new - next_to_next->theta;
      double mbal_Z_correction = mbal_correction / denominator;
      double depth_new = cum_layer_thickness_cm[layer_num] + mbal_Z_correction;

      if (!std::isfinite(depth_new)) {
        theta_new -= 1.0E-14;
        denominator = theta_new - next_to_next->theta;
        mbal_Z_correction = mbal_correction / denominator;
        depth_new = cum_layer_thickness_cm[layer_num] + mbal_Z_correction;
      }

      if (!std::isfinite(depth_new)) {
        depth_new = cum_layer_thickness_cm[layer_num] + TRUNCATION_DEPTH;
      }

      current->depth_cm = cum_layer_thickness_cm[layer_num];
      next->theta = theta_new;
      next->psi_cm = current->psi_cm;
      next->depth_cm = depth_new;

      next->layer_num = layer_num + 1;
      next->dzdt_cm_per_h = current->dzdt_cm_per_h;
      // next->is_WF_GW = TRUE;

      current->dzdt_cm_per_h = 0.0;
      current->to_bottom = TRUE;
      next->to_bottom = FALSE;
      /*
       * This new to_bottom scaffold is psi-connected to the lower-layer front.
       * For extremely dry fronts, theta -> psi -> theta can change theta
       * measurably in some soils; account for the same psi-consistent theta
       * that the boundary snap will enforce before returning the storage residual.
       */
      lgarto_update_front_from_shared_psi(current, current->psi_cm, num_layers,
                                          soil_type, frozen_factor, soil_properties);
      crossed = true;

      if (!std::isfinite(next->depth_cm)) {
        next = listDeleteFront(next->front_num, head, soil_type, soil_properties);
        if (next == NULL) {
          break;
        }
      }

      const double moving_down_residual_cm =
        mass_before_moving_down_crossing -
        lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
      if (mass_balance_flux_correction_cm != NULL &&
          std::fabs(moving_down_residual_cm) > MBAL_ITERATIVE_TOLERANCE) {
        *mass_balance_flux_correction_cm += moving_down_residual_cm;
        if (verbosity.compare("high") == 0) {
          printf("TO moving-down layer crossing left residual %.12e cm after "
                 "psi-consistent boundary remap; routing residual to "
                 "lower-boundary flux bookkeeping.\n",
                 moving_down_residual_cm);
        }
      }

      if (verbosity.compare("high") == 0) {
        printf("After layer boundary crossing for TO wetting front, moving down case ... \n");
        listPrint(*head);
      }
    }

    if (previous != NULL && current->depth_cm < cum_layer_thickness_cm[layer_num - 1] &&
        previous->to_bottom && (current->layer_num > previous->layer_num) && current->is_WF_GW) {
      if (verbosity.compare("high") == 0) {
        printf("Inside layer boundary crossing for TO wetting front, moving up case ... \n");
        listPrint(*head);
      }

      const int soil_num = soil_type[layer_num];
      const double theta_e = soil_properties[soil_num].theta_e;
      const double theta_r = soil_properties[soil_num].theta_r;
      const double vg_a = soil_properties[soil_num].vg_alpha_per_cm;
      const double vg_m = soil_properties[soil_num].vg_m;
      const double vg_n = soil_properties[soil_num].vg_n;
      const double Ksat_cm_per_h =
        soil_properties[soil_num].Ksat_cm_per_h * frozen_factor[current->layer_num];

      const double current_theta = fmin(theta_e, current->theta);
      const double overshot_depth = previous->depth_cm - current->depth_cm;
      const int soil_num_above = soil_type[layer_num - 1];

      const double above_theta_e = soil_properties[soil_num_above].theta_e;
      const double above_theta_r = soil_properties[soil_num_above].theta_r;
      const double above_vg_a = soil_properties[soil_num_above].vg_alpha_per_cm;
      const double above_vg_m = soil_properties[soil_num_above].vg_m;
      const double above_vg_n = soil_properties[soil_num_above].vg_n;

      const double Se = calc_Se_from_theta(current->theta, theta_e, theta_r);
      current->psi_cm = calc_h_from_Se(Se, vg_a, vg_m, vg_n);
      current->K_cm_per_h = calc_K_from_Se(Se, Ksat_cm_per_h, vg_m);

      const double theta_new =
        calc_theta_from_h(current->psi_cm, above_vg_a, above_vg_m, above_vg_n, above_theta_e, above_theta_r);
      const double new_boundary_theta =
        calc_theta_from_h(next->psi_cm, above_vg_a, above_vg_m, above_vg_n, above_theta_e, above_theta_r);

      const double mbal_correction = overshot_depth * (next->theta - current_theta);
      const double denominator = theta_new - new_boundary_theta;

      /* A surface/TO merge can leave a lower-layer TO support above its
         assigned layer boundary while the next TO front in that same layer has
         identical theta/psi. That segment stores no water, so the moving-up
         remap has no depth to solve for (0/0). Delete the out-of-layer support
         instead of manufacturing a non-finite depth. */
      const bool redundant_same_layer_TO_segment =
        next != NULL && !current->to_bottom && current->is_WF_GW && next->is_WF_GW &&
        current->layer_num == next->layer_num &&
        std::fabs(next->theta - current_theta) <= THRESHOLD_NO_MOISTURE_DIFF &&
        std::fabs(current->psi_cm - next->psi_cm) <=
          lgar_psi_assertion_tolerance_cm(current->psi_cm, next->psi_cm) &&
        std::fabs(mbal_correction) <= MBAL_ITERATIVE_TOLERANCE &&
        std::fabs(denominator) <= THRESHOLD_NO_MOISTURE_DIFF;
      if (redundant_same_layer_TO_segment) {
        if (verbosity.compare("high") == 0) {
          printf("TO moving-up layer crossing deleting redundant zero-storage front %d "
                 "in layer %d before depth remap.\n",
                 current->front_num, current->layer_num);
        }
        listDeleteFront(current->front_num, head, soil_type, soil_properties);
        crossed = true;
        if (verbosity.compare("high") == 0) {
          printf("After deleting redundant TO front during moving-up layer crossing ... \n");
          listPrint(*head);
        }
        continue;
      }

      const double mass_before_moving_up_crossing =
        lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
      const double mbal_Z_correction = mbal_correction / denominator;
      double depth_new = cum_layer_thickness_cm[layer_num - 1] + mbal_Z_correction;
      const double upper_layer_top_cm = cum_layer_thickness_cm[previous->layer_num - 1];
      const double upper_layer_bottom_cm = cum_layer_thickness_cm[previous->layer_num];
      const bool bounded_moving_up_remap =
        !std::isfinite(depth_new) ||
        depth_new < upper_layer_top_cm ||
        depth_new > upper_layer_bottom_cm;
      if (bounded_moving_up_remap) {
        /*
         * The algebraic moving-up remap is mass conservative only if the old
         * to_bottom scaffold can become a real TO front inside the upper layer.
         * Near saturation, the theta contrast can be so small that the finite
         * solution lands above the layer. Keep the geometry valid and let the
         * adjacent boundary-connected TO pair solve the remaining storage by
         * psi, preserving boundary continuity without editing deeper fronts.
         */
        if (std::isfinite(depth_new)) {
          depth_new = fmax(upper_layer_top_cm, fmin(depth_new, upper_layer_bottom_cm));
        }
        else {
          depth_new = upper_layer_top_cm;
        }
      }

      current->depth_cm = cum_layer_thickness_cm[layer_num - 1];
      current->theta = new_boundary_theta;
      current->psi_cm = next->psi_cm;
      previous->depth_cm = depth_new;
      current->layer_num = previous->layer_num;
      previous->dzdt_cm_per_h = current->dzdt_cm_per_h;
      current->dzdt_cm_per_h = 0.0;
      current->to_bottom = TRUE;
      previous->to_bottom = FALSE;
      crossed = true;

      double moving_up_residual_cm =
        mass_before_moving_up_crossing -
        lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
      if (std::fabs(moving_up_residual_cm) > MBAL_ITERATIVE_TOLERANCE) {
        if (bounded_moving_up_remap) {
          lgar_theta_mass_balance_correction(false, current->front_num,
                                             mass_before_moving_up_crossing,
                                             head, cum_layer_thickness_cm,
                                             soil_type, soil_properties);
          moving_up_residual_cm =
            mass_before_moving_up_crossing -
            lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
        }

        /*
         * Moving a TO front upward across a layer boundary can change storage
         * when the inherited lower-layer scaffold is hydraulically discontinuous.
         * Keep the geometry repair local and return the remaining signed storage
         * residual so the caller's lower-boundary flux bookkeeping stays closed.
         */
        if (mass_balance_flux_correction_cm != NULL &&
            std::fabs(moving_up_residual_cm) > MBAL_ITERATIVE_TOLERANCE) {
          *mass_balance_flux_correction_cm += moving_up_residual_cm;
        }
        if (verbosity.compare("high") == 0 &&
            std::fabs(moving_up_residual_cm) > MBAL_ITERATIVE_TOLERANCE) {
          printf("TO moving-up layer crossing left residual %.12e cm after local "
                 "remap/repair; routing residual to lower-boundary flux bookkeeping.\n",
                 moving_up_residual_cm);
        }
      }

      if (previous->depth_cm < 0.0 && front_num_with_negative_depth != NULL) {
        *front_num_with_negative_depth = previous->front_num;
      }

      if (verbosity.compare("high") == 0) {
        printf("After layer boundary crossing for TO wetting front, moving up case ... \n");
        listPrint(*head);
      }
    }
  }

  return crossed;
}

extern bool lgar_TO_wetting_fronts_cross_layer_boundary(int num_layers, double *cum_layer_thickness_cm,
							int *soil_type, double *frozen_factor,
							struct wetting_front **head,
							struct soil_properties_ *soil_properties,
							double *mass_balance_flux_correction_cm)
{
  int front_num_with_negative_depth = -1;
  return lgar_TO_wetting_fronts_cross_layer_boundary(&front_num_with_negative_depth, num_layers,
						       cum_layer_thickness_cm, soil_type, frozen_factor,
						       head, soil_properties,
						       mass_balance_flux_correction_cm);
}


// ############################################################################################
/*
  The function handles the basic TO-TO merge-via-depth case from old LGARTO.
  It is intentionally limited to the case where a TO wetting front overtakes the TO wetting front below it.
*/
// ############################################################################################
static bool lgarto_restore_TO_depth_merge_mass_via_bounded_depth(double target_mass,
                                                                 double column_depth,
                                                                 double *cum_layer_thickness_cm,
                                                                 struct wetting_front **head,
                                                                 struct wetting_front *front_to_adjust)
{
  if (front_to_adjust == NULL || head == NULL || *head == NULL ||
      front_to_adjust->layer_num < 1) {
    return false;
  }

  const double initial_depth_cm = front_to_adjust->depth_cm;
  double lower_depth_cm = cum_layer_thickness_cm[front_to_adjust->layer_num - 1];
  double upper_depth_cm = fmin(column_depth, cum_layer_thickness_cm[front_to_adjust->layer_num]);

  struct wetting_front *previous = NULL;
  if (front_to_adjust->front_num > 1) {
    previous = listFindFront(front_to_adjust->front_num - 1, *head, NULL);
  }

  if (previous != NULL && previous->layer_num == front_to_adjust->layer_num) {
    lower_depth_cm = fmax(lower_depth_cm, previous->depth_cm + DEPTH_AVOIDS_SAME_WF_DEPTH);
  }
  if (front_to_adjust->next != NULL &&
      front_to_adjust->next->layer_num == front_to_adjust->layer_num) {
    upper_depth_cm = fmin(upper_depth_cm,
                          front_to_adjust->next->depth_cm - DEPTH_AVOIDS_SAME_WF_DEPTH);
  }

  if (upper_depth_cm <= lower_depth_cm) {
    front_to_adjust->depth_cm = initial_depth_cm;
    return false;
  }

  const double original_depth_cm = fmax(lower_depth_cm, fmin(initial_depth_cm, upper_depth_cm));
  front_to_adjust->depth_cm = original_depth_cm;
  const double original_mass_cm = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
  if (std::fabs(target_mass - original_mass_cm) <= MBAL_ITERATIVE_TOLERANCE) {
    return true;
  }

  front_to_adjust->depth_cm = lower_depth_cm;
  const double mass_at_lower_cm = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
  front_to_adjust->depth_cm = upper_depth_cm;
  const double mass_at_upper_cm = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);

  const bool monotonic_increasing = mass_at_upper_cm > mass_at_lower_cm + MBAL_ITERATIVE_TOLERANCE;
  const bool monotonic_decreasing = mass_at_lower_cm > mass_at_upper_cm + MBAL_ITERATIVE_TOLERANCE;
  if (!monotonic_increasing && !monotonic_decreasing) {
    front_to_adjust->depth_cm = original_depth_cm;
    return false;
  }

  const double min_mass_cm = fmin(mass_at_lower_cm, mass_at_upper_cm);
  const double max_mass_cm = fmax(mass_at_lower_cm, mass_at_upper_cm);
  bool target_bracketed = true;
  if (target_mass <= min_mass_cm) {
    front_to_adjust->depth_cm = (mass_at_lower_cm <= mass_at_upper_cm) ? lower_depth_cm : upper_depth_cm;
    target_bracketed = false;
  }
  else if (target_mass >= max_mass_cm) {
    front_to_adjust->depth_cm = (mass_at_lower_cm >= mass_at_upper_cm) ? lower_depth_cm : upper_depth_cm;
    target_bracketed = false;
  }

  if (!target_bracketed) {
    if (verbosity.compare("high") == 0) {
      printf("TO depth-based merge bounded depth solve could not bracket target mass; "
             "front %d set to closest bounded depth %.17lf cm "
             "(bounds %.17lf..%.17lf cm, mass_at_lower %.17lf, mass_at_upper %.17lf, "
             "target %.17lf).\n",
             front_to_adjust->front_num,
             front_to_adjust->depth_cm,
             lower_depth_cm,
             upper_depth_cm,
             mass_at_lower_cm,
             mass_at_upper_cm,
             target_mass);
    }
    return false;
  }

  double bracket_lo_cm = lower_depth_cm;
  double bracket_hi_cm = upper_depth_cm;
  double best_depth_cm = original_depth_cm;
  double best_abs_residual_cm = std::fabs(target_mass - original_mass_cm);

  for (int iter = 0; iter < 80; iter++) {
    const double probe_depth_cm = 0.5 * (bracket_lo_cm + bracket_hi_cm);
    front_to_adjust->depth_cm = probe_depth_cm;
    const double probe_mass_cm = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
    const double probe_abs_residual_cm = std::fabs(target_mass - probe_mass_cm);

    if (probe_abs_residual_cm < best_abs_residual_cm) {
      best_abs_residual_cm = probe_abs_residual_cm;
      best_depth_cm = probe_depth_cm;
    }

    if (probe_abs_residual_cm <= MBAL_ITERATIVE_TOLERANCE) {
      best_depth_cm = probe_depth_cm;
      break;
    }

    if (monotonic_increasing) {
      if (probe_mass_cm < target_mass) {
        bracket_lo_cm = probe_depth_cm;
      }
      else {
        bracket_hi_cm = probe_depth_cm;
      }
    }
    else {
      if (probe_mass_cm > target_mass) {
        bracket_lo_cm = probe_depth_cm;
      }
      else {
        bracket_hi_cm = probe_depth_cm;
      }
    }
  }

  front_to_adjust->depth_cm = best_depth_cm;
  if (verbosity.compare("high") == 0) {
    printf("TO depth-based merge bounded depth solve updated front %d from %.17lf cm "
           "to %.17lf cm with residual %.12e cm "
           "(bounds %.17lf..%.17lf cm).\n",
           front_to_adjust->front_num,
           initial_depth_cm,
           front_to_adjust->depth_cm,
           target_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, *head),
           lower_depth_cm,
           upper_depth_cm);
  }

  return std::fabs(target_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, *head)) <=
         MBAL_ITERATIVE_TOLERANCE;
}

extern double lgarto_TO_WFs_merge_via_depth(double target_mass, double column_depth,
					    double *cum_layer_thickness_cm, struct wetting_front **head,
					    int *soil_type, struct soil_properties_ *soil_properties)
{
  bool merged = false;
  struct wetting_front *current = *head;
  struct wetting_front *next = current != NULL ? current->next : NULL;

  while (current != NULL && next != NULL) {
    if ((current->depth_cm > next->depth_cm) && current->is_WF_GW && next->is_WF_GW &&
        (next->to_bottom == FALSE) && (current->to_bottom == FALSE) && (current->theta < next->theta)) {
      if (verbosity.compare("high") == 0) {
        printf("Before TO depth-based merging ... \n");
        listPrint(*head);
      }

      next = listDeleteFront(next->front_num, head, soil_type, soil_properties);
      merged = true;

      const int merged_front_num = current->front_num;
      struct saved_front_state {
        int front_num;
        int layer_num;
        int to_bottom;
        int is_WF_GW;
        double depth_cm;
        double theta;
        double dzdt_cm_per_h;
        double K_cm_per_h;
        double psi_cm;
      };

      std::vector<saved_front_state> front_states_after_delete;
      for (struct wetting_front *front = *head; front != NULL; front = front->next) {
        front_states_after_delete.push_back({
          front->front_num,
          front->layer_num,
          front->to_bottom,
          front->is_WF_GW,
          front->depth_cm,
          front->theta,
          front->dzdt_cm_per_h,
          front->K_cm_per_h,
          front->psi_cm
        });
      }

      const double legacy_absurd_depth_limit_cm =
        column_depth > 0.0 ? 10.0 * column_depth : 1.0e4;
      bool legacy_depth_absurd = false;
      double temp_tol = MBAL_ITERATIVE_TOLERANCE;
      double factor = current->depth_cm > column_depth ? 1000.0 : 1.0;
      bool switched = false;
      int iter = 0;
      double current_mass = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);

      while (std::fabs(current_mass - target_mass) > temp_tol) {
        iter++;
        if (iter > 100000) {
          break;
        }

        if (current_mass >= target_mass) {
          current->depth_cm += 0.1 * factor;
          switched = false;
        }
        else {
          if (!switched) {
            switched = true;
            factor = factor * 0.1;
          }
          current->depth_cm -= 0.1 * factor;
        }

        current->depth_cm = fmax(0.0, current->depth_cm);
        current_mass = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);

        if (!std::isfinite(current->depth_cm) || !std::isfinite(current_mass) ||
            current->depth_cm > legacy_absurd_depth_limit_cm) {
          legacy_depth_absurd = true;
          break;
        }
      }

      legacy_depth_absurd =
        legacy_depth_absurd || !std::isfinite(current->depth_cm) ||
        !std::isfinite(current_mass) ||
        current->depth_cm > legacy_absurd_depth_limit_cm;

      bool closed_by_depth =
        std::fabs(current_mass - target_mass) <= MBAL_ITERATIVE_TOLERANCE;
      if (legacy_depth_absurd) {
        if (verbosity.compare("high") == 0) {
          printf("TO depth-based merge legacy depth solve exceeded guard "
                 "(front %d depth %.17lf cm, limit %.17lf cm, residual %.12e cm); "
                 "falling back to bounded depth solve.\n",
                 current->front_num,
                 current->depth_cm,
                 legacy_absurd_depth_limit_cm,
                 target_mass - current_mass);
        }

        for (const saved_front_state &saved : front_states_after_delete) {
          struct wetting_front *front = listFindFront(saved.front_num, *head, NULL);
          if (front == NULL) {
            continue;
          }
          front->layer_num = saved.layer_num;
          front->to_bottom = saved.to_bottom;
          front->is_WF_GW = saved.is_WF_GW;
          front->depth_cm = saved.depth_cm;
          front->theta = saved.theta;
          front->dzdt_cm_per_h = saved.dzdt_cm_per_h;
          front->K_cm_per_h = saved.K_cm_per_h;
          front->psi_cm = saved.psi_cm;
        }

        current = listFindFront(merged_front_num, *head, NULL);
        closed_by_depth =
          lgarto_restore_TO_depth_merge_mass_via_bounded_depth(target_mass, column_depth,
                                                               cum_layer_thickness_cm, head,
                                                               current);
        current_mass = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
      }
      else if (verbosity.compare("high") == 0) {
        printf("TO depth-based merge legacy depth solve accepted front %d at %.17lf cm "
               "with residual %.12e cm after %d iterations "
               "(guard %.17lf cm).\n",
               current->front_num,
               current->depth_cm,
               target_mass - current_mass,
               iter,
               legacy_absurd_depth_limit_cm);
      }

      if (legacy_depth_absurd && !closed_by_depth &&
          std::fabs(current_mass - target_mass) > MBAL_ITERATIVE_TOLERANCE) {
        struct wetting_front *psi_repair_front =
          (current->next != NULL && current->next->to_bottom) ? current->next : current;
        const int front_num_to_repair = psi_repair_front->front_num;
        std::vector<saved_front_state> front_states_before_psi_repair;
        for (struct wetting_front *front = *head; front != NULL; front = front->next) {
          front_states_before_psi_repair.push_back({
            front->front_num,
            front->layer_num,
            front->to_bottom,
            front->is_WF_GW,
            front->depth_cm,
            front->theta,
            front->dzdt_cm_per_h,
            front->K_cm_per_h,
            front->psi_cm
          });
        }
        const double residual_before_psi_repair = target_mass - current_mass;
        if (verbosity.compare("high") == 0) {
          printf("TO depth-based merge bounded depth solve left residual %.12e cm; "
                 "adjusting connected TO/to_bottom chain starting at front %d.\n",
                 residual_before_psi_repair,
                 front_num_to_repair);
        }

        lgar_theta_mass_balance_correction(false, front_num_to_repair, target_mass, head,
                                           cum_layer_thickness_cm, soil_type,
                                           soil_properties);
        current = listFindFront(merged_front_num, *head, NULL);
        current_mass = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
        const double residual_after_psi_repair = target_mass - current_mass;

        if (std::fabs(residual_after_psi_repair) >
            std::fabs(residual_before_psi_repair)) {
          for (const saved_front_state &saved : front_states_before_psi_repair) {
            struct wetting_front *front = listFindFront(saved.front_num, *head, NULL);
            if (front == NULL) {
              continue;
            }
            front->layer_num = saved.layer_num;
            front->to_bottom = saved.to_bottom;
            front->is_WF_GW = saved.is_WF_GW;
            front->depth_cm = saved.depth_cm;
            front->theta = saved.theta;
            front->dzdt_cm_per_h = saved.dzdt_cm_per_h;
            front->K_cm_per_h = saved.K_cm_per_h;
            front->psi_cm = saved.psi_cm;
          }
          current = listFindFront(merged_front_num, *head, NULL);
          current_mass = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
          if (verbosity.compare("high") == 0) {
            printf("TO depth-based merge psi-chain repair worsened residual "
                   "(before %.12e cm, after %.12e cm); restored bounded-depth state.\n",
                   residual_before_psi_repair,
                   residual_after_psi_repair);
          }
        }

        if (verbosity.compare("high") == 0) {
          printf("TO depth-based merge psi-chain repair residual %.12e cm\n",
                 target_mass - current_mass);
        }
      }

	      if (verbosity.compare("high") == 0) {
	        printf("After TO depth-based merging ... \n");
	        listPrint(*head);
	      }

      if (current == NULL) {
        break;
      }
      next = current->next;
    }

    current = current->next;
    next = current != NULL ? current->next : NULL;
  }

  if (!merged) {
    return 0.0;
  }

  return target_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
}

// ############################################################################################
/*
  The function handles the basic surface-TO merge case from old LGARTO.
  It is intentionally kept as a focused primitive rather than part of the full old correction loop.
*/
// ############################################################################################

static void lgarto_rehome_zero_depth_front_to_top_layer_from_psi(struct wetting_front *front,
                                                                 int *soil_type,
                                                                 double *frozen_factor,
                                                                 struct soil_properties_ *soil_properties)
{
  if (front == NULL || front->to_bottom == TRUE) {
    return;
  }

  front->layer_num = 1;

  const int soil_num = soil_type[front->layer_num];
  const double theta_e = soil_properties[soil_num].theta_e;
  const double theta_r = soil_properties[soil_num].theta_r;
  const double vg_a = soil_properties[soil_num].vg_alpha_per_cm;
  const double vg_m = soil_properties[soil_num].vg_m;
  const double vg_n = soil_properties[soil_num].vg_n;
  const double Ksat_cm_per_h =
    soil_properties[soil_num].Ksat_cm_per_h * frozen_factor[front->layer_num];

  front->theta = calc_theta_from_h(front->psi_cm, vg_a, vg_m, vg_n, theta_e, theta_r);
  const double Se = calc_Se_from_theta(front->theta, theta_e, theta_r);
  front->K_cm_per_h = calc_K_from_Se(Se, Ksat_cm_per_h, vg_m);
}

static double lgarto_repair_negative_depth_fronts_to_lower_boundary_flux(
  int num_layers, double *cum_layer_thickness_cm, int *soil_type,
  double *frozen_factor, struct soil_properties_ *soil_properties,
  struct wetting_front **head, const char *context)
{
  if (head == NULL || *head == NULL || cum_layer_thickness_cm == NULL) {
    return 0.0;
  }

  const double mass_before_cm = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
  const double numerical_negative_depth_tol_cm = 1.0e-6;
  bool repaired_any = false;

  for (struct wetting_front *current = *head; current != NULL;) {
    if (current->depth_cm >= 0.0) {
      current = current->next;
      continue;
    }

    repaired_any = true;
    const int front_num = current->front_num;
    const double old_depth_cm = current->depth_cm;
    const bool is_real_negative_depth =
      old_depth_cm < -numerical_negative_depth_tol_cm;

    if (current->to_bottom == TRUE) {
      if (current->layer_num >= 1 && current->layer_num <= num_layers) {
        current->depth_cm = cum_layer_thickness_cm[current->layer_num];
      }
      else {
        current->depth_cm = 0.0;
      }

      if (verbosity.compare("high") == 0) {
        printf("Snapped negative-depth to_bottom front %d from %.17lf cm to %.17lf cm "
               "during %s.\n",
               front_num, old_depth_cm, current->depth_cm,
               context != NULL ? context : "negative-depth repair");
      }
      current = current->next;
      continue;
    }

    if (is_real_negative_depth) {
      if (verbosity.compare("high") == 0) {
        printf("Deleting non-to_bottom front %d with negative depth %.17lf cm "
               "during %s; storage change will be routed to lower-boundary flux.\n",
               front_num, old_depth_cm,
               context != NULL ? context : "negative-depth repair");
      }
      current = listDeleteFront(front_num, head, soil_type, soil_properties);
      continue;
    }

    current->depth_cm = 0.0;
    lgarto_rehome_zero_depth_front_to_top_layer_from_psi(current, soil_type,
                                                         frozen_factor,
                                                         soil_properties);
    if (verbosity.compare("high") == 0) {
      printf("Clipped tiny negative-depth front %d from %.17lf cm to zero "
             "during %s.\n",
             front_num, old_depth_cm,
             context != NULL ? context : "negative-depth repair");
    }
    current = current->next;
  }

  if (!repaired_any) {
    return 0.0;
  }

  const double mass_after_cm = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
  const double lower_boundary_flux_correction_cm = mass_before_cm - mass_after_cm;
  if (verbosity.compare("high") == 0) {
    printf("Negative-depth wetting-front repair during %s routed %.17lf cm "
           "to lower-boundary flux (mass_before %.17lf, mass_after %.17lf).\n",
           context != NULL ? context : "negative-depth repair",
           lower_boundary_flux_correction_cm, mass_before_cm, mass_after_cm);
  }

  return lower_boundary_flux_correction_cm;
}

static double lgarto_surface_fronts_cross_layer_boundary_upward(
  int num_layers, double *cum_layer_thickness_cm, int *soil_type,
  double *frozen_factor, struct soil_properties_ *soil_properties,
  struct wetting_front **head, const char *context)
{
  if (head == NULL || *head == NULL || cum_layer_thickness_cm == NULL) {
    return 0.0;
  }

  const double mass_before_cm = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
  const double boundary_tol_cm = fmax(1.0e-8, 10.0 * TRUNCATION_DEPTH);
  bool repaired_any = false;

  for (struct wetting_front *current = *head; current != NULL; current = current->next) {
    if (current->is_WF_GW || current->to_bottom || current->layer_num <= 1 ||
        current->layer_num > num_layers || current->depth_cm < 0.0) {
      continue;
    }

    const double upper_boundary_cm = cum_layer_thickness_cm[current->layer_num - 1];
    if (current->depth_cm >= upper_boundary_cm - boundary_tol_cm) {
      continue;
    }

    const int old_layer = current->layer_num;
    const double old_depth_cm = current->depth_cm;
    const double old_psi_cm = current->psi_cm;
    int new_layer = old_layer;
    while (new_layer > 1 &&
           current->depth_cm < cum_layer_thickness_cm[new_layer - 1] - boundary_tol_cm) {
      new_layer--;
    }

    current->layer_num = new_layer;
    current->to_bottom = FALSE;
    current->is_WF_GW = FALSE;
    lgarto_update_front_from_shared_psi(current, old_psi_cm, num_layers, soil_type,
                                        frozen_factor, soil_properties);
    repaired_any = true;

    if (verbosity.compare("high") == 0) {
      printf("Surface upward layer crossing during %s moved front %d "
             "from layer %d to layer %d at depth %.17lf cm "
             "(psi %.17lf cm).\n",
             context != NULL ? context : "surface upward layer crossing",
             current->front_num, old_layer, current->layer_num, old_depth_cm,
             old_psi_cm);
    }

    break;
  }

  if (!repaired_any) {
    return 0.0;
  }

  listSortFrontsByDepth(*head);

  const double mass_after_cm = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
  const double lower_boundary_flux_correction_cm = mass_before_cm - mass_after_cm;
  if (verbosity.compare("high") == 0) {
    printf("Surface upward layer-crossing repair during %s routed %.17lf cm "
           "to lower-boundary flux (mass_before %.17lf, mass_after %.17lf).\n",
           context != NULL ? context : "surface upward layer crossing",
           lower_boundary_flux_correction_cm, mass_before_cm, mass_after_cm);
  }

  return lower_boundary_flux_correction_cm;
}

extern bool lgar_merge_surface_and_TO_wetting_fronts(bool merged_in_non_top_layer, int num_layers,
						     double *cum_layer_thickness_cm, int *soil_type,
						     double *frozen_factor, struct soil_properties_ *soil_properties,
						     struct wetting_front **head,
						     bool latch_surface_state_recipient_dzdt)
{
  bool merged_any = false;

  for (int wf = 1; wf != listLength(*head); wf++) {
    struct wetting_front *current = listFindFront(wf, *head, NULL);
    struct wetting_front *next = current != NULL ? current->next : NULL;

    if (current == NULL || next == NULL || next->next == NULL) {
      continue;
    }

    if ((current->depth_cm > next->depth_cm) && (current->is_WF_GW == FALSE) &&
        (next->is_WF_GW == TRUE) && (next->to_bottom == FALSE)) {
      if (!merged_any && verbosity.compare("high") == 0) {
        printf("surface-TO merging start\n");
        printf("current mass: %.17lf \n", lgar_calc_mass_bal(cum_layer_thickness_cm, *head));
      }
      merged_any = true;

      if (current->layer_num > 1) {
        merged_in_non_top_layer = true;
      }
      if (num_layers == 1) {
        merged_in_non_top_layer = false;
      }

      struct wetting_front *surface_state_recipient = NULL;

      if (current->theta > next->next->theta) {
        if (verbosity.compare("high") == 0) {
          printf("surface-TO merging, surface WF wetter case...\n");
          listPrint(*head);
        }

        const double theta_next_temp = next->theta;
        const int to_bottom_next_temp = next->to_bottom;
        const double dzdt_cm_per_h_next_temp = next->dzdt_cm_per_h;
        const double K_cm_per_h_next_temp = next->K_cm_per_h;
        const double psi_cm_next_temp = next->psi_cm;

        next->depth_cm = current->depth_cm +
                         (current->depth_cm - next->depth_cm) * (next->next->theta - next->theta) /
                           std::fabs(current->theta - next->next->theta);

        next->theta = current->theta;
        next->to_bottom = current->to_bottom;
        next->dzdt_cm_per_h = current->dzdt_cm_per_h;
        next->K_cm_per_h = current->K_cm_per_h;
        next->is_WF_GW = FALSE;
        next->psi_cm = current->psi_cm;
        surface_state_recipient = next;

        current->depth_cm = 0.0;
        current->theta = theta_next_temp;
        current->to_bottom = to_bottom_next_temp;
        current->dzdt_cm_per_h = dzdt_cm_per_h_next_temp;
        current->K_cm_per_h = K_cm_per_h_next_temp;
        current->is_WF_GW = TRUE;
        current->psi_cm = psi_cm_next_temp;
        lgarto_rehome_zero_depth_front_to_top_layer_from_psi(current, soil_type, frozen_factor,
                                                             soil_properties);
      }
      else {
        if (verbosity.compare("high") == 0) {
          printf("surface-TO merging, surface WF drier case...\n");
          listPrint(*head);
        }

        current->is_WF_GW = TRUE;
        double overshot_depth = current->depth_cm - next->depth_cm;
        double local_remap_depth_cm = next->depth_cm - overshot_depth;
        double remapped_depth_cm = local_remap_depth_cm;
        double theta_ratio = 1.0;
        bool use_mass_weighted_remap = false;

        const double theta_denominator = std::fabs(next->next->theta - current->theta);
        if (theta_denominator > 1.0e-12) {
          theta_ratio = (current->theta - next->theta) / theta_denominator;
          const double ratio_remap_depth_cm = next->depth_cm - overshot_depth * theta_ratio;
          if (std::isfinite(ratio_remap_depth_cm)) {
            remapped_depth_cm = ratio_remap_depth_cm;
            use_mass_weighted_remap = true;
          }
        }

        if (remapped_depth_cm < 0.0) {
          bool next_is_dry_zero_depth_TO_support = false;
          double next_effective_saturation = 1.0;
          if (next->depth_cm <= TRUNCATION_DEPTH &&
              next->layer_num >= 1 && next->layer_num <= num_layers) {
            const int zero_depth_soil_num = soil_type[next->layer_num];
            const double theta_r = soil_properties[zero_depth_soil_num].theta_r;
            const double theta_e = soil_properties[zero_depth_soil_num].theta_e;
            const double theta_range = theta_e - theta_r;
            if (theta_range > 0.0) {
              next_effective_saturation = (next->theta - theta_r) / theta_range;
              next_is_dry_zero_depth_TO_support =
                std::isfinite(next_effective_saturation) &&
                next_effective_saturation <= ZERO_DEPTH_TO_DRY_SUPPORT_MAX_SE;
            }
          }

          if (next_is_dry_zero_depth_TO_support) {
            const double prior_mass = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
            if (verbosity.compare("high") == 0) {
              printf("Removing zero-depth TO front %d during surface/TO drier merge "
                     "because the local mass remap would be above the surface "
                     "and the TO support is near residual saturation "
                     "(surface front %d depth %.17lf cm, remapped_depth %.17lf cm, "
                     "TO effective saturation %.17lf).\n",
                     next->front_num,
                     current->front_num,
                     current->depth_cm,
                     remapped_depth_cm,
                     next_effective_saturation);
            }

            current->is_WF_GW = FALSE;
            (void) listDeleteFront(next->front_num, head, soil_type, soil_properties);

            if (current->next != NULL &&
                current->next->is_WF_GW == FALSE &&
                current->layer_num == current->next->layer_num &&
                current->theta <= current->next->theta) {
              struct wetting_front *repair_front =
                listDeleteFront(current->front_num, head, soil_type, soil_properties);
              if (repair_front != NULL) {
                lgar_theta_mass_balance_correction(true, repair_front->front_num, prior_mass,
                                                   head, cum_layer_thickness_cm, soil_type,
                                                   soil_properties);
              }
            }
            break;
          }
          if (verbosity.compare("high") == 0) {
            printf("Event-limited surface/TO drier merge for front %d from depth %.17lf cm "
                   "to %.17lf cm to keep local remap nonnegative across TO front %d at %.17lf cm.\n",
                   current->front_num,
                   current->depth_cm,
                   next->depth_cm + TRUNCATION_DEPTH,
                   next->front_num,
                   next->depth_cm);
          }
          if (use_mass_weighted_remap && theta_ratio > 0.0) {
            const double capped_overshot_depth =
              fmax(0.0, (next->depth_cm - TRUNCATION_DEPTH) / theta_ratio);
            current->depth_cm = next->depth_cm + capped_overshot_depth;
          }
          else {
            current->depth_cm = next->depth_cm + TRUNCATION_DEPTH;
          }
          overshot_depth = current->depth_cm - next->depth_cm;
          remapped_depth_cm = use_mass_weighted_remap ?
            next->depth_cm - overshot_depth * theta_ratio :
            next->depth_cm - overshot_depth;
        }

        current->depth_cm = next->depth_cm;

        // Prefer the mass-weighted remap; it converts the geometric overshoot
        // into the water-volume correction implied by the local theta contrasts.
        next->depth_cm = remapped_depth_cm;

        const double theta_next_temp = next->theta;
        const double dzdt_cm_per_h_next_temp = next->dzdt_cm_per_h;
        const double K_cm_per_h_next_temp = next->K_cm_per_h;
        const double psi_cm_next_temp = next->psi_cm;

        next->theta = current->theta;
        next->dzdt_cm_per_h = current->dzdt_cm_per_h;
        next->K_cm_per_h = current->K_cm_per_h;
        next->psi_cm = current->psi_cm;
        // next->is_WF_GW = FALSE;
        surface_state_recipient = next;

        current->depth_cm = 0.0;
        current->theta = theta_next_temp;
        current->dzdt_cm_per_h = dzdt_cm_per_h_next_temp;
        current->K_cm_per_h = K_cm_per_h_next_temp;
        current->psi_cm = psi_cm_next_temp;
        lgarto_rehome_zero_depth_front_to_top_layer_from_psi(current, soil_type, frozen_factor,
                                                             soil_properties);
      }

      // Resolve one surface/TO inversion per correction-loop pass. The caller
      // recomputes correction_type before deciding whether another merge is needed.
      // When the pre-TO-motion surface loop converts a just-moved surface
      // state into a mobile TO front, do not let that same state move again as
      // TO in this substep. The next dzdt calculation will refresh it.
      if (latch_surface_state_recipient_dzdt &&
          surface_state_recipient != NULL &&
          surface_state_recipient->is_WF_GW == TRUE &&
          surface_state_recipient->to_bottom == FALSE) {
        surface_state_recipient->dzdt_cm_per_h = 0.0;
      }
      break;
    }
  }

  for (int wf = listLength(*head) - 1; wf != 0; wf--) {
    struct wetting_front *current = listFindFront(wf, *head, NULL);
    struct wetting_front *next = current != NULL ? current->next : NULL;
    if (current != NULL && next != NULL && current->to_bottom == TRUE &&
        next->is_WF_GW != current->is_WF_GW) {
      current->is_WF_GW = next->is_WF_GW;
    }
  }

  listSendToTop(*head);

  if (merged_any && verbosity.compare("high") == 0) {
    printf("surface-TO merging end\n");
    printf("current mass: %.17lf \n", lgar_calc_mass_bal(cum_layer_thickness_cm, *head));
  }

  return merged_in_non_top_layer;
}

// ############################################################################################
/*
  Cleanup after surface/TO merging below the top layer. This preserves layer bookkeeping and
  updates theta from psi for fronts whose identity changed during the merge.
*/
// ############################################################################################

extern void lgarto_cleanup_after_surface_TO_merging_in_layer_below_top(bool merged_in_non_top_layer,
								       int *soil_type,
								       struct soil_properties_ *soil_properties,
								       struct wetting_front **head)
{
  if (verbosity.compare("high") == 0) {
    printf("before lgarto_cleanup_after_surface_TO_merging_in_layer_below_top: \n");
    listPrint(*head);
  }

  if (merged_in_non_top_layer) {
    listSendToTop(*head);

    for (int wf = listLength(*head) - 1; wf > 0; wf--) {
      struct wetting_front *current = listFindFront(wf, *head, NULL);
      struct wetting_front *next = current != NULL ? current->next : NULL;
      if (current == NULL || next == NULL) {
        continue;
      }

      if (current->to_bottom == TRUE) {
        current->is_WF_GW = next->is_WF_GW;
      }

      if ((next->to_bottom == FALSE) && (current->to_bottom == FALSE) &&
          (current->layer_num != next->layer_num)) {
        current->layer_num = next->layer_num;
      }

      if ((current->layer_num > next->layer_num) && (next->to_bottom == TRUE)) {
        current->layer_num = next->layer_num;
      }

      const int soil_num_k = soil_type[current->layer_num];
      const double theta_e_k = soil_properties[soil_num_k].theta_e;
      const double theta_r_k = soil_properties[soil_num_k].theta_r;
      const double vg_a_k = soil_properties[soil_num_k].vg_alpha_per_cm;
      const double vg_m_k = soil_properties[soil_num_k].vg_m;
      const double vg_n_k = soil_properties[soil_num_k].vg_n;
      current->theta = calc_theta_from_h(current->psi_cm, vg_a_k, vg_m_k, vg_n_k, theta_e_k, theta_r_k);

      if ((current->to_bottom == TRUE) &&
          (current->layer_num != next->layer_num) && (current->psi_cm != next->psi_cm)) {
        current->psi_cm = next->psi_cm;
        current->theta = calc_theta_from_h(current->psi_cm, vg_a_k, vg_m_k, vg_n_k, theta_e_k, theta_r_k);
      }

      // if ((current->to_bottom == TRUE) && (next->to_bottom == TRUE) &&
      //     (current->layer_num != next->layer_num) && (current->psi_cm != next->psi_cm)) {
      //   current->psi_cm = next->psi_cm;
      //   current->theta = calc_theta_from_h(current->psi_cm, vg_a_k, vg_m_k, vg_n_k, theta_e_k, theta_r_k);
      // }

      // if ((current->to_bottom == TRUE) && (next->to_bottom == FALSE) && (current->layer_num != next->layer_num) && (current->psi_cm != next->psi_cm)){
      //   current->psi_cm = next->psi_cm;
      //   current->theta = calc_theta_from_h(current->psi_cm, vg_a_k, vg_m_k, vg_n_k, theta_e_k, theta_r_k);
      // }

    }
  }

  if (verbosity.compare("high") == 0) {
    printf("after lgarto_cleanup_after_surface_TO_merging_in_layer_below_top: \n");
    listPrint(*head);
  }
}


// ############################################################################################
/*
  the function lets wetting fronts of a sufficient depth interact with the lower boundary; called from lgar_move_wetting_fronts.
  For mobile groundwater, the current GW depth is an exchange boundary, not a hard
  surface-front clipping plane. Surface fronts that reach it are handled by the
  same merge/flux bookkeeping as lower-boundary crossing, which avoids parking a
  surface profile at the water table.
*/
// ############################################################################################

extern double lgar_wetting_front_cross_domain_boundary(double domain_depth_cm, int *soil_type,
						       double *frozen_factor, struct wetting_front** head,
						       struct soil_properties_ *soil_properties)
{
  struct wetting_front *current;
  struct wetting_front *next;
  struct wetting_front *next_to_next;
  current = *head; 
  double bottom_flux_cm = 0.0;
  int length = listLength(*head);
  
  if (verbosity.compare("high") == 0) {
    printf("Domain boundary crossing (bottom flux calc.) \n");
  }

  // local variables
  double theta_e,theta_r;
  double vg_a, vg_m, vg_n;
  double bottom_flux_cm_temp;
  int layer_num, soil_num;
    
  for (int wf=1; wf != length; wf++) {

    if (verbosity.compare("high") == 0) {
      printf("Domain boundary crossing | ***** Wetting Front = %d ****** \n", wf);
    }

    // ensure that loop iterations never exceed the total number of wetting fronts after altering the list
    if (wf >= listLength(*head))
      break;
    
    bottom_flux_cm_temp = 0.0;
    
    layer_num   = current->layer_num;
    soil_num    = soil_type[layer_num];
    
    next = current->next;
    next_to_next = current->next->next;
    
    // case : wetting front is the deepest one in the last layer (most deepested wetting front in the domain)
    /**********************************************************/
    if (next_to_next == NULL && current->depth_cm >= domain_depth_cm) {
      //  this is the water exchanged through the bottom/GW boundary
      bottom_flux_cm_temp = (current->theta - next->theta) *  (current->depth_cm - next->depth_cm);
      theta_e   = soil_properties[soil_num].theta_e;
      theta_r   = soil_properties[soil_num].theta_r;
      vg_a      = soil_properties[soil_num].vg_alpha_per_cm;
      vg_m      = soil_properties[soil_num].vg_m;
      vg_n      = soil_properties[soil_num].vg_n;
      double Ksat_cm_per_h  = soil_properties[soil_num].Ksat_cm_per_h * frozen_factor[current->layer_num];

      next->theta = current->theta;
      double Se_k = calc_Se_from_theta(current->theta,theta_e,theta_r);
      next->psi_cm = calc_h_from_Se(Se_k, vg_a, vg_m, vg_n);
      next->K_cm_per_h = calc_K_from_Se(Se_k, Ksat_cm_per_h, vg_m);
      current = listDeleteFront(current->front_num, head, soil_type, soil_properties);
      bottom_flux_cm += bottom_flux_cm_temp; 
      break;
    }
    
    current = current->next;

    if (bottom_flux_cm_temp!=0.0){//PTL: perhaps due to potential problems with precision error, this should just be replaces with a flag that gets toggled when bottom_boundary_flux_cm is changed from 0 to nonzero
      break;
    }
  }

  if (verbosity.compare("high") == 0) {
    printf("State after lowest wetting front contributes to flux through the bottom boundary...\n");
    listPrint(*head);
    printf("Bottom boundary flux = %lf \n",bottom_flux_cm);
  }

  return bottom_flux_cm;

}


// ############################################################################################
/* The function handles situation of dry over wet wetting fronts
  mainly happen when AET extracts more water from the upper wetting front
  and the front gets drier than the lower wetting front, or if free 
  drainage is enabled it can do the same thing */
// ############################################################################################
extern void lgar_fix_dry_over_wet_wetting_fronts(double *mass_change, double* cum_layer_thickness_cm, int *soil_type,
					 struct wetting_front** head, struct soil_properties_ *soil_properties)
{
  // This function will delete the wetting front that is drier than the WF below it that is in the same layer, and then it will 
  // iteratively adjust the psi and theta values of the region of the soil column that should have just 1 psi value now that a WF was deleted.
  // mass_change will return any mass balance error and should usually be 0.0
  if (verbosity.compare("high") == 0) {
    printf("Fix Dry over Wet Wetting Front (before) ... \n");
    listPrint(*head);
  }

  struct wetting_front *current;
  struct wetting_front *next;
  current = *head;
  next = current != NULL ? current->next : NULL;

  while (current != NULL && next != NULL) {
    if ((current->is_WF_GW == 0) && (next->is_WF_GW == 0)) {
      // this part fixes case of upper theta less than lower theta due to AET extraction or free drainage
      // also handles the case when the current and next wetting fronts have the same theta
      // and are within the same layer
      /***************************************************/

      if ((current->theta <= next->theta) && (current->layer_num == next->layer_num)) {
        double prior_mass = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
        current = listDeleteFront(current->front_num, head, soil_type, soil_properties); // current will be the WF directly after the one that got deleted
        if (current == NULL) {
          break;
        }

        int front_num_correction = current->front_num;
        bool use_dry_over_wet = true;
        lgar_theta_mass_balance_correction(use_dry_over_wet, front_num_correction, prior_mass, head, cum_layer_thickness_cm,
                                           soil_type, soil_properties);
        double mass_after = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
        *mass_change += (mass_after - prior_mass);

        next = current->next;
        continue;
      }
    }

    current = current->next;
    next = current != NULL ? current->next : NULL;
  }
if (verbosity.compare("high") == 0) {
  printf("Fix Dry over Wet Wetting Front (after) ... \n");
  listPrint(*head);
}
}

// ############################################################################################
/* The function handles situation of dry over wet wetting fronts
  mainly happen when AET extracts more water from the upper wetting front
  and the front gets drier than the lower wetting front */
// ############################################################################################
extern bool lgar_check_dry_over_wet_wetting_fronts(struct wetting_front* head)
{
  struct wetting_front *current = head;
  struct wetting_front *next    = current->next;
  int length = listLength(head);
  
  for (int l=1; l <= length; l++) {
    if (next != NULL) {
      
      if ( (current->theta <= next->theta) && (current->layer_num == next->layer_num) && (current->is_WF_GW==FALSE) && (next->is_WF_GW==FALSE) )
	return true;
      
      current = current->next;
      
      if (current == NULL)
	next = NULL;
      else
	next = current->next;
      
    }
    
  }
  
  return false;
}

// Estimate the positive same-substep TO/GW lower-boundary flux implied by the
// current dzdt values. This is intentionally read-only: the actual movement,
// correction, and recharge bookkeeping still happen in lgar_move_wetting_fronts.
extern double lgarto_project_TO_motion_lower_boundary_flux_cm(double timestep_h, int num_layers,
								      double *cum_layer_thickness_cm, int *soil_type,
								      struct wetting_front* head,
								      struct soil_properties_ *soil_properties,
								      double groundwater_depth_cm)
{
  if (head == NULL || timestep_h <= 0.0 || num_layers <= 0) {
    return 0.0;
  }

  const double fixed_column_depth_cm = cum_layer_thickness_cm[num_layers];
  const double column_depth_cm =
    (std::isfinite(groundwater_depth_cm) && groundwater_depth_cm > 0.0)
        ? groundwater_depth_cm
        : fixed_column_depth_cm;
  if (column_depth_cm <= 0.0) {
    return 0.0;
  }

  double projected_lower_boundary_flux_cm = 0.0;
  const int wetting_front_count = listLength(head);
  for (int wf = wetting_front_count - 1; wf >= 1; wf--) {
    struct wetting_front *current = listFindFront(wf, head, NULL);
    if (current == NULL || !current->is_WF_GW || current->to_bottom) {
      continue;
    }

    struct wetting_front *next_to_use = current->next;
    while (next_to_use != NULL && !next_to_use->is_WF_GW) {
      next_to_use = next_to_use->next;
    }

    if (next_to_use == NULL) {
      continue;
    }

    double delta_depth_cm = current->dzdt_cm_per_h * timestep_h;
    if (delta_depth_cm <= 0.0) {
      continue;
    }

    double dzdt_cap_boundary_cm = column_depth_cm;
    double allowed_boundary_overshoot_cm = 0.001 * column_depth_cm;
    if (current->layer_num >= 1 && current->layer_num <= num_layers) {
      dzdt_cap_boundary_cm = fmin(column_depth_cm, cum_layer_thickness_cm[current->layer_num]);
    }
    if (current->layer_num >= 1 && current->layer_num < num_layers) {
      const double next_layer_thickness_cm =
        cum_layer_thickness_cm[current->layer_num + 1] - cum_layer_thickness_cm[current->layer_num];
      allowed_boundary_overshoot_cm = 0.1 * next_layer_thickness_cm;
    }

    const double max_TO_dzdt_depth_cm = dzdt_cap_boundary_cm + allowed_boundary_overshoot_cm;
    if (current->depth_cm + delta_depth_cm > max_TO_dzdt_depth_cm) {
      delta_depth_cm = fmax(0.0, max_TO_dzdt_depth_cm - current->depth_cm);
    }

    if (delta_depth_cm <= 0.0) {
      continue;
    }

    double delta_theta = 0.0;
    if (current->layer_num == next_to_use->layer_num) {
      delta_theta = next_to_use->theta - current->theta;
    }
    else if (current->layer_num >= 1 && current->layer_num <= num_layers) {
      const int soil_num_current = soil_type[current->layer_num];
      const double theta_e_current = soil_properties[soil_num_current].theta_e;
      const double theta_r_current = soil_properties[soil_num_current].theta_r;
      const double vg_a_current = soil_properties[soil_num_current].vg_alpha_per_cm;
      const double vg_m_current = soil_properties[soil_num_current].vg_m;
      const double vg_n_current = soil_properties[soil_num_current].vg_n;
      const double equiv_next_theta =
        calc_theta_from_h(next_to_use->psi_cm, vg_a_current, vg_m_current, vg_n_current,
                          theta_e_current, theta_r_current);
      delta_theta = equiv_next_theta - current->theta;
    }

    const double projected_flux_cm = delta_depth_cm * delta_theta;
    if (std::isfinite(projected_flux_cm) && projected_flux_cm > 0.0) {
      projected_lower_boundary_flux_cm += projected_flux_cm;
    }
  }

  return fmax(0.0, projected_lower_boundary_flux_cm);
}

static void lgarto_set_front_hydraulic_state_from_psi(struct wetting_front *front,
						      double psi_cm,
						      int *soil_type,
						      double *frozen_factor,
						      struct soil_properties_ *soil_properties)
{
  if (front == NULL || soil_type == NULL || soil_properties == NULL ||
      front->layer_num < 1) {
    return;
  }

  if (!std::isfinite(psi_cm) || psi_cm < 0.0) {
    psi_cm = 0.0;
  }

  const int soil_num = soil_type[front->layer_num];
  const double theta_e = soil_properties[soil_num].theta_e;
  const double theta_r = soil_properties[soil_num].theta_r;
  const double vg_a = soil_properties[soil_num].vg_alpha_per_cm;
  const double vg_m = soil_properties[soil_num].vg_m;
  const double vg_n = soil_properties[soil_num].vg_n;
  const double frozen_multiplier =
    frozen_factor != NULL ? frozen_factor[front->layer_num] : 1.0;
  const double Ksat_cm_per_h = frozen_multiplier * soil_properties[soil_num].Ksat_cm_per_h;

  front->psi_cm = psi_cm;
  front->theta = calc_theta_from_h(psi_cm, vg_a, vg_m, vg_n, theta_e, theta_r);
  const double Se = calc_Se_from_theta(front->theta, theta_e, theta_r);
  front->K_cm_per_h = calc_K_from_Se(Se, Ksat_cm_per_h, vg_m);
}

extern double lgarto_submerge_wetting_fronts_below_groundwater(double groundwater_depth_cm,
							      int num_layers,
							      double *cum_layer_thickness_cm,
							      int *soil_type,
							      double *frozen_factor,
							      struct wetting_front **head,
							      struct soil_properties_ *soil_properties)
{
  if (head == NULL || *head == NULL || num_layers <= 0 ||
      cum_layer_thickness_cm == NULL || soil_type == NULL ||
      soil_properties == NULL ||
      !std::isfinite(groundwater_depth_cm) || groundwater_depth_cm < 0.0) {
    return 0.0;
  }

  const double fixed_column_depth_cm = cum_layer_thickness_cm[num_layers];
  if (!std::isfinite(fixed_column_depth_cm) || fixed_column_depth_cm <= 0.0 ||
      groundwater_depth_cm >= fixed_column_depth_cm + LOWER_BOUNDARY_FINAL_TOL_CM) {
    return 0.0;
  }

  const double submerged_depth_tol_cm =
    fmax(MOBILE_GROUNDWATER_SUBMERGENCE_TOL_CM,
         1.0e-10 * fmax(1.0, fixed_column_depth_cm));
  const double mass_before_cm = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
  bool changed = true;
  int iterations = 0;

  while (changed && iterations < MAX_NUM_WETTING_FRONTS) {
    changed = false;
    iterations++;

    for (int wf = listLength(*head) - 1; wf >= 1; wf--) {
      struct wetting_front *current = listFindFront(wf, *head, NULL);
      if (current == NULL || current->next == NULL ||
          current->is_WF_GW == FALSE ||
          current->to_bottom == TRUE ||
          current->next->is_WF_GW == FALSE ||
          !std::isfinite(current->depth_cm) ||
          current->depth_cm < groundwater_depth_cm - submerged_depth_tol_cm ||
          current->depth_cm <= ZERO_DEPTH_TO_DELETE_DEPTH_TOL_CM) {
        continue;
      }

      const int deleted_front_num = current->front_num;
      const double deleted_depth_cm = current->depth_cm;
      const double deleted_psi_cm = current->psi_cm;
      struct wetting_front *recipient = current->next;

      lgarto_set_front_hydraulic_state_from_psi(recipient, deleted_psi_cm,
                                                soil_type, frozen_factor,
                                                soil_properties);
      recipient->is_WF_GW = TRUE;

      if (verbosity.compare("high") == 0) {
        printf("Mobile groundwater submerged TO/GW front %d at %.17lf cm "
               "(GW depth %.17lf cm); merging it into front %d with psi %.17lf cm.\n",
               deleted_front_num,
               deleted_depth_cm,
               groundwater_depth_cm,
               recipient->front_num,
               deleted_psi_cm);
      }

      (void) listDeleteFront(deleted_front_num, head, soil_type, soil_properties);
      changed = true;
      break;
    }
  }

  if (iterations >= MAX_NUM_WETTING_FRONTS) {
    fprintf(stderr,
            "Error: mobile groundwater submergence repair exceeded iteration limit.\n"
            "  groundwater_depth_cm=%.17g fixed_column_depth_cm=%.17g\n"
            "  Wetting front list follows:\n",
            groundwater_depth_cm,
            fixed_column_depth_cm);
    fflush(stderr);
    listPrint(*head);
    fflush(stdout);
    abort();
  }

  const double mass_after_cm = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
  const double lower_boundary_flux_correction_cm = mass_before_cm - mass_after_cm;
  if (!std::isfinite(lower_boundary_flux_correction_cm)) {
    fprintf(stderr,
            "Error: mobile groundwater submergence repair produced non-finite flux correction.\n"
            "  mass_before_cm=%.17g mass_after_cm=%.17g groundwater_depth_cm=%.17g\n"
            "  Wetting front list follows:\n",
            mass_before_cm,
            mass_after_cm,
            groundwater_depth_cm);
    fflush(stderr);
    listPrint(*head);
    fflush(stdout);
    abort();
  }

  return lower_boundary_flux_correction_cm;
}

static double lgar_saturated_pressure_head_cm(const struct wetting_front *front,
                                              int *soil_type,
                                              struct soil_properties_ *soil_properties)
{
  (void)front;
  (void)soil_type;
  (void)soil_properties;

  /* Disabled for now: theta_e fronts in LGAR/LGARTO are often ordinary
     Green-Ampt wetted zones rather than confirmed positive-pressure saturated
     columns. Ponded head h_p remains the only explicit pressure-head addition
     until we define a conservative connected saturated-column criterion. */
  return 0.0;
}

// ############################################################################################
/* The module computes the potential infiltration capacity, fp (in the lgar manuscript),
   potential infiltration capacity = the maximum amount of water that can be inserted into
   the soil depending on the availability of water.
   this module is called when a new superficial wetting front is not created
   in the current timestep, that is precipitation in the current and previous
   timesteps was greater than zero */
// ############################################################################################
extern double lgar_insert_water(bool use_closed_form_G, int nint, double timestep_h, double AET_demand_cm, double free_drainage_subtimestep_cm, double *ponded_depth_cm,
				double *volin_this_timestep, double precip_timestep_cm, int wf_free_drainage_demand,
			        int num_layers, double ponded_depth_max_cm, int *soil_type,
				double *cum_layer_thickness_cm, double *frozen_factor,
				struct wetting_front* head, struct soil_properties_ *soil_properties,
				double *raw_fp_cm_per_h,
				double *storage_limit_fp_cm_per_h,
				double *capped_fp_cm_per_h)
{
  // note ponded_depth_cm is a pointer.   Access its value as (*ponded_depth_cm).
  int wf_that_supplies_free_drainage_demand = wf_free_drainage_demand;
  // local vars
  double theta_e, theta_r;
  double vg_a, vg_m, vg_n,Ksat_cm_per_h;
  double h_min_cm;
  struct wetting_front *current;
  struct wetting_front *current_free_drainage;
  struct wetting_front *current_free_drainage_next;
  int soil_num;
  double f_p = 0.0;
  double runoff = 0.0;

  double h_p = fmax(*ponded_depth_cm - precip_timestep_cm * timestep_h, 0.0); // water ponded on the surface

  current = head;
  current_free_drainage      = listFindFront(wf_that_supplies_free_drainage_demand, head, NULL);
  current_free_drainage_next = listFindFront(wf_that_supplies_free_drainage_demand+1, head, NULL);

  int number_of_wetting_fronts = listLength(head);

  //int last_wetting_front_index = number_of_wetting_fronts;
  int layer_num_fp = current_free_drainage->layer_num;


  double Geff;

  if (number_of_wetting_fronts == num_layers) {
    Geff = 0.0; // i.e., case of no capillary suction, dz/dt is also zero for all wetting fronts
    soil_num = soil_type[layer_num_fp];
    Ksat_cm_per_h = soil_properties[soil_num].Ksat_cm_per_h * frozen_factor[current->layer_num]; //23 feb 2024
  }
  else {

    //double theta = current_free_drainage->theta;
    double theta_below = current_free_drainage_next->theta;

    soil_num = soil_type[layer_num_fp];

    theta_e = soil_properties[soil_num].theta_e;  // rhs of the new front, assumes theta_e as per Peter
    theta_r = soil_properties[soil_num].theta_r;
    h_min_cm = soil_properties[soil_num].h_min_cm;
    vg_a     = soil_properties[soil_num].vg_alpha_per_cm;
    vg_m     = soil_properties[soil_num].vg_m;
    vg_n     = soil_properties[soil_num].vg_n;
    double lambda = soil_properties[soil_num].bc_lambda;
    double bc_psib_cm = soil_properties[soil_num].bc_psib_cm;
    Ksat_cm_per_h = soil_properties[soil_num].Ksat_cm_per_h * frozen_factor[current->layer_num];

    // Se = calc_Se_from_theta(theta,theta_e,theta_r);
    // psi_cm = calc_h_from_Se(Se, vg_a, vg_m, vg_n);

    Geff = calc_Geff(use_closed_form_G, theta_below, theta_e, theta_e, theta_r, vg_a, vg_n, vg_m, h_min_cm, Ksat_cm_per_h, nint, lambda, bc_psib_cm); 

  }

  const double saturated_pressure_head_cm =
    lgar_saturated_pressure_head_cm(current_free_drainage, soil_type, soil_properties);

  // if the free_drainage wetting front is the top most, then the potential infiltration capacity has the following simple form
  if (layer_num_fp == 1) {
      f_p = Ksat_cm_per_h * (1 + (Geff + h_p + saturated_pressure_head_cm)/current_free_drainage->depth_cm);
  }
  else {
    // see the paper "Layered Green and Ampt Infiltration With Redistribution" by La Follette et al. (https://agupubs.onlinelibrary.wiley.com/doi/pdfdirect/10.1029/2022WR033742), equations 16 or 19
    double bottom_sum = (current_free_drainage->depth_cm - cum_layer_thickness_cm[layer_num_fp-1])/Ksat_cm_per_h;

    for (int k = 1; k < layer_num_fp; k++) {
      int soil_num_k = soil_type[layer_num_fp-k];
      double Ksat_cm_per_h_k = soil_properties[soil_num_k].Ksat_cm_per_h * frozen_factor[layer_num_fp - k];

      bottom_sum += (cum_layer_thickness_cm[layer_num_fp - k] - cum_layer_thickness_cm[layer_num_fp - (k+1)])/ Ksat_cm_per_h_k;
    }

    f_p = (current_free_drainage->depth_cm / bottom_sum) + ((Geff + h_p + saturated_pressure_head_cm)*Ksat_cm_per_h/(current_free_drainage->depth_cm)); //Geff + h_p plus saturated pressure head

  }

  //this code checks if there is enough storage available for infiltrating water. That is, f_p can only be as big as there is room for water, but also considering that some water will leave via AET and free drainage. 
  double current_mass = lgar_calc_mass_bal(cum_layer_thickness_cm, head);
  double max_storage = 0.0;
  for (int k = 1; k < num_layers+1; k++) {
    int layer_num = k;
    soil_num = soil_type[layer_num];
    max_storage += soil_properties[soil_num].theta_e * (cum_layer_thickness_cm[k]-cum_layer_thickness_cm[k-1]);
  }

  const double raw_f_p_cm_per_h = f_p;
  const double storage_limit_f_p_cm_per_h =
    (max_storage + free_drainage_subtimestep_cm + AET_demand_cm - current_mass) / timestep_h;
  if (raw_fp_cm_per_h != NULL) {
    *raw_fp_cm_per_h = raw_f_p_cm_per_h;
  }
  if (storage_limit_fp_cm_per_h != NULL) {
    *storage_limit_fp_cm_per_h = storage_limit_f_p_cm_per_h;
  }

  if (f_p > storage_limit_f_p_cm_per_h){
    f_p = storage_limit_f_p_cm_per_h;
  }
  if (capped_fp_cm_per_h != NULL) {
    *capped_fp_cm_per_h = f_p;
  }

  double ponded_depth_temp = *ponded_depth_cm;

  double free_drainage_demand = 0;

  // 'if' condition is not needed ... AJ
  if ((layer_num_fp==num_layers) && (num_layers == number_of_wetting_fronts))
    ponded_depth_temp = *ponded_depth_cm - f_p * timestep_h - free_drainage_demand*0;
  else
    ponded_depth_temp = *ponded_depth_cm - f_p * timestep_h - free_drainage_demand*0;

  ponded_depth_temp   = fmax(ponded_depth_temp, 0.0);

  double fp_cm = f_p * timestep_h + free_drainage_demand/timestep_h; // infiltration in cm

  if (ponded_depth_max_cm > 0.0 ) {

    if (ponded_depth_temp < ponded_depth_max_cm) {
      runoff = 0.0;
      *volin_this_timestep = fmin(*ponded_depth_cm, fp_cm); //PTL: does this code account for the case where volin_this_timestep can not all infiltrate?
      *ponded_depth_cm     = *ponded_depth_cm - *volin_this_timestep;

      if (verbosity.compare("high") == 0){
        printf("fp_cm: %lf \n", fp_cm);
        printf("*volin_this_timestep: %lf \n", *volin_this_timestep);
        printf("runoff: %lf \n", runoff);
        printf("*ponded_depth_cm: %lf \n", *ponded_depth_cm);
      }

      return runoff;
    }
    else if (ponded_depth_temp > ponded_depth_max_cm ) {
      runoff = ponded_depth_temp - ponded_depth_max_cm;
      *ponded_depth_cm     = ponded_depth_max_cm;
      *volin_this_timestep = fp_cm;

      if (verbosity.compare("high") == 0){
        printf("fp_cm: %lf \n", fp_cm);
        printf("*volin_this_timestep: %lf \n", *volin_this_timestep);
        printf("runoff: %lf \n", runoff);
        printf("*ponded_depth_cm: %lf \n", *ponded_depth_cm);
      }

      return runoff;
    }

  }

  else {
    // if it got to this point, no ponding is allowed, either infiltrate or runoff
    // order is important here; assign zero to ponded depth once we compute volume in and runoff
    *volin_this_timestep = fmin(*ponded_depth_cm, fp_cm); //
    runoff = *ponded_depth_cm < fp_cm ? 0.0 : (*ponded_depth_cm - *volin_this_timestep);
    *ponded_depth_cm = 0.0;

  }

  if (verbosity.compare("high") == 0){
    printf("fp_cm: %lf \n", fp_cm);
    printf("*volin_this_timestep: %lf \n", *volin_this_timestep);
    printf("runoff: %lf \n", runoff);
    printf("*ponded_depth_cm: %lf \n", *ponded_depth_cm);
  }

  return runoff;
}

static struct wetting_front* lgar_insert_surficial_front(bool TO_enabled,
                                                         double dry_depth,
                                                         double theta_new,
                                                         int layer_num,
                                                         bool to_bottom,
                                                         struct wetting_front **head)
{
  if (!TO_enabled || *head == NULL) {
    listInsertFirst(dry_depth, theta_new, 1, layer_num, to_bottom, head);
    return *head;
  }

  // In TO mode, keep the corrected list geometry intact by inserting the new
  // surface front directly after the last front that is shallower than, or at,
  // the requested dry depth.
  // if ((*head)->depth_cm > dry_depth) {
  if ((*head)->depth_cm >= dry_depth) {
    listInsertFirst(dry_depth, theta_new, 1, layer_num, to_bottom, head);
    return *head;
  }

  struct wetting_front *previous = *head;
  // while (previous->next != NULL && previous->next->depth_cm <= dry_depth) {
  while (previous->next != NULL && previous->next->depth_cm < dry_depth) {
    previous = previous->next;
  }

  struct wetting_front *inserted = (struct wetting_front*) malloc(sizeof(struct wetting_front));
  if (inserted == NULL) {
    throw runtime_error("Unable to allocate a new surficial wetting front.\n");
  }

  inserted->depth_cm = dry_depth;
  inserted->theta = theta_new;
  inserted->front_num = previous->front_num + 1;
  inserted->layer_num = layer_num;
  inserted->to_bottom = to_bottom;
  inserted->dzdt_cm_per_h = 0.0;
  inserted->is_WF_GW = false;
  inserted->next = previous->next;
  previous->next = inserted;

  for (struct wetting_front *scan = inserted->next; scan != NULL; scan = scan->next) {
    scan->front_num++;
  }

  return inserted;
}

static bool lgar_reset_overtaken_groundwater_fronts_to_surface(struct wetting_front *head,
                                                               struct wetting_front *new_surface_front)
{
  bool moved_any = false;

  for (struct wetting_front *current = head; current != NULL && current != new_surface_front; current = current->next) {
    if (current->is_WF_GW && current->depth_cm > 0.0) {
      current->depth_cm = 0.0;
      moved_any = true;
    }
  }

  return moved_any;
}

static double lgar_restore_new_surface_front_mass(struct wetting_front *new_surface_front,
                                                  double target_mass,
                                                  double *cum_layer_thickness_cm,
                                                  int *soil_type,
                                                  struct wetting_front *head,
                                                  struct soil_properties_ *soil_properties)
{
  if (new_surface_front == NULL) {
    return target_mass;
  }

  if (new_surface_front->is_WF_GW) {
    return target_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, head);
  }

  const double base_depth = cum_layer_thickness_cm[new_surface_front->layer_num - 1];
  const double block_thickness_cm = new_surface_front->depth_cm - base_depth;
  if (block_thickness_cm <= 0.0) {
    return target_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, head);
  }

  const double current_mass = lgar_calc_mass_bal(cum_layer_thickness_cm, head);
  const double mass_error_cm = target_mass - current_mass;
  if (std::fabs(mass_error_cm) <= 1.0e-12) {
    return 0.0;
  }

  const double theta_before = new_surface_front->theta;
  const double psi_before = new_surface_front->psi_cm;

  const int soil_num = soil_type[new_surface_front->layer_num];
  const double theta_e = soil_properties[soil_num].theta_e;
  const double theta_r = soil_properties[soil_num].theta_r;
  const double corrected_theta = new_surface_front->theta + mass_error_cm / block_thickness_cm;

  new_surface_front->theta = fmax(theta_r, fmin(theta_e, corrected_theta));

  if (new_surface_front->theta >= theta_e) {
    new_surface_front->psi_cm = 0.0;
  }
  else {
    const double vg_alpha = soil_properties[soil_num].vg_alpha_per_cm;
    const double vg_m = soil_properties[soil_num].vg_m;
    const double vg_n = soil_properties[soil_num].vg_n;
    const double Se = calc_Se_from_theta(new_surface_front->theta, theta_e, theta_r);
    new_surface_front->psi_cm = calc_h_from_Se(Se, vg_alpha, vg_m, vg_n);
  }

  if (verbosity.compare("high") == 0) {
    printf("surface creation theta-repair updated front %d in layer %d: "
           "depth %.6f cm, theta %.15f -> %.15f, psi %.15f -> %.15f, "
           "mass_error %.12e cm\n",
           new_surface_front->front_num, new_surface_front->layer_num, new_surface_front->depth_cm,
           theta_before, new_surface_front->theta, psi_before, new_surface_front->psi_cm, mass_error_cm);
  }

  return target_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, head);
}

static double lgar_restore_new_surface_front_mass_via_depth(struct wetting_front *new_surface_front,
                                                            double target_mass,
                                                            double *cum_layer_thickness_cm,
                                                            struct wetting_front *head)
{
  if (new_surface_front == NULL || head == NULL || new_surface_front->is_WF_GW) {
    return target_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, head);
  }

  const double base_depth = cum_layer_thickness_cm[new_surface_front->layer_num - 1];
  double lower_depth = base_depth;
  double upper_depth = cum_layer_thickness_cm[new_surface_front->layer_num];

  if (new_surface_front->next != NULL && new_surface_front->next->layer_num == new_surface_front->layer_num) {
    upper_depth = new_surface_front->next->depth_cm;
  }

  if (upper_depth <= lower_depth) {
    return target_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, head);
  }

  const double original_depth = fmax(lower_depth, fmin(new_surface_front->depth_cm, upper_depth));
  new_surface_front->depth_cm = original_depth;
  const double depth_before = original_depth;

  double current_mass = lgar_calc_mass_bal(cum_layer_thickness_cm, head);
  double remaining_mass_error = target_mass - current_mass;
  if (std::fabs(remaining_mass_error) <= MBAL_ITERATIVE_TOLERANCE) {
    return remaining_mass_error;
  }

  const bool need_more_mass = remaining_mass_error > 0.0;
  double bracket_lo = need_more_mass ? original_depth : lower_depth;
  double bracket_hi = need_more_mass ? upper_depth : original_depth;

  if (bracket_hi <= bracket_lo) {
    return remaining_mass_error;
  }

  new_surface_front->depth_cm = bracket_hi;
  const double mass_at_hi = lgar_calc_mass_bal(cum_layer_thickness_cm, head);
  new_surface_front->depth_cm = bracket_lo;
  const double mass_at_lo = lgar_calc_mass_bal(cum_layer_thickness_cm, head);

  const bool monotonic_increasing = mass_at_hi > mass_at_lo + MBAL_ITERATIVE_TOLERANCE;
  if (!monotonic_increasing) {
    new_surface_front->depth_cm = original_depth;
    return target_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, head);
  }

  if (need_more_mass && mass_at_hi < target_mass) {
    new_surface_front->depth_cm = bracket_hi;
    return target_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, head);
  }

  if (!need_more_mass && mass_at_lo > target_mass) {
    new_surface_front->depth_cm = bracket_lo;
    return target_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, head);
  }

  double best_depth = original_depth;
  double best_mass_error = std::fabs(remaining_mass_error);

  for (int iter = 0; iter < 80; iter++) {
    const double probe_depth = 0.5 * (bracket_lo + bracket_hi);
    new_surface_front->depth_cm = probe_depth;
    current_mass = lgar_calc_mass_bal(cum_layer_thickness_cm, head);
    remaining_mass_error = target_mass - current_mass;

    if (std::fabs(remaining_mass_error) < best_mass_error) {
      best_mass_error = std::fabs(remaining_mass_error);
      best_depth = probe_depth;
    }

    if (std::fabs(remaining_mass_error) <= MBAL_ITERATIVE_TOLERANCE) {
      best_depth = probe_depth;
      break;
    }

    if (remaining_mass_error > 0.0) {
      bracket_lo = probe_depth;
    }
    else {
      bracket_hi = probe_depth;
    }
  }

  new_surface_front->depth_cm = best_depth;
  if (verbosity.compare("high") == 0) {
    printf("surface creation depth-repair updated front %d in layer %d: "
           "depth %.6f cm -> %.6f cm, residual %.12e cm\n",
           new_surface_front->front_num, new_surface_front->layer_num, depth_before,
           new_surface_front->depth_cm, target_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, head));
  }
  return target_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, head);
}

static void lgar_refresh_front_state_from_theta(struct wetting_front *front,
                                                int *soil_type,
                                                double *frozen_factor,
                                                struct soil_properties_ *soil_properties)
{
  if (front == NULL) {
    return;
  }

  const int soil_num = soil_type[front->layer_num];
  const double theta_e = soil_properties[soil_num].theta_e;
  const double theta_r = soil_properties[soil_num].theta_r;
  const double vg_a = soil_properties[soil_num].vg_alpha_per_cm;
  const double vg_m = soil_properties[soil_num].vg_m;
  const double vg_n = soil_properties[soil_num].vg_n;
  const double Ksat_cm_per_h =
    frozen_factor[front->layer_num] * soil_properties[soil_num].Ksat_cm_per_h;

  const double Se = calc_Se_from_theta(front->theta, theta_e, theta_r);
  front->psi_cm = calc_h_from_Se(Se, vg_a, vg_m, vg_n);
  front->K_cm_per_h = calc_K_from_Se(Se, Ksat_cm_per_h, vg_m);
}

static void lgar_refresh_all_front_states_from_theta(struct wetting_front *head,
                                                     int *soil_type,
                                                     double *frozen_factor,
                                                     struct soil_properties_ *soil_properties)
{
  for (struct wetting_front *current = head; current != NULL; current = current->next) {
    lgar_refresh_front_state_from_theta(current, soil_type, frozen_factor, soil_properties);
  }

  for (struct wetting_front *current = head; current != NULL && current->next != NULL; current = current->next) {
    struct wetting_front *next = current->next;
    if (!current->to_bottom || current->layer_num == next->layer_num) {
      continue;
    }

    const double psi_mismatch_cm = std::fabs(current->psi_cm - next->psi_cm);
    const double psi_scale_cm = fmax(1.0, fmax(std::fabs(current->psi_cm), std::fabs(next->psi_cm)));
    const double roundtrip_psi_tol_cm = fmax(1.E-3, 1.E-5 * psi_scale_cm);
    if (psi_mismatch_cm > roundtrip_psi_tol_cm) {
      continue;
    }

    current->psi_cm = next->psi_cm;
    const int soil_num = soil_type[current->layer_num];
    const double theta_e = soil_properties[soil_num].theta_e;
    const double theta_r = soil_properties[soil_num].theta_r;
    const double vg_a = soil_properties[soil_num].vg_alpha_per_cm;
    const double vg_m = soil_properties[soil_num].vg_m;
    const double vg_n = soil_properties[soil_num].vg_n;
    const double Ksat_cm_per_h =
      frozen_factor[current->layer_num] * soil_properties[soil_num].Ksat_cm_per_h;

    current->theta = calc_theta_from_h(current->psi_cm, vg_a, vg_m, vg_n, theta_e, theta_r);
    const double Se = calc_Se_from_theta(current->theta, theta_e, theta_r);
    current->K_cm_per_h = calc_K_from_Se(Se, Ksat_cm_per_h, vg_m);
  }
}

static void lgarto_cap_zero_depth_TO_psi(int num_layers,
                                         int *soil_type,
                                         double *frozen_factor,
                                         struct soil_properties_ *soil_properties,
                                         struct wetting_front *head)
{
  for (struct wetting_front *current = head; current != NULL; current = current->next) {
    if (!current->is_WF_GW ||
        current->depth_cm != 0.0 ||
        !std::isfinite(current->psi_cm) ||
        current->psi_cm <= ZERO_DEPTH_TO_PSI_CAP_CM ||
        current->layer_num < 1 ||
        current->layer_num > num_layers) {
      continue;
    }

    const double psi_before_cm = current->psi_cm;
    const int soil_num = soil_type[current->layer_num];
    const double theta_e = soil_properties[soil_num].theta_e;
    const double theta_r = soil_properties[soil_num].theta_r;
    const double vg_a = soil_properties[soil_num].vg_alpha_per_cm;
    const double vg_m = soil_properties[soil_num].vg_m;
    const double vg_n = soil_properties[soil_num].vg_n;
    const double Ksat_cm_per_h =
      frozen_factor[current->layer_num] * soil_properties[soil_num].Ksat_cm_per_h;

    current->psi_cm = ZERO_DEPTH_TO_PSI_CAP_CM;
    current->theta = calc_theta_from_h(current->psi_cm, vg_a, vg_m, vg_n, theta_e, theta_r);
    const double Se = calc_Se_from_theta(current->theta, theta_e, theta_r);
    current->K_cm_per_h = calc_K_from_Se(Se, Ksat_cm_per_h, vg_m);

    if (verbosity.compare("high") == 0) {
      printf("Capped zero-depth TO/GW front %d psi from %.17g cm to %.17g cm and refreshed theta.\n",
             current->front_num,
             psi_before_cm,
             current->psi_cm);
    }
  }
}

static struct wetting_front *lgar_find_surface_creation_repair_target(struct wetting_front *head)
{
  struct wetting_front *current = head;

  while (current != NULL && current->is_WF_GW && current->depth_cm == 0.0) {
    current = current->next;
  }

  return current != NULL ? current : head;
}

static bool lgar_should_delay_surface_creation_gw_conversion(struct wetting_front *head)
{
  struct wetting_front *repair_target = lgar_find_surface_creation_repair_target(head);
  if (repair_target == NULL || repair_target->is_WF_GW || repair_target->to_bottom ||
      repair_target->next == NULL) {
    return false;
  }

  if (!repair_target->next->is_WF_GW) {
    return false;
  }

  if (repair_target->layer_num != repair_target->next->layer_num) {
    return false;
  }

  return repair_target->theta < repair_target->next->theta;
}

static bool lgar_connected_surface_creation_structure_is_effectively_saturated(struct wetting_front *front_start,
                                                                               int *soil_type,
                                                                               struct soil_properties_ *soil_properties);

static bool lgar_created_surface_front_can_fold_into_connected_gw_chain(struct wetting_front *repair_target,
                                                                        int *soil_type,
                                                                        struct soil_properties_ *soil_properties)
{
  if (repair_target == NULL || repair_target->next == NULL) {
    return false;
  }

  struct wetting_front *chain_head = repair_target->next;
  if (!chain_head->is_WF_GW || !chain_head->to_bottom) {
    return false;
  }

  const bool downstream_gw_structure_is_saturated =
    lgar_connected_surface_creation_structure_is_effectively_saturated(chain_head, soil_type,
                                                                       soil_properties);
  const bool hydraulically_redundant_with_chain_head =
    std::fabs(repair_target->psi_cm - chain_head->psi_cm) <= 1.0e-3;

  return downstream_gw_structure_is_saturated || hydraulically_redundant_with_chain_head;
}

static struct wetting_front *lgar_fold_redundant_created_gw_front_into_connected_chain(struct wetting_front *repair_target,
                                                                                       struct wetting_front **head,
                                                                                       int *soil_type,
                                                                                       struct soil_properties_ *soil_properties)
{
  if (repair_target == NULL || head == NULL || *head == NULL || repair_target->next == NULL) {
    return repair_target;
  }

  if (repair_target->is_WF_GW &&
      !repair_target->to_bottom &&
      repair_target->next->is_WF_GW &&
      repair_target->next->to_bottom &&
      repair_target->layer_num == repair_target->next->layer_num &&
      std::fabs(repair_target->psi_cm - repair_target->next->psi_cm) <= 1.0e-10) {
    if (verbosity.compare("high") == 0) {
      printf("surface creation folding redundant GW front %d into connected chain below "
             "because it shares psi %.15f with front %d.\n",
             repair_target->front_num, repair_target->psi_cm, repair_target->next->front_num);
    }
    return listDeleteFront(repair_target->front_num, head, soil_type, soil_properties);
  }

  return repair_target;
}

static struct wetting_front *lgar_fold_redundant_created_surface_front_into_connected_chain(struct wetting_front *repair_target,
                                                                                            struct wetting_front **head,
                                                                                            int *soil_type,
                                                                                            struct soil_properties_ *soil_properties)
{
  if (repair_target == NULL || head == NULL || *head == NULL || repair_target->next == NULL) {
    return repair_target;
  }

  const bool same_layer_and_depth =
    repair_target->layer_num == repair_target->next->layer_num &&
    std::fabs(repair_target->depth_cm - repair_target->next->depth_cm) <= CREATION_COLOCATED_TOLERANCE_CM;
  const bool created_surface_not_drier =
    !repair_target->is_WF_GW &&
    !repair_target->to_bottom &&
    repair_target->next->to_bottom &&
    repair_target->theta >= repair_target->next->theta - 1.0e-10;
  const bool can_fold_into_surface_chain =
    created_surface_not_drier && !repair_target->next->is_WF_GW && same_layer_and_depth;
  const bool can_fold_into_gw_chain =
    created_surface_not_drier && repair_target->next->is_WF_GW && same_layer_and_depth &&
    lgar_created_surface_front_can_fold_into_connected_gw_chain(repair_target, soil_type,
                                                                soil_properties);

  if (can_fold_into_surface_chain || can_fold_into_gw_chain) {
    if (verbosity.compare("high") == 0) {
      printf("surface creation folding redundant surface front %d into connected chain below "
             "because it shares depth %.15f with front %d (GW=%d, to_bottom=%d, theta=%.15f, psi=%.15f).\n",
             repair_target->front_num, repair_target->depth_cm, repair_target->next->front_num,
             repair_target->next->is_WF_GW, repair_target->next->to_bottom, repair_target->next->theta,
             repair_target->next->psi_cm);
      if (can_fold_into_gw_chain) {
        printf("  fold into GW chain allowed because the downstream GW structure is effectively saturated "
               "or hydraulically redundant with the created surface front.\n");
      }
    }
    return listDeleteFront(repair_target->front_num, head, soil_type, soil_properties);
  }

  return repair_target;
}

static double lgar_restore_surface_creation_mass_on_connected_gw_chain(struct wetting_front *repair_target,
                                                                       double target_mass,
                                                                       struct wetting_front **head,
                                                                       double *cum_layer_thickness_cm,
                                                                       int *soil_type,
                                                                       struct soil_properties_ *soil_properties)
{
  if (repair_target == NULL || head == NULL || *head == NULL) {
    return target_mass;
  }

  if (verbosity.compare("high") == 0) {
    printf("surface creation chain-repair starting from front %d in layer %d "
           "(depth %.6f cm, theta %.15f, psi %.15f)\n",
           repair_target->front_num, repair_target->layer_num, repair_target->depth_cm,
           repair_target->theta, repair_target->psi_cm);
  }

  const double min_psi_cm =
    lgar_surface_creation_downstream_TO_psi_bound_cm(repair_target);
  /* Surface-creation repair may borrow storage from the adjacent TO chain, but
     it should not wet that chain past an immediately downstream TO state. */
  lgar_theta_mass_balance_correction_with_min_psi(false, repair_target->front_num,
                                                  target_mass, head,
                                                  cum_layer_thickness_cm, soil_type,
                                                  soil_properties, min_psi_cm);

  if (verbosity.compare("high") == 0) {
    printf("surface creation chain-repair completed; residual %.12e cm\n",
           target_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, *head));
    listPrint(*head);
  }

  return target_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
}

static double lgar_restore_surface_creation_mass_on_connected_surface_chain(struct wetting_front *repair_target,
                                                                            double target_mass,
                                                                            struct wetting_front **head,
                                                                            double *cum_layer_thickness_cm,
                                                                            int *soil_type,
                                                                            struct soil_properties_ *soil_properties)
{
  if (repair_target == NULL || head == NULL || *head == NULL) {
    return target_mass;
  }

  if (verbosity.compare("high") == 0) {
    printf("surface creation surface-chain repair starting from front %d in layer %d "
           "(depth %.6f cm, theta %.15f, psi %.15f)\n",
           repair_target->front_num, repair_target->layer_num, repair_target->depth_cm,
           repair_target->theta, repair_target->psi_cm);
  }

  lgar_theta_mass_balance_correction(false, repair_target->front_num, target_mass, head,
                                     cum_layer_thickness_cm, soil_type, soil_properties);

  if (verbosity.compare("high") == 0) {
    printf("surface creation surface-chain repair completed; residual %.12e cm\n",
           target_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, *head));
    listPrint(*head);
  }

  return target_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
}

static struct wetting_front *lgar_find_first_non_to_bottom_below_connected_chain(struct wetting_front *repair_target)
{
  if (repair_target == NULL) {
    return NULL;
  }

  struct wetting_front *current = repair_target;
  while (current != NULL && current->to_bottom) {
    if (current->next == NULL) {
      return NULL;
    }
    if (current->next->is_WF_GW != repair_target->is_WF_GW) {
      return NULL;
    }
    if (!current->next->to_bottom) {
      return current->next;
    }
    current = current->next;
  }

  return NULL;
}

static struct wetting_front *lgar_find_depth_target_for_connected_gw_chain(struct wetting_front *repair_target)
{
  if (repair_target == NULL || !repair_target->is_WF_GW) {
    return NULL;
  }

  if (!repair_target->to_bottom) {
    if (repair_target->next != NULL && repair_target->next->is_WF_GW) {
      return repair_target;
    }
    return NULL;
  }

  return lgar_find_first_non_to_bottom_below_connected_chain(repair_target);
}

static struct wetting_front *lgar_surface_creation_connected_TO_chain_end(struct wetting_front *repair_target)
{
  if (repair_target == NULL || !repair_target->is_WF_GW) {
    return repair_target;
  }

  struct wetting_front *chain_end = repair_target;
  while (chain_end->next != NULL &&
         chain_end->next->is_WF_GW == repair_target->is_WF_GW &&
         chain_end->next->to_bottom) {
    chain_end = chain_end->next;
  }

  if (chain_end->to_bottom &&
      chain_end->next != NULL &&
      chain_end->next->is_WF_GW == repair_target->is_WF_GW &&
      !chain_end->next->to_bottom) {
    chain_end = chain_end->next;
  }

  return chain_end;
}

static double lgar_surface_creation_downstream_TO_psi_bound_cm(struct wetting_front *repair_target)
{
  struct wetting_front *chain_end =
    lgar_surface_creation_connected_TO_chain_end(repair_target);
  if (chain_end == NULL) {
    return 0.0;
  }

  struct wetting_front *bound = chain_end->next;
  if (bound != NULL && bound->is_WF_GW &&
      std::isfinite(bound->psi_cm) && bound->psi_cm >= 0.0) {
    return bound->psi_cm;
  }

  return 0.0;
}

static bool lgar_connected_surface_creation_structure_is_effectively_saturated(struct wetting_front *front_start,
                                                                               int *soil_type,
                                                                               struct soil_properties_ *soil_properties)
{
  if (front_start == NULL) {
    return false;
  }

  for (struct wetting_front *current = front_start; current != NULL; current = current->next) {
    const int soil_num = soil_type[current->layer_num];
    const double theta_e = soil_properties[soil_num].theta_e;
    const bool near_theta_e = current->theta >= theta_e - 1.0e-8;
    const bool near_zero_psi = current->psi_cm <= 1.0e-8;

    if (!(near_theta_e || near_zero_psi)) {
      return false;
    }

    if (!current->to_bottom) {
      break;
    }
  }

  return true;
}

static double lgar_column_effective_saturation_storage_cm(int num_layers,
                                                          double *cum_layer_thickness_cm,
                                                          int *soil_type,
                                                          struct soil_properties_ *soil_properties)
{
  double max_storage_cm = 0.0;

  for (int layer_num = 1; layer_num <= num_layers; layer_num++) {
    const int soil_num = soil_type[layer_num];
    max_storage_cm +=
      soil_properties[soil_num].theta_e *
      (cum_layer_thickness_cm[layer_num] - cum_layer_thickness_cm[layer_num - 1]);
  }

  return max_storage_cm;
}

static std::vector<lgarto_surface_TO_merge_front_state>
lgarto_save_front_hydraulic_states(struct wetting_front *head)
{
  std::vector<lgarto_surface_TO_merge_front_state> states;

  for (struct wetting_front *current = head; current != NULL; current = current->next) {
    lgarto_surface_TO_merge_front_state saved = {
      current,
      current->psi_cm,
      current->theta,
      current->K_cm_per_h
    };
    states.push_back(saved);
  }

  return states;
}

static double lgar_restore_surface_creation_mass_in_available_TO_storage(
  double target_mass,
  int num_layers,
  double *cum_layer_thickness_cm,
  int *soil_type,
  struct wetting_front **head,
  struct soil_properties_ *soil_properties)
{
  if (head == NULL || *head == NULL) {
    return target_mass;
  }

  double residual_cm = target_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
  if (residual_cm <= MBAL_ITERATIVE_TOLERANCE) {
    return residual_cm;
  }
  const double initial_abs_residual_cm = fabs(residual_cm);

  const double max_storage_cm =
    lgar_column_effective_saturation_storage_cm(num_layers, cum_layer_thickness_cm,
                                                soil_type, soil_properties);
  const double storage_tol_cm = fmax(1.0e-8, MBAL_ITERATIVE_TOLERANCE);
  if (target_mass > max_storage_cm + storage_tol_cm) {
    return residual_cm;
  }

  const std::vector<lgarto_surface_TO_merge_front_state> original_states =
    lgarto_save_front_hydraulic_states(*head);
  std::vector<lgarto_surface_TO_merge_front_state> best_states = original_states;
  double best_abs_residual_cm = fabs(residual_cm);

  for (struct wetting_front *candidate = *head; candidate != NULL; candidate = candidate->next) {
    if (candidate->is_WF_GW == FALSE || candidate->to_bottom ||
        candidate->depth_cm <= ZERO_DEPTH_TO_DELETE_DEPTH_TOL_CM ||
        !std::isfinite(candidate->psi_cm) || candidate->psi_cm <= 0.0) {
      continue;
    }

    lgarto_restore_front_hydraulic_states(original_states);
    const double min_psi_cm =
      lgar_surface_creation_downstream_TO_psi_bound_cm(candidate);
    /* Keep trial storage repairs from making an upper TO chain wetter than an
       immediately downstream TO profile it drains into. */
    lgar_theta_mass_balance_correction_with_min_psi(false, candidate->front_num,
                                                    target_mass, head,
                                                    cum_layer_thickness_cm, soil_type,
                                                    soil_properties, min_psi_cm);

    const double candidate_residual_cm =
      target_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
    const double candidate_abs_residual_cm = fabs(candidate_residual_cm);
    if (candidate_abs_residual_cm < best_abs_residual_cm) {
      best_abs_residual_cm = candidate_abs_residual_cm;
      best_states = lgarto_save_front_hydraulic_states(*head);
      residual_cm = candidate_residual_cm;
    }

    if (candidate_abs_residual_cm <= MBAL_ITERATIVE_TOLERANCE) {
      break;
    }
  }

  lgarto_restore_front_hydraulic_states(best_states);
  residual_cm = target_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, *head);

  if (verbosity.compare("high") == 0 &&
      fabs(residual_cm) + 1.0e-12 < initial_abs_residual_cm) {
    printf("surface creation attempted whole-column TO storage repair; residual %.12e cm\n",
           residual_cm);
  }

  return residual_cm;
}

static double lgar_restore_surface_creation_mass_on_connected_gw_chain_via_depth(struct wetting_front *repair_target,
                                                                                 double target_mass,
                                                                                 int num_layers,
                                                                                 struct wetting_front *head,
                                                                                 double *cum_layer_thickness_cm)
{
  if (repair_target == NULL || head == NULL) {
    return target_mass;
  }

  struct wetting_front *depth_target = lgar_find_depth_target_for_connected_gw_chain(repair_target);
  if (depth_target == NULL) {
    return target_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, head);
  }

  const double original_depth = depth_target->depth_cm;
  double upper_depth = cum_layer_thickness_cm[num_layers];
  if (depth_target->next != NULL) {
    upper_depth = depth_target->next->depth_cm;
  }

  if (upper_depth <= original_depth) {
    return target_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, head);
  }

  if (verbosity.compare("high") == 0) {
    printf("surface creation connected-chain depth-repair starting from front %d in layer %d "
           "(depth %.6f cm, upper bound %.6f cm)\n",
           depth_target->front_num, depth_target->layer_num, depth_target->depth_cm, upper_depth);
  }

  const double mass_at_lo = lgar_calc_mass_bal(cum_layer_thickness_cm, head);
  depth_target->depth_cm = upper_depth;
  const double mass_at_hi = lgar_calc_mass_bal(cum_layer_thickness_cm, head);
  depth_target->depth_cm = original_depth;

  if (mass_at_hi <= mass_at_lo + MBAL_ITERATIVE_TOLERANCE) {
    return target_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, head);
  }

  if (mass_at_hi < target_mass) {
    depth_target->depth_cm = upper_depth;
    return target_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, head);
  }

  double bracket_lo = original_depth;
  double bracket_hi = upper_depth;
  double best_depth = original_depth;
  double best_mass_error = std::fabs(target_mass - mass_at_lo);

  for (int iter = 0; iter < 80; iter++) {
    const double probe_depth = 0.5 * (bracket_lo + bracket_hi);
    depth_target->depth_cm = probe_depth;
    const double current_mass = lgar_calc_mass_bal(cum_layer_thickness_cm, head);
    const double remaining_mass_error = target_mass - current_mass;

    if (std::fabs(remaining_mass_error) < best_mass_error) {
      best_mass_error = std::fabs(remaining_mass_error);
      best_depth = probe_depth;
    }

    if (std::fabs(remaining_mass_error) <= MBAL_ITERATIVE_TOLERANCE) {
      best_depth = probe_depth;
      break;
    }

    if (remaining_mass_error > 0.0) {
      bracket_lo = probe_depth;
    }
    else {
      bracket_hi = probe_depth;
    }
  }

  depth_target->depth_cm = best_depth;
  if (verbosity.compare("high") == 0) {
    printf("surface creation connected-chain depth-repair updated front %d in layer %d: "
           "depth %.6f cm -> %.6f cm, residual %.12e cm\n",
           depth_target->front_num, depth_target->layer_num, original_depth, depth_target->depth_cm,
           target_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, head));
  }

  return target_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, head);
}

static double lgar_delete_unbracketed_surface_creation_front_and_repair_TO_chain(
  struct wetting_front *repair_target,
  double target_mass,
  int num_layers,
  struct wetting_front **head,
  double *cum_layer_thickness_cm,
  int *soil_type,
  struct soil_properties_ *soil_properties,
  bool *adjusted_depth)
{
  if (repair_target == NULL || head == NULL || *head == NULL ||
      repair_target->next == NULL ||
      repair_target->is_WF_GW ||
      repair_target->to_bottom ||
      repair_target->layer_num != repair_target->next->layer_num ||
      !repair_target->next->is_WF_GW) {
    return target_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
  }

  const int deleted_front_num = repair_target->front_num;
  const int chain_front_num = repair_target->next->front_num;
  if (verbosity.compare("high") == 0) {
    printf("surface creation deleting unbracketed dry surface front %d and handing repair "
           "to adjacent TO/GW chain starting at front %d (depth %.6f cm, psi %.15f).\n",
           deleted_front_num, chain_front_num, repair_target->next->depth_cm,
           repair_target->next->psi_cm);
  }

  /*
   * If the created/local surface front has already reached residual water
   * content and the profile is still too wet, that surface front cannot
   * represent the needed local storage change. Delete the transient surface
   * discontinuity and solve on the immediately deeper connected TO/GW chain,
   * preserving boundary psi continuity through the existing chain correction.
   */
  struct wetting_front *chain_target =
    listDeleteFront(deleted_front_num, head, soil_type, soil_properties);
  if (chain_target == NULL) {
    return target_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
  }

  double remaining_mass_error =
    lgar_restore_surface_creation_mass_on_connected_gw_chain(chain_target, target_mass, head,
                                                             cum_layer_thickness_cm, soil_type,
                                                             soil_properties);
  if (remaining_mass_error > MBAL_ITERATIVE_TOLERANCE) {
    const double residual_before_depth = remaining_mass_error;
    remaining_mass_error =
      lgar_restore_surface_creation_mass_on_connected_gw_chain_via_depth(chain_target, target_mass,
                                                                         num_layers, *head,
                                                                         cum_layer_thickness_cm);
    if (adjusted_depth != NULL) {
      *adjusted_depth = *adjusted_depth ||
                        std::fabs(remaining_mass_error - residual_before_depth) > 1.0e-12;
    }
  }

  return remaining_mass_error;
}

static double lgar_route_surface_creation_residual_if_needed(double remaining_mass_error_after_depth,
                                                             int num_layers,
                                                             double *target_mass,
                                                             double *creation_excess_gw_flux_cm,
                                                             double *creation_excess_runoff_cm,
                                                             double *saturated_creation_gw_flux_capacity_cm,
                                                             struct wetting_front *routing_target,
                                                             int *soil_type,
                                                             struct soil_properties_ *soil_properties,
                                                             double *cum_layer_thickness_cm,
                                                             struct wetting_front *head,
                                                             const char *message)
{
  if (remaining_mass_error_after_depth > MBAL_ITERATIVE_TOLERANCE &&
      (creation_excess_gw_flux_cm != NULL || creation_excess_runoff_cm != NULL)) {
    const double residual_before_column_storage = remaining_mass_error_after_depth;
    remaining_mass_error_after_depth =
      lgar_restore_surface_creation_mass_in_available_TO_storage(
        *target_mass, num_layers, cum_layer_thickness_cm, soil_type, &head, soil_properties);
    if (remaining_mass_error_after_depth <= MBAL_ITERATIVE_TOLERANCE) {
      return remaining_mass_error_after_depth;
    }
    if (verbosity.compare("high") == 0 &&
        remaining_mass_error_after_depth + 1.0e-12 < residual_before_column_storage) {
      printf("surface creation whole-column storage repair reduced residual %.12e -> %.12e cm before routing.\n",
             residual_before_column_storage, remaining_mass_error_after_depth);
    }

    const bool route_to_runoff =
      lgar_connected_surface_creation_structure_is_effectively_saturated(routing_target, soil_type,
                                                                         soil_properties);

    /* If creation leaves excess water against an effectively saturated TO scaffold,
       let the lower boundary take only the hydraulic capacity that remains in this
       substep. This preserves the useful drainage behavior seen before the
       double-motion guard without allowing unlimited residual-to-recharge routing. */
    if (route_to_runoff && creation_excess_gw_flux_cm != NULL &&
        saturated_creation_gw_flux_capacity_cm != NULL &&
        *saturated_creation_gw_flux_capacity_cm > MBAL_ITERATIVE_TOLERANCE) {
      const double bounded_gw_flux_cm =
        fmin(remaining_mass_error_after_depth, *saturated_creation_gw_flux_capacity_cm);
      *creation_excess_gw_flux_cm += bounded_gw_flux_cm;
      *saturated_creation_gw_flux_capacity_cm -= bounded_gw_flux_cm;
      *target_mass -= bounded_gw_flux_cm;
      remaining_mass_error_after_depth -= bounded_gw_flux_cm;
      if (verbosity.compare("high") == 0) {
        printf("%s %.12e cm to lower-boundary flux using saturated TO capacity.\n",
               message, bounded_gw_flux_cm);
      }
    }

    if (remaining_mass_error_after_depth <= MBAL_ITERATIVE_TOLERANCE) {
      return *target_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, head);
    }

    if (route_to_runoff && creation_excess_runoff_cm != NULL) {
      *creation_excess_runoff_cm += remaining_mass_error_after_depth;
    }
    else if (creation_excess_gw_flux_cm != NULL) {
      lgarto_assert_gw_flux_mass_balance_correction_within_debug_threshold(
        remaining_mass_error_after_depth,
        "lgar_route_surface_creation_residual_if_needed",
        message,
        head);
      *creation_excess_gw_flux_cm += remaining_mass_error_after_depth;
    }
    else if (creation_excess_runoff_cm != NULL) {
      *creation_excess_runoff_cm += remaining_mass_error_after_depth;
    }

    *target_mass -= remaining_mass_error_after_depth;
    if (verbosity.compare("high") == 0) {
      printf("%s %.12e cm to %s.\n", message, remaining_mass_error_after_depth,
             route_to_runoff ? "surface runoff" : "lower-boundary flux");
    }
    remaining_mass_error_after_depth = *target_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, head);
  }

  return remaining_mass_error_after_depth;
}

extern double lgarto_apply_surface_TO_merge_creation_correction(
  const char *caller_name,
  const char *residual_description,
  int num_layers,
  double *target_mass,
  double *cum_layer_thickness_cm,
  int *soil_type,
  double *frozen_factor,
  struct wetting_front **head,
  struct soil_properties_ *soil_properties)
{
  const double mass_before_surface_TO_merge = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
  bool merged_in_non_top_layer = false;
  merged_in_non_top_layer =
    lgar_merge_surface_and_TO_wetting_fronts(merged_in_non_top_layer, num_layers,
                                             cum_layer_thickness_cm, soil_type,
                                             frozen_factor, soil_properties, head);
  const bool did_a_WF_have_negative_depth = lgarto_correct_negative_depths(head);
  lgarto_cleanup_after_surface_TO_merging_in_layer_below_top(merged_in_non_top_layer, soil_type,
                                                             soil_properties, head);
  const bool close_psis = correct_close_psis(soil_type, soil_properties, head);

  double mass_balance_flux_correction_cm =
    mass_before_surface_TO_merge - lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
  int boundary_pinned_surface_front_num = -1;
  if (fabs(mass_balance_flux_correction_cm) > MBAL_ITERATIVE_TOLERANCE) {
    mass_balance_flux_correction_cm =
      lgarto_restore_surface_TO_merge_mass_via_surface_depths(mass_before_surface_TO_merge,
                                                              num_layers,
                                                              cum_layer_thickness_cm,
                                                              *head,
                                                              &boundary_pinned_surface_front_num);
  }
  if (mass_balance_flux_correction_cm > MBAL_ITERATIVE_TOLERANCE &&
      boundary_pinned_surface_front_num > -1) {
    mass_balance_flux_correction_cm =
      lgarto_restore_surface_TO_merge_mass_via_boundary_psi(mass_before_surface_TO_merge,
                                                            boundary_pinned_surface_front_num,
                                                            num_layers,
                                                            cum_layer_thickness_cm,
                                                            soil_type,
                                                            frozen_factor,
                                                            soil_properties,
                                                            *head);
  }
  if (did_a_WF_have_negative_depth || close_psis ||
      fabs(mass_balance_flux_correction_cm) > MBAL_ITERATIVE_TOLERANCE) {
    lgarto_assert_gw_flux_mass_balance_correction_within_debug_threshold(
      mass_balance_flux_correction_cm,
      caller_name,
      residual_description,
      *head);
    if (target_mass != NULL) {
      *target_mass -= mass_balance_flux_correction_cm;
    }
  }

  return mass_balance_flux_correction_cm;
}

static double lgar_apply_surface_creation_correction(int correction_type_surf,
                                                     int num_layers,
                                                     double *target_mass,
                                                     double surface_lower_boundary_depth_cm,
                                                     double *cum_layer_thickness_cm,
                                                     int *soil_type,
                                                     double *frozen_factor,
                                                     struct wetting_front **head,
                                                     struct soil_properties_ *soil_properties)
{
  // Keep target-mass edits paired with the same signed lower-boundary flux bookkeeping.
  double lower_boundary_flux_correction_cm = 0.0;

  if (correction_type_surf == 1) {
    lgar_merge_wetting_fronts(soil_type, frozen_factor, head, soil_properties);
  }

  if (correction_type_surf == 2) {
    correct_close_psis(soil_type, soil_properties, head);
    lgar_wetting_fronts_cross_layer_boundary(num_layers, cum_layer_thickness_cm, soil_type,
                                             frozen_factor, head, soil_properties);
  }

  if (correction_type_surf == 3) {
    double flux_out_cm = lgar_wetting_front_cross_domain_boundary(surface_lower_boundary_depth_cm,
                                                                  soil_type, frozen_factor, head,
                                                                  soil_properties);
    if (!std::isfinite(flux_out_cm)) {
      flux_out_cm = 0.0;
    }
    *target_mass -= flux_out_cm;
    lower_boundary_flux_correction_cm += flux_out_cm;
  }

  if (correction_type_surf == 4) {
    double mass_change = 0.0;
    lgar_fix_dry_over_wet_wetting_fronts(&mass_change, cum_layer_thickness_cm, soil_type, head,
                                         soil_properties);
  }

  if (correction_type_surf == 5) {
    lgarto_convert_surface_fronts_drier_than_TO_below(head, "surface creation correction type 5");
  }

  if (correction_type_surf == 6) {
    (void) num_layers;
    (void) target_mass;
    (void) cum_layer_thickness_cm;
    (void) soil_type;
    (void) frozen_factor;
    (void) soil_properties;
    (void) lgarto_convert_overtaken_surface_front_above_TO_chain(
      head, "surface creation mixed surface/surface/TO correction");
  }

  if (correction_type_surf == 7) {
    lower_boundary_flux_correction_cm +=
      lgarto_apply_surface_TO_merge_creation_correction(
      "lgar_apply_surface_creation_correction",
      "pre-TO-motion surface/TO merge cleanup mass residual",
      num_layers, target_mass, cum_layer_thickness_cm, soil_type,
      frozen_factor, head, soil_properties);
  }

  if (correction_type_surf == 9) {
    const double mass_balance_flux_correction_cm =
      lgarto_surface_fronts_cross_layer_boundary_upward(
        num_layers, cum_layer_thickness_cm, soil_type, frozen_factor,
        soil_properties, head, "surface creation correction");
    *target_mass -= mass_balance_flux_correction_cm;
    lower_boundary_flux_correction_cm += mass_balance_flux_correction_cm;
  }

  return lower_boundary_flux_correction_cm;
}

static double lgar_apply_creation_general_correction(int correction_type,
                                                     int num_layers,
                                                     double *target_mass,
                                                     double column_depth,
                                                     double surface_lower_boundary_depth_cm,
                                                     double *cum_layer_thickness_cm,
                                                     int *soil_type,
                                                     double *frozen_factor,
                                                     struct wetting_front **head,
                                                     struct soil_properties_ *soil_properties)
{
  double lower_boundary_flux_correction_cm = 0.0;

  if (correction_type == 1) {
    lower_boundary_flux_correction_cm +=
      lgarto_apply_surface_TO_merge_creation_correction(
        "lgar_apply_creation_general_correction",
        "surface-creation general surface/TO merge cleanup mass residual",
        num_layers, target_mass, cum_layer_thickness_cm, soil_type,
        frozen_factor, head, soil_properties);
  }

  if (correction_type == 2) {
    const double mass_balance_flux_correction_cm =
      lgarto_TO_WFs_merge_via_depth(*target_mass, column_depth, cum_layer_thickness_cm, head,
                                    soil_type, soil_properties);
    lgarto_assert_gw_flux_mass_balance_correction_within_debug_threshold(
      mass_balance_flux_correction_cm,
      "lgar_apply_creation_general_correction",
      "surface-creation TO depth merge residual",
      *head);
    *target_mass -= mass_balance_flux_correction_cm;
    lower_boundary_flux_correction_cm += mass_balance_flux_correction_cm;
    lgar_global_psi_update(soil_type, soil_properties, head);
  }

  if (correction_type == 3) {
    const double mass_balance_flux_correction_cm =
      lgarto_TO_WFs_merge_via_theta(*target_mass, column_depth, cum_layer_thickness_cm, head,
                                    soil_type, soil_properties);
    lgarto_assert_gw_flux_mass_balance_correction_within_debug_threshold(
      mass_balance_flux_correction_cm,
      "lgar_apply_creation_general_correction",
      "surface-creation TO theta merge residual",
      *head);
    *target_mass -= mass_balance_flux_correction_cm;
    lower_boundary_flux_correction_cm += mass_balance_flux_correction_cm;
    lgar_global_psi_update(soil_type, soil_properties, head);
  }

  if (correction_type == 4) {
    int front_num_with_negative_depth = -1;
    double mass_balance_flux_correction_cm = 0.0;
    lgar_TO_wetting_fronts_cross_layer_boundary(&front_num_with_negative_depth, num_layers,
                                                cum_layer_thickness_cm, soil_type, frozen_factor, head,
                                                soil_properties,
                                                &mass_balance_flux_correction_cm);
    lgar_global_psi_update(soil_type, soil_properties, head);
    lower_boundary_flux_correction_cm += mass_balance_flux_correction_cm;

    if (front_num_with_negative_depth > -1) {
      // listDeleteFront(front_num_with_negative_depth, head, soil_type, soil_properties);
    }
    // listSendToTop(*head);
  }

  if (correction_type == 5) {
    lgar_merge_wetting_fronts(soil_type, frozen_factor, head, soil_properties);
  }

  if (correction_type == 6) {
    correct_close_psis(soil_type, soil_properties, head);
    lgar_wetting_fronts_cross_layer_boundary(num_layers, cum_layer_thickness_cm, soil_type,
                                             frozen_factor, head, soil_properties);
  }

  if (correction_type == 7) {
    double flux_out_cm = lgar_wetting_front_cross_domain_boundary(surface_lower_boundary_depth_cm,
                                                                  soil_type, frozen_factor, head,
                                                                  soil_properties);
    if (!std::isfinite(flux_out_cm)) {
      flux_out_cm = 0.0;
    }
    *target_mass -= flux_out_cm;
    lower_boundary_flux_correction_cm += flux_out_cm;
  }

  if (correction_type == 8) {
    lgarto_convert_surface_fronts_drier_than_TO_below(head, "surface creation correction type 8");
  }

  if (correction_type == 9) {
    const double flux_out_cm =
      lgarto_truncate_last_layer_GW_overshoot(cum_layer_thickness_cm, num_layers,
                                              head, soil_type, soil_properties);
    lgarto_assert_gw_flux_mass_balance_correction_within_debug_threshold(
      flux_out_cm,
      "lgar_apply_creation_general_correction",
      "surface-creation last-layer TO/GW overshoot truncation",
      *head);
    *target_mass -= flux_out_cm;
    lower_boundary_flux_correction_cm = flux_out_cm;
  }

  if (correction_type == 10) {
    const double mass_balance_flux_correction_cm =
      lgarto_surface_fronts_cross_layer_boundary_upward(
        num_layers, cum_layer_thickness_cm, soil_type, frozen_factor,
        soil_properties, head, "surface-creation general correction");
    *target_mass -= mass_balance_flux_correction_cm;
    lower_boundary_flux_correction_cm += mass_balance_flux_correction_cm;
  }

  return lower_boundary_flux_correction_cm;
}

static void lgar_normalize_after_surface_front_creation(int num_layers,
                                                        double target_mass,
                                                        double surface_lower_boundary_depth_cm,
                                                        double *cum_layer_thickness_cm,
                                                        int *soil_type,
                                                        double *frozen_factor,
                                                        struct wetting_front *new_surface_front,
                                                        struct wetting_front **head,
                                                        struct soil_properties_ *soil_properties,
                                                        double *creation_excess_gw_flux_cm,
                                                        double *creation_excess_runoff_cm,
                                                        double saturated_creation_gw_flux_capacity_cm)
{
  if (*head == NULL) {
    return;
  }

  const double column_depth = cum_layer_thickness_cm[num_layers];
  (void) lgar_reset_overtaken_groundwater_fronts_to_surface(*head, new_surface_front);
  bool delayed_surface_creation_gw_conversion = false;

  for (int phase = 0; phase < 8; phase++) {
    bool delayed_gw_conversion_this_phase = false;
    int iteration = 0;
    while (iteration < 200) {
      iteration++;

      const int correction_type_surf =
        lgarto_correction_type_surf(num_layers, cum_layer_thickness_cm, head,
                                    surface_lower_boundary_depth_cm);
      if (correction_type_surf != 0) {
        if (correction_type_surf == 5 && !delayed_surface_creation_gw_conversion &&
            lgar_should_delay_surface_creation_gw_conversion(*head)) {
          delayed_surface_creation_gw_conversion = true;
          delayed_gw_conversion_this_phase = true;
          if (verbosity.compare("high") == 0) {
            struct wetting_front *repair_target = lgar_find_surface_creation_repair_target(*head);
            if (repair_target != NULL && repair_target->next != NULL) {
              printf("surface creation delaying surface correction type 5 for front %d in layer %d "
                     "until after surface repair (theta %.15f, psi %.15f; GW below front %d "
                     "theta %.15f, psi %.15f).\n",
                     repair_target->front_num, repair_target->layer_num, repair_target->theta,
                     repair_target->psi_cm, repair_target->next->front_num, repair_target->next->theta,
                     repair_target->next->psi_cm);
            }
          }
          break;
        }

        const double lower_boundary_flux_correction_cm =
          lgar_apply_surface_creation_correction(correction_type_surf, num_layers, &target_mass,
                                                 surface_lower_boundary_depth_cm, cum_layer_thickness_cm,
                                                 soil_type, frozen_factor, head,
                                                 soil_properties);
        if (creation_excess_gw_flux_cm != NULL) {
          *creation_excess_gw_flux_cm += lower_boundary_flux_correction_cm;
        }
        continue;
      }

      const int correction_type =
        lgarto_correction_type(num_layers, cum_layer_thickness_cm, head,
                               surface_lower_boundary_depth_cm);
      if (correction_type != 0) {
        if (correction_type == 8 && !delayed_surface_creation_gw_conversion &&
            lgar_should_delay_surface_creation_gw_conversion(*head)) {
          delayed_surface_creation_gw_conversion = true;
          delayed_gw_conversion_this_phase = true;
          if (verbosity.compare("high") == 0) {
            struct wetting_front *repair_target = lgar_find_surface_creation_repair_target(*head);
            if (repair_target != NULL && repair_target->next != NULL) {
              printf("surface creation delaying correction type 8 for front %d in layer %d "
                     "until after surface repair (theta %.15f, psi %.15f; GW below front %d "
                     "theta %.15f, psi %.15f).\n",
                     repair_target->front_num, repair_target->layer_num, repair_target->theta,
                     repair_target->psi_cm, repair_target->next->front_num, repair_target->next->theta,
                     repair_target->next->psi_cm);
            }
          }
          break;
        }

        const double lower_boundary_flux_correction_cm =
          lgar_apply_creation_general_correction(correction_type, num_layers, &target_mass, column_depth,
                                                 surface_lower_boundary_depth_cm, cum_layer_thickness_cm,
                                                 soil_type, frozen_factor, head,
                                                 soil_properties);
        if (creation_excess_gw_flux_cm != NULL) {
          *creation_excess_gw_flux_cm += lower_boundary_flux_correction_cm;
        }
        continue;
      }

      break;
    }

    struct wetting_front *repair_target = lgar_find_surface_creation_repair_target(*head);
    if (repair_target == NULL) {
      break;
    }

    repair_target = lgar_fold_redundant_created_gw_front_into_connected_chain(repair_target, head,
                                                                              soil_type, soil_properties);
    repair_target = lgar_fold_redundant_created_surface_front_into_connected_chain(repair_target, head,
                                                                                   soil_type, soil_properties);
    if (repair_target == NULL) {
      break;
    }

    if (verbosity.compare("high") == 0) {
      printf("surface creation repair target after phase %d is front %d: "
             "depth %.6f cm, theta %.15f, psi %.15f, layer %d, GW=%d, to_bottom=%d\n",
             phase + 1, repair_target->front_num, repair_target->depth_cm, repair_target->theta,
             repair_target->psi_cm, repair_target->layer_num, repair_target->is_WF_GW,
             repair_target->to_bottom);
    }

    double remaining_mass_error = 0.0;
    double remaining_mass_error_after_depth = 0.0;
    bool adjusted_depth = false;

    if (repair_target->is_WF_GW && !repair_target->to_bottom) {
      remaining_mass_error_after_depth =
        lgar_restore_surface_creation_mass_on_connected_gw_chain(repair_target, target_mass, head,
                                                                 cum_layer_thickness_cm, soil_type,
                                                                 soil_properties);
      if (remaining_mass_error_after_depth > MBAL_ITERATIVE_TOLERANCE) {
        const double residual_before_depth = remaining_mass_error_after_depth;
        remaining_mass_error_after_depth =
          lgar_restore_surface_creation_mass_on_connected_gw_chain_via_depth(repair_target, target_mass,
                                                                             num_layers, *head,
                                                                             cum_layer_thickness_cm);
        adjusted_depth =
          std::fabs(remaining_mass_error_after_depth - residual_before_depth) > 1.0e-12;
      }
      if (remaining_mass_error_after_depth > MBAL_ITERATIVE_TOLERANCE &&
          (creation_excess_gw_flux_cm != NULL || creation_excess_runoff_cm != NULL)) {
        remaining_mass_error_after_depth = lgar_route_surface_creation_residual_if_needed(
          remaining_mass_error_after_depth, num_layers, &target_mass, creation_excess_gw_flux_cm,
          creation_excess_runoff_cm, &saturated_creation_gw_flux_capacity_cm, repair_target, soil_type, soil_properties,
          cum_layer_thickness_cm, *head,
          "surface creation connected TO-chain repair could not fully close mass; routing");
      }
    }
    else if (repair_target->is_WF_GW && repair_target->to_bottom) {
      remaining_mass_error_after_depth =
        lgar_restore_surface_creation_mass_on_connected_gw_chain(repair_target, target_mass, head,
                                                                 cum_layer_thickness_cm, soil_type,
                                                                 soil_properties);
      if (remaining_mass_error_after_depth > MBAL_ITERATIVE_TOLERANCE) {
        const double residual_before_depth = remaining_mass_error_after_depth;
        remaining_mass_error_after_depth =
          lgar_restore_surface_creation_mass_on_connected_gw_chain_via_depth(repair_target, target_mass,
                                                                             num_layers, *head,
                                                                             cum_layer_thickness_cm);
        adjusted_depth =
          std::fabs(remaining_mass_error_after_depth - residual_before_depth) > 1.0e-12;
      }
      if (remaining_mass_error_after_depth > MBAL_ITERATIVE_TOLERANCE &&
          (creation_excess_gw_flux_cm != NULL || creation_excess_runoff_cm != NULL)) {
        remaining_mass_error_after_depth = lgar_route_surface_creation_residual_if_needed(
          remaining_mass_error_after_depth, num_layers, &target_mass, creation_excess_gw_flux_cm,
          creation_excess_runoff_cm, &saturated_creation_gw_flux_capacity_cm, repair_target, soil_type, soil_properties,
          cum_layer_thickness_cm, *head,
          "surface creation connected-chain repair could not fully close mass; routing");
      }
    }
    else if (!repair_target->is_WF_GW && repair_target->to_bottom) {
      remaining_mass_error_after_depth =
        lgar_restore_surface_creation_mass_on_connected_surface_chain(repair_target, target_mass, head,
                                                                      cum_layer_thickness_cm, soil_type,
                                                                      soil_properties);
      if (remaining_mass_error_after_depth > MBAL_ITERATIVE_TOLERANCE) {
        const double residual_before_depth = remaining_mass_error_after_depth;
        remaining_mass_error_after_depth =
          lgar_restore_surface_creation_mass_on_connected_gw_chain_via_depth(repair_target, target_mass,
                                                                             num_layers, *head,
                                                                             cum_layer_thickness_cm);
        adjusted_depth =
          std::fabs(remaining_mass_error_after_depth - residual_before_depth) > 1.0e-12;
      }
      if (remaining_mass_error_after_depth > MBAL_ITERATIVE_TOLERANCE &&
          (creation_excess_gw_flux_cm != NULL || creation_excess_runoff_cm != NULL)) {
        remaining_mass_error_after_depth = lgar_route_surface_creation_residual_if_needed(
          remaining_mass_error_after_depth, num_layers, &target_mass, creation_excess_gw_flux_cm,
          creation_excess_runoff_cm, &saturated_creation_gw_flux_capacity_cm, repair_target, soil_type, soil_properties,
          cum_layer_thickness_cm, *head,
          "surface creation connected surface-chain repair could not fully close mass; routing");
      }
    }
    else {
      remaining_mass_error =
        lgar_restore_new_surface_front_mass(repair_target, target_mass, cum_layer_thickness_cm, soil_type,
                                            *head, soil_properties);

      remaining_mass_error_after_depth = remaining_mass_error;
      if (std::fabs(remaining_mass_error) > MBAL_ITERATIVE_TOLERANCE) {
        remaining_mass_error_after_depth =
          lgar_restore_new_surface_front_mass_via_depth(repair_target, target_mass, cum_layer_thickness_cm,
                                                        *head);
        adjusted_depth = std::fabs(remaining_mass_error_after_depth - remaining_mass_error) > 1.0e-12;
      }

      const int repair_soil_num = soil_type[repair_target->layer_num];
      const double repair_theta_e = soil_properties[repair_soil_num].theta_e;
      const double repair_theta_r = soil_properties[repair_soil_num].theta_r;
      if (remaining_mass_error_after_depth > MBAL_ITERATIVE_TOLERANCE &&
          repair_target->theta >= repair_theta_e - 1.0e-10 &&
          repair_target->next != NULL &&
          repair_target->layer_num == repair_target->next->layer_num &&
          (repair_target->next->is_WF_GW || repair_target->next->to_bottom)) {
        struct wetting_front *chain_target = repair_target->next;
        if (verbosity.compare("high") == 0) {
          printf("surface creation handing residual %.12e cm from saturated front %d to connected chain "
                 "starting at front %d (GW=%d, depth %.6f cm).\n",
                 remaining_mass_error_after_depth, repair_target->front_num, chain_target->front_num,
                 chain_target->is_WF_GW, chain_target->depth_cm);
        }

        if (chain_target->is_WF_GW) {
          remaining_mass_error_after_depth =
            lgar_restore_surface_creation_mass_on_connected_gw_chain(chain_target, target_mass, head,
                                                                     cum_layer_thickness_cm, soil_type,
                                                                     soil_properties);
        }
        else {
          remaining_mass_error_after_depth =
            lgar_restore_surface_creation_mass_on_connected_surface_chain(chain_target, target_mass, head,
                                                                          cum_layer_thickness_cm, soil_type,
                                                                          soil_properties);
        }

        if (remaining_mass_error_after_depth > MBAL_ITERATIVE_TOLERANCE) {
          const double residual_before_depth = remaining_mass_error_after_depth;
          remaining_mass_error_after_depth =
            lgar_restore_surface_creation_mass_on_connected_gw_chain_via_depth(chain_target, target_mass,
                                                                               num_layers, *head,
                                                                               cum_layer_thickness_cm);
          adjusted_depth = adjusted_depth ||
                           std::fabs(remaining_mass_error_after_depth - residual_before_depth) > 1.0e-12;
        }

        if (remaining_mass_error_after_depth > MBAL_ITERATIVE_TOLERANCE &&
            (creation_excess_gw_flux_cm != NULL || creation_excess_runoff_cm != NULL)) {
          remaining_mass_error_after_depth = lgar_route_surface_creation_residual_if_needed(
            remaining_mass_error_after_depth, num_layers, &target_mass, creation_excess_gw_flux_cm,
            creation_excess_runoff_cm, &saturated_creation_gw_flux_capacity_cm, chain_target, soil_type, soil_properties,
            cum_layer_thickness_cm, *head,
            "surface creation connected-chain handoff could not fully close mass; routing");
        }
      }

      bool deleted_unbracketed_surface_front = false;
      if (remaining_mass_error_after_depth < -MBAL_ITERATIVE_TOLERANCE &&
          repair_target->theta <= repair_theta_r + 1.0e-10 &&
          repair_target->next != NULL &&
          repair_target->layer_num == repair_target->next->layer_num &&
          repair_target->next->is_WF_GW) {
        remaining_mass_error_after_depth =
          lgar_delete_unbracketed_surface_creation_front_and_repair_TO_chain(
            repair_target, target_mass, num_layers, head, cum_layer_thickness_cm,
            soil_type, soil_properties, &adjusted_depth);
        adjusted_depth = true;
        deleted_unbracketed_surface_front = true;
      }

      if (!deleted_unbracketed_surface_front &&
          remaining_mass_error_after_depth > MBAL_ITERATIVE_TOLERANCE &&
          repair_target->next != NULL &&
          repair_target->layer_num == repair_target->next->layer_num &&
          std::fabs(repair_target->depth_cm - repair_target->next->depth_cm) <= CREATION_COLOCATED_TOLERANCE_CM) {
        if (verbosity.compare("high") == 0) {
          printf("surface creation deleting zero-thickness front %d and handing repair to colocated front %d "
                 "(GW=%d, to_bottom=%d, depth %.6f cm).\n",
                 repair_target->front_num, repair_target->next->front_num, repair_target->next->is_WF_GW,
                 repair_target->next->to_bottom, repair_target->next->depth_cm);
        }

        struct wetting_front *coincident_target =
          listDeleteFront(repair_target->front_num, head, soil_type, soil_properties);
        if (coincident_target != NULL) {
          if (coincident_target->is_WF_GW) {
            remaining_mass_error_after_depth =
              lgar_restore_surface_creation_mass_on_connected_gw_chain(coincident_target, target_mass, head,
                                                                       cum_layer_thickness_cm, soil_type,
                                                                       soil_properties);
            if (remaining_mass_error_after_depth > MBAL_ITERATIVE_TOLERANCE) {
              const double residual_before_depth = remaining_mass_error_after_depth;
              remaining_mass_error_after_depth =
                lgar_restore_surface_creation_mass_on_connected_gw_chain_via_depth(coincident_target, target_mass,
                                                                                   num_layers, *head,
                                                                                   cum_layer_thickness_cm);
              adjusted_depth = adjusted_depth ||
                               std::fabs(remaining_mass_error_after_depth - residual_before_depth) > 1.0e-12;
            }
          }
          else if (!coincident_target->is_WF_GW && coincident_target->to_bottom) {
            remaining_mass_error_after_depth =
              lgar_restore_surface_creation_mass_on_connected_surface_chain(coincident_target, target_mass, head,
                                                                            cum_layer_thickness_cm, soil_type,
                                                                            soil_properties);
            if (remaining_mass_error_after_depth > MBAL_ITERATIVE_TOLERANCE) {
              const double residual_before_depth = remaining_mass_error_after_depth;
              remaining_mass_error_after_depth =
                lgar_restore_surface_creation_mass_on_connected_gw_chain_via_depth(coincident_target, target_mass,
                                                                                   num_layers, *head,
                                                                                   cum_layer_thickness_cm);
              adjusted_depth = adjusted_depth ||
                               std::fabs(remaining_mass_error_after_depth - residual_before_depth) > 1.0e-12;
            }
          }
          else {
            remaining_mass_error_after_depth =
              lgar_restore_new_surface_front_mass(coincident_target, target_mass, cum_layer_thickness_cm,
                                                  soil_type, *head, soil_properties);
            if (std::fabs(remaining_mass_error_after_depth) > MBAL_ITERATIVE_TOLERANCE) {
              const double residual_before_depth = remaining_mass_error_after_depth;
              remaining_mass_error_after_depth =
                lgar_restore_new_surface_front_mass_via_depth(coincident_target, target_mass,
                                                              cum_layer_thickness_cm, *head);
              adjusted_depth = adjusted_depth ||
                               std::fabs(remaining_mass_error_after_depth - residual_before_depth) > 1.0e-12;
            }
          }

          remaining_mass_error_after_depth = lgar_route_surface_creation_residual_if_needed(
            remaining_mass_error_after_depth, num_layers, &target_mass, creation_excess_gw_flux_cm,
            creation_excess_runoff_cm, &saturated_creation_gw_flux_capacity_cm, coincident_target, soil_type, soil_properties,
            cum_layer_thickness_cm, *head,
            "surface creation colocated-front handoff could not fully close mass; routing");
        }
      }
    }

    if (verbosity.compare("high") == 0 && std::fabs(remaining_mass_error_after_depth) > 1.0e-10) {
      printf("surface creation normalization left residual mass error of %.12e cm after phase %d\n",
             remaining_mass_error_after_depth, phase + 1);
    }

    if (std::fabs(remaining_mass_error_after_depth) <= MBAL_ITERATIVE_TOLERANCE &&
        !adjusted_depth && !delayed_gw_conversion_this_phase) {
      const int post_repair_correction_type_surf =
        lgarto_correction_type_surf(num_layers, cum_layer_thickness_cm, head,
                                    surface_lower_boundary_depth_cm);
      if (post_repair_correction_type_surf != 0) {
        const double lower_boundary_flux_correction_cm =
          lgar_apply_surface_creation_correction(post_repair_correction_type_surf, num_layers,
                                                 &target_mass, surface_lower_boundary_depth_cm,
                                                 cum_layer_thickness_cm, soil_type,
                                                 frozen_factor, head, soil_properties);
        if (creation_excess_gw_flux_cm != NULL) {
          *creation_excess_gw_flux_cm += lower_boundary_flux_correction_cm;
        }
        continue;
      }

      const int post_repair_correction_type =
        lgarto_correction_type(num_layers, cum_layer_thickness_cm, head,
                               surface_lower_boundary_depth_cm);
      if (post_repair_correction_type != 0) {
        const double lower_boundary_flux_correction_cm =
          lgar_apply_creation_general_correction(post_repair_correction_type, num_layers, &target_mass,
                                                 column_depth, surface_lower_boundary_depth_cm,
                                                 cum_layer_thickness_cm, soil_type, frozen_factor,
                                                 head, soil_properties);
        if (creation_excess_gw_flux_cm != NULL) {
          *creation_excess_gw_flux_cm += lower_boundary_flux_correction_cm;
        }
        continue;
      }

      break;
    }
  }

  lgar_refresh_all_front_states_from_theta(*head, soil_type, frozen_factor, soil_properties);
}

// ######################################################################################
/* This subroutine is called iff there is no surfacial front, it creates a new front and
   inserts ponded depth, and will return some amount if can't fit all water */
// ######################################################################################
extern void lgar_create_surficial_front(bool TO_enabled, int num_layers, double *ponded_depth_cm, double *volin, double dry_depth,
						double theta1, int *soil_type, double *cum_layer_thickness_cm,
						double *frozen_factor, struct wetting_front** head, struct soil_properties_ *soil_properties,
						double *creation_excess_gw_flux_cm,
						double *creation_excess_runoff_cm,
						double saturated_creation_gw_flux_capacity_cm,
						double groundwater_depth_cm)
{
  // into the soil.  Note ponded_depth_cm is a pointer.   Access its value as (*ponded_depth_cm).

  // local vars
  double theta_e,Se,theta_r;
  double delta_theta;
  double vg_alpha_per_cm, vg_m, vg_n, Ksat_cm_per_h;

  bool to_bottom = FALSE;
  struct wetting_front *current;
  int layer_num,soil_num;
  const double prior_mass = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
  const double column_depth = cum_layer_thickness_cm[num_layers];
  double surface_lower_boundary_depth_cm = column_depth;
  if (std::isfinite(groundwater_depth_cm) && groundwater_depth_cm > 0.0 &&
      groundwater_depth_cm < column_depth) {
    surface_lower_boundary_depth_cm = groundwater_depth_cm;
  }

  layer_num = 1;   // we only create new surfacial fronts in the first layer
  soil_num = soil_type[layer_num];

  theta_e = soil_properties[soil_num].theta_e;  // rhs of the new front, assumes theta_e as per Peter
  theta_r = soil_properties[soil_num].theta_r;
  delta_theta =  theta_e - theta1;

  double theta_new;

  if(dry_depth * delta_theta > (*ponded_depth_cm))  // all the ponded depth enters the soil
    {
      *volin = *ponded_depth_cm;
      theta_new = fmin(theta1 + (*ponded_depth_cm) /dry_depth, theta_e);
      current = lgar_insert_surficial_front(TO_enabled, dry_depth, theta_new, layer_num, to_bottom, head);
      *ponded_depth_cm = 0.0;
      //hp_cm =0.0;
    }
  else  // not all ponded depth fits in
    {
      *volin = dry_depth * delta_theta;
      *ponded_depth_cm -= dry_depth * delta_theta;
      theta_new = theta_e; //fmin(theta1 + (*ponded_depth_cm) /dry_depth, theta_e);
      current = lgar_insert_surficial_front(TO_enabled, dry_depth, theta_e, layer_num, to_bottom, head);
      //hp_cm = *ponded_depth_cm;
    }

  vg_alpha_per_cm    = soil_properties[soil_num].vg_alpha_per_cm;
  vg_m               = soil_properties[soil_num].vg_m;
  vg_n               = soil_properties[soil_num].vg_n;
  Ksat_cm_per_h      = soil_properties[soil_num].Ksat_cm_per_h * frozen_factor[layer_num];

  Se = calc_Se_from_theta(theta_new,theta_e,theta_r);
  current->psi_cm = calc_h_from_Se(Se, vg_alpha_per_cm , vg_m, vg_n);

  current->K_cm_per_h = calc_K_from_Se(Se, Ksat_cm_per_h, vg_m) * frozen_factor[layer_num]; // AJ - K_temp in python version for 1st layer

  current->dzdt_cm_per_h = 0.0; //for now assign 0 to dzdt as it will be computed/updated in lgar_dzdt_calc function

  if (verbosity.compare("high") == 0) {
    std::cerr<<"Print after WF creation, no correction yet: \n";
    listPrint(*head);
  }

	  if (TO_enabled && current != NULL) {
	    const double target_mass = prior_mass + *volin;
		    lgar_normalize_after_surface_front_creation(num_layers, target_mass,
		                                                surface_lower_boundary_depth_cm,
		                                                cum_layer_thickness_cm, soil_type,
		                                                frozen_factor, current, head,
		                                                soil_properties, creation_excess_gw_flux_cm,
		                                                creation_excess_runoff_cm,
		                                                saturated_creation_gw_flux_capacity_cm);
	  }

  if (!TO_enabled && current->next!=NULL){// sometimes a new WF immediately has to merge with another WF
    bool had_to_merge = false;
    while ( (current->depth_cm > current->next->depth_cm) && (current->layer_num == current->next->layer_num) && !(current->next->to_bottom)){
      // Technically this should be replaced with the function that iteratively checks for merging, layer crossing, lower boundary crossing, and dry over wet, 
      // but because all we are doing is adding a new WF onto a linked list that will not need correction because correction was just done on it, it could be that all we need to do here is merge
      // because the resulting depths should all be in the top layer
      lgar_merge_wetting_fronts(soil_type, frozen_factor, head, soil_properties);
      had_to_merge = true;
    }
    if (had_to_merge){ //pretty sure this is not necessary but keeping it in
      lgar_wetting_fronts_cross_layer_boundary(num_layers, cum_layer_thickness_cm, soil_type, frozen_factor,
          head, soil_properties);
    }
  }

  return;

}

// ############################################################################################
/* This routine calculates the "dry depth" of a newly created wetting front in the top soil layer after
   a non-rainy period or a big increase in rainrate  on an unsaturated first layer.
   Note: Calculation of the initial depth of a new wetting front in the first layer uses the concept of "dry depth",
   described in the 2015 GARTO paper (Lai et al., An efficient and guaranteed stable numerical method ffor
   continuous modeling of infiltration and redistribution with a shallow dynamic water table). */
// ############################################################################################
extern double lgar_calc_dry_depth(bool TO_enabled, bool use_closed_form_G, int nint, double timestep_h, double *delta_theta, int *soil_type,
				  double *cum_layer_thickness_cm, double *frozen_factor,
				  struct wetting_front* head, struct soil_properties_ *soil_properties)
{

  // local variables
  struct wetting_front *current;
  double theta1,theta2,theta_e,theta_r;
  double vg_alpha_per_cm,vg_n,vg_m,Ksat_cm_per_h,h_min_cm;
  double tau;
  double Geff;
  double dry_depth;
  int    soil_num;
  int    layer_num;

	current=head;

	if (listLength_surface(head) > 0) {
	  while (current != NULL && current->is_WF_GW) {
	    current = current->next;
	  }
	}

	if (current != NULL && current->depth_cm == 0.0) {
	  while (current != NULL && current->depth_cm == 0.0) {
	    current = current->next;
	  }
	}

	if (current == NULL) {
	  current = head;
	}

	layer_num  = current->layer_num;
  soil_num   = soil_type[layer_num];

  // copy values of soil properties into shorter variable names to improve readability
  theta_r         = soil_properties[soil_num].theta_r;
  vg_alpha_per_cm = soil_properties[soil_num].vg_alpha_per_cm;
  vg_m            = soil_properties[soil_num].vg_m;
  vg_n            = soil_properties[soil_num].vg_n;
  Ksat_cm_per_h   = soil_properties[soil_num].Ksat_cm_per_h * frozen_factor[layer_num];
  h_min_cm        = soil_properties[soil_num].h_min_cm;
  double lambda = soil_properties[soil_num].bc_lambda;
  double bc_psib_cm = soil_properties[soil_num].bc_psib_cm;

  // these are the limits of integration
  theta1   = current->theta;                 // water content of the first (most surficial) existing wetting front
  theta_e  = soil_properties[soil_num].theta_e;
  theta2 = theta_e;

  *delta_theta = theta_e - current->theta;  // return the delta_theta value to the calling function

  tau  = timestep_h * Ksat_cm_per_h/(theta_e-current->theta); //3600

  Geff = calc_Geff(use_closed_form_G, theta1, theta2, theta_e, theta_r, vg_alpha_per_cm, vg_n, vg_m, h_min_cm, Ksat_cm_per_h, nint, lambda, bc_psib_cm); 

  // note that dry depth originally has a factor of 0.5 in front
  dry_depth = 0.5 * (tau + sqrt( tau*tau + 4.0*tau*Geff) );

  //when dry depth greater than layer 1 thickness, set dry depth to layer 1 thickness
  dry_depth = fmin(cum_layer_thickness_cm[layer_num], dry_depth);

  if (TO_enabled){
    // dry_depth = dry_depth * 0.99; //0.2
    // dry_depth = dry_depth - CREATION_COLOCATED_TOLERANCE_CM;
  }
  
  return dry_depth;

}

// ###########################################################################
/* function to calculate the amount of soil moisture (total mass of water)
   in the profile (cm) */
// ###########################################################################
double lgar_calc_mass_bal(double *cum_layer_thickness, struct wetting_front* head)
{

  struct wetting_front* current;
  struct wetting_front* next;     // Beware, might get a little confusing, because there exists a next->next

  double sum=0.0;
  double base_depth;
  int layer;

  if(head == NULL) return 0.0;

  current=head;

  do
    {
      layer=current->layer_num;
      base_depth=cum_layer_thickness[layer-1];   // note cum_layer_thickness[0]=0.0;

      if(current->next != NULL) {            // this is not the last entry in the list
	next=current->next;
	if(next->layer_num == current->layer_num)
	  sum += (current->depth_cm - base_depth) * (current->theta - next->theta); // note no need for fabs() here otherwise we get more mass for the case dry-over-wet front within a layer
	else
	  sum += (current->depth_cm - base_depth) * current->theta;
      }
      else { // this is the last entry in the list.  This must be the deepest front in the final layer
	layer=current->layer_num;
	sum+=current->theta * (current->depth_cm - base_depth);
      }

      current=current->next;

    } while(current != NULL );   // putting conditional at end of do looop makes sure it executes once.

  return sum;
}

// ############################################################################################
/* The module reads the soil parameters.
   Open file to read in the van Genuchten parameters for standard soil types*/
// ############################################################################################
extern int lgar_read_vG_param_file(char const* vG_param_file_name, int num_soil_types, double wilting_point_psi_cm,
				    struct soil_properties_ *soil_properties, bool log_mode)
{

  if (verbosity.compare("high") == 0) {
    std::cerr<<"Reading van Genuchten parameters files...\n";
  }

  // local vars
  FILE *in_vG_params_fptr = NULL;
  char jstr[256];
  char soil_name[30];
  // bool error;
  int length;
  int num_soils_in_file = 0;             // soil counter
  int soil = 1;
  double theta_e,theta_r,vg_n,vg_m,vg_alpha_per_cm,Ksat_cm_per_h;  // shorthand variable names
  double m,p,lambda;

  // open the file
  if((in_vG_params_fptr=fopen(vG_param_file_name,"r"))==NULL) {
    printf("Can't open input file named %s. Program stopped.\n",vG_param_file_name); exit(-1);
  }

  fgets(jstr,255,in_vG_params_fptr);   // read the header line and ignore

  //for(soil=1;soil<=num_soil_types;soil++) {// read the num_soil_types lines with data
  while (fgets(jstr,255,in_vG_params_fptr) != NULL) {

    sscanf(jstr,"%s %lf %lf %lf %lf %lf",soil_name,&theta_r,&theta_e,&vg_alpha_per_cm,&vg_n,&Ksat_cm_per_h);
    length=strlen(soil_name);

    if(length>MAX_SOIL_NAME_CHARS) {
      printf("While reading vG soil parameter file: %s, soil name longer than allowed.  Increase MAX_SOIL_NAME_CHARS\n",
	     vG_param_file_name);
      printf("Program stopped.\n");
      exit(0);
    }

    strcpy(soil_properties[soil].soil_name,soil_name);
    if (log_mode){
      vg_alpha_per_cm = pow(10.0, vg_alpha_per_cm);
      Ksat_cm_per_h   = pow(10.0, Ksat_cm_per_h);
    }
    vg_m = 1-1/vg_n;
    soil_properties[soil].theta_r         = theta_r;
    soil_properties[soil].theta_e         = theta_e;
    soil_properties[soil].vg_alpha_per_cm = vg_alpha_per_cm; // cm^(-1)
    soil_properties[soil].vg_n            = vg_n;
    soil_properties[soil].vg_m            = vg_m;
    soil_properties[soil].Ksat_cm_per_h   = Ksat_cm_per_h;
    soil_properties[soil].theta_wp = calc_theta_from_h(wilting_point_psi_cm, vg_alpha_per_cm,
						       vg_m, vg_n, theta_e, theta_r);

    // Given van Genuchten parameters calculate estimates of Brooks & Corey bc_lambda and bc_psib
    if (1.0 < vg_n) {
      m = 1.0 - 1.0 / vg_n;
      p = 1.0 + 2.0 / m;
      soil_properties[soil].bc_lambda  = 2.0 / (p - 3.0);
      soil_properties[soil].bc_psib_cm = (p + 3.0) * (147.8 + 8.1 * p + 0.092 * p * p) /
	(2.0 * vg_alpha_per_cm * p * (p - 1.0) * (55.6 + 7.4 * p + p * p));
      assert(0.0 < soil_properties[soil].bc_psib_cm);
    }
    else {
      fprintf(stderr, "ERROR: van Genuchten parameter n must be greater than 1\n");
      //error = TRUE;  // TODO FIXME - what todo in this ccase?
    }

    /* this is the effective capillary drive after */
    /* Morel-Seytoux et al. (1996) eqn. 13 or 15 */
    /* psi should not be less than this value.  */
    lambda=soil_properties[soil].bc_lambda;
    soil_properties[soil].h_min_cm = soil_properties[soil].bc_psib_cm*(2.0+3.0/lambda)/(1.0+3.0/lambda);
    num_soils_in_file++;
    soil++;

    // the soil_properties array has allocation of num_soil_types, so break if condition is satisfied
    if (num_soil_types == num_soils_in_file)
      break;
  }

  fclose(in_vG_params_fptr);      // close the file, we're done with it

  return num_soils_in_file;
}

// ############################################################################################
/* code to calculate velocity of fronts
   equations with full description are provided in the lgar paper (currently under review) */
// ############################################################################################
extern void lgar_dzdt_calc(bool use_closed_form_G, int nint, int num_layers, double h_p, double subtimestep_h, int *soil_type, double *cum_layer_thickness_cm,
				   double *frozen_factor, struct wetting_front* head, struct soil_properties_ *soil_properties, bool switch_caching, int cache_count, int new_front,
				   double groundwater_depth_cm)
{
  if (verbosity.compare("high") == 0) {
    std::cerr<<"Calculating dz/dt .... \n";
  }

  struct wetting_front* current;
  struct wetting_front* next;
  struct wetting_front* previous;
  struct wetting_front* next_to_use;

  double vg_alpha_per_cm,vg_n,vg_m,Ksat_cm_per_h,theta_e,theta_r;  // local variables to make things clearer
  double delta_theta;
  double Geff;
  double depth_cm;    // the absolute depth down to a wetting front from the surface
  double h_min_cm;
  double K_cm_per_h;  // unsaturated hydraulic conductivity K(theta) at the RHS of the current wetting front (cm/h)
  double theta1, theta2;  // limits of integration on Geff from theta1 to theta2
  double bottom_sum;  // store a running sum of L_n/K(theta_n) n increasing from 1 to N-1, as we go down in layers N
  double dzdt;
  int    soil_num, layer_num;
  int count_correct_for_crossing = 0;


  if(head == NULL) {
    stringstream errMsg;
    errMsg << "lgar derivative function called for empty list (no wetting front exists) \n";
    throw runtime_error(errMsg.str());
  }

  // make sure to use previous state values as current state is updated during the timestep (that's how it is done is Peter's python version)

  current = head;

  do {  // loop through the wetting fronts
    dzdt = 0.0;
    bool apply_layer_crossing_limit = true;

    // copy structure elements into shorter variables names to increase readability
    // WETTING FRONT PROPERTIES
    layer_num    = current->layer_num;    // what layer the front is in
    K_cm_per_h   = current->K_cm_per_h;   // K(theta)

    if (K_cm_per_h < 0) {
      printf("K is negative (layer_num, wf_num, K): %d %d %lf \n", layer_num, current->front_num, K_cm_per_h);
      listPrint(head);
      printf("Is your n value very close to 1? Very small n values can cause K to become 0. \n");
      //The parameter n must physically attain a value greater than 1. However, when n is small, and apparently less than 1.02, sometimes n can make K evaluate to 0, for larger values of psi.
      //So, checking for K_cm_per_h <= 0 has been replaced by checking if K_cm_per_h is negative. K_cm_per_h should never be negative (although perhaps machine precision could make this occur, although we haven't seen it yet), but mathematically can be 0 in some rare cases. 
      abort();
    }

    depth_cm = current->depth_cm;     // absolute Z to this wetting front measured down from land surface

    // SOIL PROPERTIES
    soil_num        = soil_type[layer_num];
    vg_alpha_per_cm = soil_properties[soil_num].vg_alpha_per_cm;
    vg_n            = soil_properties[soil_num].vg_n;
    vg_m            = soil_properties[soil_num].vg_m;
    theta_e         = soil_properties[soil_num].theta_e;
    theta_r         = soil_properties[soil_num].theta_r;
    h_min_cm        = soil_properties[soil_num].h_min_cm;
    Ksat_cm_per_h   = soil_properties[soil_num].Ksat_cm_per_h * frozen_factor[current->layer_num];
    double lambda = soil_properties[soil_num].bc_lambda;
    double bc_psib_cm = soil_properties[soil_num].bc_psib_cm;

    next = current->next;    // the next element in the linked list
    if (next == NULL) break; // we're done calculating dZ/dt's because we're at the end of the list

    if (current->is_WF_GW) {
      apply_layer_crossing_limit = false;

      if (current->to_bottom == TRUE) {
        dzdt = 0.0;
      }
      else {
        next_to_use = next;
        previous = listFindFront(current->front_num - 1, head, NULL);
        const double fixed_column_depth_cm = cum_layer_thickness_cm[num_layers];
        double D =
          (std::isfinite(groundwater_depth_cm) && groundwater_depth_cm > 0.0)
              ? groundwater_depth_cm
              : fixed_column_depth_cm;
        // Mobile groundwater is a hydraulic boundary; the to_bottom scaffold
        // still stays on fixed soil-layer boundaries in this implementation.
        if (D <= current->depth_cm) {
          D = current->depth_cm + fmax(1.0e-9, 1.0e-9 * fixed_column_depth_cm);
        }

        if (next_to_use->is_WF_GW == FALSE) {
          while (next_to_use != NULL && next_to_use->is_WF_GW == FALSE) { //because there technically can be surface WFs between TO WFs just so long as the higher TO WFs have a depth of 0
            next_to_use = next_to_use->next;
          }
          while (next_to_use != NULL && next_to_use->to_bottom == TRUE && next_to_use->next != NULL) {
            next_to_use = next_to_use->next;
          }
        }

        if (next_to_use == NULL) {
          dzdt = 0.0;
        }
        else if (current->layer_num == num_layers) {
          //Because the WF is in the lowest layer, the single layer form of dZ/dt for the TO WF can be used. 
          //note that TO WFs extend from the bottom of the model domain up, whereas surface WFs extend from the top of the model domain down. Therefore, a TO WF spans just 1 layer if it is only in the bottom layer, but a TO WF that is in more than the bottom layer spans more than 1 layer. 
          double avoid_div_by_zero_factor = 1.0E-9;
          if ((next_to_use->theta - current->theta) == 1.0E-9) {
            avoid_div_by_zero_factor = 1.0E-10;
          }

          dzdt = (next_to_use->K_cm_per_h - current->K_cm_per_h) /
                 (next_to_use->theta - current->theta + avoid_div_by_zero_factor) *
                 (1.0 - (next_to_use->psi_cm + 1.0E-6) /
                  (D - current->depth_cm + avoid_div_by_zero_factor));

          if ((dzdt * subtimestep_h + current->depth_cm) > D) {
            dzdt = (D - current->depth_cm - 1.0E-6) / subtimestep_h;
          }

          if (isnan(dzdt)) {
            dzdt = 0.0;
          }
        }
        else { //the TO WF is in a layer that is not the deepest layer, so multilayer composite dzdt has to be calculated
          double K_composite = 0.0;
          double K_composite_left = 0.0;
          double denominator = 0.0;
          double denominator_left = 0.0;

          double avoid_div_by_zero_factor = 1.0E-16;

          for (int k = num_layers; k > (current->layer_num - 1); k--) {
            int soil_num_loc = soil_type[k];
            double theta_prev_loc = calc_theta_from_h(next_to_use->psi_cm, soil_properties[soil_num_loc].vg_alpha_per_cm,
                                                      soil_properties[soil_num_loc].vg_m,
                                                      soil_properties[soil_num_loc].vg_n, soil_properties[soil_num_loc].theta_e,
                                                      soil_properties[soil_num_loc].theta_r);

            double Se_prev_loc = calc_Se_from_theta(theta_prev_loc, soil_properties[soil_num_loc].theta_e,
                                                    soil_properties[soil_num_loc].theta_r);

            double K_cm_per_h_prev_loc = calc_K_from_Se(Se_prev_loc,
                                                        soil_properties[soil_num_loc].Ksat_cm_per_h * frozen_factor[k],
                                                        soil_properties[soil_num_loc].vg_m);
            if (K_cm_per_h_prev_loc==0.0){
              K_cm_per_h_prev_loc = K_cm_per_h_prev_loc + avoid_div_by_zero_factor;
            }
            if (k == current->layer_num) {
              denominator += (cum_layer_thickness_cm[k] - current->depth_cm) / K_cm_per_h_prev_loc;
            }
            else {
              denominator += (cum_layer_thickness_cm[k] - cum_layer_thickness_cm[k - 1]) / K_cm_per_h_prev_loc;
            }
          }

          double temp_psi = current->psi_cm;
          double temp_theta = current->theta;
          bool identical_thetas = false;
          if (current->psi_cm == next_to_use->psi_cm) {
            temp_psi += 0.1;
            identical_thetas = true;
            int soil_num_loc = soil_type[current->layer_num];
            temp_theta = calc_theta_from_h(temp_psi, soil_properties[soil_num_loc].vg_alpha_per_cm,
                                           soil_properties[soil_num_loc].vg_m,
                                           soil_properties[soil_num_loc].vg_n, soil_properties[soil_num_loc].theta_e,
                                           soil_properties[soil_num_loc].theta_r);
          }

          for (int k = num_layers; k > (current->layer_num - 1); k--) {
            int soil_num_loc = soil_type[k];
            double theta_prev_loc = calc_theta_from_h(temp_psi, soil_properties[soil_num_loc].vg_alpha_per_cm,
                                                      soil_properties[soil_num_loc].vg_m,
                                                      soil_properties[soil_num_loc].vg_n, soil_properties[soil_num_loc].theta_e,
                                                      soil_properties[soil_num_loc].theta_r);

            double Se_prev_loc = calc_Se_from_theta(theta_prev_loc, soil_properties[soil_num_loc].theta_e,
                                                    soil_properties[soil_num_loc].theta_r);

            double K_cm_per_h_prev_loc = calc_K_from_Se(Se_prev_loc,
                                                        soil_properties[soil_num_loc].Ksat_cm_per_h * frozen_factor[k],
                                                        soil_properties[soil_num_loc].vg_m);
            if (K_cm_per_h_prev_loc==0.0){
              K_cm_per_h_prev_loc = K_cm_per_h_prev_loc + avoid_div_by_zero_factor;
            }
            if (k == layer_num) {
              if (previous != NULL) {
                denominator_left += (cum_layer_thickness_cm[k] - current->depth_cm) / K_cm_per_h_prev_loc;
              }
              else {
                denominator_left += cum_layer_thickness_cm[k] / K_cm_per_h_prev_loc;
              }
            }
            else {
              denominator_left += (cum_layer_thickness_cm[k] - cum_layer_thickness_cm[k - 1]) / K_cm_per_h_prev_loc;
            }
          }

          if (denominator==0.0){
            denominator = denominator + avoid_div_by_zero_factor;
          }
          if (denominator_left==0.0){
            denominator_left = denominator_left + avoid_div_by_zero_factor;
          }

          K_composite = (D - current->depth_cm) / denominator;
          if (previous == NULL) {
            K_composite_left = D / denominator_left;
          }
          else {
            K_composite_left = (D - current->depth_cm) / denominator_left;
          }

          if (current->layer_num == next_to_use->layer_num) {
            if (!identical_thetas) {
              dzdt = (K_composite - K_composite_left) / (next_to_use->theta - current->theta + 1.0E-9) *
                     (1.0 - next_to_use->psi_cm / (D - current->depth_cm));
            }
            else {
              dzdt = (K_composite - K_composite_left) / (next_to_use->theta - temp_theta + 1.0E-9) *
                     (1.0 - next_to_use->psi_cm / (D - current->depth_cm));
            }
          }
          else {
            int soil_num_loc = soil_type[current->layer_num];
            double theta_temp = calc_theta_from_h(next_to_use->psi_cm, soil_properties[soil_num_loc].vg_alpha_per_cm,
                                                  soil_properties[soil_num_loc].vg_m,
                                                  soil_properties[soil_num_loc].vg_n, soil_properties[soil_num_loc].theta_e,
                                                  soil_properties[soil_num_loc].theta_r);
            if (!identical_thetas) {
              dzdt = (K_composite - K_composite_left) / (theta_temp - current->theta + 1.0E-9) *
                     (1.0 - next_to_use->psi_cm / (D - current->depth_cm));
            }
            else {
              dzdt = (K_composite - K_composite_left) / (next_to_use->theta - theta_temp + 1.0E-9) *
                     (1.0 - next_to_use->psi_cm / (D - current->depth_cm));
            }
          }
        }

        if ((current->depth_cm + dzdt * subtimestep_h) < 0.0) {
          dzdt = current->depth_cm / subtimestep_h;
        }

        if (current->next != NULL && current->next->theta == current->theta) {
          dzdt = 0.0;
        }

        if (next_to_use != NULL) {
          const double hydrostatic_target_depth_cm = D - next_to_use->psi_cm;
          dzdt = lgar_limit_TO_dzdt_near_hydrostatic_target(
              dzdt, current->depth_cm, hydrostatic_target_depth_cm, D, subtimestep_h);
        }
      }
    }
    else {
      theta1 = next->theta;
      theta2 = current->theta;

      bottom_sum = 0.0;  // needed ffor multi-layered dz/dt equation.  Equal to sum from n=1 to N-1 of (L_n/K_n(theta_n))

      if(current->to_bottom == TRUE) {
        if(layer_num > 1)
	  current->dzdt_cm_per_h = 0.0;
        else
	  current->dzdt_cm_per_h = 0.0;

        current = current->next;  // point to the next link
        continue;                 // go to next front, this one fully penetrates the layer
      }
      else if(layer_num > 1) {
        bottom_sum += (current->depth_cm-cum_layer_thickness_cm[layer_num-1])/K_cm_per_h;
      }

      if(theta1 > theta2) {
        printf("Calculating dzdt : theta1 > theta2 = (%lf, %lf) ... aborting \n", theta1, theta2);  // this should never happen
        exit(0);
      }

      Geff = calc_Geff(use_closed_form_G, theta1, theta2, theta_e, theta_r, vg_alpha_per_cm, vg_n, vg_m, h_min_cm, Ksat_cm_per_h, nint, lambda, bc_psib_cm); 
      delta_theta = current->theta - next->theta;
      const double saturated_pressure_head_cm =
        lgar_saturated_pressure_head_cm(current, soil_type, soil_properties);

      if(current->layer_num == 1) { // this front is in the upper layer
        if (delta_theta > 0){
	  dzdt = 1.0/delta_theta*(Ksat_cm_per_h*(Geff+h_p+saturated_pressure_head_cm)/current->depth_cm+current->K_cm_per_h);}
        else{
	  dzdt = 0.0;}
      }
      else {  // we are in the second or greater layer
        double denominator = bottom_sum;

        for (int k = 1; k < layer_num; k++) {
	  int soil_num_loc = soil_type[layer_num-k]; // _loc denotes the soil_num is local to this loop
	  double theta_prev_loc = calc_theta_from_h(current->psi_cm, soil_properties[soil_num_loc].vg_alpha_per_cm,
						    soil_properties[soil_num_loc].vg_m,
						    soil_properties[soil_num_loc].vg_n,soil_properties[soil_num_loc].theta_e,
						    soil_properties[soil_num_loc].theta_r);


	  double Se_prev_loc = calc_Se_from_theta(theta_prev_loc,soil_properties[soil_num_loc].theta_e,soil_properties[soil_num_loc].theta_r);

	  double K_cm_per_h_prev_loc = calc_K_from_Se(Se_prev_loc,soil_properties[soil_num_loc].Ksat_cm_per_h * frozen_factor[layer_num-k],
						      soil_properties[soil_num_loc].vg_m);

	  denominator += (cum_layer_thickness_cm[k] - cum_layer_thickness_cm[k-1])/ K_cm_per_h_prev_loc;

        }

        double numerator = depth_cm;// + (Geff +h_p)* Ksat_cm_per_h / K_cm_per_h;

        if (delta_theta > 0)
	  dzdt = (1.0/delta_theta) * ( (numerator / denominator) + Ksat_cm_per_h*(Geff+h_p+saturated_pressure_head_cm)/depth_cm );
        else
	  dzdt = 0.0;
      }

      if ((dzdt == 0.0) && (current->to_bottom==FALSE)){
        //in lgar_move, we have: "if (current->dzdt_cm_per_h == 0.0 && current->to_bottom == FALSE) { // a new front was just created, so don't update it."
        //the issue here is that when theta approaches theta_r, then dzdt can in some cases numerically evaluate to 0, even if the wetting front has to_bottom==FALSE.
        //so, there are cases where a WF should be moving very slowly, but not being completely still. 
        dzdt = 1e-9;
      }
    }

    if (dzdt>1e4){//insanity check
      dzdt = 1e4;
    }
    
    double largest_K_s = 0.0;
    for (int ii=1; ii<=num_layers; ii++) {
      int soil_num = soil_type[ii];
      double temp_K_s = soil_properties[soil_num].Ksat_cm_per_h;
      largest_K_s = fmax(temp_K_s, largest_K_s);
    }

    if (dzdt>100*largest_K_s){//insanity check; was 1E4 but now addtionally defining based on K_s
      dzdt = 100*largest_K_s;
    }

    if (dzdt<-100*largest_K_s){//insanity check; was -1E4 (LGARTO thing)
      dzdt = -100*largest_K_s;
    }

    bool is_next_to_bottom = false;
    if (next){
      is_next_to_bottom = next->to_bottom;
    }

    if (apply_layer_crossing_limit && is_next_to_bottom && dzdt*subtimestep_h + current->depth_cm > FACTOR_LIMITS_LAYER_CROSSING_SPEED*(cum_layer_thickness_cm[current->layer_num])){
      count_correct_for_crossing ++;
      // idea is that we want to make sure that layer boundary crossing does not go crazy, should never cross more than
      // some factor of the current layer depth into the next layer. In most situations it won't matter, but in for example sand
      // over clay, can have the case that layer boundary crossing immediately goes super deep because dzdt is large for sand,
      // and this very sudden expansion into clay for a single time step is unrealistic 
      // this will probably not be desirable for the fracture domain and is more of a soil matrix property 
      dzdt = fabs(current->depth_cm - FACTOR_LIMITS_LAYER_CROSSING_SPEED*(cum_layer_thickness_cm[current->layer_num]))/subtimestep_h + count_correct_for_crossing*DEPTH_AVOIDS_SAME_WF_DEPTH; //count_correct_for_crossing is just so that a bunch of WFs don't end up at exactly the same depth 
    }

    if (switch_caching){
      if (current->front_num!=new_front){
        dzdt = dzdt*(cache_count);
      }
    }

    current->dzdt_cm_per_h = dzdt;

    current = current->next;  // point to the next link

    if (verbosity.compare("high") == 0) {
      printf("dzdt: %lf \n", dzdt);
    }

  } while(current != NULL );   // putting conditional at end of do looop makes sure it executes at least once

}

// ############################################################################################
/* The function does mass balance for a wetting front to get an updated theta.
   The head (psi) value is iteratively altered until the error between prior mass and new mass
   is within a tolerance. This is only used updating WF theta after it moves.*/
// ############################################################################################
extern double lgar_theta_mass_balance(int layer_num, int soil_num, double psi_cm, double new_mass,
				      double prior_mass, double precip_mass_to_add, double *AET_demand_cm, double *delta_theta, double *delta_thickness,
				      int *soil_type, struct soil_properties_ *soil_properties,
                      bool allow_legacy_aet_bookkeeping_adjustment,
                      double psi_upper_limit_cm)
{

  const double psi_upper_limit_for_solve =
    (std::isfinite(psi_upper_limit_cm) && psi_upper_limit_cm > PSI_UPPER_LIM)
      ? psi_upper_limit_cm
      : PSI_UPPER_LIM;
  double psi_cm_loc = psi_cm; // location psi
  double delta_mass = fabs(new_mass - prior_mass); // mass different between the new and prior
  // double original_delta_mass = delta_mass; 

  double factor = fmax(1.0,psi_cm/10.0);//was 1.0 previously. This code is far faster and seems to avoid loops with >10000 iterations.
  if (psi_cm>1.E4){ // in very dry cases, mass conservation will take longer to achieve if using a small factor
    factor = psi_cm;
  }
  bool switched = false; // flag that determines capillary head to be incremented or decremented

  double theta             = 0; // this will be updated and returned
  double psi_cm_loc_prev   = psi_cm_loc;
  double delta_mass_prev   = delta_mass;
  int count_no_mass_change = 0;
  int break_no_mass_change = 5;
  bool wanted_to_saturate_flag = FALSE;
  
  // check if the difference is less than the tolerance
  if (delta_mass <= MBAL_ITERATIVE_TOLERANCE) {
    theta = calc_theta_from_h(psi_cm_loc, soil_properties[soil_num].vg_alpha_per_cm,
			      soil_properties[soil_num].vg_m, soil_properties[soil_num].vg_n,
			      soil_properties[soil_num].theta_e,soil_properties[soil_num].theta_r);
    return theta;
  }

  // the loop increments/decrements the capillary head until mass difference between
  // the new and prior is within the tolerance
  int iter = 0;
  bool iter_aug_flag = FALSE;
  bool iter_aug_flag_extreme = FALSE;
  int first_speedup_thresh = MAX_ITER_MBAL_LOOP/100;
  int second_speedup_thresh = MAX_ITER_MBAL_LOOP/10;

  while (delta_mass > MBAL_ITERATIVE_TOLERANCE) {
    iter++;

    if (iter>first_speedup_thresh && iter_aug_flag==FALSE){
      factor = factor*100;
      iter_aug_flag = TRUE;
    }

    if (iter>second_speedup_thresh && iter_aug_flag_extreme==FALSE){
      factor = factor*100;
      iter_aug_flag_extreme = TRUE;
    }

    if (new_mass > prior_mass) {
      psi_cm_loc += 0.1 * factor;
      switched = false;
    }
    else {
      if (!switched) {
	switched = true;
	factor = factor * 0.1;
      }
      
      psi_cm_loc_prev = psi_cm_loc;
      psi_cm_loc -= 0.1 * factor;

      if (psi_cm_loc<0.0){
        wanted_to_saturate_flag = TRUE;
      }
      
      if (psi_cm_loc < 0 && psi_cm_loc_prev != 0) {
	/* this is for the extremely rare case when iterative psi_cm_loc calculation temporarily
	   yields a negative value and the actual answer for psi_cm_loc is nonzero. For example
	   when a completely saturated wetting front with a tiny amount of ET should yield a resulting
	   theta that is slightly below saturation. */
        psi_cm_loc = psi_cm_loc_prev * 0.1;
      }
      
    }

    if (psi_cm_loc<0.0){ // addresses the case where psi_cm_loc was 0 at the start so psi_cm_loc_prev is 0
      psi_cm_loc = 0.0;
      factor = factor * 0.5; // allow factor to smoothly approach 0 in case the solution is slightly above psi = 0
    }

    double theta_layer;
    double mass_layers= 0.0;

    theta = calc_theta_from_h(psi_cm_loc, soil_properties[soil_num].vg_alpha_per_cm, soil_properties[soil_num].vg_m,
			      soil_properties[soil_num].vg_n,soil_properties[soil_num].theta_e,
			      soil_properties[soil_num].theta_r);

    mass_layers += delta_thickness[layer_num] * (theta - delta_theta[layer_num]);

    for (int k=1; k<layer_num; k++) {
      int soil_num_loc =  soil_type[k]; // _loc denotes the variable is local to the loop

      theta_layer = calc_theta_from_h(psi_cm_loc, soil_properties[soil_num_loc].vg_alpha_per_cm,
				      soil_properties[soil_num_loc].vg_m, soil_properties[soil_num_loc].vg_n,
				      soil_properties[soil_num_loc].theta_e, soil_properties[soil_num_loc].theta_r);

      mass_layers += delta_thickness[k] * (theta_layer - delta_theta[k]);
    }

    new_mass = mass_layers;
    delta_mass = fabs(new_mass - prior_mass);

    // stop the loop if the error between the current and previous psi is less than 10^-15
    // 1. enough accuracy, 2. the algorithm can't improve the error further,
    // 3. avoid infinite loop, 4. handles case where theta is very close to theta_r and convergence might be possible but would be extremely slow
    // 5. handles a corner case when prior mass is tiny (e.g., <1.E-5)
    // printf("A1 = %.20f, %.18f %.18f %.18f %.18f \n ",fabs(psi_cm_loc - psi_cm_loc_prev) , psi_cm_loc, psi_cm_loc_prev, factor, delta_mass);
    
    if (fabs(psi_cm_loc - psi_cm_loc_prev) < 1E-15 && factor < 1E-13) break;

    // another condition to avoid infinite loop when the error does not improve
    if (fabs(delta_mass - delta_mass_prev) < 1E-15)
      count_no_mass_change++;
    else
      count_no_mass_change = 0;

    // break the loop if the mass does not change in the five consecutive iterations.
    if (count_no_mass_change == break_no_mass_change && precip_mass_to_add < 1.E-12){ // made it so that this check only occurs if there is no infiltration, because there is a case where precip on a very dry wetting front will need more than 5 iterations to have its mass change at all
      break;
    }

    if (iter>MAX_ITER_MBAL_LOOP){ //limit the total number of iterations 
      break;
    }

    if ( (psi_cm_loc > psi_upper_limit_for_solve) && (iter > first_speedup_thresh) ){ //unrealistic pressures, but there are some cases where convergence is possible even at large psi values, and there is a case where AET, free drainage, or WF movement can bring psi above PSI_UPPER_LIM, so we do want to allow a few iterations
      psi_cm_loc = psi_upper_limit_for_solve;
      theta = calc_theta_from_h(psi_cm_loc, soil_properties[soil_num].vg_alpha_per_cm, soil_properties[soil_num].vg_m,
				soil_properties[soil_num].vg_n,soil_properties[soil_num].theta_e,
				soil_properties[soil_num].theta_r);
      break;
    }

    // -ve pressure will return NAN, so terminate the loop if previous psi is way small and current psi is zero
    // the wetting front is almost saturated
    if (psi_cm_loc <= 0 && psi_cm_loc_prev < 0) break;

    delta_mass_prev = delta_mass;

  }

  //There is a rare case where mass balance closure would require that theta<theta_r. 
  //However, the above loop can never increase psi to the point where theta<theta_r, because theta must always be between theta_r and theta_r, because of the van Genuchten model (calc_theta_from_h).
  //If we get to the case where theta<theta_r would be necessary for mass balance closure, then the above loop will break before delta_mass <= tolerance.
  //In this case, the remaining mass balance error is put into AET. This should usually be acceptable, because it will often be the AET flux that would have been the trigger for theta < theta_r for mass conservation, so reducing AET works in this case.
  if (allow_legacy_aet_bookkeeping_adjustment &&
      (delta_mass > MBAL_ITERATIVE_TOLERANCE) && (!wanted_to_saturate_flag)){//the second condition is necessary because count_no_mass_change == break_no_mass_change in the loop above will trigger when the model approaches saturation; in this event the extra water should go into runoff or recharge (handled elsewhere), because the soil saturates, rather than AET
    *AET_demand_cm = *AET_demand_cm - fabs(delta_mass - MBAL_ITERATIVE_TOLERANCE);
  }

  if (allow_legacy_aet_bookkeeping_adjustment &&
      (theta>=soil_properties[soil_num].theta_e) && (psi_cm_loc!=0.0) ){
    //addresses a very rare case. Sometimes when psi gets very close to 0 but is not 0, calc_theta_from_h will actually yield theta_e for a very small nonzero psi value (for example psi=1e-3 or something like that).
    //This can happen for example when the model domain is very close to saturation, and the number of WFs == the number of layers, but there is a little bit of AET so the resulting model state should have just slightly less water than complete saturation.
    //However, layers above the current one might not have the property that this small zonzero psi value yields theta = theta_e.
    //This leads to the case where the mass balance correctly closes, but with theta=theta_e for the current layer and not with theta = theta_e for some higher layer(s).
    //Then, later, the psi value for the current layer is set to 0.0 because its theta value is theta_e, and then the layers above will have their psi values set to 0.0, and then theta = theta_e everywhere in the model domain, whereas it should account for the recent small AET that happened.
    //Just setting the AET to 0 in these cases closes mass balance; another option would be to augment the recharge. Error is very rare and seems to happen once every 100k parameter sets or so, using yearlong simulations.
    *AET_demand_cm = 0.0;
  }
  
  return theta;

}

static bool lgar_restore_theta_merge_mass_via_depth(double target_mass, double column_depth,
					       double *cum_layer_thickness_cm, struct wetting_front **head,
					       struct wetting_front *front_to_adjust)
{
  if (front_to_adjust == NULL) {
    return false;
  }

  if (front_to_adjust->to_bottom) {
    /*
     * to_bottom fronts are structural layer-boundary markers.  They may share
     * theta/psi with a connected TO chain, but their depth is not a valid mass
     * correction degree of freedom; return the residual to the caller instead.
     */
    return false;
  }

  const double initial_depth_cm = front_to_adjust->depth_cm;
  double lower_depth_cm = cum_layer_thickness_cm[front_to_adjust->layer_num - 1];
  double upper_depth_cm = fmin(column_depth, cum_layer_thickness_cm[front_to_adjust->layer_num]);

  if (front_to_adjust->front_num > 1) {
    struct wetting_front *previous = listFindFront(front_to_adjust->front_num - 1, *head, NULL);
    if (previous != NULL) {
      lower_depth_cm = fmax(lower_depth_cm, previous->depth_cm + DEPTH_AVOIDS_SAME_WF_DEPTH);
    }
  }

  if (front_to_adjust->next != NULL && front_to_adjust->next->layer_num == front_to_adjust->layer_num) {
    upper_depth_cm = fmin(upper_depth_cm, front_to_adjust->next->depth_cm - DEPTH_AVOIDS_SAME_WF_DEPTH);
  }

  if (upper_depth_cm <= lower_depth_cm) {
    front_to_adjust->depth_cm = initial_depth_cm;
    return false;
  }

  const double original_depth_cm = fmax(lower_depth_cm, fmin(initial_depth_cm, upper_depth_cm));
  front_to_adjust->depth_cm = original_depth_cm;
  const double original_mass_cm = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
  if (std::fabs(target_mass - original_mass_cm) <= MBAL_ITERATIVE_TOLERANCE) {
    return true;
  }

  front_to_adjust->depth_cm = lower_depth_cm;
  const double mass_at_lower_cm = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
  front_to_adjust->depth_cm = upper_depth_cm;
  const double mass_at_upper_cm = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);

  const bool monotonic_increasing = mass_at_upper_cm > mass_at_lower_cm + MBAL_ITERATIVE_TOLERANCE;
  const bool monotonic_decreasing = mass_at_lower_cm > mass_at_upper_cm + MBAL_ITERATIVE_TOLERANCE;
  if (!monotonic_increasing && !monotonic_decreasing) {
    front_to_adjust->depth_cm = original_depth_cm;
    return false;
  }

  double best_depth_cm = original_depth_cm;
  double best_mass_error_cm = std::fabs(target_mass - original_mass_cm);

  auto maybe_update_best = [&](double probe_depth_cm) {
    front_to_adjust->depth_cm = probe_depth_cm;
    const double probe_mass_cm = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
    const double probe_mass_error_cm = std::fabs(target_mass - probe_mass_cm);
    if (probe_mass_error_cm < best_mass_error_cm) {
      best_mass_error_cm = probe_mass_error_cm;
      best_depth_cm = probe_depth_cm;
    }
    return probe_mass_cm;
  };

  const double min_mass_cm = fmin(mass_at_lower_cm, mass_at_upper_cm);
  const double max_mass_cm = fmax(mass_at_lower_cm, mass_at_upper_cm);
  if (target_mass <= min_mass_cm) {
    best_depth_cm = (mass_at_lower_cm <= mass_at_upper_cm) ? lower_depth_cm : upper_depth_cm;
    front_to_adjust->depth_cm = best_depth_cm;
    return true;
  }
  if (target_mass >= max_mass_cm) {
    best_depth_cm = (mass_at_lower_cm >= mass_at_upper_cm) ? lower_depth_cm : upper_depth_cm;
    front_to_adjust->depth_cm = best_depth_cm;
    return true;
  }

  double bracket_lo_cm = lower_depth_cm;
  double bracket_hi_cm = upper_depth_cm;
  for (int iter = 0; iter < 80; iter++) {
    const double probe_depth_cm = 0.5 * (bracket_lo_cm + bracket_hi_cm);
    const double probe_mass_cm = maybe_update_best(probe_depth_cm);

    if (std::fabs(target_mass - probe_mass_cm) <= MBAL_ITERATIVE_TOLERANCE) {
      best_depth_cm = probe_depth_cm;
      break;
    }

    if (monotonic_increasing) {
      if (probe_mass_cm < target_mass) {
        bracket_lo_cm = probe_depth_cm;
      }
      else {
        bracket_hi_cm = probe_depth_cm;
      }
    }
    else {
      if (probe_mass_cm > target_mass) {
        bracket_lo_cm = probe_depth_cm;
      }
      else {
        bracket_hi_cm = probe_depth_cm;
      }
    }
  }

  front_to_adjust->depth_cm = best_depth_cm;
  return true;
}

static void lgar_restore_upstream_to_bottom_chain_psi(struct wetting_front *surviving_front,
					       double psi_cm,
					       int *soil_type,
					       struct soil_properties_ *soil_properties,
					       struct wetting_front *head)
{
  if (surviving_front == NULL || surviving_front->front_num <= 1 || head == NULL) {
    return;
  }

  struct wetting_front *current = listFindFront(surviving_front->front_num - 1, head, NULL);
  while (current != NULL && current->to_bottom) {
    current->psi_cm = psi_cm;
    const int soil_num = soil_type[current->layer_num];
    current->theta = calc_theta_from_h(psi_cm,
                                       soil_properties[soil_num].vg_alpha_per_cm,
                                       soil_properties[soil_num].vg_m,
                                       soil_properties[soil_num].vg_n,
                                       soil_properties[soil_num].theta_e,
                                       soil_properties[soil_num].theta_r);
    const double Se = calc_Se_from_theta(current->theta,
                                         soil_properties[soil_num].theta_e,
                                         soil_properties[soil_num].theta_r);
    current->K_cm_per_h = calc_K_from_Se(Se,
                                         soil_properties[soil_num].Ksat_cm_per_h,
                                         soil_properties[soil_num].vg_m);

    if (current->front_num == 1) {
      break;
    }
    current = listFindFront(current->front_num - 1, head, NULL);
  }
}

extern double lgarto_TO_WFs_merge_via_theta(double target_mass, double column_depth,
					    double *cum_layer_thickness_cm, struct wetting_front **head,
					    int *soil_type, struct soil_properties_ *soil_properties)
{
  if (verbosity.compare("high") == 0) {
    printf("before lgarto_TO_WFs_merge_via_theta ... \n");
    listPrint(*head);
  }

  struct wetting_front *current = *head;
  struct wetting_front *next = current != NULL ? current->next : NULL;
  struct wetting_front *previous = *head;
  double mass_diff = 0.0;

  for (int wf = 1; wf != listLength(*head); wf++) {
    if (current == NULL || current->next == NULL) {
      break;
    }

    if ((current->theta > next->theta) && current->is_WF_GW && next->is_WF_GW &&
        (next->layer_num == current->layer_num)) {
      if (current->depth_cm == 0.0) {
        current = listDeleteFront(current->front_num, head, soil_type, soil_properties);
        mass_diff = target_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
      }
      else {
        if (next->to_bottom == FALSE && previous->to_bottom == FALSE) {
          current = listDeleteFront(current->front_num, head, soil_type, soil_properties);

          if (current != NULL) {
            /*
             * If the theta merge cannot close mass by moving the surviving
             * front within its own layer, return the residual as a signed flux
             * correction.  The old unbounded loop could satisfy local mass by
             * pushing a non-to_bottom scaffold through lower layers, leaving an
             * impossible TO geometry for later corrections to repair.
             */
            (void) lgar_restore_theta_merge_mass_via_depth(target_mass, column_depth,
                                                           cum_layer_thickness_cm, head,
                                                           current);
            mass_diff = target_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
          }
	    }
	        else {
	          const double deleted_psi_cm = current->psi_cm;
	          const bool upstream_to_bottom_chain = (previous != NULL && previous->to_bottom);
	          current = listDeleteFront(current->front_num, head, soil_type, soil_properties);
	          if (current != NULL) {
	            if (upstream_to_bottom_chain) {
	              lgar_restore_upstream_to_bottom_chain_psi(current, deleted_psi_cm, soil_type,
	                                                       soil_properties, *head);
	            }
	            lgar_theta_mass_balance_correction(false, current->front_num, target_mass, head,
	                                              cum_layer_thickness_cm, soil_type, soil_properties);
	            if (std::fabs(target_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, *head)) >
	                MBAL_ITERATIVE_TOLERANCE) {
	              (void) lgar_restore_theta_merge_mass_via_depth(target_mass, column_depth,
	                                                            cum_layer_thickness_cm, head, current);
	            }
	          }
	          mass_diff = target_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
	        }
	      }
	    }

    previous = current;
    if (previous == NULL || previous->front_num == listLength(*head)) {
      break;
    }
    current = current->next;
    next = current != NULL ? current->next : NULL;
  }

  if (verbosity.compare("high") == 0) {
    printf("after lgarto_TO_WFs_merge_via_theta ... \n");
    listPrint(*head);
  }

  return mass_diff;
}

extern bool correct_close_psis(int *soil_type, struct soil_properties_ *soil_properties,
			       struct wetting_front **head)
{
  if (verbosity.compare("high") == 0) {
    printf("before correct_close_psis: \n");
    listPrint(*head);
  }

  bool close_psis = false;
  struct wetting_front *current = *head;
  struct wetting_front *next = current != NULL ? current->next : NULL;

  for (int wf = 1; wf != listLength(*head); wf++) {
    if (current == NULL || next == NULL) {
      break;
    }

    if ((current->layer_num == next->layer_num) && (current->is_WF_GW == FALSE) &&
        (next->is_WF_GW == TRUE) && (std::fabs(current->psi_cm - next->psi_cm) < 1.E-3)) {
      current = listDeleteFront(current->front_num, head, soil_type, soil_properties);
      close_psis = true;
      break;
    }

    current = next;
    next = current->next;
  }

  if (verbosity.compare("high") == 0) {
    printf("after correct_close_psis: \n");
    listPrint(*head);
  }

  return close_psis;
}

// this function handles wetting front merging and layer boundary crossing iteratively, so we no longer need to rely on a predetermined order of these events. Theoretically this is slightly faster and more stable.
// the function has the name "lgarto" rather than "lgar" because future work will further incorporate lgarto
extern int lgarto_correction_type(int num_layers, double* cum_layer_thickness_cm, struct wetting_front** head,
                                  double vadose_lower_boundary_depth_cm){
  if (*head == NULL || (*head)->next == NULL) {
    return 0;
  }

  int correction_type = 0;
  struct wetting_front *current = *head;
  struct wetting_front *next = current->next;
  struct wetting_front *next_to_next = NULL;
  bool TO_layer_cross = FALSE;
  const double fixed_lower_boundary_depth_cm = cum_layer_thickness_cm[num_layers];
  const double effective_lower_boundary_depth_cm =
    (std::isfinite(vadose_lower_boundary_depth_cm) && vadose_lower_boundary_depth_cm > 0.0)
      ? fmin(fixed_lower_boundary_depth_cm, vadose_lower_boundary_depth_cm)
      : fixed_lower_boundary_depth_cm;

  if (next != NULL) {
    next_to_next = current->next->next;
  }

  for (int wf = 1; wf != listLength(*head); wf++) {
    struct wetting_front *previous = listFindFront(wf - 1, *head, NULL);

    const double boundary_tol_cm = fmax(1.0e-8, 10.0 * TRUNCATION_DEPTH);
    if ((current->is_WF_GW == FALSE) && (!current->to_bottom) &&
        (current->layer_num > 1) &&
        (current->depth_cm >= 0.0) &&
        (current->depth_cm < cum_layer_thickness_cm[current->layer_num - 1] - boundary_tol_cm)) {
      correction_type = 10; // rare upward surface WF layer-boundary crossing after remap/merge
      break;
    }

    if (next != NULL) {
      if (lgarto_is_last_layer_GW_overshoot_pair(current, num_layers)) {
        correction_type = 9;
        break;
      }
      if ((current->depth_cm > next->depth_cm) && (current->is_WF_GW == FALSE) &&
          (next->is_WF_GW == TRUE) && (next->to_bottom == FALSE)) {
        correction_type = 1;
      }
      if ((current->depth_cm > next->depth_cm) && (current->is_WF_GW == TRUE) &&
          (next->is_WF_GW == TRUE) && (next->to_bottom == FALSE) &&
          (current->to_bottom == FALSE) && (current->theta < next->theta)) {
        correction_type = 2;
      }
      if ((current->theta > next->theta) && (current->is_WF_GW == TRUE) &&
          (next->is_WF_GW == TRUE) && (next->layer_num == current->layer_num)) {
        correction_type = 3;
      }
    }

    if (previous != NULL) {
      if ((current->depth_cm < cum_layer_thickness_cm[current->layer_num - 1]) &&
          (previous->to_bottom == TRUE) && (current->layer_num > previous->layer_num) &&
          (current->is_WF_GW == TRUE)) {
        correction_type = 4;
        TO_layer_cross = TRUE;
      }
    }

    if (next != NULL) {
      if ((current->depth_cm > cum_layer_thickness_cm[current->layer_num]) &&
          (next->to_bottom == TRUE) && (current->is_WF_GW == TRUE) &&
          (current->layer_num != num_layers)) {
        correction_type = 4;
        TO_layer_cross = TRUE;
      }
    }

    if (next != NULL) {
      if ((current->is_WF_GW == FALSE) && (next->is_WF_GW == FALSE) &&
          (current->theta > next->theta) && (current->depth_cm > next->depth_cm) &&
          (current->layer_num == next->layer_num) && (!next->to_bottom)) {
        correction_type = 5;
      }
      if ((current->is_WF_GW == FALSE) &&
          (current->depth_cm > cum_layer_thickness_cm[current->layer_num]) &&
          (next->depth_cm == cum_layer_thickness_cm[current->layer_num]) &&
          (current->theta > next->theta) && (current->layer_num != num_layers)) {
        correction_type = 6;
      }
    }

    const double current_layer_lower_boundary_cm =
      fmin(cum_layer_thickness_cm[current->layer_num], effective_lower_boundary_depth_cm);
    if ((next_to_next == NULL) && (current->depth_cm > current_layer_lower_boundary_cm) &&
        (current->is_WF_GW == FALSE)) {
      correction_type = 7;
    }

    if (next != NULL && current->is_WF_GW == FALSE && next->is_WF_GW == TRUE &&
        current->theta < next->theta && current->layer_num == next->layer_num &&
        correction_type == 0) {
      correction_type = 8;
    }

    if (next != NULL && current->to_bottom == TRUE && current->is_WF_GW == FALSE &&
        current->next->is_WF_GW == TRUE) {
      correction_type = 8;
    }

    current = next;
    if (current == NULL) {
      break;
    }
    next = current->next;
    if (next != NULL) {
      next_to_next = current->next->next;
    }
  }

  if (correction_type == 9 || correction_type == 10) {
    if (verbosity.compare("high") == 0) {
      printf("computed correction type: %d \n", correction_type);
    }
    return correction_type;
  }

  if (TO_layer_cross) {
    correction_type = 4;
  }
  if (verbosity.compare("high") == 0) {
    printf("computed correction type: %d \n", correction_type);
  }
  return correction_type;
}

extern int lgarto_correction_type_surf(int num_layers, double* cum_layer_thickness_cm, struct wetting_front** head,
                                       double vadose_lower_boundary_depth_cm){
  int correction_type_surf = 0;
  struct wetting_front *current = *head;
  struct wetting_front *next = current->next;
  struct wetting_front *next_to_next = NULL;
  bool lgarto_active = false;
  const double fixed_lower_boundary_depth_cm = cum_layer_thickness_cm[num_layers];
  const double effective_lower_boundary_depth_cm =
    (std::isfinite(vadose_lower_boundary_depth_cm) && vadose_lower_boundary_depth_cm > 0.0)
      ? fmin(fixed_lower_boundary_depth_cm, vadose_lower_boundary_depth_cm)
      : fixed_lower_boundary_depth_cm;
  for (struct wetting_front *front = *head; front != NULL; front = front->next) {
    if (front->is_WF_GW) {
      lgarto_active = true;
      break;
    }
  }
  // // this will be necessary for lgarto
  // struct wetting_front *top_most_TO_front_below_surfs = NULL;
  // double top_most_TO_front_below_surfs_psi_cm = 1.E16;
  // if (listLength(*head)>(listLength_surface(*head)+listLength_TO_WFs_above_surface_WFs(*head))){
  //   top_most_TO_front_below_surfs = listFindFront(listLength_surface(*head)+listLength_TO_WFs_above_surface_WFs(*head) + 1, *head, NULL);
  //   top_most_TO_front_below_surfs_psi_cm = top_most_TO_front_below_surfs->psi_cm;
  // }

  if (next!=NULL){
    next_to_next = current->next->next;
  }

  for (int wf = 1; wf != (listLength(*head)); wf++) {

    if ( (current->depth_cm > 0.0) && (current->is_WF_GW) ){
      break;
    }

    const double boundary_tol_cm = fmax(1.0e-8, 10.0 * TRUNCATION_DEPTH);
    if ((current->is_WF_GW == FALSE) && (!current->to_bottom) &&
        (current->layer_num > 1) &&
        (current->depth_cm >= 0.0) &&
        (current->depth_cm < cum_layer_thickness_cm[current->layer_num - 1] - boundary_tol_cm)) {
      correction_type_surf = 9; // rare upward surface WF layer-boundary crossing after remap/merge
      break;
    }

    if (next!=NULL){
	      if (lgarto_active && lgarto_surface_front_overtook_surface_front_above_TO_chain(current)) {
	        correction_type_surf = 6; // mixed surface/surface/TO overtake correction
	        break;
	      }
		      if ( (current->is_WF_GW==0) && (next->is_WF_GW==1) && (current->depth_cm > next->depth_cm) && (!next->to_bottom) && (next->next != NULL) ){
		        correction_type_surf = 7; // surface/TO merge before TO/GW motion reopens the crossing gap
		        break;
		      }
	      if ( (current->is_WF_GW==0) && (next->is_WF_GW==0) && (current->theta>next->theta) && (current->depth_cm > next->depth_cm) && (current->layer_num == next->layer_num) && (!next->to_bottom) ){
	      // if ( (current->theta>next->theta) && (current->depth_cm > next->depth_cm) && (current->layer_num == next->layer_num) && (!next->to_bottom) ){
		        correction_type_surf = 1; //this is surface-surface WF merging
        break;
      }
      // if ( (current->is_WF_GW==0) && (current->depth_cm > cum_layer_thickness_cm[current->layer_num]) && (next->depth_cm == cum_layer_thickness_cm[current->layer_num]) && (current->theta>next->theta) && (current->layer_num!=num_layers) && (current->psi_cm<top_most_TO_front_below_surfs_psi_cm) ){
      if ( (current->is_WF_GW==0) && (current->depth_cm > cum_layer_thickness_cm[current->layer_num]) && (next->depth_cm == cum_layer_thickness_cm[current->layer_num]) && (!lgarto_active || next->to_bottom) && (current->theta>next->theta) && (current->layer_num!=num_layers) ){
      // if ( (current->depth_cm > cum_layer_thickness_cm[current->layer_num]) && (next->depth_cm == cum_layer_thickness_cm[current->layer_num] && (next->to_bottom)) && (current->theta>next->theta) && (current->layer_num!=num_layers) ){
        correction_type_surf = 2; //this is surface WF layer bdy crossing
        break;
      }
    }
    const double current_layer_lower_boundary_cm =
      fmin(cum_layer_thickness_cm[current->layer_num], effective_lower_boundary_depth_cm);
    if ( (next_to_next == NULL) && (current->depth_cm > current_layer_lower_boundary_cm) && (current->is_WF_GW==0) ){
    // if ( (next_to_next == NULL) && (current->depth_cm > cum_layer_thickness_cm[current->layer_num]) ){
      correction_type_surf = 3; //this is a surface WF crossing the model lower bdy
      break;
    }
    if (lgar_check_dry_over_wet_wetting_fronts(*head)){
      correction_type_surf = 4;
      break;
    }

    /*
     * Leave surface-to-TO conversion for the TO/general correction loop
     * (correction type 8), after TO motion and AET have had a chance to run
     * on the pre-conversion scaffold.
     */

    current = next;
    next = current->next;
    if (next!=NULL){
      next_to_next = current->next->next;
    }
  }

  if (verbosity.compare("high") == 0){
    printf("computed correction type for surface WFs: %d \n", correction_type_surf);
  }
  return correction_type_surf;
}

static int lgar_count_leading_zero_depth_groundwater_fronts(const struct wetting_front *head)
{
  int count = 0;
  const struct wetting_front *current = head;
  while (current != NULL && current->is_WF_GW == TRUE && current->layer_num == 1 &&
         fabs(current->depth_cm) <= CREATION_COLOCATED_TOLERANCE_CM) {
    count++;
    current = current->next;
  }
  return count;
}

static bool lgar_are_leading_zero_depth_groundwater_fronts_redundant(const struct wetting_front *front_a,
                                                                     const struct wetting_front *front_b)
{
  if (front_a == NULL || front_b == NULL) {
    return false;
  }

  const bool both_leading_zero_depth_groundwater =
    front_a->is_WF_GW == TRUE && front_b->is_WF_GW == TRUE &&
    front_a->layer_num == 1 && front_b->layer_num == 1 &&
    fabs(front_a->depth_cm) <= CREATION_COLOCATED_TOLERANCE_CM &&
    fabs(front_b->depth_cm) <= CREATION_COLOCATED_TOLERANCE_CM;

  if (!both_leading_zero_depth_groundwater) {
    return false;
  }

  const bool same_theta = fabs(front_a->theta - front_b->theta) < THRESHOLD_NO_MOISTURE_DIFF;
  const bool same_psi = fabs(front_a->psi_cm - front_b->psi_cm) < 1.0e-3;
  const bool same_to_bottom = (front_a->to_bottom == front_b->to_bottom);

  return same_to_bottom && (same_theta || same_psi);
}

static void lgar_trim_excess_leading_zero_depth_groundwater_fronts(struct wetting_front **head,
                                                                   int *soil_type,
                                                                   struct soil_properties_ *soil_properties)
{
  if (head == NULL || *head == NULL) {
    return;
  }

  while (*head != NULL && (*head)->next != NULL &&
         lgar_are_leading_zero_depth_groundwater_fronts_redundant(*head, (*head)->next)) {
    if (verbosity.compare("high") == 0) {
      printf("Trimming redundant leading zero-depth TO/GW front with nearly identical state.\n");
    }

    listDeleteFront(1, head, soil_type, soil_properties);
  }
}

static void lgar_cap_leading_zero_depth_groundwater_fronts(struct wetting_front **head,
                                                           int *soil_type,
                                                           struct soil_properties_ *soil_properties)
{
  if (head == NULL || *head == NULL) {
    return;
  }

  const int max_leading_zero_depth_groundwater_fronts = 6; //controls num of 0 depth GW fronts
  int leading_zero_depth_groundwater_fronts =
    lgar_count_leading_zero_depth_groundwater_fronts(*head);

  while (leading_zero_depth_groundwater_fronts > max_leading_zero_depth_groundwater_fronts) {
    if (verbosity.compare("high") == 0) {
      printf("Trimming excess leading zero-depth TO/GW front. Count %d exceeds cap %d.\n",
             leading_zero_depth_groundwater_fronts,
             max_leading_zero_depth_groundwater_fronts);
    }

    listDeleteFront(1, head, soil_type, soil_properties);
    leading_zero_depth_groundwater_fronts =
      lgar_count_leading_zero_depth_groundwater_fronts(*head);
  }
}

static void lgar_dry_zero_depth_groundwater_fronts_to_surface_profile(
  struct wetting_front **head,
  int *soil_type,
  struct soil_properties_ *soil_properties)
{
  if (head == NULL || *head == NULL || listLength_surface(*head) == 0) {
    return;
  }

  double driest_surface_psi_cm = -1.0;
  for (struct wetting_front *current = *head; current != NULL; current = current->next) {
    if (current->is_WF_GW == FALSE && std::isfinite(current->psi_cm)) {
      driest_surface_psi_cm = fmax(driest_surface_psi_cm, current->psi_cm);
    }
  }

  if (driest_surface_psi_cm < 0.0) {
    return;
  }

  double target_zero_depth_TO_psi_cm = std::nextafter(driest_surface_psi_cm, HUGE_VAL);
  if (!std::isfinite(target_zero_depth_TO_psi_cm)) {
    target_zero_depth_TO_psi_cm = driest_surface_psi_cm;
  }

  for (struct wetting_front *current = *head; current != NULL; current = current->next) {
    if (current->is_WF_GW == FALSE || fabs(current->depth_cm) > 1.0e-12) {
      continue;
    }

    if (std::isfinite(current->psi_cm) &&
        current->psi_cm >= target_zero_depth_TO_psi_cm) {
      continue;
    }

    /* A zero-depth TO/GW support is metadata for the top of the TO chain,
       not a finite storage volume. Once surface fronts exist at the end of
       a correction pass, those zero-depth supports must be at least as dry as
       the driest active surface front; otherwise hidden wet TO metadata can
       sit above a drier surface profile without changing storage. */
    current->psi_cm = target_zero_depth_TO_psi_cm;
    const int soil_num = soil_type[current->layer_num];
    const double theta_e = soil_properties[soil_num].theta_e;
    const double theta_r = soil_properties[soil_num].theta_r;
    const double vg_a = soil_properties[soil_num].vg_alpha_per_cm;
    const double vg_m = soil_properties[soil_num].vg_m;
    const double vg_n = soil_properties[soil_num].vg_n;
    current->theta = calc_theta_from_h(current->psi_cm, vg_a, vg_m, vg_n, theta_e, theta_r);
    const double Se = calc_Se_from_theta(current->theta, theta_e, theta_r);
    current->K_cm_per_h = calc_K_from_Se(Se, soil_properties[soil_num].Ksat_cm_per_h, vg_m);
  }
}

static void lgar_dry_zero_depth_groundwater_fronts_to_surface_TO_chain(
  struct wetting_front **head,
  int *soil_type,
  struct soil_properties_ *soil_properties,
  double *cum_layer_thickness_cm)
{
  if (head == NULL || *head == NULL || soil_type == NULL ||
      soil_properties == NULL || cum_layer_thickness_cm == NULL) {
    return;
  }

  const lgarto_surface_TO_support_ordering_state state =
    lgarto_find_surface_TO_support_ordering_state(*head);
  if (lgarto_zero_depth_TO_support_ordering_is_valid(state) ||
      state.first_TO_below_surface == NULL) {
    return;
  }

  const double target_zero_depth_TO_psi_cm = state.first_TO_below_surface->psi_cm;
  if (!std::isfinite(target_zero_depth_TO_psi_cm) ||
      target_zero_depth_TO_psi_cm < 0.0) {
    return;
  }

  const double mass_before_cm = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);

  for (struct wetting_front *current = *head; current != NULL; current = current->next) {
    if (current->is_WF_GW == FALSE) {
      break;
    }

    if (current->is_WF_GW == FALSE || current->to_bottom == TRUE ||
        current->layer_num != 1 ||
        std::fabs(current->depth_cm) > ZERO_DEPTH_TO_DELETE_DEPTH_TOL_CM ||
        (std::isfinite(current->psi_cm) &&
         current->psi_cm >= target_zero_depth_TO_psi_cm)) {
      continue;
    }

    /* Zero-depth TO/GW supports carry no finite storage. When surface fronts
       exist, they should be metadata for the finite TO chain below the surface
       stack, not a hidden wetter state above it. Dry only those metadata fronts
       and leave the finite first-below-surface TO front unchanged. */
    current->psi_cm = target_zero_depth_TO_psi_cm;
    const int soil_num = soil_type[current->layer_num];
    const double theta_e = soil_properties[soil_num].theta_e;
    const double theta_r = soil_properties[soil_num].theta_r;
    const double vg_a = soil_properties[soil_num].vg_alpha_per_cm;
    const double vg_m = soil_properties[soil_num].vg_m;
    const double vg_n = soil_properties[soil_num].vg_n;
    current->theta = calc_theta_from_h(current->psi_cm, vg_a, vg_m, vg_n, theta_e, theta_r);
    const double Se = calc_Se_from_theta(current->theta, theta_e, theta_r);
    current->K_cm_per_h = calc_K_from_Se(Se, soil_properties[soil_num].Ksat_cm_per_h, vg_m);
  }

  const double mass_after_cm = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
  const double mass_change_cm = mass_after_cm - mass_before_cm;
  if (!std::isfinite(mass_change_cm) ||
      std::fabs(mass_change_cm) > ZERO_DEPTH_TO_DELETE_MASS_TOL_CM) {
    fprintf(stderr,
            "Error: drying zero-depth TO/GW supports to the first below-surface TO front changed profile storage by %.17g cm.\n"
            "  tolerance=%.17g cm target_zero_depth_TO_psi_cm=%.17g\n"
            "  Wetting front list follows:\n",
            mass_change_cm,
            ZERO_DEPTH_TO_DELETE_MASS_TOL_CM,
            target_zero_depth_TO_psi_cm);
    fflush(stderr);
    listPrint(*head);
    fflush(stdout);
    abort();
  }

  if (verbosity.compare("high") == 0) {
    printf("Dried leading zero-depth TO/GW supports to psi %.17g cm to match "
           "the first TO front below active surface fronts (mass change %.17g cm).\n",
           target_zero_depth_TO_psi_cm,
           mass_change_cm);
  }
}

static bool lgar_is_deletable_zero_depth_groundwater_support(const struct wetting_front *front)
{
  return front != NULL &&
         front->is_WF_GW == TRUE &&
         front->to_bottom == FALSE &&
         front->layer_num == 1 &&
         fabs(front->depth_cm) <= ZERO_DEPTH_TO_DELETE_DEPTH_TOL_CM;
}

static struct wetting_front *lgar_unlink_zero_depth_groundwater_support(
  struct wetting_front *front,
  struct wetting_front **head)
{
  if (front == NULL || head == NULL || *head == NULL) {
    return front;
  }

  struct wetting_front *previous = NULL;
  struct wetting_front *current = *head;
  while (current != NULL && current != front) {
    previous = current;
    current = current->next;
  }

  if (current == NULL) {
    return front;
  }

  struct wetting_front *after_delete = current->next;
  if (previous == NULL) {
    *head = after_delete;
  }
  else {
    previous->next = after_delete;
  }

  free(current);
  for (struct wetting_front *renumber = after_delete; renumber != NULL;
       renumber = renumber->next) {
    renumber->front_num--;
  }

  return after_delete;
}

static struct wetting_front *lgar_delete_zero_depth_groundwater_support_if_mass_neutral(
  struct wetting_front *front,
  struct wetting_front **head,
  int *soil_type,
  struct soil_properties_ *soil_properties,
  double *cum_layer_thickness_cm,
  const char *reason)
{
  if (front == NULL || head == NULL || *head == NULL || cum_layer_thickness_cm == NULL) {
    return front;
  }
  (void)soil_type;
  (void)soil_properties;

  const int deleted_front_num = front->front_num;
  const double mass_before_cm = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
  struct wetting_front *after_delete =
    lgar_unlink_zero_depth_groundwater_support(front, head);
  const double mass_after_cm = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
  const double mass_change_cm = mass_after_cm - mass_before_cm;

  if (verbosity.compare("high") == 0) {
    printf("Deleted zero-depth TO/GW front %d: %s (mass change %.17g cm).\n",
           deleted_front_num, reason, mass_change_cm);
  }

  if (std::fabs(mass_change_cm) > ZERO_DEPTH_TO_DELETE_MASS_TOL_CM) {
    fprintf(stderr,
            "Error: zero-depth TO/GW support deletion changed profile storage by %.17g cm.\n"
            "  deleted_front_num=%d reason=%s tolerance=%.17g cm\n"
            "  Wetting front list after deletion follows:\n",
            mass_change_cm,
            deleted_front_num,
            reason,
            ZERO_DEPTH_TO_DELETE_MASS_TOL_CM);
    fflush(stderr);
    listPrint(*head);
    fflush(stdout);
    abort();
  }

  return after_delete;
}

static bool lgar_zero_depth_groundwater_support_duplicates_next_TO_state(
  const struct wetting_front *front)
{
  if (!lgar_is_deletable_zero_depth_groundwater_support(front) ||
      front->next == NULL ||
      front->next->is_WF_GW == FALSE ||
      !std::isfinite(front->theta) ||
      !std::isfinite(front->next->theta) ||
      !std::isfinite(front->psi_cm) ||
      !std::isfinite(front->next->psi_cm)) {
    return false;
  }

  const double theta_tol = 1.0e-12;
  const double psi_scale_cm =
    fmax(1.0, fmax(std::fabs(front->psi_cm), std::fabs(front->next->psi_cm)));
  const double psi_tol_cm = fmax(1.0e-8, 1.0e-12 * psi_scale_cm);

  return std::fabs(front->theta - front->next->theta) <= theta_tol &&
         std::fabs(front->psi_cm - front->next->psi_cm) <= psi_tol_cm;
}

static void lgar_delete_zero_depth_groundwater_supports_matching_next_TO_state(
  struct wetting_front **head,
  int *soil_type,
  struct soil_properties_ *soil_properties,
  double *cum_layer_thickness_cm)
{
  if (head == NULL || *head == NULL || cum_layer_thickness_cm == NULL) {
    return;
  }

  bool deleted_front = true;
  while (deleted_front) {
    deleted_front = false;
    for (struct wetting_front *current = *head; current != NULL && current->next != NULL;
         current = current->next) {
      if (!lgar_zero_depth_groundwater_support_duplicates_next_TO_state(current)) {
        continue;
      }

      /* A zero-depth non-to_bottom TO/GW front has no finite storage interval.
         If it duplicates the hydraulic state of the next TO/GW front, keeping
         it only gives later TO-only routines an artificial support front to
         select. Delete it only after verifying the removal is mass-neutral. */
      lgar_delete_zero_depth_groundwater_support_if_mass_neutral(
        current, head, soil_type, soil_properties, cum_layer_thickness_cm,
        "zero-depth TO/GW support duplicated the next TO/GW front state");
      deleted_front = true;
      break;
    }
  }
}

static void lgar_delete_zero_depth_groundwater_supports_causing_same_layer_TO_psi_reversal(
  struct wetting_front **head,
  int *soil_type,
  struct soil_properties_ *soil_properties,
  double *cum_layer_thickness_cm)
{
  if (head == NULL || *head == NULL || cum_layer_thickness_cm == NULL) {
    return;
  }

  bool deleted_front = true;
  while (deleted_front) {
    deleted_front = false;
    for (struct wetting_front *current = *head; current != NULL && current->next != NULL;
         current = current->next) {
      struct wetting_front *next = current->next;
      if (!lgar_is_deletable_zero_depth_groundwater_support(current) ||
          next->is_WF_GW == FALSE ||
          next->layer_num != current->layer_num ||
          next->depth_cm <= current->depth_cm + ZERO_DEPTH_TO_DELETE_DEPTH_TOL_CM ||
          !std::isfinite(current->psi_cm) ||
          !std::isfinite(next->psi_cm)) {
        continue;
      }

      const double psi_tol_cm = lgar_psi_assertion_tolerance_cm(current->psi_cm, next->psi_cm);
      if (next->psi_cm <= current->psi_cm + psi_tol_cm) {
        continue;
      }

      /* This is the same structural case caught by lgar_assert_to_psi_monotonic_with_depth,
         but the offending shallower front is a zero-depth non-storage TO/GW support.
         Deleting that support is a local metadata cleanup; widening the TO ordering
         assertion would hide real dry-over-wet errors in finite intervals. */
      lgar_delete_zero_depth_groundwater_support_if_mass_neutral(
        current, head, soil_type, soil_properties, cum_layer_thickness_cm,
        "zero-depth TO/GW support caused same-layer TO psi ordering reversal");
      deleted_front = true;
      break;
    }
  }
}

static void lgar_delete_redundant_finite_same_layer_TO_fronts(
  struct wetting_front **head,
  int *soil_type,
  struct soil_properties_ *soil_properties,
  double *cum_layer_thickness_cm)
{
  if (head == NULL || *head == NULL || cum_layer_thickness_cm == NULL) {
    return;
  }

  bool deleted_front = true;
  while (deleted_front) {
    deleted_front = false;
    for (struct wetting_front *current = *head; current != NULL && current->next != NULL;
         current = current->next) {
      struct wetting_front *next = current->next;
      if (current->is_WF_GW == FALSE || next->is_WF_GW == FALSE ||
          current->to_bottom == TRUE || next->to_bottom == TRUE ||
          current->layer_num != next->layer_num ||
          current->depth_cm >= next->depth_cm - ZERO_DEPTH_TO_DELETE_DEPTH_TOL_CM ||
          current->depth_cm <= cum_layer_thickness_cm[current->layer_num - 1] +
                               ZERO_DEPTH_TO_DELETE_DEPTH_TOL_CM ||
          !std::isfinite(current->theta) || !std::isfinite(next->theta) ||
          !std::isfinite(current->psi_cm) || !std::isfinite(next->psi_cm)) {
        continue;
      }

      const bool same_theta =
        std::fabs(current->theta - next->theta) <= 1.0e-12;
      const bool same_psi =
        std::fabs(current->psi_cm - next->psi_cm) <=
        lgar_psi_assertion_tolerance_cm(current->psi_cm, next->psi_cm);
      if (!same_theta || !same_psi) {
        continue;
      }

      /* Two adjacent finite TO/GW fronts in the same layer with the same
         hydraulic state do not define a real discontinuity. Delete the
         shallower front only, and verify that the operation is mass-neutral. */
      const int deleted_front_num = current->front_num;
      const double mass_before_cm = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
      listDeleteFront(current->front_num, head, soil_type, soil_properties);
      const double mass_after_cm = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
      const double mass_change_cm = mass_after_cm - mass_before_cm;

      if (verbosity.compare("high") == 0) {
        printf("Deleted redundant finite same-layer TO/GW front %d "
               "(shallower duplicate; mass change %.17g cm).\n",
               deleted_front_num, mass_change_cm);
      }

      if (!std::isfinite(mass_change_cm) ||
          std::fabs(mass_change_cm) > FINITE_SAME_LAYER_TO_DELETE_MASS_TOL_CM) {
        fprintf(stderr,
                "Error: finite same-layer TO/GW duplicate deletion changed profile storage by %.17g cm.\n"
                "  deleted_front_num=%d tolerance=%.17g cm\n"
                "  Wetting front list after deletion follows:\n",
                mass_change_cm,
                deleted_front_num,
                FINITE_SAME_LAYER_TO_DELETE_MASS_TOL_CM);
        fflush(stderr);
        listPrint(*head);
        fflush(stdout);
        abort();
      }

      deleted_front = true;
      break;
    }
  }
}

extern void lgar_clean_redundant_fronts(struct wetting_front** head, int *soil_type,
                                        struct soil_properties_ *soil_properties,
                                        bool apply_zero_depth_groundwater_cap,
                                        double *cum_layer_thickness_cm,
                                        double domain_depth_cm){
  (void)domain_depth_cm;

  if (verbosity.compare("high") == 0) {
    printf("before lgar_clean_redundant_fronts: \n");
    listPrint(*head);
  }

  lgar_trim_excess_leading_zero_depth_groundwater_fronts(head, soil_type, soil_properties);
  if (apply_zero_depth_groundwater_cap) {
    lgar_cap_leading_zero_depth_groundwater_fronts(head, soil_type, soil_properties);
  }
  lgar_dry_zero_depth_groundwater_fronts_to_surface_profile(head, soil_type, soil_properties);
  lgar_dry_zero_depth_groundwater_fronts_to_surface_TO_chain(
    head, soil_type, soil_properties, cum_layer_thickness_cm);
  lgar_delete_zero_depth_groundwater_supports_matching_next_TO_state(
    head, soil_type, soil_properties, cum_layer_thickness_cm);
  /* Do not delete very dry zero-depth supports just because surface fronts exist:
     surface cleanup can later fold those fronts into TO/GW, and the zero-depth
     scaffold is needed to keep rainfall entry from seeing a finite-depth head.
     Only delete a zero-depth support when it is the local cause of a same-layer
     TO psi-ordering reversal. */
  lgar_delete_zero_depth_groundwater_supports_causing_same_layer_TO_psi_reversal(
    head, soil_type, soil_properties, cum_layer_thickness_cm);
  lgar_delete_redundant_finite_same_layer_TO_fronts(
    head, soil_type, soil_properties, cum_layer_thickness_cm);

  if (*head == NULL || (*head)->next == NULL) {
    if (verbosity.compare("high") == 0) {
      printf("after lgar_clean_redundant_fronts: \n");
      listPrint(*head);
    }
    return;
  }

  struct wetting_front *current;
  struct wetting_front *next;
  current = *head;
  next = current->next;
  for (int wf = 1; wf != (listLength(*head)); wf++) {
    const bool same_layer = (current->layer_num == next->layer_num);
    const bool same_theta = (fabs(current->theta - next->theta) < THRESHOLD_NO_MOISTURE_DIFF);
    const bool colocated_in_depth = (fabs(current->depth_cm - next->depth_cm) <= CREATION_COLOCATED_TOLERANCE_CM);

    if (same_layer && same_theta && colocated_in_depth &&
        current->to_bottom == FALSE && next->to_bottom == FALSE) {
      // Only delete fronts when they are both moisture-redundant and
      // essentially co-located in depth. Deleting separated fronts can alter
      // connected TO/to_bottom chains through later psi-continuity repair.
      current = listDeleteFront(current->front_num, head, soil_type, soil_properties); 
      break;
    }

    current = next;
    next = current->next;
  }

  if (verbosity.compare("high") == 0) {
    printf("after lgar_clean_redundant_fronts: \n");
    listPrint(*head);
  }
}

extern double calc_min_water_possible_for_free_drainage_wetting_front(int wf_free_drainage, struct wetting_front** head, int *soil_type, struct soil_properties_ *soil_properties){
	/* The region of the soil column from which AET and free drainage are extracted is equal to the most surficial region sharing a single psi value, which 
     can span multiple layers. This function calculates the minimum amount of water that this region can hold, and AET and free drainage will be augmented such that
     they can not yield a theta value below the threshold for maximum psi. Initially this just checked such that storage would not go below theta_r, but for consistency with
     the mass balance loops that augment theta and psi values, this should check that we do not exceed our maximum psi value. */
  double min_storage = 0.0;
  double previous_depth = 0.0;

  struct wetting_front *current;
  current = *head;
  int layer_num;
  int soil_num;
  double theta_r;
  double theta_e;
  double vg_n;
  double vg_m;
  double vg_a;
  while (current!=NULL){
    layer_num  = current->layer_num;
    soil_num   = soil_type[layer_num];
    theta_r    = soil_properties[soil_num].theta_r;
    theta_e    = soil_properties[soil_num].theta_e;
    vg_a       = soil_properties[soil_num].vg_alpha_per_cm;
    vg_m       = soil_properties[soil_num].vg_m;
    vg_n       = soil_properties[soil_num].vg_n;
    double min_theta = calc_theta_from_h(PSI_UPPER_LIM, vg_a, vg_m, vg_n, theta_e, theta_r);
    min_storage += min_theta*(current->depth_cm - previous_depth);
    if (current->front_num==wf_free_drainage){
      break;
    }
    previous_depth = current->depth_cm;
    current = current->next;
  }
  return min_storage;

}

extern double calc_storage_in_free_drainage_wetting_front(int wf_free_drainage, struct wetting_front** head){
  /* The region of the soil column from which AET and free drainage are extracted is equal to the most surficial region sharing a single psi value, which 
     can span multiple layers. This function calculates the amount of water that this region currently holds, and AET and free drainage will be augmented such that
     they can not yield a theta value below theta_r. */
  double storage = 0.0;
  double previous_depth = 0.0;

  struct wetting_front *current;
  current = *head;
  while (current!=NULL){
    storage += current->theta*(current->depth_cm - previous_depth);
    if (current->front_num==wf_free_drainage){
      break;
    }
    previous_depth = current->depth_cm;
    current = current->next;
  }
  return storage;
}

// ############################################################################################
/* The function does mass balance for a wetting front to get an updated theta.
   The head (psi) value is iteratively altered until the error between prior mass and new mass
   is within a tolerance. This is only called to correct psi and theta in the case of dry
   over wet wetting fronts, or after layer boundary crossing. It is simpler than 
   lgar_theta_mass_balance because it does not need information about old WFs or external fluxes
   and is called far less often.*/
// ############################################################################################
static double lgar_apply_theta_mass_balance_correction_psi(
  bool use_dry_over_wet, double psi_cm_loc, struct wetting_front *current,
  struct wetting_front **head, double *cum_layer_thickness_cm, int *soil_type,
  struct soil_properties_ *soil_properties)
{
  psi_cm_loc = fmax(0.0, fmin(PSI_UPPER_LIM, psi_cm_loc));

  int layer_num = current->layer_num;
  int soil_num = soil_type[layer_num];
  current->psi_cm = psi_cm_loc;
  current->theta = calc_theta_from_h(psi_cm_loc,
                                     soil_properties[soil_num].vg_alpha_per_cm,
                                     soil_properties[soil_num].vg_m,
                                     soil_properties[soil_num].vg_n,
                                     soil_properties[soil_num].theta_e,
                                     soil_properties[soil_num].theta_r);

  struct wetting_front *next = current->next;
  struct wetting_front *before_next = current;
  bool skip_bottom_chain_below = false;
  if (next){
    skip_bottom_chain_below = use_dry_over_wet && next->to_bottom && !current->to_bottom;
  }

  if (next && !skip_bottom_chain_below){
    double theta_e;
    double theta_r;
    double vg_a;
    double vg_m;
    double vg_n;
    while (next->to_bottom){
      next->psi_cm = before_next->psi_cm;
      layer_num = next->layer_num;
      soil_num = soil_type[layer_num];
      theta_e = soil_properties[soil_num].theta_e;
      theta_r = soil_properties[soil_num].theta_r;
      vg_a = soil_properties[soil_num].vg_alpha_per_cm;
      vg_m = soil_properties[soil_num].vg_m;
      vg_n = soil_properties[soil_num].vg_n;

      next->theta = calc_theta_from_h(next->psi_cm, vg_a, vg_m, vg_n, theta_e, theta_r);

      before_next = next;
      next = next->next;
      if (!next){
        break;
      }
    }
    if (next){
      if (!next->to_bottom && before_next->to_bottom){
        next->psi_cm = before_next->psi_cm;
        layer_num = next->layer_num;
        soil_num = soil_type[layer_num];
        theta_e = soil_properties[soil_num].theta_e;
        theta_r = soil_properties[soil_num].theta_r;
        vg_a = soil_properties[soil_num].vg_alpha_per_cm;
        vg_m = soil_properties[soil_num].vg_m;
        vg_n = soil_properties[soil_num].vg_n;

        next->theta = calc_theta_from_h(next->psi_cm, vg_a, vg_m, vg_n, theta_e, theta_r);
      }
    }
  }

  return lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
}

static double lgar_theta_mass_balance_correction_with_min_psi(
  bool use_dry_over_wet, int front_num, double prior_mass, struct wetting_front** head,
  double *cum_layer_thickness_cm, int *soil_type,
  struct soil_properties_ *soil_properties, double min_psi_cm)
{
  if (head == NULL || *head == NULL) {
    return prior_mass;
  }

  struct wetting_front *current = listFindFront(front_num, *head, NULL);
  if (current == NULL) {
    return prior_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
  }

  if (current->front_num > 1) {
    struct wetting_front *front = listFindFront(current->front_num - 1, *head, NULL);
    if (front != NULL && front->to_bottom) {
      current = front;
    }
  }
  while (current->to_bottom) {
    if (current->front_num == 1) {
      break;
    }
    current = listFindFront(current->front_num - 1, *head, NULL);
    if (current == NULL) {
      return prior_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
    }
    if (!current->to_bottom) {
      current = current->next;
      break;
    }
  }

  const double mass_at_original = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
  if (fabs(mass_at_original - prior_mass) <= MBAL_ITERATIVE_TOLERANCE) {
    return prior_mass - mass_at_original;
  }

  const double psi_wet_bound_cm =
    fmax(0.0, fmin(PSI_UPPER_LIM, std::isfinite(min_psi_cm) ? min_psi_cm : 0.0));
  const double original_psi_cm = current->psi_cm;
  const double mass_at_wet_bound =
    lgar_apply_theta_mass_balance_correction_psi(use_dry_over_wet, psi_wet_bound_cm,
                                                 current, head, cum_layer_thickness_cm,
                                                 soil_type, soil_properties);
  const double mass_at_dry_limit =
    lgar_apply_theta_mass_balance_correction_psi(use_dry_over_wet, PSI_UPPER_LIM,
                                                 current, head, cum_layer_thickness_cm,
                                                 soil_type, soil_properties);

  double best_psi_cm = psi_wet_bound_cm;
  double best_abs_error_cm = fabs(mass_at_wet_bound - prior_mass);
  if (fabs(mass_at_dry_limit - prior_mass) < best_abs_error_cm) {
    best_psi_cm = PSI_UPPER_LIM;
    best_abs_error_cm = fabs(mass_at_dry_limit - prior_mass);
  }

  if (std::isfinite(original_psi_cm) &&
      original_psi_cm >= psi_wet_bound_cm &&
      original_psi_cm <= PSI_UPPER_LIM) {
    const double mass_at_bounded_original =
      lgar_apply_theta_mass_balance_correction_psi(use_dry_over_wet, original_psi_cm,
                                                   current, head, cum_layer_thickness_cm,
                                                   soil_type, soil_properties);
    if (fabs(mass_at_bounded_original - prior_mass) < best_abs_error_cm) {
      best_psi_cm = original_psi_cm;
      best_abs_error_cm = fabs(mass_at_bounded_original - prior_mass);
    }
  }

  const double min_mass_cm = fmin(mass_at_wet_bound, mass_at_dry_limit);
  const double max_mass_cm = fmax(mass_at_wet_bound, mass_at_dry_limit);
  const bool target_bracketed =
    prior_mass >= min_mass_cm - MBAL_ITERATIVE_TOLERANCE &&
    prior_mass <= max_mass_cm + MBAL_ITERATIVE_TOLERANCE &&
    fabs(mass_at_dry_limit - mass_at_wet_bound) > MBAL_ITERATIVE_TOLERANCE;

  if (target_bracketed) {
    double bracket_lo_psi_cm = psi_wet_bound_cm;
    double bracket_hi_psi_cm = PSI_UPPER_LIM;
    const bool mass_increases_with_psi = mass_at_dry_limit > mass_at_wet_bound;

    for (int iter = 0; iter < MAX_ITER_MBAL_LOOP; iter++) {
      const double probe_psi_cm = 0.5 * (bracket_lo_psi_cm + bracket_hi_psi_cm);
      const double probe_mass_cm =
        lgar_apply_theta_mass_balance_correction_psi(use_dry_over_wet, probe_psi_cm,
                                                     current, head,
                                                     cum_layer_thickness_cm,
                                                     soil_type, soil_properties);
      const double probe_abs_error_cm = fabs(probe_mass_cm - prior_mass);
      if (probe_abs_error_cm < best_abs_error_cm) {
        best_psi_cm = probe_psi_cm;
        best_abs_error_cm = probe_abs_error_cm;
      }
      if (probe_abs_error_cm <= MBAL_ITERATIVE_TOLERANCE) {
        break;
      }

      if (mass_increases_with_psi) {
        if (probe_mass_cm < prior_mass) {
          bracket_lo_psi_cm = probe_psi_cm;
        }
        else {
          bracket_hi_psi_cm = probe_psi_cm;
        }
      }
      else {
        if (probe_mass_cm > prior_mass) {
          bracket_lo_psi_cm = probe_psi_cm;
        }
        else {
          bracket_hi_psi_cm = probe_psi_cm;
        }
      }
    }
  }

  lgar_apply_theta_mass_balance_correction_psi(use_dry_over_wet, best_psi_cm,
                                               current, head, cum_layer_thickness_cm,
                                               soil_type, soil_properties);
  const double residual_cm = prior_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
  if (verbosity.compare("high") == 0 &&
      best_abs_error_cm > MBAL_ITERATIVE_TOLERANCE) {
    printf("Bounded theta mass-balance correction left residual %.12e cm "
           "(mass_at_min_psi %.17lf, mass_at_psi_upper %.17lf, target %.17lf, "
           "min_psi %.17lf).\n",
           residual_cm,
           mass_at_wet_bound, mass_at_dry_limit, prior_mass,
           psi_wet_bound_cm);
  }

  return residual_cm;
}

extern void lgar_theta_mass_balance_correction(bool use_dry_over_wet, int front_num, double prior_mass, struct wetting_front** head, double *cum_layer_thickness_cm, int *soil_type, struct soil_properties_ *soil_properties){
  struct wetting_front *current;
  current = listFindFront(front_num, *head, NULL);

  if (current->front_num > 1){ //check if any to_bottom WFs above need to be included in the mass balance iterations. That is, if the WF above is to_bottom, include all to_bottom WFs directly above
    struct wetting_front *front = listFindFront(current->front_num-1, *head, NULL);
    if (front->to_bottom){
      current = front;
    }
  }
  while (current->to_bottom){
    if (current->front_num==1){
      break;
    }
    current = listFindFront(current->front_num-1, *head, NULL);
    if (!current->to_bottom){
      current = current->next;
      break;
    }
  }

  const double mass_at_original = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
  if (fabs(mass_at_original - prior_mass) <= MBAL_ITERATIVE_TOLERANCE) {
    return;
  }

  const double original_psi_cm = current->psi_cm;
  const double mass_at_saturation =
    lgar_apply_theta_mass_balance_correction_psi(use_dry_over_wet, 0.0, current,
                                                 head, cum_layer_thickness_cm,
                                                 soil_type, soil_properties);
  const double mass_at_dry_limit =
    lgar_apply_theta_mass_balance_correction_psi(use_dry_over_wet, PSI_UPPER_LIM,
                                                 current, head, cum_layer_thickness_cm,
                                                 soil_type, soil_properties);

  double best_psi_cm = 0.0;
  double best_abs_error_cm = fabs(mass_at_saturation - prior_mass);
  if (fabs(mass_at_dry_limit - prior_mass) < best_abs_error_cm) {
    best_psi_cm = PSI_UPPER_LIM;
    best_abs_error_cm = fabs(mass_at_dry_limit - prior_mass);
  }

  if (std::isfinite(original_psi_cm) &&
      original_psi_cm >= 0.0 && original_psi_cm <= PSI_UPPER_LIM) {
    const double mass_at_bounded_original =
      lgar_apply_theta_mass_balance_correction_psi(use_dry_over_wet, original_psi_cm,
                                                   current, head, cum_layer_thickness_cm,
                                                   soil_type, soil_properties);
    if (fabs(mass_at_bounded_original - prior_mass) < best_abs_error_cm) {
      best_psi_cm = original_psi_cm;
      best_abs_error_cm = fabs(mass_at_bounded_original - prior_mass);
    }
  }

  const double min_mass_cm = fmin(mass_at_saturation, mass_at_dry_limit);
  const double max_mass_cm = fmax(mass_at_saturation, mass_at_dry_limit);
  const bool target_bracketed =
    prior_mass >= min_mass_cm - MBAL_ITERATIVE_TOLERANCE &&
    prior_mass <= max_mass_cm + MBAL_ITERATIVE_TOLERANCE &&
    fabs(mass_at_dry_limit - mass_at_saturation) > MBAL_ITERATIVE_TOLERANCE;

  if (target_bracketed) {
    double bracket_lo_psi_cm = 0.0;
    double bracket_hi_psi_cm = PSI_UPPER_LIM;
    const bool mass_increases_with_psi = mass_at_dry_limit > mass_at_saturation;

    /*
     * The legacy loop assumed increasing psi always reduced storage. That is
     * false when a prior correction leaves a same-layer depth inversion, so
     * infer the actual response direction and solve within physical psi bounds.
     */
    for (int iter = 0; iter < MAX_ITER_MBAL_LOOP; iter++) {
      const double probe_psi_cm = 0.5 * (bracket_lo_psi_cm + bracket_hi_psi_cm);
      const double probe_mass_cm =
        lgar_apply_theta_mass_balance_correction_psi(use_dry_over_wet, probe_psi_cm,
                                                     current, head,
                                                     cum_layer_thickness_cm,
                                                     soil_type, soil_properties);
      const double probe_abs_error_cm = fabs(probe_mass_cm - prior_mass);
      if (probe_abs_error_cm < best_abs_error_cm) {
        best_psi_cm = probe_psi_cm;
        best_abs_error_cm = probe_abs_error_cm;
      }
      if (probe_abs_error_cm <= MBAL_ITERATIVE_TOLERANCE) {
        break;
      }

      if (mass_increases_with_psi) {
        if (probe_mass_cm < prior_mass) {
          bracket_lo_psi_cm = probe_psi_cm;
        }
        else {
          bracket_hi_psi_cm = probe_psi_cm;
        }
      }
      else {
        if (probe_mass_cm > prior_mass) {
          bracket_lo_psi_cm = probe_psi_cm;
        }
        else {
          bracket_hi_psi_cm = probe_psi_cm;
        }
      }
    }
  }

  lgar_apply_theta_mass_balance_correction_psi(use_dry_over_wet, best_psi_cm,
                                               current, head, cum_layer_thickness_cm,
                                               soil_type, soil_properties);
  if (verbosity.compare("high") == 0 &&
      best_abs_error_cm > MBAL_ITERATIVE_TOLERANCE) {
    printf("Bounded theta mass-balance correction left residual %.12e cm "
           "(mass_at_psi0 %.17lf, mass_at_psi_upper %.17lf, target %.17lf).\n",
           prior_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, *head),
           mass_at_saturation, mass_at_dry_limit, prior_mass);
  }
}

#endif
