#ifndef AET_CXX_INCLUDE
#define AET_CXX_INCLUDE

#include "../include/all.hxx"

static constexpr double ROOT_ZONE_TO_POPULATION_MASS_TOLERANCE_CM = 1.0e-10;

//################################################################################
/* authors : Fred Ogden and Ahmad Jan and Peter La Follette
   year    : 2022
   the code computes actual evapotranspiration given PET.
   It uses an S-shaped function used in HYDRUS-1D (Simunek & Sejna, 2018).
   AET = PET * 1/(1 + (h/h_50) )^3
   h is the capillary head at the surface and
   h_50 is the capillary head at which AET = 0.5 * PET. */
//################################################################################

static double calc_aet_stress_demand_cm(double psi_cm,
                                        double PET_timestep_cm,
                                        double time_step_h,
                                        double wilting_point_psi_cm,
                                        double field_capacity_psi_cm,
                                        int soil_num,
                                        struct soil_properties_ *soil_properties,
                                        bool apply_extreme_dry_cutoff)
{
  if (apply_extreme_dry_cutoff && psi_cm >= 1.0e6) {
    return 0.0;
  }

  const double theta_e = soil_properties[soil_num].theta_e;
  const double theta_r = soil_properties[soil_num].theta_r;
  const double vg_a = soil_properties[soil_num].vg_alpha_per_cm;
  const double vg_m = soil_properties[soil_num].vg_m;
  const double vg_n = soil_properties[soil_num].vg_n;

  const double theta_fc =
    calc_theta_from_h(field_capacity_psi_cm, vg_a, vg_m, vg_n, theta_e, theta_r);
  const double theta_wp =
    calc_theta_from_h(wilting_point_psi_cm, vg_a, vg_m, vg_n, theta_e, theta_r);
  const double theta_50 = 0.5 * (theta_fc - theta_wp) + theta_wp;

  const double Se = calc_Se_from_theta(theta_50, theta_e, theta_r);
  const double psi_50_cm = calc_h_from_Se(Se, vg_a, vg_m, vg_n);
  const double h_ratio = 1.0 + pow(psi_cm / psi_50_cm, 3.0);

  double actual_ET_demand = PET_timestep_cm * (1.0 / h_ratio) * time_step_h;
  if (actual_ET_demand < 0.0) {
    actual_ET_demand = 0.0;
  }
  else if (actual_ET_demand > (PET_timestep_cm * time_step_h)) {
    actual_ET_demand = PET_timestep_cm * time_step_h;
  }

  return actual_ET_demand;
}

static bool is_to_aet_eligible_front(const wetting_front *front)
{
  return front != NULL && front->is_WF_GW == TRUE && front->to_bottom == FALSE &&
         front->depth_cm > 0.0 && front->next != NULL && front->next->is_WF_GW == TRUE;
}

static bool is_to_aet_eligible_front(const wetting_front *front, bool allow_zero_depth)
{
  return front != NULL && front->is_WF_GW == TRUE && front->to_bottom == FALSE &&
         ((allow_zero_depth && front->depth_cm >= 0.0) || front->depth_cm > 0.0) &&
         front->next != NULL && front->next->is_WF_GW == TRUE;
}

static bool is_zero_depth_to_support_front(const wetting_front *front)
{
  return is_to_aet_eligible_front(front, true) && fabs(front->depth_cm) <= 1.0e-12;
}

static int count_no_surface_to_aet_population_fronts(const wetting_front *head)
{
  int count = 0;
  for (const wetting_front *current = head; current != NULL; current = current->next) {
    if (is_zero_depth_to_support_front(current)) {
      count++;
    }
  }
  return count;
}

static wetting_front *find_active_zero_depth_to_support_front(wetting_front *head)
{
  wetting_front *active_front = NULL;
  for (wetting_front *current = head; current != NULL; current = current->next) {
    if (is_zero_depth_to_support_front(current)) {
      active_front = current;
      continue;
    }

    if (current->depth_cm > 1.0e-12) {
      break;
    }
  }

  return active_front;
}

static double find_zero_depth_support_interval_bottom_cm(const wetting_front *active_front,
                                                         double root_zone_depth_cm)
{
  if (active_front == NULL) {
    return 0.0;
  }

  for (const wetting_front *current = active_front->next; current != NULL; current = current->next) {
    if (current->depth_cm <= 1.0e-12) {
      continue;
    }

    if (current->is_WF_GW == TRUE && current->to_bottom == TRUE) {
      const wetting_front *connected_chain_front = current;
      while (connected_chain_front != NULL && connected_chain_front->is_WF_GW == TRUE &&
             connected_chain_front->to_bottom == TRUE) {
        if (connected_chain_front->next == NULL) {
          return root_zone_depth_cm;
        }

        if (connected_chain_front->next->depth_cm >= root_zone_depth_cm - 1.0e-12) {
          return root_zone_depth_cm;
        }

        connected_chain_front = connected_chain_front->next;
      }

      if (connected_chain_front == NULL) {
        return root_zone_depth_cm;
      }

      return fmin(connected_chain_front->depth_cm, root_zone_depth_cm);
    }

    return fmin(current->depth_cm, root_zone_depth_cm);
  }

  return root_zone_depth_cm;
}

static bool insert_zero_depth_to_population_front(double root_zone_depth_cm,
                                                  int *soil_type,
                                                  double *frozen_factor,
                                                  soil_properties_ *soil_properties,
                                                  wetting_front **head)
{
  if (head == NULL || *head == NULL) {
    return false;
  }

  wetting_front *current = *head;
  if (current->psi_cm >= 1.0e6) {
    return false;
  }

  const double psi_increment_cm = 0.25 * root_zone_depth_cm;
  const double new_psi_cm = current->psi_cm + psi_increment_cm;
  const int layer_num = 1;
  const int soil_num = soil_type[layer_num];
  const double theta_e = soil_properties[soil_num].theta_e;
  const double theta_r = soil_properties[soil_num].theta_r;
  const double vg_a = soil_properties[soil_num].vg_alpha_per_cm;
  const double vg_m = soil_properties[soil_num].vg_m;
  const double vg_n = soil_properties[soil_num].vg_n;
  const double theta_new = calc_theta_from_h(new_psi_cm, vg_a, vg_m, vg_n, theta_e, theta_r);

  listInsertFirst(0.0, theta_new, 1, layer_num, FALSE, head);
  wetting_front *inserted = *head;
  inserted->is_WF_GW = TRUE;
  inserted->to_bottom = FALSE;
  inserted->psi_cm = new_psi_cm;

  const double Se = calc_Se_from_theta(theta_new, theta_e, theta_r);
  const double Ksat_cm_per_h = soil_properties[soil_num].Ksat_cm_per_h * frozen_factor[layer_num];
  inserted->K_cm_per_h = calc_K_from_Se(Se, Ksat_cm_per_h, vg_m);
  inserted->dzdt_cm_per_h = 0.0;
  return true;
}

struct RootZoneTOPopulationCandidate {
  wetting_front *lower_front;
  double top_depth_cm;
  double bottom_depth_cm;
  double psi_new_cm;
  double upper_theta;
};

static bool is_mobile_to_front_in_root_zone(const wetting_front *front,
                                            double top_of_to_root_zone_cm,
                                            double root_zone_depth_cm)
{
  return is_to_aet_eligible_front(front) &&
         front->depth_cm > top_of_to_root_zone_cm + 1.0e-8 &&
         front->depth_cm < root_zone_depth_cm - 1.0e-8;
}

static int count_mobile_to_fronts_in_root_zone(const wetting_front *head,
                                               double top_of_to_root_zone_cm,
                                               double root_zone_depth_cm)
{
  int count = 0;
  for (const wetting_front *current = head; current != NULL; current = current->next) {
    if (is_mobile_to_front_in_root_zone(current, top_of_to_root_zone_cm, root_zone_depth_cm)) {
      count++;
    }
  }
  return count;
}

static bool has_sparse_no_surface_to_population(const wetting_front *head,
                                                double root_zone_depth_cm)
{
  return count_mobile_to_fronts_in_root_zone(head, 0.0, root_zone_depth_cm) < 2;
}

static bool find_root_zone_to_population_candidate(const wetting_front *head,
                                                   double top_of_to_root_zone_cm,
                                                   double root_zone_depth_cm,
                                                   double field_capacity_psi_cm,
                                                   RootZoneTOPopulationCandidate *candidate)
{
  if (head == NULL || candidate == NULL || top_of_to_root_zone_cm >= root_zone_depth_cm) {
    return false;
  }

  const double seed_increment_cm = fmax(1.0, fmin(field_capacity_psi_cm, 0.25 * root_zone_depth_cm));

  const wetting_front *first_relevant_gw = NULL;
  for (const wetting_front *current = head; current != NULL; current = current->next) {
    if (!current->is_WF_GW) {
      continue;
    }
    if (current->depth_cm > top_of_to_root_zone_cm + 1.0e-8) {
      first_relevant_gw = current;
      break;
    }
  }

  if (first_relevant_gw != NULL && first_relevant_gw->to_bottom &&
      first_relevant_gw->depth_cm > top_of_to_root_zone_cm + 1.0e-8) {
    candidate->lower_front = const_cast<wetting_front *>(first_relevant_gw);
    candidate->top_depth_cm = top_of_to_root_zone_cm;
    candidate->bottom_depth_cm = fmin(first_relevant_gw->depth_cm, root_zone_depth_cm);
    candidate->psi_new_cm = first_relevant_gw->psi_cm + seed_increment_cm;
    candidate->upper_theta = first_relevant_gw->theta;
    return candidate->bottom_depth_cm > candidate->top_depth_cm + 1.0e-8;
  }

  double interval_top_cm = top_of_to_root_zone_cm;
  const wetting_front *previous_relevant = NULL;

  bool found = false;
  double best_interval_thickness_cm = -1.0;
  double best_interval_top_cm = 0.0;
  double best_interval_bottom_cm = 0.0;
  double best_psi_new_cm = 0.0;
  double best_upper_theta = 0.0;
  wetting_front *best_lower_front = NULL;

  for (const wetting_front *current = head; current != NULL; current = current->next) {
    if (!current->is_WF_GW) {
      continue;
    }

    if (current->depth_cm <= top_of_to_root_zone_cm + 1.0e-8) {
      if (top_of_to_root_zone_cm <= 1.0e-8) {
        previous_relevant = current;
      }
      continue;
    }

    const double interval_bottom_cm = fmin(current->depth_cm, root_zone_depth_cm);
    const double interval_thickness_cm = interval_bottom_cm - interval_top_cm;

    if (previous_relevant != NULL && previous_relevant->to_bottom) {
      if (current->depth_cm >= root_zone_depth_cm - 1.0e-8) {
        break;
      }
      previous_relevant = current;
      interval_top_cm = current->depth_cm;
      continue;
    }

    if (interval_thickness_cm > best_interval_thickness_cm + 1.0e-8) {
      double psi_new_cm = current->psi_cm + seed_increment_cm;
      if (previous_relevant != NULL && previous_relevant->psi_cm > current->psi_cm + 1.0e-3) {
        psi_new_cm = 0.5 * (previous_relevant->psi_cm + current->psi_cm);
      }

      best_interval_thickness_cm = interval_thickness_cm;
      best_interval_top_cm = interval_top_cm;
      best_interval_bottom_cm = interval_bottom_cm;
      best_psi_new_cm = psi_new_cm;
      best_upper_theta = (previous_relevant != NULL) ? previous_relevant->theta : current->theta;
      best_lower_front = const_cast<wetting_front *>(current);
      found = true;
    }

    if (current->depth_cm >= root_zone_depth_cm - 1.0e-8) {
      break;
    }

    previous_relevant = current;
    interval_top_cm = current->depth_cm;
  }

  if (!found || best_lower_front == NULL || best_interval_bottom_cm <= best_interval_top_cm + 1.0e-8) {
    return false;
  }

  candidate->lower_front = best_lower_front;
  candidate->top_depth_cm = best_interval_top_cm;
  candidate->bottom_depth_cm = best_interval_bottom_cm;
  candidate->psi_new_cm = best_psi_new_cm;
  candidate->upper_theta = best_upper_theta;
  return true;
}

static int layer_num_from_depth_cm(int num_layers,
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

static wetting_front *insert_root_zone_to_front(double depth_cm,
                                                double psi_new_cm,
                                                wetting_front *lower_front,
                                                int num_layers,
                                                int *soil_type,
                                                double *cum_layer_thickness_cm,
                                                double *frozen_factor,
                                                soil_properties_ *soil_properties,
                                                wetting_front **head)
{
  if (lower_front == NULL || head == NULL || *head == NULL) {
    return NULL;
  }

  const int layer_num = layer_num_from_depth_cm(num_layers, cum_layer_thickness_cm, depth_cm);
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

  wetting_front *inserted = NULL;
  if (lower_front == *head) {
    listInsertFirst(depth_cm, theta_new, 1, layer_num, FALSE, head);
    inserted = *head;
  }
  else {
    wetting_front *previous = *head;
    while (previous != NULL && previous->next != lower_front) {
      previous = previous->next;
    }

    if (previous == NULL) {
      return NULL;
    }

    inserted = (wetting_front *) malloc(sizeof(wetting_front));
    inserted->depth_cm = depth_cm;
    inserted->theta = theta_new;
    inserted->front_num = lower_front->front_num;
    inserted->layer_num = layer_num;
    inserted->to_bottom = FALSE;
    inserted->dzdt_cm_per_h = 0.0;
    inserted->is_WF_GW = TRUE;
    inserted->next = lower_front;
    previous->next = inserted;

    for (wetting_front *current = lower_front; current != NULL; current = current->next) {
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

static wetting_front *find_depth_repair_target_below_connected_to_chain(wetting_front *lower_front)
{
  if (lower_front == NULL) {
    return NULL;
  }

  if (!lower_front->to_bottom) {
    return lower_front;
  }

  wetting_front *current = lower_front;
  while (current != NULL && current->to_bottom) {
    if (current->next == NULL || current->next->is_WF_GW != lower_front->is_WF_GW) {
      return NULL;
    }
    if (!current->next->to_bottom) {
      return current->next;
    }
    current = current->next;
  }

  return NULL;
}

static bool restore_root_zone_to_population_mass_locally(wetting_front *inserted,
                                                         wetting_front *lower_front,
                                                         double upper_theta)
{
  if (inserted == NULL || lower_front == NULL || lower_front->to_bottom || lower_front->next == NULL) {
    return false;
  }

  const double next_depth_cm = lower_front->next->depth_cm;
  if (next_depth_cm <= inserted->depth_cm + 1.0e-8) {
    return false;
  }

  const double theta_new = inserted->theta;
  const double theta_lower = lower_front->theta;
  const double denominator = theta_new - theta_lower;
  if (fabs(denominator) <= 1.0e-12) {
    return false;
  }

  const double z_inserted = inserted->depth_cm;
  const double z_lower = lower_front->depth_cm;
  const double solved_lower_depth_cm =
    (((upper_theta - theta_lower) * z_lower) + ((theta_new - upper_theta) * z_inserted)) /
    denominator;

  if (solved_lower_depth_cm <= z_inserted + 1.0e-8 ||
      solved_lower_depth_cm >= next_depth_cm - 1.0e-8) {
    return false;
  }

  lower_front->depth_cm = solved_lower_depth_cm;
  return true;
}

static double restore_root_zone_to_population_mass_via_depth(wetting_front *lower_front,
                                                             double target_mass,
                                                             int num_layers,
                                                             wetting_front *head,
                                                             double *cum_layer_thickness_cm)
{
  wetting_front *depth_target = find_depth_repair_target_below_connected_to_chain(lower_front);
  if (depth_target == NULL) {
    return target_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, head);
  }

  const double original_depth_cm = depth_target->depth_cm;
  double upper_depth_cm = cum_layer_thickness_cm[num_layers];
  if (depth_target->next != NULL) {
    upper_depth_cm = depth_target->next->depth_cm;
  }

  if (upper_depth_cm <= original_depth_cm + 1.0e-12) {
    return target_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, head);
  }

  const double mass_at_lo = lgar_calc_mass_bal(cum_layer_thickness_cm, head);
  depth_target->depth_cm = upper_depth_cm;
  const double mass_at_hi = lgar_calc_mass_bal(cum_layer_thickness_cm, head);
  depth_target->depth_cm = original_depth_cm;

  if (mass_at_hi <= mass_at_lo + ROOT_ZONE_TO_POPULATION_MASS_TOLERANCE_CM) {
    return target_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, head);
  }

  if (mass_at_hi < target_mass) {
    depth_target->depth_cm = upper_depth_cm;
    return target_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, head);
  }

  double bracket_lo = original_depth_cm;
  double bracket_hi = upper_depth_cm;
  double best_depth_cm = original_depth_cm;
  double best_mass_error_cm = fabs(target_mass - mass_at_lo);

  for (int iter = 0; iter < 80; iter++) {
    const double probe_depth_cm = 0.5 * (bracket_lo + bracket_hi);
    depth_target->depth_cm = probe_depth_cm;
    const double current_mass = lgar_calc_mass_bal(cum_layer_thickness_cm, head);
    const double remaining_mass_error_cm = target_mass - current_mass;

    if (fabs(remaining_mass_error_cm) < best_mass_error_cm) {
      best_mass_error_cm = fabs(remaining_mass_error_cm);
      best_depth_cm = probe_depth_cm;
    }

    if (fabs(remaining_mass_error_cm) <= ROOT_ZONE_TO_POPULATION_MASS_TOLERANCE_CM) {
      best_depth_cm = probe_depth_cm;
      break;
    }

    if (remaining_mass_error_cm > 0.0) {
      bracket_lo = probe_depth_cm;
    }
    else {
      bracket_hi = probe_depth_cm;
    }
  }

  depth_target->depth_cm = best_depth_cm;
  return target_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, head);
}

static bool lgarto_try_insert_root_zone_to_population_front(int num_layers,
                                                            double top_of_to_root_zone_cm,
                                                            double root_zone_depth_cm,
                                                            double field_capacity_psi_cm,
                                                            int *soil_type,
                                                            double *cum_layer_thickness_cm,
                                                            double *frozen_factor,
                                                            soil_properties_ *soil_properties,
                                                            wetting_front **head)
{
  RootZoneTOPopulationCandidate candidate = {};
  if (!find_root_zone_to_population_candidate(*head, top_of_to_root_zone_cm, root_zone_depth_cm,
                                              field_capacity_psi_cm, &candidate)) {
    return false;
  }

  const double prior_mass = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
  wetting_front *snapshot = listCopy(*head, NULL);
  const double insert_depth_cm = 0.5 * (candidate.top_depth_cm + candidate.bottom_depth_cm);

  wetting_front *inserted =
    insert_root_zone_to_front(insert_depth_cm, candidate.psi_new_cm, candidate.lower_front, num_layers,
                              soil_type, cum_layer_thickness_cm, frozen_factor, soil_properties, head);
  if (inserted == NULL || inserted->next == NULL || inserted->next->is_WF_GW == FALSE) {
    if (*head != NULL) {
      listDelete(*head);
    }
    *head = snapshot;
    return false;
  }

  wetting_front *lower_front = inserted->next;
  if (verbosity.compare("high") == 0) {
    printf("root-zone TO population candidate: insert depth %.6f cm, psi %.6f cm, "
           "upper_theta %.6f, lower front %d at depth %.6f cm (to_bottom=%d)\n",
           insert_depth_cm, candidate.psi_new_cm, candidate.upper_theta, lower_front->front_num,
           lower_front->depth_cm, lower_front->to_bottom ? 1 : 0);
  }
  double residual_mass_error_cm = prior_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
  if (residual_mass_error_cm > ROOT_ZONE_TO_POPULATION_MASS_TOLERANCE_CM &&
      restore_root_zone_to_population_mass_locally(inserted, lower_front, candidate.upper_theta)) {
    residual_mass_error_cm = prior_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
    if (verbosity.compare("high") == 0) {
      printf("root-zone TO population local mass repair residual = %.17lf cm\n",
             residual_mass_error_cm);
    }
  }

  if (residual_mass_error_cm > ROOT_ZONE_TO_POPULATION_MASS_TOLERANCE_CM) {
    if (verbosity.compare("high") == 0) {
      printf("root-zone TO population falling back to iterative theta repair with residual %.17lf cm\n",
             residual_mass_error_cm);
    }
    lgar_theta_mass_balance_correction(false, lower_front->front_num, prior_mass, head,
                                       cum_layer_thickness_cm, soil_type, soil_properties);
    residual_mass_error_cm = prior_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
  }

  if (residual_mass_error_cm > ROOT_ZONE_TO_POPULATION_MASS_TOLERANCE_CM) {
    if (verbosity.compare("high") == 0) {
      printf("root-zone TO population falling back to depth repair with residual %.17lf cm\n",
             residual_mass_error_cm);
    }
    residual_mass_error_cm =
      restore_root_zone_to_population_mass_via_depth(lower_front, prior_mass, num_layers, *head,
                                                     cum_layer_thickness_cm);
  }

  if (fabs(prior_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, *head)) >
      100.0 * ROOT_ZONE_TO_POPULATION_MASS_TOLERANCE_CM) {
    listDelete(*head);
    *head = snapshot;
    return false;
  }

  if (verbosity.compare("high") == 0) {
    printf("root-zone TO population inserted front %d at depth %.6f cm with psi %.6f cm; "
           "current mobile root-zone TO count = %d\n",
           inserted->front_num, inserted->depth_cm, inserted->psi_cm,
           count_mobile_to_fronts_in_root_zone(*head, top_of_to_root_zone_cm, root_zone_depth_cm));
    listPrint(*head);
  }

  listDelete(snapshot);
  return true;
}

static void lgarto_ensure_root_zone_to_population(int num_layers,
                                                  double deepest_surface_depth_cm,
                                                  double root_zone_depth_cm,
                                                  double PET_timestep_cm,
                                                  double field_capacity_psi_cm,
                                                  int *soil_type,
                                                  double *cum_layer_thickness_cm,
                                                  double *frozen_factor,
                                                  soil_properties_ *soil_properties,
                                                  wetting_front **head,
                                                  bool allow_sparse_existing_mobile_fronts)
{
  (void) num_layers;
  (void) deepest_surface_depth_cm;

  if (head == NULL || *head == NULL || root_zone_depth_cm <= 0.0 || PET_timestep_cm <= 0.0) {
    return;
  }

  if (listLength_surface(*head) != 0) {
    return;
  }

  const int mobile_to_front_count =
    count_mobile_to_fronts_in_root_zone(*head, 0.0, root_zone_depth_cm);

  if (!allow_sparse_existing_mobile_fronts && mobile_to_front_count > 0) {
    return;
  }

  if (allow_sparse_existing_mobile_fronts && mobile_to_front_count >= 2) {
    return;
  }

  const int min_mobile_to_fronts_in_root_zone = 3;
  for (int population_iter = 0; population_iter < min_mobile_to_fronts_in_root_zone; population_iter++) {
    if (count_no_surface_to_aet_population_fronts(*head) >= min_mobile_to_fronts_in_root_zone) {
      break;
    }

    if (!insert_zero_depth_to_population_front(root_zone_depth_cm, soil_type, frozen_factor,
                                               soil_properties, head)) {
      break;
    }

    if (verbosity.compare("high") == 0) {
      printf("states after zero-depth TO root-zone population insertion:\n");
      listPrint(*head);
    }
  }
}


extern double calc_aet(double PET_timestep_cm, double time_step_h, double wilting_point_psi_cm, double field_capacity_psi_cm,
		       int *soil_type, double AET_thresh_Theta, double AET_expon,
		       struct wetting_front* head, struct soil_properties_ *soil_properties)
{

  if (verbosity.compare("high") == 0) {
    printf("Computing AET... \n");
    printf("Note: AET_thresh_theta = %lf and AET_expon = %lf are not used in the computation of the current AET model. \n", AET_thresh_Theta, AET_expon);
  }
  
  double actual_ET_demand = 0.0;
  struct wetting_front *current;
  current = head;

  if (current != NULL) {
    const int layer_num = current->layer_num;
    const int soil_num = soil_type[layer_num];
    actual_ET_demand = calc_aet_stress_demand_cm(current->psi_cm, PET_timestep_cm, time_step_h,
                                                 wilting_point_psi_cm, field_capacity_psi_cm,
                                                 soil_num, soil_properties, true);
  }

  if (verbosity.compare("high") == 0) {
    printf("AET =  %14.10f \n",actual_ET_demand);
  }
  
  return actual_ET_demand;

}

extern double calc_aet(bool TO_enabled, double PET_timestep_cm, double time_step_h, double wilting_point_psi_cm,
		       double field_capacity_psi_cm, double root_zone_depth_cm, double *surf_frac_rz,
		       int *soil_type, double AET_thresh_Theta, double AET_expon, struct wetting_front* head,
		       struct soil_properties_ *soil_properties, double *surf_AET_vec)
{
  if (verbosity.compare("high") == 0) {
    printf("Computing AET... \n");
    printf("Note: AET_thresh_theta = %lf and AET_expon = %lf are not used in the computation of the current AET model. \n",
           AET_thresh_Theta, AET_expon);
  }

  const int front_count = listLength(head);
  for (int front_num = 0; front_num <= front_count; front_num++) {
    surf_AET_vec[front_num] = 0.0;
  }
  *surf_frac_rz = 0.0;

  if (!TO_enabled) {
    const double actual_ET_demand =
      calc_aet(PET_timestep_cm, time_step_h, wilting_point_psi_cm, field_capacity_psi_cm,
               soil_type, AET_thresh_Theta, AET_expon, head, soil_properties);
    if (head != NULL) {
      const int wetting_front_free_drainage_temp = wetting_front_free_drainage(head);
      surf_AET_vec[wetting_front_free_drainage_temp] = actual_ET_demand;
    }
    return actual_ET_demand;
  }

  if (head == NULL || root_zone_depth_cm <= 0.0 || PET_timestep_cm <= 0.0 || listLength_surface(head) == 0) {
    return 0.0;
  }

  std::vector<double> surface_thickness_by_front_cm(front_count + 1, 0.0);
  wetting_front *current = head;
  while (current != NULL && current->is_WF_GW) {
    current = current->next;
  }

  double prior_depth_cm = 0.0;
  double total_surface_thickness_in_root_zone_cm = 0.0;
  double pending_surface_thickness_cm = 0.0;
  wetting_front *deepest_surface_front = nullptr;
  while (current != NULL && !current->is_WF_GW) {
    deepest_surface_front = current;

    if (total_surface_thickness_in_root_zone_cm < root_zone_depth_cm) {
      const double interval_top_cm = fmin(prior_depth_cm, root_zone_depth_cm);
      const double interval_bottom_cm = fmin(current->depth_cm, root_zone_depth_cm);
      const double thickness_cm = fmax(0.0, interval_bottom_cm - interval_top_cm);
      pending_surface_thickness_cm += thickness_cm;
      total_surface_thickness_in_root_zone_cm += thickness_cm;
    }

    if (current->to_bottom == FALSE) {
      surface_thickness_by_front_cm[current->front_num] += pending_surface_thickness_cm;
      pending_surface_thickness_cm = 0.0;

      if (total_surface_thickness_in_root_zone_cm >= root_zone_depth_cm ||
          current->depth_cm >= root_zone_depth_cm) {
        break;
      }
    }

    prior_depth_cm = current->depth_cm;

    if (total_surface_thickness_in_root_zone_cm >= root_zone_depth_cm &&
        pending_surface_thickness_cm <= 0.0) {
      break;
    }

    current = current->next;
  }

  if (pending_surface_thickness_cm > 0.0 && deepest_surface_front != NULL) {
    // Preserve the old LGARTO behavior for a root-zone segment that ends inside a
    // chain of to_bottom surface fronts: keep carrying that thickness to the
    // deepest available surface front rather than allocating it to the interface
    // fronts themselves.
    surface_thickness_by_front_cm[deepest_surface_front->front_num] += pending_surface_thickness_cm;
  }

  *surf_frac_rz = fmin(1.0, total_surface_thickness_in_root_zone_cm / root_zone_depth_cm);

  const int soil_num = soil_type[1];
  double cumulative_surface_aet_cm = 0.0;
  for (current = head; current != NULL; current = current->next) {
    if (current->is_WF_GW) {
      continue;
    }

    const double thickness_cm = surface_thickness_by_front_cm[current->front_num];
    if (thickness_cm <= 0.0) {
      continue;
    }

    const double frac_of_root_zone = fmin(1.0, thickness_cm / root_zone_depth_cm);
    double actual_ET_demand =
      calc_aet_stress_demand_cm(current->psi_cm, PET_timestep_cm, time_step_h, wilting_point_psi_cm,
                                field_capacity_psi_cm, soil_num, soil_properties, false);
    actual_ET_demand *= frac_of_root_zone;
    actual_ET_demand = fmax(0.0, fmin(actual_ET_demand, PET_timestep_cm * time_step_h * frac_of_root_zone));

    surf_AET_vec[current->front_num] = actual_ET_demand;
    cumulative_surface_aet_cm += actual_ET_demand;

    if (verbosity.compare("high") == 0) {
      printf("\n ");
      printf("AET =  %14.10f \n", actual_ET_demand);
      printf("PET_timestep_cm*time_step_h =  %14.10f \n", PET_timestep_cm * time_step_h);
      printf("front_num = %d, frac_of_root_zone = %lf \n", current->front_num, frac_of_root_zone);
      printf("\n ");
    }
  }

  return cumulative_surface_aet_cm;
}

static double calc_aet_for_individual_TO_WFs(int WF_num,
                                             double WF_thickness_cm,
                                             double rooting_zone_depth_cm,
                                             double PET_timestep_cm,
                                             double time_step_h,
                                             double wilting_point_psi_cm,
                                             double field_capacity_psi_cm,
                                             int *soil_type,
                                             struct soil_properties_ *soil_properties,
                                             struct wetting_front **head)
{
  wetting_front *current = listFindFront(WF_num, *head, NULL);
  if (current == NULL || rooting_zone_depth_cm <= 0.0) {
    return 0.0;
  }

  const double frac = fmin(WF_thickness_cm / rooting_zone_depth_cm, 1.0);
  if (verbosity.compare("high") == 0) {
    printf("frac, used in calc_aet_for_individual_TO_WFs: %lf \n", frac);
  }

  if (current->is_WF_GW == FALSE || frac <= 0.0) {
    return 0.0;
  }

  // Intentionally preserve the old LGARTO behavior: TO-root-zone stress uses layer-1
  // soil properties even for deeper TO fronts.
  const int soil_num = soil_type[1];
  double actual_ET_demand =
    calc_aet_stress_demand_cm(current->psi_cm, PET_timestep_cm, time_step_h,
                              wilting_point_psi_cm, field_capacity_psi_cm,
                              soil_num, soil_properties, false);
  actual_ET_demand *= frac;

  if (verbosity.compare("high") == 0) {
    printf("current->front_num: %d \n", current->front_num);
    printf("AET component =  %14.17f \n", actual_ET_demand);
  }

  return actual_ET_demand;
}

static double cap_to_aet_extraction_to_local_interval_cm(const wetting_front *front,
                                                         double requested_aet_cm)
{
  if (front == NULL || front->next == NULL) {
    return 0.0;
  }

  const double delta_theta = front->next->theta - front->theta;
  if (delta_theta <= 1.0e-12) {
    return 0.0;
  }

  const double available_depth_cm = front->next->depth_cm - front->depth_cm - 1.0e-8;
  if (available_depth_cm <= 0.0) {
    return 0.0;
  }

  const double max_extractable_aet_cm = available_depth_cm * delta_theta;
  return fmax(0.0, fmin(requested_aet_cm, max_extractable_aet_cm));
}

static wetting_front *find_surface_supported_to_reference_front(wetting_front *head)
{
  wetting_front *first_groundwater_front = NULL;
  bool seen_surface_front = false;

  for (wetting_front *current = head; current != NULL; current = current->next) {
    if (current->is_WF_GW == TRUE) {
      if (first_groundwater_front == NULL) {
        first_groundwater_front = current;
      }
      if (seen_surface_front) {
        return current;
      }
    }
    else {
      seen_surface_front = true;
    }
  }

  return first_groundwater_front;
}

static wetting_front *find_surface_supported_to_aet_target_front(wetting_front *head)
{
  for (wetting_front *current = find_surface_supported_to_reference_front(head);
       current != NULL;
       current = current->next) {
    if (current->is_WF_GW == TRUE && current->to_bottom == FALSE &&
        current->next != NULL && current->next->is_WF_GW == TRUE) {
      return current;
    }
  }

  return NULL;
}

static wetting_front *find_deepest_surface_front(wetting_front *head)
{
  wetting_front *deepest_surface_front = NULL;
  for (wetting_front *current = head; current != NULL; current = current->next) {
    if (current->is_WF_GW == FALSE) {
      deepest_surface_front = current;
    }
  }
  return deepest_surface_front;
}

static void refresh_front_state_from_theta(wetting_front *front,
                                           int *soil_type,
                                           double *frozen_factor,
                                           soil_properties_ *soil_properties)
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
  const double Se = calc_Se_from_theta(front->theta, theta_e, theta_r);
  const double Ksat_cm_per_h = soil_properties[soil_num].Ksat_cm_per_h * frozen_factor[front->layer_num];

  front->psi_cm = calc_h_from_Se(Se, vg_a, vg_m, vg_n);
  front->K_cm_per_h = calc_K_from_Se(Se, Ksat_cm_per_h, vg_m);
}

static void refresh_front_state_from_psi(wetting_front *front,
                                         int *soil_type,
                                         double *frozen_factor,
                                         soil_properties_ *soil_properties)
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
  const double Ksat_cm_per_h = soil_properties[soil_num].Ksat_cm_per_h * frozen_factor[front->layer_num];

  front->theta = calc_theta_from_h(front->psi_cm, vg_a, vg_m, vg_n, theta_e, theta_r);
  const double Se = calc_Se_from_theta(front->theta, theta_e, theta_r);
  front->K_cm_per_h = calc_K_from_Se(Se, Ksat_cm_per_h, vg_m);
}

static void propagate_surface_chain_psi_from_below(int *soil_type,
                                                   double *frozen_factor,
                                                   soil_properties_ *soil_properties,
                                                   wetting_front **head)
{
  if (head == NULL || *head == NULL) {
    return;
  }

  const int wetting_front_count = listLength(*head);
  for (int wf = wetting_front_count - 1; wf >= 1; wf--) {
    wetting_front *current = listFindFront(wf, *head, NULL);
    if (current == NULL || current->next == NULL) {
      continue;
    }

    if (current->is_WF_GW == FALSE && current->next->is_WF_GW == FALSE && current->to_bottom == TRUE) {
      current->psi_cm = current->next->psi_cm;
      refresh_front_state_from_psi(current, soil_type, frozen_factor, soil_properties);
    }
  }
}

static void propagate_groundwater_chain_psi_from_below(int *soil_type,
                                                       double *frozen_factor,
                                                       soil_properties_ *soil_properties,
                                                       wetting_front **head)
{
  if (head == NULL || *head == NULL) {
    return;
  }

  const int wetting_front_count = listLength(*head);
  for (int wf = wetting_front_count - 1; wf >= 1; wf--) {
    wetting_front *current = listFindFront(wf, *head, NULL);
    if (current == NULL || current->next == NULL) {
      continue;
    }

    if (current->is_WF_GW == TRUE && current->next->is_WF_GW == TRUE && current->to_bottom == TRUE) {
      current->psi_cm = current->next->psi_cm;
      refresh_front_state_from_psi(current, soil_type, frozen_factor, soil_properties);
    }
  }
}

static double lgarto_extract_missing_to_aet_from_surface_WFs(double missing_to_aet_cm,
                                                             double *cum_layer_thickness_cm,
                                                             int *soil_type,
                                                             double *frozen_factor,
                                                             soil_properties_ *soil_properties,
                                                             wetting_front **head)
{
  if (head == NULL || *head == NULL || listLength_surface(*head) == 0 || missing_to_aet_cm <= 0.0) {
    return 0.0;
  }

  wetting_front *current = find_deepest_surface_front(*head);
  if (current == NULL) {
    return 0.0;
  }

  const double mass_before = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);

  if (current->layer_num == 1) {
    if (current->depth_cm > 0.0) {
      const int soil_num = soil_type[current->layer_num];
      const double theta_e = soil_properties[soil_num].theta_e;
      const double theta_r = soil_properties[soil_num].theta_r;
      const double theta_reduction = missing_to_aet_cm / current->depth_cm;
      const double theta_new =
        fmax(theta_r + 1.0e-12, fmin(current->theta - theta_reduction, theta_e));
      current->theta = theta_new;
      refresh_front_state_from_theta(current, soil_type, frozen_factor, soil_properties);
    }
  }
  else if (current->next != NULL) {
    const int layer_num = current->layer_num;
    const int soil_num = soil_type[layer_num];
    double *delta_thetas = (double *)malloc(sizeof(double) * (layer_num + 1));
    double *delta_thickness = (double *)malloc(sizeof(double) * (layer_num + 1));

    const double psi_cm = current->psi_cm;
    const double psi_cm_below = current->next->psi_cm;
    double new_mass =
      (current->depth_cm - cum_layer_thickness_cm[layer_num - 1]) *
      (current->theta - current->next->theta);
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

    delta_thetas[layer_num] = current->next->theta;
    delta_thickness[layer_num] = current->depth_cm - cum_layer_thickness_cm[layer_num - 1];
    prior_mass -= missing_to_aet_cm;

    double aet_demand_cm = missing_to_aet_cm;
    const double theta_new =
      lgar_theta_mass_balance(layer_num, soil_num, current->psi_cm, new_mass, prior_mass, 0.0,
                              &aet_demand_cm, delta_thetas, delta_thickness, soil_type,
                              soil_properties);
    current->theta =
      fmax(soil_properties[soil_num].theta_r, fmin(theta_new, soil_properties[soil_num].theta_e));
    refresh_front_state_from_theta(current, soil_type, frozen_factor, soil_properties);

    free(delta_thetas);
    free(delta_thickness);
  }

  propagate_surface_chain_psi_from_below(soil_type, frozen_factor, soil_properties, head);
  return fmax(0.0, mass_before - lgar_calc_mass_bal(cum_layer_thickness_cm, *head));
}

static double lgarto_extract_missing_to_aet_from_to_chain(wetting_front *target_front,
                                                          double missing_to_aet_cm,
                                                          double *cum_layer_thickness_cm,
                                                          int *soil_type,
                                                          double *frozen_factor,
                                                          soil_properties_ *soil_properties,
                                                          wetting_front **head)
{
  if (target_front == NULL || head == NULL || *head == NULL || missing_to_aet_cm <= 0.0) {
    return 0.0;
  }

  const double mass_before = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
  const double target_mass = mass_before - missing_to_aet_cm;
  if (target_mass >= mass_before) {
    return 0.0;
  }

  const double original_psi_cm = target_front->psi_cm;
  double psi_lo_cm = fmax(0.0, original_psi_cm);
  double mass_lo = mass_before;

  auto apply_trial_psi = [&](double psi_trial_cm) {
    target_front->psi_cm = psi_trial_cm;
    refresh_front_state_from_psi(target_front, soil_type, frozen_factor, soil_properties);
    propagate_groundwater_chain_psi_from_below(soil_type, frozen_factor, soil_properties, head);
    return lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
  };

  double psi_hi_cm = fmax(psi_lo_cm + 1.0, psi_lo_cm * 2.0);
  double mass_hi = apply_trial_psi(psi_hi_cm);
  int expand_iter = 0;
  while (mass_hi > target_mass && psi_hi_cm < 1.0e7 && expand_iter < 80) {
    psi_hi_cm = fmax(psi_hi_cm + 1.0, psi_hi_cm * 2.0);
    mass_hi = apply_trial_psi(psi_hi_cm);
    expand_iter++;
  }

  double best_psi_cm = psi_lo_cm;
  double best_mass_error_cm = fabs(mass_lo - target_mass);
  if (fabs(mass_hi - target_mass) < best_mass_error_cm) {
    best_mass_error_cm = fabs(mass_hi - target_mass);
    best_psi_cm = psi_hi_cm;
  }

  if (mass_hi <= target_mass) {
    for (int iter = 0; iter < 80; iter++) {
      const double psi_mid_cm = 0.5 * (psi_lo_cm + psi_hi_cm);
      const double mass_mid = apply_trial_psi(psi_mid_cm);
      const double mass_error_cm = fabs(mass_mid - target_mass);

      if (mass_error_cm < best_mass_error_cm) {
        best_mass_error_cm = mass_error_cm;
        best_psi_cm = psi_mid_cm;
      }

      if (mass_error_cm <= 1.0e-10) {
        best_psi_cm = psi_mid_cm;
        break;
      }

      if (mass_mid > target_mass) {
        psi_lo_cm = psi_mid_cm;
      }
      else {
        psi_hi_cm = psi_mid_cm;
      }
    }
  }

  apply_trial_psi(best_psi_cm);
  const double extracted_cm = fmax(0.0, mass_before - lgar_calc_mass_bal(cum_layer_thickness_cm, *head));

  if (verbosity.compare("high") == 0 && extracted_cm > 0.0) {
    printf("TO AET direct fallback extracted %.17lf cm from TO front %d via psi increase to %.17lf cm.\n",
           extracted_cm, target_front->front_num, target_front->psi_cm);
  }

  if (extracted_cm <= 0.0) {
    target_front->psi_cm = original_psi_cm;
    refresh_front_state_from_psi(target_front, soil_type, frozen_factor, soil_properties);
    propagate_groundwater_chain_psi_from_below(soil_type, frozen_factor, soil_properties, head);
  }

  return extracted_cm;
}

extern double lgarto_calc_aet_from_TO_WFs(int num_layers,
					  double deepest_surf_depth_at_start,
					  double root_zone_depth_cm,
					  double PET_timestep_cm,
					  double timestep_h,
					  double surf_frac_rz,
					  bool allow_root_zone_to_population,
					  double wilting_point_psi_cm,
					  double field_capacity_psi_cm,
					  int *soil_type,
					  double *cum_layer_thickness_cm,
					  double *frozen_factor,
					  struct soil_properties_ *soil_properties,
					  struct wetting_front **head)
{
  (void) deepest_surf_depth_at_start;

  if (head == NULL || *head == NULL || root_zone_depth_cm <= 0.0 || PET_timestep_cm <= 0.0) {
    return 0.0;
  }

  for (wetting_front *current = *head; current != NULL && current->next != NULL; ) {
    wetting_front *next = current->next;
    if (current->to_bottom == FALSE && current->is_WF_GW == TRUE && next->is_WF_GW == TRUE &&
        fabs(current->psi_cm - next->psi_cm) < 1.0e-3) {
      current = listDeleteFront(current->front_num, head, soil_type, soil_properties);
      if (current == NULL) {
        break;
      }
      continue;
    }
    current = next;
  }

  if (verbosity.compare("high") == 0) {
    printf("calculated mass before AET extraction: %.17lf \n",
           lgar_calc_mass_bal(cum_layer_thickness_cm, *head));
    listPrint(*head);
  }

  const int surface_front_count = listLength_surface(*head);
  double deepest_surface_depth_cm = 0.0;
  for (wetting_front *current = *head; current != NULL; current = current->next) {
    if (!current->is_WF_GW) {
      deepest_surface_depth_cm = current->depth_cm;
    }
  }

  const double top_of_to_root_zone_cm = fmin(root_zone_depth_cm, deepest_surface_depth_cm);
  surf_frac_rz = (root_zone_depth_cm > 0.0) ? fmin(1.0, top_of_to_root_zone_cm / root_zone_depth_cm) : 0.0;
  if (surf_frac_rz >= (1.0 - 1.0e-5)) {
    return 0.0;
  }

  const bool no_surface_needs_zero_depth_support =
    (surface_front_count == 0 &&
     count_mobile_to_fronts_in_root_zone(*head, top_of_to_root_zone_cm, root_zone_depth_cm) == 0);

  if (allow_root_zone_to_population && no_surface_needs_zero_depth_support) {
    lgarto_ensure_root_zone_to_population(num_layers, deepest_surface_depth_cm, root_zone_depth_cm,
                                          PET_timestep_cm, field_capacity_psi_cm, soil_type,
                                          cum_layer_thickness_cm, frozen_factor, soil_properties, head,
                                          false);
  }

  struct TOAETAllocation {
    wetting_front *front;
    double thickness_cm;
  };

  std::vector<TOAETAllocation> allocations;
  double segment_top_cm = top_of_to_root_zone_cm;
  wetting_front *active_zero_depth_support_front = NULL;
  if (no_surface_needs_zero_depth_support) {
    active_zero_depth_support_front = find_active_zero_depth_to_support_front(*head);
    if (active_zero_depth_support_front != NULL) {
      const double support_interval_bottom_cm =
        find_zero_depth_support_interval_bottom_cm(active_zero_depth_support_front, root_zone_depth_cm);
      const double support_thickness_cm = fmax(0.0, support_interval_bottom_cm - top_of_to_root_zone_cm);
      if (support_thickness_cm > 0.0) {
        allocations.push_back({active_zero_depth_support_front, support_thickness_cm});
        segment_top_cm = fmax(segment_top_cm, support_interval_bottom_cm);
      }
    }
  }

  for (wetting_front *current = *head; current != NULL; current = current->next) {
    if (!is_to_aet_eligible_front(current)) {
      continue;
    }

    if (current == active_zero_depth_support_front) {
      continue;
    }

    if (current->depth_cm >= root_zone_depth_cm) {
      break;
    }

    const double thickness_cm = fmax(0.0, current->depth_cm - segment_top_cm);
    allocations.push_back({current, thickness_cm});
    segment_top_cm = fmax(segment_top_cm, current->depth_cm);
  }

  if (!allocations.empty() && allocations.back().front != active_zero_depth_support_front) {
    allocations.back().thickness_cm +=
      fmax(0.0, root_zone_depth_cm - allocations.back().front->depth_cm);
  }

  wetting_front *surface_supported_interval_limited_front = NULL;
  if (surface_front_count > 0) {
    surface_supported_interval_limited_front = find_surface_supported_to_aet_target_front(*head);
  }

  double cumulative_ET_from_TO_WFs_cm = 0.0;
  double unmet_interval_limited_to_aet_cm = 0.0;
  for (const TOAETAllocation &allocation : allocations) {
    wetting_front *current = allocation.front;
    if (current == NULL || current->next == NULL || allocation.thickness_cm <= 0.0) {
      continue;
    }

    const double temp_AET_value =
      calc_aet_for_individual_TO_WFs(current->front_num, allocation.thickness_cm, root_zone_depth_cm,
                                     PET_timestep_cm, timestep_h, wilting_point_psi_cm,
                                     field_capacity_psi_cm, soil_type, soil_properties, head);
    if (temp_AET_value <= 0.0) {
      continue;
    }

    const double delta_theta = current->next->theta - current->theta;
    if (fabs(delta_theta) <= 1.0e-12) {
      if (verbosity.compare("high") == 0) {
        printf("Skipping TO AET extraction for front %d due to near-zero delta_theta.\n",
               current->front_num);
      }
      continue;
    }

    const bool limit_this_to_aet_extraction_to_local_interval =
      (surface_supported_interval_limited_front != NULL &&
       current == surface_supported_interval_limited_front &&
       current->next != NULL &&
       current->next->to_bottom == TRUE);
    const double capped_AET_value =
      limit_this_to_aet_extraction_to_local_interval ?
      cap_to_aet_extraction_to_local_interval_cm(current, temp_AET_value) :
      temp_AET_value;
    if (capped_AET_value <= 0.0) {
      if (verbosity.compare("high") == 0) {
        printf("Skipping TO AET extraction for front %d because it would overtake the next GW front.\n",
               current->front_num);
      }
      continue;
    }

    if (verbosity.compare("high") == 0 && capped_AET_value + 1.0e-12 < temp_AET_value) {
      printf("Capping TO AET extraction for front %d from %.17lf cm to %.17lf cm to keep it within its local interval.\n",
             current->front_num, temp_AET_value, capped_AET_value);
    }

    unmet_interval_limited_to_aet_cm += fmax(0.0, temp_AET_value - capped_AET_value);
    current->depth_cm += capped_AET_value / delta_theta;
    cumulative_ET_from_TO_WFs_cm += capped_AET_value;
  }

  if (surface_front_count > 0) {
    wetting_front *reference_front = find_surface_supported_to_aet_target_front(*head);
    const double missing_to_thickness_cm =
      fmax(0.0, root_zone_depth_cm - deepest_surface_depth_cm);

    if (reference_front != NULL && missing_to_thickness_cm > 0.0) {
      double missing_to_aet_cm =
        calc_aet_for_individual_TO_WFs(reference_front->front_num, missing_to_thickness_cm,
                                       root_zone_depth_cm, PET_timestep_cm, timestep_h,
                                       wilting_point_psi_cm, field_capacity_psi_cm, soil_type,
                                       soil_properties, head);
      missing_to_aet_cm = fmax(0.0, missing_to_aet_cm - cumulative_ET_from_TO_WFs_cm);
      missing_to_aet_cm += unmet_interval_limited_to_aet_cm;

      if (missing_to_aet_cm > 1.0e-12) {
        const double extracted_to_aet_from_surface_cm =
          lgarto_extract_missing_to_aet_from_surface_WFs(missing_to_aet_cm, cum_layer_thickness_cm,
                                                         soil_type, frozen_factor, soil_properties, head);
        cumulative_ET_from_TO_WFs_cm += extracted_to_aet_from_surface_cm;
        missing_to_aet_cm = fmax(0.0, missing_to_aet_cm - extracted_to_aet_from_surface_cm);

        if (verbosity.compare("high") == 0 && extracted_to_aet_from_surface_cm > 0.0) {
          printf("TO AET fallback extracted %.17lf cm from surface WFs using reference front %d.\n",
                 extracted_to_aet_from_surface_cm, reference_front->front_num);
        }
      }

      if (missing_to_aet_cm > 1.0e-12) {
        cumulative_ET_from_TO_WFs_cm +=
          lgarto_extract_missing_to_aet_from_to_chain(reference_front, missing_to_aet_cm,
                                                      cum_layer_thickness_cm, soil_type, frozen_factor,
                                                      soil_properties, head);
      }
    }
  }
  else if (unmet_interval_limited_to_aet_cm > 1.0e-12 && !allocations.empty()) {
    wetting_front *fallback_front = allocations.back().front;
    cumulative_ET_from_TO_WFs_cm +=
      lgarto_extract_missing_to_aet_from_to_chain(fallback_front, unmet_interval_limited_to_aet_cm,
                                                  cum_layer_thickness_cm, soil_type, frozen_factor,
                                                  soil_properties, head);
  }

  const bool no_surface_needs_post_aet_sparse_support =
    (surface_front_count == 0 && PET_timestep_cm > 0.0 &&
     has_sparse_no_surface_to_population(*head, root_zone_depth_cm));

  if (allow_root_zone_to_population && no_surface_needs_post_aet_sparse_support &&
      cumulative_ET_from_TO_WFs_cm > 0.0) {
    lgarto_ensure_root_zone_to_population(num_layers, deepest_surface_depth_cm, root_zone_depth_cm,
                                          PET_timestep_cm, field_capacity_psi_cm, soil_type,
                                          cum_layer_thickness_cm, frozen_factor, soil_properties, head,
                                          true);
  }

  if (verbosity.compare("high") == 0 && cumulative_ET_from_TO_WFs_cm > 0.0) {
    printf("WFs after TO WF depths updated via AET extraction: \n");
    listPrint(*head);
  }

  return cumulative_ET_from_TO_WFs_cm;
}

#endif
