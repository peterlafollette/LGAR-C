#ifndef AET_CXX_INCLUDE
#define AET_CXX_INCLUDE

#include "../include/all.hxx"

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

  if (current->depth_cm == 0.0 || current->is_WF_GW == FALSE || frac <= 0.0) {
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

extern double lgarto_calc_aet_from_TO_WFs(int num_layers,
					  double deepest_surf_depth_at_start,
					  double root_zone_depth_cm,
					  double PET_timestep_cm,
					  double timestep_h,
					  double surf_frac_rz,
					  double wilting_point_psi_cm,
					  double field_capacity_psi_cm,
					  int *soil_type,
					  double *cum_layer_thickness_cm,
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

  struct TOAETAllocation {
    wetting_front *front;
    double thickness_cm;
  };

  std::vector<TOAETAllocation> allocations;
  double segment_top_cm = top_of_to_root_zone_cm;
  for (wetting_front *current = *head; current != NULL; current = current->next) {
    if (!is_to_aet_eligible_front(current)) {
      continue;
    }

    if (current->depth_cm >= root_zone_depth_cm) {
      break;
    }

    const double thickness_cm = fmax(0.0, current->depth_cm - segment_top_cm);
    allocations.push_back({current, thickness_cm});
    segment_top_cm = fmax(segment_top_cm, current->depth_cm);
  }

  if (allocations.empty()) {
    return 0.0;
  }

  allocations.back().thickness_cm += fmax(0.0, root_zone_depth_cm - allocations.back().front->depth_cm);

  double cumulative_ET_from_TO_WFs_cm = 0.0;
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

    current->depth_cm += temp_AET_value / delta_theta;
    cumulative_ET_from_TO_WFs_cm += temp_AET_value;
  }

  if (verbosity.compare("high") == 0 && cumulative_ET_from_TO_WFs_cm > 0.0) {
    printf("WFs after TO WF depths updated via AET extraction: \n");
    listPrint(*head);
  }

  return cumulative_ET_from_TO_WFs_cm;
}

#endif
