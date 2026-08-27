#ifndef CR_CXX_INCLUDE
#define CR_CXX_INCLUDE

#include "../include/all.hxx"
#include <cmath>

static double calc_single_CR_Q(double subtimestep_h,
                               double a,
                               double b,
                               double input_cm_per_h,
                               double discharge_threshold_cm,
                               double *CR_storage_cm)
{
    if (CR_storage_cm == NULL) {
        return 0.0;
    }

    double storage_cm = fmax(*CR_storage_cm, 0.0);
    if (subtimestep_h <= 0.0) {
        *CR_storage_cm = storage_cm;
        return 0.0;
    }

    const double threshold_cm = fmax(discharge_threshold_cm, 0.0);
    const double input_cm = subtimestep_h * input_cm_per_h;
    const double active_storage_start_cm = fmax(storage_cm - threshold_cm, 0.0);

    double Q_cm = 0.0;
    if (a > 0.0 && b > 0.0 && active_storage_start_cm > 0.0) {
        Q_cm = subtimestep_h * (a * pow(active_storage_start_cm, b));
    }

    const double storage_after_input_cm = fmax(storage_cm + input_cm, 0.0);
    const double streamflow_available_cm =
        fmax(storage_after_input_cm - threshold_cm, 0.0);
    Q_cm = fmin(fmax(Q_cm, 0.0), streamflow_available_cm);

    storage_cm = storage_after_input_cm - Q_cm;
    if (storage_cm < 1.0e-12) {
        storage_cm = 0.0;
    }
    *CR_storage_cm = storage_cm;

    return Q_cm;
}

static double extract_from_CR_storage(double demand_cm, double *CR_storage_cm)
{
    if (CR_storage_cm == NULL || demand_cm <= 0.0 || *CR_storage_cm <= 0.0) {
        return 0.0;
    }

    const double extraction_cm = fmin(fmax(*CR_storage_cm, 0.0), demand_cm);
    *CR_storage_cm -= extraction_cm;
    if (*CR_storage_cm < 1.0e-12) {
        *CR_storage_cm = 0.0;
    }
    return extraction_cm;
}

extern double calc_CR_Q(
    double subtimestep_h,
    double a_fast, double a_slow,
    double b_fast, double b_slow,
    double fast_discharge_threshold_cm,
    double slow_discharge_threshold_cm,
    double frac_slow,  // fraction (0 - 1) of recharge going to slow reservoir
    double precip_for_CR_subtimestep_cm_per_h,
    double *CR_fast_storage_cm,
    double *CR_slow_storage_cm)
{
    // Partition recharge between fast and slow reservoirs
    double input_slow = precip_for_CR_subtimestep_cm_per_h * frac_slow;
    double input_fast = precip_for_CR_subtimestep_cm_per_h - input_slow; // implicit (1 - frac_slow)

    return calc_CR_Q_explicit_inputs(
        subtimestep_h, a_fast, a_slow, b_fast, b_slow,
        fast_discharge_threshold_cm, slow_discharge_threshold_cm,
        input_fast, input_slow, CR_fast_storage_cm, CR_slow_storage_cm);
}

extern double calc_CR_Q_explicit_inputs(
    double subtimestep_h,
    double a_fast, double a_slow,
    double b_fast, double b_slow,
    double fast_discharge_threshold_cm,
    double slow_discharge_threshold_cm,
    double fast_input_cm_per_h,
    double slow_input_cm_per_h,
    double *CR_fast_storage_cm,
    double *CR_slow_storage_cm)
{

    // === FAST reservoir outflow ===
    double Q_fast = calc_single_CR_Q(subtimestep_h, a_fast, b_fast,
                                     fast_input_cm_per_h,
                                     fast_discharge_threshold_cm,
                                     CR_fast_storage_cm);

    // === SLOW reservoir outflow ===
    double Q_slow = calc_single_CR_Q(subtimestep_h, a_slow, b_slow,
                                     slow_input_cm_per_h,
                                     slow_discharge_threshold_cm,
                                     CR_slow_storage_cm);

    return Q_fast + Q_slow;
}

extern double lgar_cap_CR_input_by_available_storage(
    double raw_CR_input_cm,
    double available_storage_cm,
    double same_substep_drainage_cm,
    double *accepted_CR_input_cm)
{
    if (!std::isfinite(raw_CR_input_cm) || raw_CR_input_cm <= 0.0) {
        if (accepted_CR_input_cm != NULL) {
            *accepted_CR_input_cm = 0.0;
        }
        return 0.0;
    }

    const double available_cm =
        std::isfinite(available_storage_cm) ? fmax(available_storage_cm, 0.0) : 0.0;
    const double drainage_allowance_cm =
        std::isfinite(same_substep_drainage_cm) ? fmax(same_substep_drainage_cm, 0.0) : 0.0;
    const double accepted_cm =
        fmin(raw_CR_input_cm, available_cm + drainage_allowance_cm);
    const double rejected_cm = fmax(0.0, raw_CR_input_cm - accepted_cm);

    if (accepted_CR_input_cm != NULL) {
        *accepted_CR_input_cm = accepted_cm;
    }
    return rejected_cm;
}

extern void lgar_partition_lower_boundary_flux_for_CR(
    bool route_positive_lower_boundary_flux_to_CR,
    double lower_boundary_flux_cm,
    double *percolation_cm,
    double *CR_input_cm,
    double *CR_fast_storage_cm,
    double *CR_slow_storage_cm,
    double *CR_storage_exchange_cm)
{
    if (percolation_cm == NULL || CR_input_cm == NULL) {
        return;
    }

    if (CR_storage_exchange_cm != NULL) {
        *CR_storage_exchange_cm = 0.0;
    }

    if (!route_positive_lower_boundary_flux_to_CR) {
        *percolation_cm += lower_boundary_flux_cm;
        return;
    }

    if (lower_boundary_flux_cm > 0.0) {
        *CR_input_cm += lower_boundary_flux_cm;
        if (CR_storage_exchange_cm != NULL) {
            *CR_storage_exchange_cm += lower_boundary_flux_cm;
        }
        return;
    }

    if (lower_boundary_flux_cm < 0.0) {
        double negative_recharge_demand_cm = -lower_boundary_flux_cm;
        double reservoir_extraction_cm =
            extract_from_CR_storage(negative_recharge_demand_cm,
                                    CR_fast_storage_cm);
        negative_recharge_demand_cm -= reservoir_extraction_cm;

        if (negative_recharge_demand_cm > 0.0) {
            const double slow_reservoir_extraction_cm =
                extract_from_CR_storage(negative_recharge_demand_cm,
                                        CR_slow_storage_cm);
            reservoir_extraction_cm += slow_reservoir_extraction_cm;
            negative_recharge_demand_cm -= slow_reservoir_extraction_cm;
        }

        if (CR_storage_exchange_cm != NULL) {
            *CR_storage_exchange_cm -= reservoir_extraction_cm;
        }

        if (negative_recharge_demand_cm > 0.0) {
            *percolation_cm -= negative_recharge_demand_cm;
        }
    }
}

extern void lgar_partition_dual_lower_boundary_fluxes_for_CR(
    bool route_positive_lower_boundary_flux_to_CR,
    double matrix_lower_boundary_flux_cm,
    double fracture_lower_boundary_flux_cm,
    double frac_slow,
    double *percolation_cm,
    double *fast_CR_input_cm,
    double *slow_CR_input_cm,
    double *CR_fast_storage_cm,
    double *CR_slow_storage_cm,
    double *CR_storage_exchange_cm)
{
    if (percolation_cm == NULL || fast_CR_input_cm == NULL ||
        slow_CR_input_cm == NULL) {
        return;
    }

    *percolation_cm = 0.0;
    *fast_CR_input_cm = 0.0;
    *slow_CR_input_cm = 0.0;
    if (CR_storage_exchange_cm != NULL) {
        *CR_storage_exchange_cm = 0.0;
    }

    if (!route_positive_lower_boundary_flux_to_CR) {
        *percolation_cm = matrix_lower_boundary_flux_cm +
                          fracture_lower_boundary_flux_cm;
        return;
    }

    // Combine both upward demands before touching shared storage so matrix
    // and fracture call order cannot affect which reservoir supplies them.
    const double upward_demand_cm =
        -fmin(matrix_lower_boundary_flux_cm, 0.0) -
        fmin(fracture_lower_boundary_flux_cm, 0.0);
    if (upward_demand_cm > 0.0) {
        double ignored_CR_input_cm = 0.0;
        lgar_partition_lower_boundary_flux_for_CR(
            true, -upward_demand_cm, percolation_cm, &ignored_CR_input_cm,
            CR_fast_storage_cm, CR_slow_storage_cm, CR_storage_exchange_cm);
    }

    const double matrix_recharge_cm =
        fmax(matrix_lower_boundary_flux_cm, 0.0);
    const double fracture_recharge_cm =
        fmax(fracture_lower_boundary_flux_cm, 0.0);
    const double slow_fraction = fmax(0.0, fmin(1.0, frac_slow));
    *slow_CR_input_cm = slow_fraction * matrix_recharge_cm;
    *fast_CR_input_cm = matrix_recharge_cm - *slow_CR_input_cm +
                        fracture_recharge_cm;
    if (CR_storage_exchange_cm != NULL) {
        *CR_storage_exchange_cm += *fast_CR_input_cm + *slow_CR_input_cm;
    }
}

#endif
