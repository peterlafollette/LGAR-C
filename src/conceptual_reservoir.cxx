#ifndef CR_CXX_INCLUDE
#define CR_CXX_INCLUDE

#include "../include/all.hxx"

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

    // === FAST reservoir outflow ===
    double Q_fast = calc_single_CR_Q(subtimestep_h, a_fast, b_fast, input_fast,
                                     fast_discharge_threshold_cm,
                                     CR_fast_storage_cm);

    // === SLOW reservoir outflow ===
    double Q_slow = calc_single_CR_Q(subtimestep_h, a_slow, b_slow, input_slow,
                                     slow_discharge_threshold_cm,
                                     CR_slow_storage_cm);

    return Q_fast + Q_slow;
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

#endif
