#ifndef CQM_CXX_INCLUDE
#define CQM_CXX_INCLUDE

#include "../include/all.hxx"

extern double calc_QF_Q(double subtimestep_h, double a, double b, double precip_for_QF_subtimestep_cm, double *QF_storage_cm){
    double QF = subtimestep_h * (a * pow(*QF_storage_cm, b));
    if (*QF_storage_cm < 0.01){ // idea here is that we do want the GW contribution to actually become 0 in arid or semi arid environments
        QF = 0.0;
    }
    if (*QF_storage_cm + (subtimestep_h * precip_for_QF_subtimestep_cm - QF) > 0.0){
        *QF_storage_cm += (subtimestep_h * precip_for_QF_subtimestep_cm - QF);
    }
    else {
        QF = *QF_storage_cm + precip_for_QF_subtimestep_cm;
        *QF_storage_cm = 0.0;
    }
    return(QF);
}

#endif

