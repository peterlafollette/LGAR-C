#include "../include/all.hxx"
#include <iostream>
#include <fstream>
#include <string.h>
#include <sstream>
#include <stddef.h>
#include <cstddef>

using namespace std;

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

typedef struct {
  double depth_cm;        // INTERNAL units (cm) after dividing by 10
  double theta;
  int layer_num;
  int front_num;
  double psi_cm;          // INTERNAL units (cm) after dividing by 10
  double dzdt_cm_per_h;   // INTERNAL units (cm/h), stored directly in file
  bool is_WF_GW;          // true for TO/GW wetting fronts, false for surface wetting fronts
} WFRecord;

static int cmp_record_by_front_num(const void *a, const void *b)
{
  const WFRecord *ra = (const WFRecord*)a;
  const WFRecord *rb = (const WFRecord*)b;
  return (ra->front_num - rb->front_num);
}

// Reads the next non-empty, non-comment line into buf. Returns true on success.
static bool read_next_data_line(FILE *fp, char *buf, size_t buflen)
{
  while (fgets(buf, buflen, fp)) {
    // skip comment lines and empty lines
    if (buf[0] == '#') continue;

    // skip whitespace-only lines
    bool any = false;
    for (size_t i = 0; buf[i] != '\0'; i++) {
      if (buf[i] != ' ' && buf[i] != '\t' && buf[i] != '\r' && buf[i] != '\n') {
        any = true;
        break;
      }
    }
    if (!any) continue;

    return true;
  }
  return false;
}

/*
  Parses a state line like:
    [(180.000000,0.197519,1,1,20000.000000,0.000010,1)|(...)]
  Using FILE units:
    depth and psi are written as *10 in the file and converted back to internal cm by /10.
    dzdt is written/read directly in cm/h.
    is_WF_GW is optional for backward compatibility: 1 means TO/GW, 0 means surface.
  This parser REQUIRES dzdt to be present. It will fail on old 5-field files.
*/
static int parse_state_line_to_records(const char *line_in, WFRecord **recs_out, int *n_out)
{
  *recs_out = NULL;
  *n_out = 0;

  char *line = strdup(line_in);
  if (!line) return 1;

  char *lbr = strchr(line, '[');
  char *rbr = strrchr(line, ']');
  if (!lbr || !rbr || rbr <= lbr) {
    free(line);
    return 2;
  }
  *rbr = '\0';
  char *p = lbr + 1;

  int cap = 16;
  int n = 0;
  WFRecord *recs = (WFRecord*)malloc(sizeof(WFRecord) * cap);
  if (!recs) {
    free(line);
    return 3;
  }

  char *saveptr = NULL;
  for (char *tok = strtok_r(p, "|", &saveptr); tok != NULL; tok = strtok_r(NULL, "|", &saveptr)) {
    while (*tok == ' ' || *tok == '\t' || *tok == ',') tok++;

    double depth_file, theta, psi_file, dzdt_cm_per_h;
    int layer_num, front_num;
    int is_WF_GW_int = 0;

    int matched = sscanf(tok, " ( %lf , %lf , %d , %d , %lf , %lf , %d ) ",
                         &depth_file, &theta, &layer_num, &front_num,
                         &psi_file, &dzdt_cm_per_h, &is_WF_GW_int);

    if (matched != 7) {
      matched = sscanf(tok, " ( %lf , %lf , %d , %d , %lf , %lf ) ",
                       &depth_file, &theta, &layer_num, &front_num,
                       &psi_file, &dzdt_cm_per_h);
      is_WF_GW_int = 0;
    }

    if (matched != 6 && matched != 7) {
      free(recs);
      free(line);
      return 4;
    }

    if (n == cap) {
      cap *= 2;
      WFRecord *tmp = (WFRecord*)realloc(recs, sizeof(WFRecord) * cap);
      if (!tmp) {
        free(recs);
        free(line);
        return 5;
      }
      recs = tmp;
    }

    recs[n].depth_cm       = depth_file / 10.0;   // undo writer scaling
    recs[n].theta          = theta;
    recs[n].layer_num      = layer_num;
    recs[n].front_num      = front_num;
    recs[n].psi_cm         = psi_file / 10.0;     // undo writer scaling
    recs[n].dzdt_cm_per_h  = dzdt_cm_per_h;       // stored directly
    recs[n].is_WF_GW       = is_WF_GW_int != 0;
    n++;
  }

  free(line);
  *recs_out = recs;
  *n_out = n;
  return 0;
}


extern void InitializeWettingFrontsFromCSV(
    int num_layers,
    const char *data_layers_csv_path,
    int *layer_soil_type,
    double *cum_layer_thickness_cm,
    double *frozen_factor,
    struct wetting_front **head,
    struct soil_properties_ *soil_properties)
{
  (void)cum_layer_thickness_cm;  // currently unused in this loader

  *head = NULL;

  FILE *fp = fopen(data_layers_csv_path, "r");
  if (!fp) {
    fprintf(stderr, "ERROR: could not open data_layers file: %s\n", data_layers_csv_path);
    exit(1);
  }

  // 1) Find the first timestep header line "# Timestep = ..."
  char line[65536];
  bool found_header = false;
  while (fgets(line, sizeof(line), fp)) {
    if (strncmp(line, "# Timestep", 10) == 0) {
      found_header = true;
      break;
    }
  }
  if (!found_header) {
    fclose(fp);
    fprintf(stderr, "ERROR: no '# Timestep' header found in %s\n", data_layers_csv_path);
    exit(1);
  }

  // 2) Next non-comment line should be the state line "[...]"
  char state_line[65536];
  if (!read_next_data_line(fp, state_line, sizeof(state_line))) {
    fclose(fp);
    fprintf(stderr, "ERROR: found timestep header but no state line after it in %s\n", data_layers_csv_path);
    exit(1);
  }
  fclose(fp);

  WFRecord *recs = NULL;
  int nrecs = 0;
  int perr = parse_state_line_to_records(state_line, &recs, &nrecs);
  if (perr != 0 || !recs || nrecs <= 0) {
    fprintf(stderr, "ERROR: failed to parse state line in %s (err=%d). "
                    "Expected tuples of the form "
                    "(depth_x10,theta,layer_num,front_num,psi_x10,dzdt_cm_per_h[,is_WF_GW])\n",
            data_layers_csv_path, perr);
    exit(1);
  }

  // 3) Insert in front_num order so listInsertFront() works (and does not renumber)
  qsort(recs, nrecs, sizeof(WFRecord), cmp_record_by_front_num);

  // Sanity: ensure we start at 1 and are contiguous
  if (recs[0].front_num != 1) {
    fprintf(stderr, "ERROR: first front_num in file is %d, but listInsertFront requires starting at 1\n",
            recs[0].front_num);
    free(recs);
    exit(1);
  }

  for (int i = 1; i < nrecs; i++) {
    if (recs[i].front_num != recs[i-1].front_num + 1) {
      fprintf(stderr, "ERROR: front_num values are not contiguous at %d -> %d\n",
              recs[i-1].front_num, recs[i].front_num);
      free(recs);
      exit(1);
    }
  }

  for (int i = 0; i < nrecs; i++) {
    WFRecord r = recs[i];

    if (r.layer_num < 1 || r.layer_num > num_layers) {
      fprintf(stderr, "ERROR: record layer_num=%d out of range 1..%d\n", r.layer_num, num_layers);
      free(recs);
      exit(1);
    }

    // Temporarily set bottom_flag false; we will set to_bottom correctly in a post-pass
    bool bottom_flag = false;

    struct wetting_front *current =
        listInsertFront(r.depth_cm, r.theta, r.front_num, r.layer_num, bottom_flag, head);

    if (current == NULL) {
      fprintf(stderr, "ERROR: listInsertFront returned NULL inserting front_num=%d\n", r.front_num);
      free(recs);
      exit(1);
    }

    // Trust file theta, psi, and dzdt
    current->psi_cm = r.psi_cm;
    current->dzdt_cm_per_h = r.dzdt_cm_per_h;
    current->is_WF_GW = r.is_WF_GW;

    // Compute K from theta
    int soil = layer_soil_type[r.layer_num];
    double Se = calc_Se_from_theta(current->theta,
                                   soil_properties[soil].theta_e,
                                   soil_properties[soil].theta_r);

    double Ksat_cm_per_h = frozen_factor[r.layer_num] * soil_properties[soil].Ksat_cm_per_h;
    current->K_cm_per_h = calc_K_from_Se(Se, Ksat_cm_per_h, soil_properties[soil].vg_m);
  }

  // 4) Post-pass: set to_bottom correctly:
  //    true when this is the deepest front in its layer (next front has greater layer_num, or next is NULL)
  {
    struct wetting_front *cur = *head;
    while (cur != NULL) {
      if (cur->next == NULL) {
        cur->to_bottom = true;
      } else {
        cur->to_bottom = (cur->next->layer_num > cur->layer_num);
      }
      cur = cur->next;
    }
  }

  free(recs);
}


typedef struct {
  double CR_fast_storage_cm;
  double CR_slow_storage_cm;
  double volon_timestep_cm;
  bool runoff_in_prev_step;
  double precip_previous_timestep_cm;
  bool has_groundwater_depth_cm;
  double groundwater_depth_cm;
} nonvadoseRestartState;

static bool parse_bool_01(const char *s, bool *out)
{
  if (strcmp(s, "1") == 0) {
    *out = true;
    return true;
  }
  if (strcmp(s, "0") == 0) {
    *out = false;
    return true;
  }
  return false;
}

static int parse_non_vadose_state_kv_line(const char *line_in, nonvadoseRestartState *rst)
{
  char *line = strdup(line_in);
  if (!line) return 1;

  bool got_fast = false;
  bool got_slow = false;
  bool got_volon = false;
  bool got_runoff = false;
  bool got_precip_prev = false;
  rst->has_groundwater_depth_cm = false;
  rst->groundwater_depth_cm = 0.0;

  char *saveptr = NULL;
  for (char *tok = strtok_r(line, ",\r\n", &saveptr);
       tok != NULL;
       tok = strtok_r(NULL, ",\r\n", &saveptr)) {

    while (*tok == ' ' || *tok == '\t') tok++;

    char *eq = strchr(tok, '=');
    if (!eq) {
      free(line);
      return 2;
    }

    *eq = '\0';
    const char *key = tok;
    const char *val = eq + 1;

    if (strcmp(key, "CR_fast_storage_cm") == 0) {
      rst->CR_fast_storage_cm = strtod(val, NULL);
      got_fast = true;
    }
    else if (strcmp(key, "CR_slow_storage_cm") == 0) {
      rst->CR_slow_storage_cm = strtod(val, NULL);
      got_slow = true;
    }
    else if (strcmp(key, "volon_timestep_cm") == 0) {
      rst->volon_timestep_cm = strtod(val, NULL);
      got_volon = true;
    }
    else if (strcmp(key, "runoff_in_prev_step") == 0) {
      if (!parse_bool_01(val, &rst->runoff_in_prev_step)) {
        free(line);
        return 3;
      }
      got_runoff = true;
    }
    else if (strcmp(key, "precip_previous_timestep_cm") == 0) {
      rst->precip_previous_timestep_cm = strtod(val, NULL);
      got_precip_prev = true;
    }
    else if (strcmp(key, "groundwater_depth_cm") == 0) {
      rst->groundwater_depth_cm = strtod(val, NULL);
      rst->has_groundwater_depth_cm = true;
    }
  }

  free(line);

  if (!(got_fast && got_slow && got_volon && got_runoff && got_precip_prev)) {
    return 4;
  }

  return 0;
}

extern void InitializenonvadoseStateFromCSV(
    const char *non_vadose_state_csv_path,
    struct model_state *state)
{
  FILE *fp = fopen(non_vadose_state_csv_path, "r");
  if (!fp) {
    fprintf(stderr, "ERROR: could not open non vadose state file: %s\n",
            non_vadose_state_csv_path);
    exit(1);
  }

  char line[65536];

  if (!read_next_data_line(fp, line, sizeof(line))) {
    fclose(fp);
    fprintf(stderr, "ERROR: no data lines found in non vadose state file: %s\n",
            non_vadose_state_csv_path);
    exit(1);
  }

  nonvadoseRestartState rst;
  int perr = 0;

  // Allow either:
  //   Time,CR_fast_storage_cm,...
  // followed by:
  //   2024-...,CR_fast_storage_cm=...
  //
  // or simply:
  //   CR_fast_storage_cm=...,CR_slow_storage_cm=...
  if (strncmp(line, "Time,", 5) == 0) {
    if (!read_next_data_line(fp, line, sizeof(line))) {
      fclose(fp);
      fprintf(stderr, "ERROR: header found but no restart row in %s\n",
              non_vadose_state_csv_path);
      exit(1);
    }

    char *comma = strchr(line, ',');
    if (!comma || *(comma + 1) == '\0') {
      fclose(fp);
      fprintf(stderr, "ERROR: malformed restart row in %s\n",
              non_vadose_state_csv_path);
      exit(1);
    }

    perr = parse_non_vadose_state_kv_line(comma + 1, &rst);
  }
  else {
    perr = parse_non_vadose_state_kv_line(line, &rst);
  }

  fclose(fp);

  if (perr != 0) {
    fprintf(stderr, "ERROR: failed to parse non vadose restart state in %s (err=%d)\n",
            non_vadose_state_csv_path, perr);
    exit(1);
  }

  state->lgar_mass_balance.CR_fast_storage_cm = rst.CR_fast_storage_cm;
  state->lgar_mass_balance.CR_slow_storage_cm = rst.CR_slow_storage_cm;
  state->lgar_mass_balance.volon_timestep_cm = rst.volon_timestep_cm;
  state->lgar_bmi_params.runoff_in_prev_step = rst.runoff_in_prev_step;
  state->lgar_bmi_params.precip_previous_timestep_cm = rst.precip_previous_timestep_cm;
  if (rst.has_groundwater_depth_cm) {
    state->lgar_bmi_params.groundwater_depth_cm = rst.groundwater_depth_cm;
  }

  // Keep the "current total CR volume" field consistent with the restored storages.
  state->lgar_mass_balance.volCRend_cm =
      rst.CR_fast_storage_cm + rst.CR_slow_storage_cm;

  state->lgar_mass_balance.volCRend_timestep_cm =
      rst.CR_fast_storage_cm + rst.CR_slow_storage_cm;
}


typedef struct {
  int num_giuh_ordinates;
  double *queue_vals;
} GIUHRestartState;

static int parse_giuh_state_line(
    const char *line_in,
    GIUHRestartState *rst)
{
  rst->num_giuh_ordinates = -1;
  rst->queue_vals = NULL;

  char *line = strdup(line_in);
  if (!line) return 1;

  char *num_ptr = strstr(line, "num_giuh_ordinates=");
  char *queue_ptr = strstr(line, "queue=[");
  if (!num_ptr || !queue_ptr) {
    free(line);
    return 2;
  }

  if (sscanf(num_ptr, "num_giuh_ordinates=%d", &rst->num_giuh_ordinates) != 1) {
    free(line);
    return 3;
  }

  char *lbr = strchr(queue_ptr, '[');
  char *rbr = strrchr(queue_ptr, ']');
  if (!lbr || !rbr || rbr <= lbr) {
    free(line);
    return 4;
  }

  *rbr = '\0';
  char *vals = lbr + 1;

  int expected_n = rst->num_giuh_ordinates + 1;
  rst->queue_vals = (double*)malloc(sizeof(double) * expected_n);
  if (!rst->queue_vals) {
    free(line);
    return 5;
  }

  int nread = 0;
  char *saveptr = NULL;
  for (char *tok = strtok_r(vals, ",", &saveptr);
       tok != NULL;
       tok = strtok_r(NULL, ",", &saveptr)) {

    while (*tok == ' ' || *tok == '\t') tok++;

    if (nread >= expected_n) {
      free(rst->queue_vals);
      free(line);
      return 6;
    }

    rst->queue_vals[nread++] = strtod(tok, NULL);
  }

  free(line);

  if (nread != expected_n) {
    free(rst->queue_vals);
    rst->queue_vals = NULL;
    return 7;
  }

  return 0;
}

extern void InitializeGIUHRunoffQueueFromCSV(
    const char *giuh_state_csv_path,
    double *giuh_runoff_queue,
    int num_giuh_ordinates)
{
  FILE *fp = fopen(giuh_state_csv_path, "r");
  if (!fp) {
    fprintf(stderr, "ERROR: could not open GIUH state file: %s\n", giuh_state_csv_path);
    exit(1);
  }

  char line[65536];
  if (!read_next_data_line(fp, line, sizeof(line))) {
    fclose(fp);
    fprintf(stderr, "ERROR: no data lines found in GIUH state file: %s\n", giuh_state_csv_path);
    exit(1);
  }

  GIUHRestartState rst;
  int perr = 0;

  if (strncmp(line, "Time,", 5) == 0) {
    if (!read_next_data_line(fp, line, sizeof(line))) {
      fclose(fp);
      fprintf(stderr, "ERROR: header found but no GIUH restart row in %s\n", giuh_state_csv_path);
      exit(1);
    }

    char *comma = strchr(line, ',');
    if (!comma || *(comma + 1) == '\0') {
      fclose(fp);
      fprintf(stderr, "ERROR: malformed GIUH restart row in %s\n", giuh_state_csv_path);
      exit(1);
    }

    perr = parse_giuh_state_line(comma + 1, &rst);
  }
  else {
    perr = parse_giuh_state_line(line, &rst);
  }

  fclose(fp);

  if (perr != 0) {
    fprintf(stderr, "ERROR: failed to parse GIUH restart state in %s (err=%d)\n",
            giuh_state_csv_path, perr);
    exit(1);
  }

  if (rst.num_giuh_ordinates != num_giuh_ordinates) {
    fprintf(stderr,
            "ERROR: GIUH restart file has num_giuh_ordinates=%d but model expects %d\n",
            rst.num_giuh_ordinates, num_giuh_ordinates);
    free(rst.queue_vals);
    exit(1);
  }

  for (int i = 0; i <= num_giuh_ordinates; i++) {
    giuh_runoff_queue[i] = rst.queue_vals[i];
  }

  free(rst.queue_vals);
}
