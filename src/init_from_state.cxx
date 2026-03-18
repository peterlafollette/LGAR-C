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
    [(180.000000,0.197519,1,1,20000.000000,0.000010)|(...)]
  Using FILE units:
    depth and psi are written as *10 in the file and converted back to internal cm by /10.
    dzdt is written/read directly in cm/h.
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

    int matched = sscanf(tok, " ( %lf , %lf , %d , %d , %lf , %lf ) ",
                         &depth_file, &theta, &layer_num, &front_num,
                         &psi_file, &dzdt_cm_per_h);

    if (matched != 6) {
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
                    "(depth_x10,theta,layer_num,front_num,psi_x10,dzdt_cm_per_h)\n",
            data_layers_csv_path, perr);
    exit(1);
  }

  // 3) Insert in front_num order so listInsertFront() works (and doesn’t renumber)
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