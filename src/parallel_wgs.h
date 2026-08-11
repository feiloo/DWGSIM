#ifndef PARALLEL_WGS_H
#define PARALLEL_WGS_H

#include "dwgsim_opt.h"

int dwgsim_parallel_wgs_supported(const dwgsim_opt_t *opt,
                                  const char *reference_path);
int dwgsim_parallel_wgs_run(const dwgsim_opt_t *opt,
                            const char *reference_path,
                            const char *output_prefix);

#endif
