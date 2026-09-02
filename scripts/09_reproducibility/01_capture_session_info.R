# Capture the exact software environment after a successful validated rerun.
source(file.path('config','project_paths.R'))
source(file.path('R','project_io.R'))
out <- ensure_dir(file.path(RESULTS_DIR,'reproducibility'))
sink(file.path(out,'sessionInfo.txt')); print(sessionInfo()); sink()
if(requireNamespace('renv',quietly=TRUE)) message('renv is available. After validating outputs, run renv::snapshot() to create renv.lock.')
