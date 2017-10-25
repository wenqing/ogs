AddTest(
    NAME LARGE_2D_ThermalConvection_constviscosity
    PATH Parabolic/HT/ConstViscosity
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_5500x5500.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    ABSTOL 1e-10 RELTOL 1e-16
    DIFF_DATA
    square_5500x5500.vtu ConstViscosityThermalConvection_pcs_0_ts_149_t_50000000000.000000.vtu T_ref T
    VIS ConstViscosityThermalConvection_pcs_0_ts_149_t_50000000000.000000.vtu
)

AddTest(
    NAME LARGE_2D_ThermalConvection_constviscosityStaggeredScheme
    PATH Parabolic/HT/StaggeredCoupling/ConstViscosity
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_5500x5500_staggered_scheme.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    ABSTOL 1e-10 RELTOL 1.e-9
    DIFF_DATA
    square_5500x5500.vtu ConstViscosityThermalConvection_pcs_1_ts_149_t_50000000000.000000.vtu T_ref T
    square_5500x5500.vtu ConstViscosityThermalConvection_pcs_1_ts_149_t_50000000000.000000.vtu p_ref p
    VIS ConstViscosityThermalConvection_pcs_1_ts_149_t_50000000000.000000.vtu
)

# MPI/PETSc tests
AddTest(
    NAME Parallel_LARGE_2D_ThermalConvection_constviscosity
    PATH Parabolic/HT/ConstViscosity
    EXECUTABLE_ARGS square_5500x5500.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 4
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    ABSTOL 1e-15 RELTOL 1e-14
    DIFF_DATA
    ConstViscosityThermalConvection_pcs_0_ts_149_t_50000000000_000000_0.vtu ConstViscosityThermalConvection_pcs_0_ts_149_t_50000000000_000000_0.vtu p p
    ConstViscosityThermalConvection_pcs_0_ts_149_t_50000000000_000000_1.vtu ConstViscosityThermalConvection_pcs_0_ts_149_t_50000000000_000000_1.vtu p p
    ConstViscosityThermalConvection_pcs_0_ts_149_t_50000000000_000000_2.vtu ConstViscosityThermalConvection_pcs_0_ts_149_t_50000000000_000000_2.vtu p p
    ConstViscosityThermalConvection_pcs_0_ts_149_t_50000000000_000000_3.vtu ConstViscosityThermalConvection_pcs_0_ts_149_t_50000000000_000000_3.vtu p p
    ConstViscosityThermalConvection_pcs_0_ts_149_t_50000000000_000000_0.vtu ConstViscosityThermalConvection_pcs_0_ts_149_t_50000000000_000000_0.vtu T T
    ConstViscosityThermalConvection_pcs_0_ts_149_t_50000000000_000000_1.vtu ConstViscosityThermalConvection_pcs_0_ts_149_t_50000000000_000000_1.vtu T T
    ConstViscosityThermalConvection_pcs_0_ts_149_t_50000000000_000000_2.vtu ConstViscosityThermalConvection_pcs_0_ts_149_t_50000000000_000000_2.vtu T T
    ConstViscosityThermalConvection_pcs_0_ts_149_t_50000000000_000000_3.vtu ConstViscosityThermalConvection_pcs_0_ts_149_t_50000000000_000000_3.vtu T T
)

AddTest(
    NAME LARGE_2D_ThermalConvection_constviscosityStaggeredScheme
    PATH Parabolic/HT/StaggeredCoupling/ConstViscosity
    EXECUTABLE_ARGS square_5500x5500_staggered_scheme.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 1
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    ABSTOL 1e-10 RELTOL 1.e-9
    DIFF_DATA
    square_5500x5500.vtu ConstViscosityThermalConvection_pcs_1_ts_149_t_50000000000_000000_0.vtu T_ref T
    square_5500x5500.vtu ConstViscosityThermalConvection_pcs_1_ts_149_t_50000000000_000000_0.vtu p_ref p
    VIS ConstViscosityThermalConvection_pcs_1_ts_149_t_50000000000.000000_0.vtu
)
