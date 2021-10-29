# to run:     mpiexec -n 6 ../fish-opt -i boussinesq32.i --n-threads=12
# peacock:    ~/projects/moose/python/peacock/peacock -r boussinesq32_out.e


[Mesh]
  file = 2d_inner_fluid_quad_corse.e
[]

mu = 0.001241                   # The viscosity [Pa*s]
rho = 1611                      # Density [kg/m³]
k = 1.10                        # Thermal conductivity [W/m*K]
cp = 2097.8                     # Specific heat capacity [J/kg*K]
alpha = 0.0003035               # Thermal expansion
temp_ref = 873.15               # Reference temperature for the thermal expansion law

vel = 'velocity'
velocity_interp_method = 'rc'
advected_interp_method = 'upwind'
cold_temp = 660
hot_temp = 740

[Variables]
  [u]
    type = INSFVVelocityVariable
    initial_condition = 1e-8
  []
  [v]
    type = INSFVVelocityVariable
    initial_condition = 1e-8
    scaling = 1e-2
  []
  [pressure]
    type = INSFVPressureVariable
  []
  [T]
    type = INSFVEnergyVariable
    scaling = 1e-6
  []
  [lambda]
    family = SCALAR
    order = FIRST
  []
[]

[AuxVariables]
  [U]
    order = CONSTANT
    family = MONOMIAL
    fv = true
  []
  [vel_x]
    order = FIRST
    family = MONOMIAL
  []
  [vel_y]
    order = FIRST
    family = MONOMIAL
  []
  [viz_T]
    order = FIRST
    family = MONOMIAL
  []
[]

[ICs]
  [temperature_gradient]
    type = FunctionIC
    variable = T
    function = 'gradient_temp'
  []
[]

[Functions]
  # linear temperature profile for initialization
  [gradient_temp]
    type = ParsedFunction
    value = 'if(x<0.01, hot_temp, if(x>0.59, cold_temp, hot_temp - (hot_temp-cold_temp) * x / 0.59))'
    vars = 'cold_temp hot_temp'
    vals = '${cold_temp} ${hot_temp}'
  []
[]

[AuxKernels]
  [mag]
    type = VectorMagnitudeAux
    variable = U
    x = u
    y = v
    execute_on = 'initial timestep_end'
  []
  [vel_x]
    type = ParsedAux
    variable = vel_x
    function = 'u'
    execute_on = 'initial timestep_end'
    args = 'u'
  []
  [vel_y]
    type = ParsedAux
    variable = vel_y
    function = 'v'
    execute_on = 'initial timestep_end'
    args = 'v'
  []
  [viz_T]
    type = ParsedAux
    variable = viz_T
    function = 'T'
    execute_on = 'initial timestep_end'
    args = 'T'
  []
[]

[FVKernels]
  [mass]
    type = INSFVMassAdvection
    variable = pressure
    vel = ${vel}
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    u = u
    v = v
    pressure = pressure
    mu = 'mu'
    rho = ${rho}
  []
  [mean_zero_pressure]
    type = FVScalarLagrangeMultiplier
    variable = pressure
    lambda = lambda
  []

  [u_time]
    type = INSFVMomentumTimeDerivative
    variable = 'u'
    rho = ${rho}
  []
  [u_advection]
    type = INSFVMomentumAdvection
    variable = u
    advected_quantity = 'rhou'
    vel = ${vel}
    velocity_interp_method = ${velocity_interp_method}
    advected_interp_method = ${advected_interp_method}
    pressure = pressure
    u = u
    v = v
    mu = 'mu'
    rho = ${rho}
  []
  [u_viscosity]
    type = FVDiffusion
    variable = u
    coeff = ${mu}
  []
  [u_pressure]
    type = INSFVMomentumPressure
    variable = u
    momentum_component = 'x'
    pressure = pressure
  []

  [v_time]
    type = INSFVMomentumTimeDerivative
    variable = v
    rho = ${rho}
  []
  [v_advection]
    type = INSFVMomentumAdvection
    variable = v
    advected_quantity = 'rhov'
    vel = ${vel}
    velocity_interp_method = ${velocity_interp_method}
    advected_interp_method = ${advected_interp_method}
    pressure = pressure
    u = u
    v = v
    mu = 'mu'
    rho = ${rho}
  []
  [v_viscosity]
    type = FVDiffusion
    variable = v
    coeff = ${mu}
  []
  [v_pressure]
    type = INSFVMomentumPressure
    variable = v
    momentum_component = 'y'
    pressure = pressure
  []
  [v_buoyancy]
    type = INSFVMomentumBoussinesq
    variable = v
    T_fluid = T
    gravity = '0 -9.82 0'
    rho = ${rho}
    ref_temperature = ${temp_ref}
    momentum_component = 'y'
  []
  [v_gravity]
    type = INSFVMomentumGravity
    variable = v
    gravity = '0 -9.82 0'
    rho = ${rho}
    momentum_component = 'y'
  []

  [temp_time]
    type = INSFVEnergyTimeDerivative
    variable = T
    rho = ${rho}
    cp_name = 'cp'
  []
  [temp_conduction]
    type = FVDiffusion
    coeff = 'k'
    variable = T
  []
  [temp_advection]
    type = INSFVEnergyAdvection
    variable = T
    vel = ${vel}
    velocity_interp_method = ${velocity_interp_method}
    advected_interp_method = ${advected_interp_method}
    pressure = pressure
    u = u
    v = v
    mu = 'mu'
    rho = ${rho}
  []
[]

[FVBCs]
  [free_slip_x]
    type = INSFVNaturalFreeSlipBC
    variable = u
    boundary = 'walls_top walls_right walls_bottom walls_left'
  []
  [free_slip_y]
    type = INSFVNaturalFreeSlipBC
    variable = v
    boundary = 'walls_top walls_right walls_bottom walls_left'
  []

  [T_hot]
    type = FVDirichletBC
    variable = 'T'
    boundary = 'walls_left'
    value = ${hot_temp}
  []

  [T_cold]
    type = FVDirichletBC
    variable = 'T'
    boundary = 'walls_right'
    value = ${cold_temp}
  []
[]

[Functions]
  [mu_rampdown]
    type = PiecewiseLinear
    x = '0.1 1 2 3'
    y = '${fparse 1e3*mu} ${fparse 1e2*mu} ${fparse 1e1*mu} ${mu}'
  []
[]

[Postprocessors]
  # To watch the rampdown
  [mu]
    type = FunctionValuePostprocessor
    function = mu_rampdown
  []
  # To evaluate convergence (visual inspection still necessary!)
  [average_temperature]
    type = ElementAverageValue
    variable = T
  []
[]

[Materials]
  [const_functor]
    type = ADGenericConstantFunctorMaterial
    prop_names = 'cp k alpha_b'
    prop_values = '${cp} ${k} ${alpha}'
  []
  [mu_functor]
    type = ADGenericFunctionFunctorMaterial
    prop_names = 'mu'
    prop_values = 'mu_rampdown'
  []
  [ins_fv]
    type = INSFVMaterial
    u = 'u'
    v = 'v'
    pressure = 'pressure'
    temperature = 'T'
    rho = ${rho}
  []
[]
[Executioner]
  type = Transient
  dt = 0.01
  dtmax = 1
  # 1s after the end of the mu_rampdown
  # going further could be required to relax the system to steady state further
  end_time = 1
  timestep_tolerance = 1e-3
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_pc_type -sub_pc_factor_shift_type'
  petsc_options_value = 'asm      1200                lu           NONZERO'
  nl_abs_tol = 1e-3 nl_rel_tol = 1e-3
  line_search = 'none'
  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
[]
[Debug]
  show_var_residual_norms = true
[]
[Outputs]
  exodus = true
[]
