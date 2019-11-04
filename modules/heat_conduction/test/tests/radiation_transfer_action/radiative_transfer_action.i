[Problem]
  kernel_coverage_check = false
[]

[Mesh]
  type = MeshGeneratorMesh

  [./cmg]
    type = CartesianMeshGenerator
    dim = 2
    dx = '1 1.3 1.9'
    ix = '3 3 3'
    dy = '2 1.2 0.9'
    iy = '3 3 3'
    subdomain_id = '0 1 0
                    4 5 2
                    0 3 0'
  [../]

  [./inner_bottom]
    type = SideSetsBetweenSubdomainsGenerator
    input = cmg
    master_block = 1
    paired_block = 5
    new_boundary = 'inner_bottom'
  [../]

  [./inner_left]
    type = SideSetsBetweenSubdomainsGenerator
    input = inner_bottom
    master_block = 4
    paired_block = 5
    new_boundary = 'inner_left'
  [../]

  [./inner_right]
    type = SideSetsBetweenSubdomainsGenerator
    input = inner_left
    master_block = 2
    paired_block = 5
    new_boundary = 'inner_right'
  [../]

  [./inner_top]
    type = SideSetsBetweenSubdomainsGenerator
    input = inner_right
    master_block = 3
    paired_block = 5
    new_boundary = 'inner_top'
  [../]

  [./rename]
    type = RenameBlockGenerator
    old_block_id = '1 2 3 4'
    new_block_id = '0 0 0 0'
    input = inner_top
  [../]
[]

[Variables]
  [./temperature]
    block = 0
  [../]
[]

[Kernels]
  [./heat_conduction]
    type = HeatConduction
    variable = temperature
    block = 0
    diffusion_coefficient = 5
  [../]
[]

[GrayDiffuseRadiation]
  [./cavity]
    sidesets = '4 5 6 7'
    emissivity = '0.9 0.8 0.4 1'
    n_patches = '2 2 2 3'
    #partitioners = 'metis metis metis metis'
    final_mesh_generator = rename
    temperature = temperature
    adiabatic_sidesets = '7'
    fixed_temperature_sidesets = '4'
    fixed_boundary_temperatures = '1200'
  [../]
[]

[BCs]
  [./left]
    type = DirichletBC
    variable = temperature
    boundary = left
    value = 600
  [../]

  [./right]
    type = DirichletBC
    variable = temperature
    boundary = right
    value = 300
  [../]
[]

[Postprocessors]
  [./average_T_inner_right]
    type = SideAverageValue
    variable = temperature
    boundary = inner_right
  [../]
[]

[Executioner]
  type = Steady
[]

[Outputs]
  exodus = true
[]