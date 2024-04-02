list(APPEND copy_dir_examples_files
    examples/cbm18_dens8.grid_nx68ny64_profiles.nc
    )

foreach(fn ${copy_dir_examples_files})
   list(APPEND copy_dir_examples_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_examples_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_examples_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/)
list(APPEND copy_dir_examples_1D-conduction_files
    examples/1D-conduction/BOUT.inp
    examples/1D-conduction/README.md
    )

foreach(fn ${copy_dir_examples_1D-conduction_files})
   list(APPEND copy_dir_examples_1D-conduction_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_examples_1D-conduction_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/1D-conduction/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_examples_1D-conduction_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/1D-conduction/)
list(APPEND copy_dir_examples_1D-hydrogen_files
    examples/1D-hydrogen/BOUT.inp
    examples/1D-hydrogen/README.md
    examples/1D-hydrogen/plot_convergence.py
    )

foreach(fn ${copy_dir_examples_1D-hydrogen_files})
   list(APPEND copy_dir_examples_1D-hydrogen_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_examples_1D-hydrogen_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/1D-hydrogen/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_examples_1D-hydrogen_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/1D-hydrogen/)
list(APPEND copy_dir_examples_1D-hydrogen-equalT_files
    examples/1D-hydrogen-equalT/BOUT.inp
    )

foreach(fn ${copy_dir_examples_1D-hydrogen-equalT_files})
   list(APPEND copy_dir_examples_1D-hydrogen-equalT_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_examples_1D-hydrogen-equalT_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/1D-hydrogen-equalT/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_examples_1D-hydrogen-equalT_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/1D-hydrogen-equalT/)
list(APPEND copy_dir_examples_1D-hydrogen_qn_files
    examples/1D-hydrogen/qn/BOUT.inp
    )

foreach(fn ${copy_dir_examples_1D-hydrogen_qn_files})
   list(APPEND copy_dir_examples_1D-hydrogen_qn_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_examples_1D-hydrogen_qn_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/1D-hydrogen/qn/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_examples_1D-hydrogen_qn_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/1D-hydrogen/qn/)
list(APPEND copy_dir_examples_1D-neon_files
    examples/1D-neon/BOUT.inp
    examples/1D-neon/README.md
    examples/1D-neon/add_vars_to_restarts.py
    examples/1D-neon/makeplot.py
    )

foreach(fn ${copy_dir_examples_1D-neon_files})
   list(APPEND copy_dir_examples_1D-neon_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_examples_1D-neon_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/1D-neon/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_examples_1D-neon_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/1D-neon/)
list(APPEND copy_dir_examples_1D-neon-source_files
    examples/1D-neon-source/BOUT.inp
    )

foreach(fn ${copy_dir_examples_1D-neon-source_files})
   list(APPEND copy_dir_examples_1D-neon-source_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_examples_1D-neon-source_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/1D-neon-source/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_examples_1D-neon-source_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/1D-neon-source/)
list(APPEND copy_dir_examples_1D-periodic_files
    examples/1D-periodic/BOUT.inp
    examples/1D-periodic/README.md
    examples/1D-periodic/makeplot.py
    examples/1D-periodic/plotstuff.py
    )

foreach(fn ${copy_dir_examples_1D-periodic_files})
   list(APPEND copy_dir_examples_1D-periodic_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_examples_1D-periodic_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/1D-periodic/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_examples_1D-periodic_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/1D-periodic/)
list(APPEND copy_dir_examples_1D-recycling_files
    examples/1D-recycling/1d_recycling.gif
    examples/1D-recycling/1d_recycling_20.png
    examples/1D-recycling/BOUT.inp
    examples/1D-recycling/README.md
    examples/1D-recycling/makeplot.py
    examples/1D-recycling/plot_convergence.py
    examples/1D-recycling/plotstuff.py
    )

foreach(fn ${copy_dir_examples_1D-recycling_files})
   list(APPEND copy_dir_examples_1D-recycling_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_examples_1D-recycling_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/1D-recycling/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_examples_1D-recycling_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/1D-recycling/)
list(APPEND copy_dir_examples_1D-recycling_cvode_files
    examples/1D-recycling/cvode/BOUT.inp
    )

foreach(fn ${copy_dir_examples_1D-recycling_cvode_files})
   list(APPEND copy_dir_examples_1D-recycling_cvode_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_examples_1D-recycling_cvode_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/1D-recycling/cvode/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_examples_1D-recycling_cvode_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/1D-recycling/cvode/)
list(APPEND copy_dir_examples_1D-recycling-with-Tt-control_files
    examples/1D-recycling-with-Tt-control/BOUT.inp
    examples/1D-recycling-with-Tt-control/README.md
    examples/1D-recycling-with-Tt-control/control.png
    examples/1D-recycling-with-Tt-control/plot_control.py
    )

foreach(fn ${copy_dir_examples_1D-recycling-with-Tt-control_files})
   list(APPEND copy_dir_examples_1D-recycling-with-Tt-control_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_examples_1D-recycling-with-Tt-control_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/1D-recycling-with-Tt-control/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_examples_1D-recycling-with-Tt-control_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/1D-recycling-with-Tt-control/)
list(APPEND copy_dir_examples_1D-sheath_files
    examples/1D-sheath/BOUT.inp
    examples/1D-sheath/README.md
    examples/1D-sheath/plotstuff.py
    )

foreach(fn ${copy_dir_examples_1D-sheath_files})
   list(APPEND copy_dir_examples_1D-sheath_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_examples_1D-sheath_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/1D-sheath/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_examples_1D-sheath_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/1D-sheath/)
list(APPEND copy_dir_examples_1D-sheath-conduction_files
    examples/1D-sheath-conduction/BOUT.inp
    examples/1D-sheath-conduction/README.md
    examples/1D-sheath-conduction/plotstuff.py
    )

foreach(fn ${copy_dir_examples_1D-sheath-conduction_files})
   list(APPEND copy_dir_examples_1D-sheath-conduction_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_examples_1D-sheath-conduction_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/1D-sheath-conduction/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_examples_1D-sheath-conduction_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/1D-sheath-conduction/)
list(APPEND copy_dir_examples_1D-te-ti_files
    examples/1D-te-ti/1d_te_ti.gif
    examples/1D-te-ti/1d_te_ti.png
    examples/1D-te-ti/BOUT.inp
    examples/1D-te-ti/README.md
    examples/1D-te-ti/makeplot.py
    examples/1D-te-ti/plotstuff.py
    )

foreach(fn ${copy_dir_examples_1D-te-ti_files})
   list(APPEND copy_dir_examples_1D-te-ti_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_examples_1D-te-ti_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/1D-te-ti/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_examples_1D-te-ti_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/1D-te-ti/)
list(APPEND copy_dir_examples_2D-drift-plane-turbulence-te-ti_files
    examples/2D-drift-plane-turbulence-te-ti/BOUT.inp
    examples/2D-drift-plane-turbulence-te-ti/plotstuff.py
    )

foreach(fn ${copy_dir_examples_2D-drift-plane-turbulence-te-ti_files})
   list(APPEND copy_dir_examples_2D-drift-plane-turbulence-te-ti_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_examples_2D-drift-plane-turbulence-te-ti_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/2D-drift-plane-turbulence-te-ti/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_examples_2D-drift-plane-turbulence-te-ti_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/2D-drift-plane-turbulence-te-ti/)
list(APPEND copy_dir_examples_blob2d_files
    examples/blob2d/BOUT.inp
    examples/blob2d/README.md
    examples/blob2d/blob2d.png
    examples/blob2d/blob_size_scan.py
    examples/blob2d/blob_velocity.py
    examples/blob2d/makeplot.py
    examples/blob2d/plotstuff.py
    )

foreach(fn ${copy_dir_examples_blob2d_files})
   list(APPEND copy_dir_examples_blob2d_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_examples_blob2d_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/blob2d/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_examples_blob2d_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/blob2d/)
list(APPEND copy_dir_examples_blob2d-te_files
    examples/blob2d-te/BOUT.inp
    examples/blob2d-te/plotstuff.py
    )

foreach(fn ${copy_dir_examples_blob2d-te_files})
   list(APPEND copy_dir_examples_blob2d-te_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_examples_blob2d-te_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/blob2d-te/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_examples_blob2d-te_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/blob2d-te/)
list(APPEND copy_dir_examples_blob2d-te-ti_files
    examples/blob2d-te-ti/BOUT.inp
    examples/blob2d-te-ti/README.md
    examples/blob2d-te-ti/blob2d-te-ti.png
    examples/blob2d-te-ti/makeplot.py
    examples/blob2d-te-ti/plotstuff.py
    )

foreach(fn ${copy_dir_examples_blob2d-te-ti_files})
   list(APPEND copy_dir_examples_blob2d-te-ti_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_examples_blob2d-te-ti_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/blob2d-te-ti/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_examples_blob2d-te-ti_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/blob2d-te-ti/)
list(APPEND copy_dir_examples_blob2d-vpol_files
    examples/blob2d-vpol/BOUT.inp
    examples/blob2d-vpol/README.md
    examples/blob2d-vpol/blob2d.png
    examples/blob2d-vpol/makeplot.py
    )

foreach(fn ${copy_dir_examples_blob2d-vpol_files})
   list(APPEND copy_dir_examples_blob2d-vpol_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_examples_blob2d-vpol_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/blob2d-vpol/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_examples_blob2d-vpol_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/blob2d-vpol/)
list(APPEND copy_dir_examples_equilibriate_files
    examples/equilibriate/BOUT.inp
    examples/equilibriate/equilibriate.png
    examples/equilibriate/makeplot.py
    )

foreach(fn ${copy_dir_examples_equilibriate_files})
   list(APPEND copy_dir_examples_equilibriate_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_examples_equilibriate_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/equilibriate/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_examples_equilibriate_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/equilibriate/)
list(APPEND copy_dir_examples_linear_files
    examples/linear/README.md
    examples/linear/compare_profiles.py
    examples/linear/compare_psd.py
    examples/linear/compare_timeseries.py
    examples/linear/plot_crossection.py
    examples/linear/plot_psd.py
    examples/linear/plot_timeseries.py
    examples/linear/psd-comparison-loglin.png
    )

foreach(fn ${copy_dir_examples_linear_files})
   list(APPEND copy_dir_examples_linear_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_examples_linear_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/linear/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_examples_linear_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/linear/)
list(APPEND copy_dir_examples_linear_annulus-isothermal-d_files
    examples/linear/annulus-isothermal-d/BOUT.inp
    )

foreach(fn ${copy_dir_examples_linear_annulus-isothermal-d_files})
   list(APPEND copy_dir_examples_linear_annulus-isothermal-d_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_examples_linear_annulus-isothermal-d_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/linear/annulus-isothermal-d/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_examples_linear_annulus-isothermal-d_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/linear/annulus-isothermal-d/)
list(APPEND copy_dir_examples_linear_annulus-isothermal-d-2_files
    examples/linear/annulus-isothermal-d-2/BOUT.inp
    examples/linear/annulus-isothermal-d-2/Ne_2_8_0_timeseries.png
    examples/linear/annulus-isothermal-d-2/Ne_crossection.png
    )

foreach(fn ${copy_dir_examples_linear_annulus-isothermal-d-2_files})
   list(APPEND copy_dir_examples_linear_annulus-isothermal-d-2_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_examples_linear_annulus-isothermal-d-2_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/linear/annulus-isothermal-d-2/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_examples_linear_annulus-isothermal-d-2_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/linear/annulus-isothermal-d-2/)
list(APPEND copy_dir_examples_linear_annulus-isothermal-he_files
    examples/linear/annulus-isothermal-he/BOUT.inp
    examples/linear/annulus-isothermal-he/Ne_2_8_0_timeseries.png
    examples/linear/annulus-isothermal-he/Ne_crossection.png
    examples/linear/annulus-isothermal-he/README.md
    )

foreach(fn ${copy_dir_examples_linear_annulus-isothermal-he_files})
   list(APPEND copy_dir_examples_linear_annulus-isothermal-he_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_examples_linear_annulus-isothermal-he_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/linear/annulus-isothermal-he/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_examples_linear_annulus-isothermal-he_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/linear/annulus-isothermal-he/)
list(APPEND copy_dir_examples_linear_annulus-isothermal-he-emag_files
    examples/linear/annulus-isothermal-he-emag/BOUT.inp
    )

foreach(fn ${copy_dir_examples_linear_annulus-isothermal-he-emag_files})
   list(APPEND copy_dir_examples_linear_annulus-isothermal-he-emag_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_examples_linear_annulus-isothermal-he-emag_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/linear/annulus-isothermal-he-emag/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_examples_linear_annulus-isothermal-he-emag_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/linear/annulus-isothermal-he-emag/)
list(APPEND copy_dir_examples_linear_annulus-te-he-fixedneutrals_files
    examples/linear/annulus-te-he-fixedneutrals/BOUT.inp
    )

foreach(fn ${copy_dir_examples_linear_annulus-te-he-fixedneutrals_files})
   list(APPEND copy_dir_examples_linear_annulus-te-he-fixedneutrals_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_examples_linear_annulus-te-he-fixedneutrals_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/linear/annulus-te-he-fixedneutrals/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_examples_linear_annulus-te-he-fixedneutrals_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/linear/annulus-te-he-fixedneutrals/)
list(APPEND copy_dir_examples_solkit-comparison_files
    examples/solkit-comparison/makeplots.py
    )

foreach(fn ${copy_dir_examples_solkit-comparison_files})
   list(APPEND copy_dir_examples_solkit-comparison_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_examples_solkit-comparison_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/solkit-comparison/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_examples_solkit-comparison_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/solkit-comparison/)
list(APPEND copy_dir_examples_solkit-comparison_single-temperature_files
    examples/solkit-comparison/single-temperature/BOUT.inp
    )

foreach(fn ${copy_dir_examples_solkit-comparison_single-temperature_files})
   list(APPEND copy_dir_examples_solkit-comparison_single-temperature_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_examples_solkit-comparison_single-temperature_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/solkit-comparison/single-temperature/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_examples_solkit-comparison_single-temperature_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/solkit-comparison/single-temperature/)
list(APPEND copy_dir_examples_solkit-comparison_two-temperatures_files
    examples/solkit-comparison/two-temperatures/BOUT.inp
    )

foreach(fn ${copy_dir_examples_solkit-comparison_two-temperatures_files})
   list(APPEND copy_dir_examples_solkit-comparison_two-temperatures_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_examples_solkit-comparison_two-temperatures_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/solkit-comparison/two-temperatures/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_examples_solkit-comparison_two-temperatures_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/solkit-comparison/two-temperatures/)
list(APPEND copy_dir_examples_stellarator-2pt-model_files
    examples/stellarator-2pt-model/BOUT.inp
    )

foreach(fn ${copy_dir_examples_stellarator-2pt-model_files})
   list(APPEND copy_dir_examples_stellarator-2pt-model_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_examples_stellarator-2pt-model_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/stellarator-2pt-model/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_examples_stellarator-2pt-model_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/stellarator-2pt-model/)
list(APPEND copy_dir_examples_tcv-x21_files
    examples/tcv-x21/convert_to_tcvx21.py
    examples/tcv-x21/gather_data.py
    examples/tcv-x21/make_tcvx21_plots.py
    )

foreach(fn ${copy_dir_examples_tcv-x21_files})
   list(APPEND copy_dir_examples_tcv-x21_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_examples_tcv-x21_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/tcv-x21/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_examples_tcv-x21_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/tcv-x21/)
list(APPEND copy_dir_examples_tcv-x21_data_files
    examples/tcv-x21/data/BOUT.inp
    )

foreach(fn ${copy_dir_examples_tcv-x21_data_files})
   list(APPEND copy_dir_examples_tcv-x21_data_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_examples_tcv-x21_data_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/tcv-x21/data/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_examples_tcv-x21_data_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/tcv-x21/data/)
list(APPEND copy_dir_examples_tokamak_files
    examples/tokamak/README.md
    examples/tokamak/adjust_curvature.py
    examples/tokamak/compass-36x48.grd.nc
    examples/tokamak/tokamak.nc
    )

foreach(fn ${copy_dir_examples_tokamak_files})
   list(APPEND copy_dir_examples_tokamak_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_examples_tokamak_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/tokamak/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_examples_tokamak_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/tokamak/)
list(APPEND copy_dir_examples_tokamak_diffusion_files
    examples/tokamak/diffusion/BOUT.inp
    examples/tokamak/diffusion/README.md
    examples/tokamak/diffusion/plotstuff.py
    )

foreach(fn ${copy_dir_examples_tokamak_diffusion_files})
   list(APPEND copy_dir_examples_tokamak_diffusion_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_examples_tokamak_diffusion_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/tokamak/diffusion/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_examples_tokamak_diffusion_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/tokamak/diffusion/)
list(APPEND copy_dir_examples_tokamak_diffusion-conduction_files
    examples/tokamak/diffusion-conduction/BOUT.inp
    examples/tokamak/diffusion-conduction/README.md
    examples/tokamak/diffusion-conduction/plotstuff.py
    )

foreach(fn ${copy_dir_examples_tokamak_diffusion-conduction_files})
   list(APPEND copy_dir_examples_tokamak_diffusion-conduction_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_examples_tokamak_diffusion-conduction_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/tokamak/diffusion-conduction/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_examples_tokamak_diffusion-conduction_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/tokamak/diffusion-conduction/)
list(APPEND copy_dir_examples_tokamak_diffusion-flow-evolveT_files
    examples/tokamak/diffusion-flow-evolveT/BOUT.inp
    examples/tokamak/diffusion-flow-evolveT/README.md
    examples/tokamak/diffusion-flow-evolveT/plotstuff.py
    )

foreach(fn ${copy_dir_examples_tokamak_diffusion-flow-evolveT_files})
   list(APPEND copy_dir_examples_tokamak_diffusion-flow-evolveT_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_examples_tokamak_diffusion-flow-evolveT_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/tokamak/diffusion-flow-evolveT/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_examples_tokamak_diffusion-flow-evolveT_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/tokamak/diffusion-flow-evolveT/)
list(APPEND copy_dir_examples_tokamak_diffusion-transport_files
    examples/tokamak/diffusion-transport/BOUT.inp
    examples/tokamak/diffusion-transport/README.md
    examples/tokamak/diffusion-transport/plotstuff.py
    )

foreach(fn ${copy_dir_examples_tokamak_diffusion-transport_files})
   list(APPEND copy_dir_examples_tokamak_diffusion-transport_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_examples_tokamak_diffusion-transport_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/tokamak/diffusion-transport/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_examples_tokamak_diffusion-transport_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/tokamak/diffusion-transport/)
list(APPEND copy_dir_examples_tokamak_heat-transport_files
    examples/tokamak/heat-transport/BOUT.inp
    examples/tokamak/heat-transport/README.md
    examples/tokamak/heat-transport/analyse_fluxes.py
    )

foreach(fn ${copy_dir_examples_tokamak_heat-transport_files})
   list(APPEND copy_dir_examples_tokamak_heat-transport_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_examples_tokamak_heat-transport_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/tokamak/heat-transport/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_examples_tokamak_heat-transport_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/tokamak/heat-transport/)
list(APPEND copy_dir_examples_tokamak_isothermal_files
    examples/tokamak/isothermal/BOUT.inp
    examples/tokamak/isothermal/makeplots.py
    )

foreach(fn ${copy_dir_examples_tokamak_isothermal_files})
   list(APPEND copy_dir_examples_tokamak_isothermal_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_examples_tokamak_isothermal_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/tokamak/isothermal/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_examples_tokamak_isothermal_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/tokamak/isothermal/)
list(APPEND copy_dir_examples_tokamak_linear-transport_files
    examples/tokamak/linear-transport/BOUT.inp
    examples/tokamak/linear-transport/plotstuff.py
    )

foreach(fn ${copy_dir_examples_tokamak_linear-transport_files})
   list(APPEND copy_dir_examples_tokamak_linear-transport_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_examples_tokamak_linear-transport_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/tokamak/linear-transport/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_examples_tokamak_linear-transport_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/tokamak/linear-transport/)
list(APPEND copy_dir_examples_tokamak_recycling_files
    examples/tokamak/recycling/BOUT.inp
    examples/tokamak/recycling/README.md
    examples/tokamak/recycling/analyse_fluxes.py
    examples/tokamak/recycling/plot_ddts.py
    examples/tokamak/recycling/plot_history.py
    examples/tokamak/recycling/plotstuff.py
    )

foreach(fn ${copy_dir_examples_tokamak_recycling_files})
   list(APPEND copy_dir_examples_tokamak_recycling_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_examples_tokamak_recycling_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/tokamak/recycling/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_examples_tokamak_recycling_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/tokamak/recycling/)
list(APPEND copy_dir_examples_tokamak_recycling-dthe_files
    examples/tokamak/recycling-dthe/BOUT.inp
    examples/tokamak/recycling-dthe/README.md
    examples/tokamak/recycling-dthe/makeplots.py
    examples/tokamak/recycling-dthe/plotstuff.py
    )

foreach(fn ${copy_dir_examples_tokamak_recycling-dthe_files})
   list(APPEND copy_dir_examples_tokamak_recycling-dthe_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_examples_tokamak_recycling-dthe_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/tokamak/recycling-dthe/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_examples_tokamak_recycling-dthe_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/tokamak/recycling-dthe/)
list(APPEND copy_dir_examples_tokamak_recycling-dthe-drifts_files
    examples/tokamak/recycling-dthe-drifts/BOUT.inp
    examples/tokamak/recycling-dthe-drifts/plotstuff.py
    )

foreach(fn ${copy_dir_examples_tokamak_recycling-dthe-drifts_files})
   list(APPEND copy_dir_examples_tokamak_recycling-dthe-drifts_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_examples_tokamak_recycling-dthe-drifts_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/tokamak/recycling-dthe-drifts/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_examples_tokamak_recycling-dthe-drifts_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/tokamak/recycling-dthe-drifts/)
list(APPEND copy_dir_examples_tokamak_recycling-dthene_files
    examples/tokamak/recycling-dthene/BOUT.inp
    examples/tokamak/recycling-dthene/README.md
    examples/tokamak/recycling-dthene/pe_nvt_nne_2d.png
    examples/tokamak/recycling-dthene/plotstuff.py
    )

foreach(fn ${copy_dir_examples_tokamak_recycling-dthene_files})
   list(APPEND copy_dir_examples_tokamak_recycling-dthene_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_examples_tokamak_recycling-dthene_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/tokamak/recycling-dthene/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_examples_tokamak_recycling-dthene_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/tokamak/recycling-dthene/)
list(APPEND copy_dir_examples_tokamak_turbulence_files
    examples/tokamak/turbulence/BOUT.inp
    )

foreach(fn ${copy_dir_examples_tokamak_turbulence_files})
   list(APPEND copy_dir_examples_tokamak_turbulence_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_examples_tokamak_turbulence_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/tokamak/turbulence/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_examples_tokamak_turbulence_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples/tokamak/turbulence/)
list(APPEND copy_dir_json_database_files
    json_database/acd96_c.json
    json_database/acd96_li.json
    json_database/acd96_n.json
    json_database/acd96_ne.json
    json_database/ccd89_li.json
    json_database/ccd89_ne.json
    json_database/ccd96_c.json
    json_database/ecd89_li.json
    json_database/ecd96_c.json
    json_database/ecd96_li.json
    json_database/ecd96_n.json
    json_database/ecd96_ne.json
    json_database/plt96_c.json
    json_database/plt96_li.json
    json_database/plt96_n.json
    json_database/plt96_ne.json
    json_database/prb96_c.json
    json_database/prb96_li.json
    json_database/prb96_n.json
    json_database/prb96_ne.json
    json_database/prc89_li.json
    json_database/prc96_c.json
    json_database/scd96_c.json
    json_database/scd96_li.json
    json_database/scd96_n.json
    json_database/scd96_ne.json
    )

foreach(fn ${copy_dir_json_database_files})
   list(APPEND copy_dir_json_database_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_json_database_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/json_database/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_json_database_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/json_database/)
list(APPEND copy_dir_tests_integrated_files
    tests/integrated/test_suite
    )

foreach(fn ${copy_dir_tests_integrated_files})
   list(APPEND copy_dir_tests_integrated_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_tests_integrated_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_tests_integrated_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/)
list(APPEND copy_dir_tests_integrated_1D-fluid_files
    tests/integrated/1D-fluid/README.md
    tests/integrated/1D-fluid/fluid_norm.png
    tests/integrated/1D-fluid/mms.py
    tests/integrated/1D-fluid/runtest
    )

foreach(fn ${copy_dir_tests_integrated_1D-fluid_files})
   list(APPEND copy_dir_tests_integrated_1D-fluid_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_tests_integrated_1D-fluid_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/1D-fluid/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_tests_integrated_1D-fluid_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/1D-fluid/)
list(APPEND copy_dir_tests_integrated_1D-fluid_data_files
    tests/integrated/1D-fluid/data/BOUT.inp
    )

foreach(fn ${copy_dir_tests_integrated_1D-fluid_data_files})
   list(APPEND copy_dir_tests_integrated_1D-fluid_data_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_tests_integrated_1D-fluid_data_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/1D-fluid/data/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_tests_integrated_1D-fluid_data_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/1D-fluid/data/)
list(APPEND copy_dir_tests_integrated_1D-recycling_files
    tests/integrated/1D-recycling/README.md
    tests/integrated/1D-recycling/runtest
    )

foreach(fn ${copy_dir_tests_integrated_1D-recycling_files})
   list(APPEND copy_dir_tests_integrated_1D-recycling_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_tests_integrated_1D-recycling_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/1D-recycling/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_tests_integrated_1D-recycling_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/1D-recycling/)
list(APPEND copy_dir_tests_integrated_1D-recycling_data_files
    tests/integrated/1D-recycling/data/BOUT.inp
    )

foreach(fn ${copy_dir_tests_integrated_1D-recycling_data_files})
   list(APPEND copy_dir_tests_integrated_1D-recycling_data_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_tests_integrated_1D-recycling_data_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/1D-recycling/data/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_tests_integrated_1D-recycling_data_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/1D-recycling/data/)
list(APPEND copy_dir_tests_integrated_diffusion_files
    tests/integrated/diffusion/README.md
    tests/integrated/diffusion/runtest
    )

foreach(fn ${copy_dir_tests_integrated_diffusion_files})
   list(APPEND copy_dir_tests_integrated_diffusion_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_tests_integrated_diffusion_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/diffusion/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_tests_integrated_diffusion_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/diffusion/)
list(APPEND copy_dir_tests_integrated_diffusion_data_files
    tests/integrated/diffusion/data/BOUT.inp
    )

foreach(fn ${copy_dir_tests_integrated_diffusion_data_files})
   list(APPEND copy_dir_tests_integrated_diffusion_data_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_tests_integrated_diffusion_data_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/diffusion/data/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_tests_integrated_diffusion_data_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/diffusion/data/)
list(APPEND copy_dir_tests_integrated_drift-wave_files
    tests/integrated/drift-wave/README.md
    tests/integrated/drift-wave/analysis.py
    tests/integrated/drift-wave/drift-wave.png
    tests/integrated/drift-wave/runtest
    )

foreach(fn ${copy_dir_tests_integrated_drift-wave_files})
   list(APPEND copy_dir_tests_integrated_drift-wave_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_tests_integrated_drift-wave_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/drift-wave/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_tests_integrated_drift-wave_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/drift-wave/)
list(APPEND copy_dir_tests_integrated_drift-wave_data_files
    tests/integrated/drift-wave/data/BOUT.inp
    )

foreach(fn ${copy_dir_tests_integrated_drift-wave_data_files})
   list(APPEND copy_dir_tests_integrated_drift-wave_data_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_tests_integrated_drift-wave_data_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/drift-wave/data/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_tests_integrated_drift-wave_data_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/drift-wave/data/)
list(APPEND copy_dir_tests_integrated_evolve_density_files
    tests/integrated/evolve_density/runtest
    )

foreach(fn ${copy_dir_tests_integrated_evolve_density_files})
   list(APPEND copy_dir_tests_integrated_evolve_density_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_tests_integrated_evolve_density_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/evolve_density/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_tests_integrated_evolve_density_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/evolve_density/)
list(APPEND copy_dir_tests_integrated_evolve_density_data_files
    tests/integrated/evolve_density/data/BOUT.inp
    )

foreach(fn ${copy_dir_tests_integrated_evolve_density_data_files})
   list(APPEND copy_dir_tests_integrated_evolve_density_data_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_tests_integrated_evolve_density_data_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/evolve_density/data/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_tests_integrated_evolve_density_data_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/evolve_density/data/)
list(APPEND copy_dir_tests_integrated_neutral_mixed_files
    tests/integrated/neutral_mixed/runtest
    )

foreach(fn ${copy_dir_tests_integrated_neutral_mixed_files})
   list(APPEND copy_dir_tests_integrated_neutral_mixed_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_tests_integrated_neutral_mixed_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/neutral_mixed/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_tests_integrated_neutral_mixed_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/neutral_mixed/)
list(APPEND copy_dir_tests_integrated_neutral_mixed_data_files
    tests/integrated/neutral_mixed/data/BOUT.inp
    )

foreach(fn ${copy_dir_tests_integrated_neutral_mixed_data_files})
   list(APPEND copy_dir_tests_integrated_neutral_mixed_data_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_tests_integrated_neutral_mixed_data_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/neutral_mixed/data/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_tests_integrated_neutral_mixed_data_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/neutral_mixed/data/)
list(APPEND copy_dir_tests_integrated_neutral_parallel_diffusion_files
    tests/integrated/neutral_parallel_diffusion/runtest
    )

foreach(fn ${copy_dir_tests_integrated_neutral_parallel_diffusion_files})
   list(APPEND copy_dir_tests_integrated_neutral_parallel_diffusion_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_tests_integrated_neutral_parallel_diffusion_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/neutral_parallel_diffusion/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_tests_integrated_neutral_parallel_diffusion_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/neutral_parallel_diffusion/)
list(APPEND copy_dir_tests_integrated_neutral_parallel_diffusion_data_files
    tests/integrated/neutral_parallel_diffusion/data/BOUT.inp
    )

foreach(fn ${copy_dir_tests_integrated_neutral_parallel_diffusion_data_files})
   list(APPEND copy_dir_tests_integrated_neutral_parallel_diffusion_data_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_tests_integrated_neutral_parallel_diffusion_data_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/neutral_parallel_diffusion/data/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_tests_integrated_neutral_parallel_diffusion_data_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/neutral_parallel_diffusion/data/)
list(APPEND copy_dir_tests_integrated_snb_files
    tests/integrated/snb/README.md
    tests/integrated/snb/runtest
    )

foreach(fn ${copy_dir_tests_integrated_snb_files})
   list(APPEND copy_dir_tests_integrated_snb_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_tests_integrated_snb_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/snb/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_tests_integrated_snb_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/snb/)
list(APPEND copy_dir_tests_integrated_snb_nonuniform_files
    tests/integrated/snb/nonuniform/BOUT.inp
    )

foreach(fn ${copy_dir_tests_integrated_snb_nonuniform_files})
   list(APPEND copy_dir_tests_integrated_snb_nonuniform_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_tests_integrated_snb_nonuniform_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/snb/nonuniform/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_tests_integrated_snb_nonuniform_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/snb/nonuniform/)
list(APPEND copy_dir_tests_integrated_snb_nonuniform_area_files
    tests/integrated/snb/nonuniform_area/BOUT.inp
    )

foreach(fn ${copy_dir_tests_integrated_snb_nonuniform_area_files})
   list(APPEND copy_dir_tests_integrated_snb_nonuniform_area_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_tests_integrated_snb_nonuniform_area_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/snb/nonuniform_area/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_tests_integrated_snb_nonuniform_area_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/snb/nonuniform_area/)
list(APPEND copy_dir_tests_integrated_snb_uniform_files
    tests/integrated/snb/uniform/BOUT.inp
    )

foreach(fn ${copy_dir_tests_integrated_snb_uniform_files})
   list(APPEND copy_dir_tests_integrated_snb_uniform_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_tests_integrated_snb_uniform_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/snb/uniform/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_tests_integrated_snb_uniform_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/snb/uniform/)
list(APPEND copy_dir_tests_integrated_sod-shock_files
    tests/integrated/sod-shock/README.md
    tests/integrated/sod-shock/runtest
    tests/integrated/sod-shock/sod_shock.png
    )

foreach(fn ${copy_dir_tests_integrated_sod-shock_files})
   list(APPEND copy_dir_tests_integrated_sod-shock_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_tests_integrated_sod-shock_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/sod-shock/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_tests_integrated_sod-shock_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/sod-shock/)
list(APPEND copy_dir_tests_integrated_sod-shock_data_files
    tests/integrated/sod-shock/data/BOUT.inp
    )

foreach(fn ${copy_dir_tests_integrated_sod-shock_data_files})
   list(APPEND copy_dir_tests_integrated_sod-shock_data_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_tests_integrated_sod-shock_data_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/sod-shock/data/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_tests_integrated_sod-shock_data_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/sod-shock/data/)
list(APPEND copy_dir_tests_integrated_sod-shock-energy_files
    tests/integrated/sod-shock-energy/README.md
    tests/integrated/sod-shock-energy/runtest
    tests/integrated/sod-shock-energy/sod_shock_energy.png
    )

foreach(fn ${copy_dir_tests_integrated_sod-shock-energy_files})
   list(APPEND copy_dir_tests_integrated_sod-shock-energy_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_tests_integrated_sod-shock-energy_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/sod-shock-energy/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_tests_integrated_sod-shock-energy_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/sod-shock-energy/)
list(APPEND copy_dir_tests_integrated_sod-shock-energy_data_files
    tests/integrated/sod-shock-energy/data/BOUT.inp
    )

foreach(fn ${copy_dir_tests_integrated_sod-shock-energy_data_files})
   list(APPEND copy_dir_tests_integrated_sod-shock-energy_data_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_tests_integrated_sod-shock-energy_data_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/sod-shock-energy/data/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_tests_integrated_sod-shock-energy_data_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/sod-shock-energy/data/)
list(APPEND copy_dir_tests_integrated_toro-1_files
    tests/integrated/toro-1/README.md
    tests/integrated/toro-1/makeplot.py
    )

foreach(fn ${copy_dir_tests_integrated_toro-1_files})
   list(APPEND copy_dir_tests_integrated_toro-1_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_tests_integrated_toro-1_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/toro-1/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_tests_integrated_toro-1_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/toro-1/)
list(APPEND copy_dir_tests_integrated_toro-1_data_files
    tests/integrated/toro-1/data/BOUT.inp
    )

foreach(fn ${copy_dir_tests_integrated_toro-1_data_files})
   list(APPEND copy_dir_tests_integrated_toro-1_data_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_tests_integrated_toro-1_data_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/toro-1/data/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_tests_integrated_toro-1_data_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/toro-1/data/)
list(APPEND copy_dir_tests_integrated_toro-2_files
    tests/integrated/toro-2/README.md
    tests/integrated/toro-2/makeplot.py
    tests/integrated/toro-2/toro-2.png
    )

foreach(fn ${copy_dir_tests_integrated_toro-2_files})
   list(APPEND copy_dir_tests_integrated_toro-2_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_tests_integrated_toro-2_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/toro-2/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_tests_integrated_toro-2_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/toro-2/)
list(APPEND copy_dir_tests_integrated_toro-2_data_files
    tests/integrated/toro-2/data/BOUT.inp
    )

foreach(fn ${copy_dir_tests_integrated_toro-2_data_files})
   list(APPEND copy_dir_tests_integrated_toro-2_data_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_tests_integrated_toro-2_data_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/toro-2/data/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_tests_integrated_toro-2_data_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/toro-2/data/)
list(APPEND copy_dir_tests_integrated_toro-2-energy_files
    tests/integrated/toro-2-energy/README.md
    tests/integrated/toro-2-energy/makeplot.py
    )

foreach(fn ${copy_dir_tests_integrated_toro-2-energy_files})
   list(APPEND copy_dir_tests_integrated_toro-2-energy_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_tests_integrated_toro-2-energy_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/toro-2-energy/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_tests_integrated_toro-2-energy_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/toro-2-energy/)
list(APPEND copy_dir_tests_integrated_toro-2-energy_data_files
    tests/integrated/toro-2-energy/data/BOUT.inp
    )

foreach(fn ${copy_dir_tests_integrated_toro-2-energy_data_files})
   list(APPEND copy_dir_tests_integrated_toro-2-energy_data_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_tests_integrated_toro-2-energy_data_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/toro-2-energy/data/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_tests_integrated_toro-2-energy_data_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/toro-2-energy/data/)
list(APPEND copy_dir_tests_integrated_toro-3_files
    tests/integrated/toro-3/README.md
    tests/integrated/toro-3/runtest
    tests/integrated/toro-3/toro-3.png
    )

foreach(fn ${copy_dir_tests_integrated_toro-3_files})
   list(APPEND copy_dir_tests_integrated_toro-3_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_tests_integrated_toro-3_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/toro-3/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_tests_integrated_toro-3_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/toro-3/)
list(APPEND copy_dir_tests_integrated_toro-3_data_files
    tests/integrated/toro-3/data/BOUT.inp
    )

foreach(fn ${copy_dir_tests_integrated_toro-3_data_files})
   list(APPEND copy_dir_tests_integrated_toro-3_data_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_tests_integrated_toro-3_data_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/toro-3/data/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_tests_integrated_toro-3_data_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/toro-3/data/)
list(APPEND copy_dir_tests_integrated_toro-3-energy_files
    tests/integrated/toro-3-energy/README.md
    tests/integrated/toro-3-energy/runtest
    tests/integrated/toro-3-energy/toro-3-energy.png
    )

foreach(fn ${copy_dir_tests_integrated_toro-3-energy_files})
   list(APPEND copy_dir_tests_integrated_toro-3-energy_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_tests_integrated_toro-3-energy_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/toro-3-energy/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_tests_integrated_toro-3-energy_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/toro-3-energy/)
list(APPEND copy_dir_tests_integrated_toro-3-energy_data_files
    tests/integrated/toro-3-energy/data/BOUT.inp
    )

foreach(fn ${copy_dir_tests_integrated_toro-3-energy_data_files})
   list(APPEND copy_dir_tests_integrated_toro-3-energy_data_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_tests_integrated_toro-3-energy_data_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/toro-3-energy/data/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_tests_integrated_toro-3-energy_data_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/toro-3-energy/data/)
list(APPEND copy_dir_tests_integrated_toro-4_files
    tests/integrated/toro-4/README.md
    tests/integrated/toro-4/makeplot.py
    )

foreach(fn ${copy_dir_tests_integrated_toro-4_files})
   list(APPEND copy_dir_tests_integrated_toro-4_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_tests_integrated_toro-4_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/toro-4/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_tests_integrated_toro-4_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/toro-4/)
list(APPEND copy_dir_tests_integrated_toro-4_data_files
    tests/integrated/toro-4/data/BOUT.inp
    )

foreach(fn ${copy_dir_tests_integrated_toro-4_data_files})
   list(APPEND copy_dir_tests_integrated_toro-4_data_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_tests_integrated_toro-4_data_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/toro-4/data/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_tests_integrated_toro-4_data_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/toro-4/data/)
list(APPEND copy_dir_tests_integrated_toro-4-energy_files
    tests/integrated/toro-4-energy/README.md
    tests/integrated/toro-4-energy/makeplot.py
    )

foreach(fn ${copy_dir_tests_integrated_toro-4-energy_files})
   list(APPEND copy_dir_tests_integrated_toro-4-energy_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_tests_integrated_toro-4-energy_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/toro-4-energy/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_tests_integrated_toro-4-energy_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/toro-4-energy/)
list(APPEND copy_dir_tests_integrated_toro-4-energy_data_files
    tests/integrated/toro-4-energy/data/BOUT.inp
    )

foreach(fn ${copy_dir_tests_integrated_toro-4-energy_data_files})
   list(APPEND copy_dir_tests_integrated_toro-4-energy_data_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_tests_integrated_toro-4-energy_data_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/toro-4-energy/data/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_tests_integrated_toro-4-energy_data_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/toro-4-energy/data/)
list(APPEND copy_dir_tests_integrated_toro-5_files
    tests/integrated/toro-5/README.md
    tests/integrated/toro-5/runtest
    )

foreach(fn ${copy_dir_tests_integrated_toro-5_files})
   list(APPEND copy_dir_tests_integrated_toro-5_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_tests_integrated_toro-5_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/toro-5/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_tests_integrated_toro-5_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/toro-5/)
list(APPEND copy_dir_tests_integrated_toro-5_data_files
    tests/integrated/toro-5/data/BOUT.inp
    )

foreach(fn ${copy_dir_tests_integrated_toro-5_data_files})
   list(APPEND copy_dir_tests_integrated_toro-5_data_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_tests_integrated_toro-5_data_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/toro-5/data/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_tests_integrated_toro-5_data_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/toro-5/data/)
list(APPEND copy_dir_tests_integrated_toro-5-energy_files
    tests/integrated/toro-5-energy/README.md
    tests/integrated/toro-5-energy/runtest
    )

foreach(fn ${copy_dir_tests_integrated_toro-5-energy_files})
   list(APPEND copy_dir_tests_integrated_toro-5-energy_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_tests_integrated_toro-5-energy_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/toro-5-energy/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_tests_integrated_toro-5-energy_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/toro-5-energy/)
list(APPEND copy_dir_tests_integrated_toro-5-energy_data_files
    tests/integrated/toro-5-energy/data/BOUT.inp
    )

foreach(fn ${copy_dir_tests_integrated_toro-5-energy_data_files})
   list(APPEND copy_dir_tests_integrated_toro-5-energy_data_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_tests_integrated_toro-5-energy_data_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/toro-5-energy/data/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_tests_integrated_toro-5-energy_data_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/toro-5-energy/data/)
list(APPEND copy_dir_tests_integrated_vorticity_files
    tests/integrated/vorticity/runtest
    )

foreach(fn ${copy_dir_tests_integrated_vorticity_files})
   list(APPEND copy_dir_tests_integrated_vorticity_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_tests_integrated_vorticity_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/vorticity/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_tests_integrated_vorticity_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/vorticity/)
list(APPEND copy_dir_tests_integrated_vorticity_data_files
    tests/integrated/vorticity/data/BOUT.inp
    )

foreach(fn ${copy_dir_tests_integrated_vorticity_data_files})
   list(APPEND copy_dir_tests_integrated_vorticity_data_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_tests_integrated_vorticity_data_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/vorticity/data/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_tests_integrated_vorticity_data_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/integrated/vorticity/data/)
list(APPEND copy_dir_tests_unit_files
    tests/unit/hermes_test_main.cxx
    tests/unit/main.cxx
    tests/unit/test_amjuel_hyd_recombination.cxx
    tests/unit/test_anomalous_diffusion.cxx
    tests/unit/test_collisions.cxx
    tests/unit/test_component.cxx
    tests/unit/test_component_scheduler.cxx
    tests/unit/test_diamagnetic_drift.cxx
    tests/unit/test_electron_force_balance.cxx
    tests/unit/test_extras.cxx
    tests/unit/test_extras.hxx
    tests/unit/test_fixed_density.cxx
    tests/unit/test_fixed_fraction_ions.cxx
    tests/unit/test_fixed_velocity.cxx
    tests/unit/test_hydrogen_charge_exchange.cxx
    tests/unit/test_integrate.cxx
    tests/unit/test_ionisation.cxx
    tests/unit/test_isothermal.cxx
    tests/unit/test_noflow_boundary.cxx
    tests/unit/test_recycling.cxx
    tests/unit/test_sheath_boundary.cxx
    tests/unit/test_sheath_closure.cxx
    tests/unit/test_snb_conduction.cxx
    tests/unit/test_sound_speed.cxx
    tests/unit/test_zero_current.cxx
    )

foreach(fn ${copy_dir_tests_unit_files})
   list(APPEND copy_dir_tests_unit_files_abs ${CMAKE_SOURCE_DIR}/${fn})
endforeach()
message(INFO ${copy_dir_tests_unit_files_abs})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E make_directory
		   $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/unit/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_if_different
                   ${copy_dir_tests_unit_files_abs} $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests/unit/)
