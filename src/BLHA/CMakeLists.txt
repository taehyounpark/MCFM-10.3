if(${with_library})
  target_sources(libmcfm PRIVATE
    CXX_Interface.cxx
    Flavor_Map.cxx
    )
  if (DEFINED OLDIR)
     target_sources(oltest PRIVATE
       oltest.cxx
       )
     target_sources(oltest2 PRIVATE
       oltest2.cxx
       )
  endif()
  if (DEFINED RCLDIR)
     target_sources(rcltest PRIVATE
       rcltest.cxx
       )
     target_sources(rcltest2 PRIVATE
       rcltest2.cxx
       )
  endif()
  if (DEFINED RCLHEFTDIR)
     target_sources(rcltest_heft PRIVATE
       rcltest.cxx
       )
  endif()
  if (DEFINED MLDIR)
     target_sources(mltest PRIVATE
       CXX_Interface.cxx
       Flavor_Map.cxx
       mltest.cxx
       )
     target_sources(mltest2 PRIVATE
       CXX_Interface.cxx
       Flavor_Map.cxx
       mltest2.cxx
       )
  endif()
endif()
