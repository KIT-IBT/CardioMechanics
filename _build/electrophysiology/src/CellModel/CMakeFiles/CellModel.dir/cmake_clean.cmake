file(REMOVE_RECURSE
  "../../../lib/macosx/libCellModel.a"
  "../../../lib/macosx/libCellModel.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang CXX)
  include(CMakeFiles/CellModel.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
