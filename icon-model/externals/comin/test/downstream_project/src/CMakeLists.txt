add_library(dummy_plugin MODULE dummy.c)
target_link_libraries(dummy_plugin PRIVATE ComIn::ComIn)
