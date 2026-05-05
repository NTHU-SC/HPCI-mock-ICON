# Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#

# ~~~
# autoconf_configure(INPUTS <inputs>,
#                    OUTPUTS <outputs>,
#                    VARIABLES <variables>,
#                    VALUES <values>)
# ~~~
# Iterates over the Autoconf templates from the <inputs> list, replaces the
# variable references from the <variables> list with the <values> (only the
# Autoconf-like @VAR@ forms) and stores the result to the respective files from
# the <outputs> list. The <inputs> and the <outputs> must have the same length.
# The <variables> and the <values> must have the same length. The function helps
# in keeping the scope of the main script free from variables needed for
# substitution only.
#
function(autoconf_configure)
  cmake_parse_arguments(
    PARSE_ARGV 0 ARG "" "" "INPUTS;OUTPUTS;VARIABLES;VALUES"
  )

  foreach(variable value IN ZIP_LISTS ARG_VARIABLES ARG_VALUES)
    set(${variable} "${value}")
  endforeach()

  foreach(input output IN ZIP_LISTS ARG_INPUTS ARG_OUTPUTS)
    configure_file(${input} ${output} @ONLY)
  endforeach()
endfunction()
