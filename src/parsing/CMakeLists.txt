# This file is part of crest.
# SPDX-Identifier: LGPL-3.0-or-later
#
# crest is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# crest is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with crest.  If not, see <https://www.gnu.org/licenses/>.

set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND srcs
   "${dir}/confparse2.f90"
   "${dir}/parse_block.f90"
   "${dir}/parse_calcdata.f90"
   "${dir}/parse_maindata.f90"
   "${dir}/parse_datastruct.f90"
   "${dir}/parse_inputfile.f90"
   "${dir}/parse_keyvalue.f90"
)

set(srcs ${srcs} PARENT_SCOPE)






