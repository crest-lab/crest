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

if get_option('man_pages')
  asciidoctor = find_program('asciidoctor', required: false)
  if not asciidoctor.found()
    error('Could not find asciidoctor.')
  endif
  foreach man_name : ['crest']
    man_page = custom_target(
      'man_@0@'.format(man_name),
      input: '@0@.adoc'.format(man_name),
      output: '@0@.1'.format(man_name),
      command: [
        asciidoctor,
	'--attribute=release-version=@0@'.format(meson.project_version()),
	'--backend=manpage',
	'--doctype=manpage',
	'--out-file=@OUTPUT@',
	'@INPUT@',
      ],
    )
  endforeach
endif
