#!/usr/bin/env bash
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

if [ -z "$MESON_DIST_ROOT" ]; then
    echo "Error: This script must be run via meson dist."
    exit 1
fi

echo "Setting up symlinks for nested subprojects in distribution tarball..."

cd "$MESON_DIST_ROOT/subprojects" || exit 1

if [ -d "test-drive" ] && [ -d "toml-f" ]; then
    echo "  -> Symlinking test-drive into toml-f/subprojects"
    mkdir -p toml-f/subprojects
    cd toml-f/subprojects
    ln -sfn ../../test-drive test-drive
    cd ../../
fi

if [ -d "mstore" ] && [ -d "dftd4" ]; then
    echo "  -> Symlinking mstore into dftd4/subprojects"
    mkdir -p dftd4/subprojects
    cd dftd4/subprojects
    ln -sfn ../../mstore mstore
    cd ../../
fi
if [ -d "multicharge" ] && [ -d "dftd4" ]; then
    echo "  -> Symlinking multicharge into dftd4/subprojects"
    mkdir -p dftd4/subprojects
    cd dftd4/subprojects
    ln -sfn ../../multicharge multicharge
    cd ../../
fi

if [ -d "dftd4" ] && [ -d "tblite" ]; then
    echo "  -> Symlinking dftd4 into tblite/subprojects"
    mkdir -p tblite/subprojects
    cd tblite/subprojects
    ln -sfn ../../dftd4 dftd4
    cd ../../
fi
if [ -d "mctc-lib" ] && [ -d "tblite" ]; then
    echo "  -> Symlinking mctc-lib into tblite/subprojects"
    mkdir -p tblite/subprojects
    cd tblite/subprojects
    ln -sfn ../../mctc-lib mctc-lib
    cd ../../
fi
if [ -d "s-dftd3" ] && [ -d "tblite" ]; then
    echo "  -> Symlinking s-dftd3 into tblite/subprojects"
    mkdir -p tblite/subprojects
    cd tblite/subprojects
    ln -sfn ../../s-dftd3 s-dftd3
    cd ../../
fi

echo "Symlink setup complete."
