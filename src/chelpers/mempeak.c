/*
! This file is part of crest.
!
! Copyright (C) 2025 Philipp Pracht
!
! crest is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! crest is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with crest.  If not, see <https://www.gnu.org/licenses/>.

*/
#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <psapi.h>
#include <windows.h>
__declspec(dllexport) long long get_peak_rss_kb(void) {
  PROCESS_MEMORY_COUNTERS pmc;
  if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc))) {
    // bytes -> kilobytes
    return (long long)(pmc.PeakWorkingSetSize / 1024);
  }
  return -1;
}
#else
#include <sys/resource.h>
#include <sys/time.h>
long long get_peak_rss_kb(void) {
  struct rusage ru;
  if (getrusage(RUSAGE_SELF, &ru) == 0) {
// On Linux: ru_maxrss is in kilobytes
// On macOS/BSD: ru_maxrss is in bytes — convert to KB
#ifdef __APPLE__
    return (long long)(ru.ru_maxrss / 1024);
#else
    return (long long)ru.ru_maxrss;
#endif
  }
  return -1;
}
#endif
