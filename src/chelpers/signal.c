/*
 *  --- taken from the internet ---
 *
 *  It overrides the FORTRAN intrinsic signal routine
 *  and gives the possibility to use custom signal
 *  handlers instead of the standard ones, pretty simple
 *  code but actually impossible to do with FORTRAN only.
 */
#include <signal.h>

typedef void (*sighandler_t)(int);

void signal_( int* signum, sighandler_t handler)
{
   signal(*signum, handler);
}

/* Called from Fortran via ISO_C_BINDING for ifx builds.
 * handler is a BIND(C) Fortran subroutine with no dummy arguments;
 * the int signum passed by the OS is silently ignored. */
void crest_install_signal(int signum, void (*handler)(void))
{
    signal(signum, (sighandler_t)handler);
}
