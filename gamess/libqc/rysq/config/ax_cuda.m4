
AC_DEFUN([AX_CUDA],
[
AC_ARG_WITH([cuda],
            [AS_HELP_STRING([--with-cuda@<:@=DIR@:>@],
                            [use cuda (default is yes) - it is possible to specify
 			     the root directory for cuda (optional)])],
            [if test "$withval" = "no"; then
                 CUDA_ROOT="/bin/false"
             elif test "$withval" != "yes"; then
                 CUDA_ROOT="$withval"
             fi])

# ## Cuda default device 
# AC_ARG_WITH(cuda-default,
#             [AS_HELP_STRING([--with-cuda-default[=0]],[default device])],
#             [cuda_default_device=$withval],[cuda_default_device=0])
# # AC_DEFINE_UNQUOTED(CUDA_DEFAULT_DEVICE, $cuda_default_device, [default device])

## CUDA device emulation
AC_ARG_WITH(cuda-emulation,
            AS_HELP_STRING([--with-cuda-emulation],[CUDA device emulation]),
            [],[withval="no"])
if test "$withval" = "yes"; then
   CUDAFLAGS="$CUDAFLAGS --device-emulation --device-compilation='C++'"
fi

## Cuda architecture
AC_ARG_WITH(cuda-arch,
            AS_HELP_STRING([--with-cuda-arch[=arch]],[CUDA device arch]),
            [with_cuda_arch=$withval],[with_cuda_arch=sm_13])
# Set compiler flags
CUDAFLAGS="$CUDAFLAGS -arch=$with_cuda_arch"

## verbos if debugging enable
if test "$enable_debug" = "yes" ; then
    CUDAFLAGS="$CUDAFLAGS -g -O0 --ptxas-options='-v'"
fi

])

AC_DEFUN([AX_HAVE_CUDA],
[
HAVE_CUDA=""
if test -n "$CUDA"; then
    HAVE_CUDA="yes"
    AC_DEFINE(HAVE_CUDA, 1, [have CUDA])
fi
AM_CONDITIONAL(HAVE_CUDA, test -n "$HAVE_CUDA")
AC_SUBST([HAVE_CUDA], `expr "$HAVE_CUDA" != ""`)
])

 #fi  # HAVE_CUDA