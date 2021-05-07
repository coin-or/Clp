
# CLP_ABOCA([default_buildtype])

#  where default_buildtype ($1) is one of:
#  1 - build Abc serial but no inherit code (default when built)
#  2 - build Abc serial and inherit code
#  3 - build Abc cilk parallel but no inherit code
#  4 - build Abc cilk parallel and inherit code

AC_DEFUN([AC_CLP_ABOCA],
[
  aboca_dflt_buildtype=m4_default([$1],1)
  AC_ARG_ENABLE([aboca],
    [AS_HELP_STRING([--enable-aboca],[Enables build of the Aboca solver;
       1 - build Abc serial but no inherit code (default when built);
       2 - build Abc serial and inherit code;
       3 - build Abc cilk parallel but no inherit code;
       4 - build Abc cilk parallel and inherit code])],
      [use_aboca=$enableval],
      [use_aboca=no])
  if test "$use_aboca" = yes ; then
    use_aboca=$aboca_dflt_buildtype
  fi
  case "$use_aboca" in 
    1 | 2 | 3 | 4)
      AC_DEFINE_UNQUOTED([CLP_HAS_ABC],[$use_aboca],
        [Define to 1, 2, 3, or 4 if Aboca should be built.])
      use_aboca=yes
      ;;
    no)
      ;;
    *)
      AC_MSG_WARN([invalid argument for --enable-aboca: $use_aboca])
      use_aboca=no
      ;;
  esac
  AM_CONDITIONAL(COIN_HAS_ABC,[test "$use_aboca" = yes])
])    # AC_CLP_ABOCA

