#!/bin/sh
set -e
if test "$CONFIGURATION" = "Debug"; then :
  cd /Users/jd827/Software/CardioMechanics/_buildXC
  make -f /Users/jd827/Software/CardioMechanics/_buildXC/CMakeScripts/ReRunCMake.make
fi
if test "$CONFIGURATION" = "Release"; then :
  cd /Users/jd827/Software/CardioMechanics/_buildXC
  make -f /Users/jd827/Software/CardioMechanics/_buildXC/CMakeScripts/ReRunCMake.make
fi
if test "$CONFIGURATION" = "MinSizeRel"; then :
  cd /Users/jd827/Software/CardioMechanics/_buildXC
  make -f /Users/jd827/Software/CardioMechanics/_buildXC/CMakeScripts/ReRunCMake.make
fi
if test "$CONFIGURATION" = "RelWithDebInfo"; then :
  cd /Users/jd827/Software/CardioMechanics/_buildXC
  make -f /Users/jd827/Software/CardioMechanics/_buildXC/CMakeScripts/ReRunCMake.make
fi

