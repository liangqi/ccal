/*
   As part of the ccal program by Zhuo Meng, this version is
   distributed under the terms of the GNU Lesser General Public License as
   published by the Free Software Foundation; either version 2 of the
   License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

   Adapted from novas.h from the NOVAS-C package
   The whole package can be obtained from
   http://aa.usno.navy.mil/AA/software/novas/novas_c/novasc_info.html

   The original header comments:

   NOVAS-C Version 2.0 (1 Nov 98)
   Header file for novas.c

   Naval Observatory Vector Astrometry Subroutines
   C Version

   U. S. Naval Observatory
   Astronomical Applications Dept.
   3450 Massachusetts Ave., NW
   Washington, DC  20392-5420
*/

#ifndef _NOVAS_H
   #define _NOVAS_H

   #ifndef __STDIO__
      #include <stdio.h>
   #endif

   #ifndef __MATH__
      #include <math.h>
   #endif

   #ifndef __STRING__
      #include <string.h>
   #endif

   #ifndef __STDLIB__
      #include <stdlib.h>
   #endif

   #ifndef __CTYPE__
      #include <ctype.h>
   #endif

   #ifndef _CONSTS_
      #include "novascon.h"
   #endif

   #ifndef _SOLSYS_
      #include "solarsystem.h"
   #endif

/*
   Define "origin" constants.
*/

   #define BARYC  0
   #define HELIOC 1

/*
   Function prototypes
*/

   void earthtilt (double tjd,
                   double *mobl, double *tobl, double *eqeq,
                   double *psi, double *eps);

   short int aberration (double *pos, double *vel, double lighttime,
                         double *pos2);

   void precession (double tjd1, double *pos, double tjd2,
                    double *pos2);

   short int nutate (double tjd, short int fn1, double *pos,
                     double *pos2);

   short int nutation_angles (double tdbtime,
                              double *longnutation,
                              double *obliqnutation);

   void fund_args (double t,
                   double a[5]);

   void tdb2tdt (double tdb,
                 double *tdtjd, double *secdiff);

   void radec2vector (double ra, double dec, double dist,
                      double *vector);

   short int solarsystem (double tjd, short int body, short int origin, 
                          double *pos, double *vel);

   double julian_date (short int year, short int month, short int day,
                       double hour);

   void cal_date (double tjd,
                  short int *year, short int *month, short int *day,
                  double *hour);

#endif
