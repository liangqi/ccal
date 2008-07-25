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

   Adapted from novascon.c from the NOVAS-C package
   The whole package can be obtained from
   http://aa.usno.navy.mil/AA/software/novas/novas_c/novasc_info.html

   The original header comments:

   NOVAS-C Version 2.0 (1 Nov 98)
   Constants file

   Naval Observatory Vector Astrometry Subroutines
   C Version

   U. S. Naval Observatory
   Astronomical Applications Dept.
   3450 Massachusetts Ave., NW
   Washington, DC  20392-5420
*/

#ifndef _CONSTS_H
   #include "novascon.h"
#endif

const short int FN0 = 0;

/*
   TDB Julian date of epoch J2000.0.
*/

const double T0 = 2451545.00000000;

/*
   Speed of light in AU/Day.
*/

const double C = 173.14463348;

/*
   Value of pi in radians.
*/

const double TWOPI = 6.28318530717958647692;

/*
   Angle conversion constants.
*/

const double RAD2SEC = 206264.806247096355;
const double DEG2RAD = 0.017453292519943296;
const double RAD2DEG = 57.295779513082321;

