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

   The whole package can be obtained from
   http://aa.usno.navy.mil/AA/software/novas/novas_c/novasc_info.html

   NOVAS-C Version 2.0 (1 Nov 98)
   Header file for all source files containing versions of
   NOVAS-C function 'solarsystem'

   Naval Observatory Vector Astrometry Subroutines
   C Version

   U. S. Naval Observatory
   Astronomical Applications Dept.
   3450 Massachusetts Ave., NW
   Washington, DC  20392-5420
*/

#ifndef _SOLSYS_
   #define _SOLSYS_

/*
   Function prototypes
*/

   short int solarsystem (double tjd, short int body, short int origin,

                          double *pos, double *vel);

#endif
