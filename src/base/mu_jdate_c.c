/*
  !---------------------------------- LICENCE BEGIN -------------------------------
  ! GEM - Library of kernel routines for the GEM numerical atmospheric model
  ! Copyright (C) 1990-2010 - Division de Recherche en Prevision Numerique
  !                       Environnement Canada
  ! This library is free software; you can redistribute it and/or modify it 
  ! under the terms of the GNU Lesser General Public License as published by
  ! the Free Software Foundation, version 2.1 of the License. This library is
  ! distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
  ! without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
  ! PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
  ! You should have received a copy of the GNU Lesser General Public License
  ! along with this library; if not, write to the Free Software Foundation, Inc.,
  ! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
  !---------------------------------- LICENCE END ---------------------------------
*/

//TODO: EPOCH should be -47120101 not -47131125(24?)
/*
 * Functions in this modules are meant to manipulate date time
 * precise to the seconds and valid 
 * from YEAR=-4713, MONTH=11, DAY=25
 * for a periode of 7980 YEARS
 * 
 * Leap seconds are ignored.
 * Leap year is computed (cannot be ignored)
 * Time is considered to be GTM, thus:
 * - time zone is ignored
 * - daylight saving time is ignored
 *
 * The date-time stamp (js) is in seconds equivalent to julian day + time of day
 * js type is long long (64bits or more))
 *
 * Author: S. Chamberland (Jan 2016)
 * Author: FLIEGEL, FLANDERN, C. THIBEAULT (NOV 79) for Julian day convert to and from functions
 */

#include <stdbool.h>

static const long long MU_JDATE_MAX = 464269103999;
static const long long MU_JDATE_GREG0 = 148699584000;
// static const long long MU_JDATE_GREG0_NOLEAP = 148699584000;

static const int MU_JDATE_EPOCH_YY_NOLEAP = -4712;
static const int MU_JDATE_DAYS_YEAR_NOLEAP = 365;

static bool ignore_leap_year = (0 != 0);

void mu_set_leap_year(int i) {
  ignore_leap_year = (i != 0);
  return;
}

/*
  Return True if leap year
 */
bool mu_is_leapyear(int yy) {
    if (ignore_leap_year) return (1==0);
    return  yy % 4 == 0 && (yy % 100 != 0 || yy % 400 == 0);
}

/*
 * Return numbers of days in month (ignore leap year))
 */
int mu_days_in_month_noleap(int mo) {
    int a[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    return a[mo-1];
}

/*
 * Return numbers of days in month for specify year
 */
int mu_days_in_month(int yy, int mo) {
  if (mu_is_leapyear(yy)) {
        return mo != 2 ? mu_days_in_month_noleap(mo) : 29;
    } else {
        return mu_days_in_month_noleap(mo);
    }
}


/*
 * Return day of year corresponding to month, day, ignoring leap year
 */
int mu_day_of_year_noleap(int mo, int dd) {
    int month;
    int doy = dd;
    for (month = 1; month < mo; month ++) {
      doy += mu_days_in_month_noleap(month);
    }
    return doy;
}

/*
 * Return day of year corresponding to month, day for specified year
 */
int mu_day_of_year(int yy, int mo, int dd) {
    int month;
    int doy = dd;
    for (month = 1; month < mo; month ++) {
      doy += mu_days_in_month(yy,month);
    }
    return doy;
}

/*
 * Return  month, day corresponding to day of year, ignoring leap year
 */
void mu_doy2md_noleap(int doy, int* mo, int* dd) {
  static const int month_len[] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
  int day_of_month = doy;
  int month;
  for (month = 0; month < 12; month ++) {
    int mlen = month_len[month];
    if (day_of_month <= mlen)
      break;
    day_of_month -= mlen;
  }
  *mo = month+1;
  *dd = day_of_month;
  return;  
}

/*
 * Return  month, day corresponding to day of year for specified year
 */
void mu_doy2md_leap(int yy, int doy, int* mo, int* dd) {
  static const int month_len[] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
  int leap = (yy % 4 == 0) && (yy % 100 != 0 || yy % 400 == 0);
  int day_of_month = doy;
  int month;
  for (month = 0; month < 12; month ++) {
    int mlen = month_len[month];
    if (leap && month == 1) mlen ++;
    if (day_of_month <= mlen)
      break;
    day_of_month -= mlen;
  }
  *mo = month+1;
  *dd = day_of_month;
  return;  
}

/*
 *AUTHOR: FLIEGEL AND FLANDERN
 *REVISION 001  C. THIBEAULT   -  NOV 79  DOCUMENTATION
 *REVISION 002  C. THIBEAULT   -  MAR 83  DATEC (ANCIENNEMENT DATE)
 *                                CHANGEMENT DE NOM A CAUSE DE LA ROUTINE
 *                                DU SYSTEME QUI PORTE LE MEME NOM
 *REVISION 003  S. Chamebrland -  Jan 2016 - C version, leap year ignore
 *OBJECT: CONVERTS A JULIAN DAY TO THREE INTEGERS, YEAR(I), MONTH(J) AND DAY(K).
 *ARGUMENTS
 *   IN    - JD - JULIAN DAY, A UNIQUE INTEGER WHICH MAPS ONE-TO-ONE ONTO
 *                TRIPLES OF INTEGERS REPRESENTING YEAR, MONTH,
 *                DAY OF MONTH.
 *   OUT   - I  - YEAR
 *         - J  - MONTH
 *         - K  - DAY OF THE MONTH
 *NOTES    - COPIED FROM "COMMUNICATIONS OF THE ACM" (1968), PAGE 657.
 *         - IT COVERS A PERIOD OF 7980 YEARS WITH DAY 1 STARTING
 *           AT YEAR=-4713, MONTH=11, DAY=25
 */

void mu_jd2ymd_leap(long jd, int* yy, int* mo, int* dd) {
  long L = jd+68569;
  long N = 4*L/146097;
  long I, J, K;
  L = L-(146097*N+3)/4;
  I = 4000*(L+1)/1461001;
  L = L-1461*I/4+31;
  J = 80*L/2447;
  K = L-2447*J/80;
  L = J/11;
  J = J+2-12*L;
  I = 100*(N-49)+I+L;
  *yy = (int)I;
  *mo = (int)J;
  *dd = (int)K;
  return;
}

/*
 *AUTHOR:  C. THIBEAULT  -  JAN 1980
 *REVISION 001 C. THIBEAULT   - JAN 80  DOCUMENTATION
 *REVISION 002 S. Chamebrland -  Jan 2016 - C version, leap year ignore
 *OBJECT: COMPUTES THE JULIAN CALENDAR DAY, GIVEN A YEAR, A MONTH AND A DAY.
 *ARGUMENTS
 *         - JD - JULIAN DAY, A UNIQUE INTEGER WHICH MAPS ONE-TO-ONE
 *                ONTO TRIPLES OF INTEGERS REPRESENTING YEAR, MONTH,
 *                DAY OF THE MONTH.
 *         - I  - YEAR
 *         - J  - MONTH
 *         - K  - DAY OF THE MONTH.
 *NOTES    - COPIED FROM "COMMUNICATIONS OF THE ACM" (1968), PAGE 657.
 *         - IT COVERS A PERIOD OF 7980 YEARS WITH DAY 1 STARTING
 *           AT YEAR=-4713, MONTH=11, DAY=25.
 */
void mu_ymd2jd_leap(long *jd, int yy, int mo, int dd) {
  long I = (long)yy;
  long J = (long)mo;
  long K = (long)dd;
  *jd = K-32075+1461*(I+4800+(J-14)/12)/4
    +367*(J-2-(J-14)/12*12)/12-3
    *((I+4900+(J-14)/12)/100)/4;
  return;
}


/*
 * Julian day back and fort from gregorian date
 * while ignoring leap years
 */
void mu_jd2ymd_noleap(long jd, int* yy, int* mo, int* dd) {
  *yy = jd / MU_JDATE_DAYS_YEAR_NOLEAP;
  long doy = jd  - *yy * MU_JDATE_DAYS_YEAR_NOLEAP + 1;
  *yy += MU_JDATE_EPOCH_YY_NOLEAP;
  mu_doy2md_noleap(doy, mo, dd);
  return;
}

void mu_ymd2jd_noleap(long *jd, int yy, int mo, int dd) {
  *jd = (yy-MU_JDATE_EPOCH_YY_NOLEAP) * MU_JDATE_DAYS_YEAR_NOLEAP 
    + mu_day_of_year_noleap(mo,dd) - 1;
  return;
}

void mu_jd2ymd(long jd, int* yy, int* mo, int* dd) {
  if (ignore_leap_year) {
    mu_jd2ymd_noleap(jd, yy, mo, dd);
  } else {
    mu_jd2ymd_leap(jd, yy, mo, dd);
  }
  return;
}

void mu_ymd2jd(long *jd, int yy, int mo, int dd) {
  if (ignore_leap_year) {
    mu_ymd2jd_noleap(jd, yy, mo, dd);
  } else {
    mu_ymd2jd_leap(jd, yy, mo, dd);
  }
  return;
}


/*
 * Compute number of seconds equivalent of julian day
 * Ignores leap seconds
 */
long long mu_jd2js(long* jd) {
  return (long long)*jd * 86400;
}

/*
 * Compute julian day equivalent of number of seconds (round up to 00Z)
 * Ignores leap seconds
 */
long mu_js2jd(long long js) {
  return (long)(js / 86400);
}

/*
 * Compute number of julian seconds equivalent of julian day with time
 * Ignores leap seconds
 */
void mu_jdhms2js(long long* js, long jd, int hh, int mn, int ss) {
  *js = (long long)jd * 86400 + hh * 3600 + mn * 60 + ss;
  return;
}

/*
 * Compute julian day with time equivalent of number of julian seconds
 * Ignores leap seconds
 */
void mu_js2jdhms(long long js, long* jd, int* hh, int* mn, int* ss) {
  *jd = (long)(js / 86400);
  *hh = (int)(js % 86400)/3600;
  *mn = (int)(js - *jd*86400 - *hh*3600)/60;
  *ss = (int)(js - *jd*86400 - *hh*3600 - *mn*60);
  return;
}

/*
 * Compute number of julian seconds equivalent of date time
 * Ignores leap seconds
 * if ignore_leap_year: js = js0 - jsmax
 */
void mu_ymdhms2js_leap(long long* js, 
                  int yy, int mo, int dd, int hh, int mn, int ss) {
  long jd;
  mu_ymd2jd_leap(&jd, yy, mo, dd);
  mu_jdhms2js(js, jd, hh, mn, ss);
  return;
}

void mu_ymdhms2js_noleap(long long* js, 
                  int yy, int mo, int dd, int hh, int mn, int ss) {
  long jd;
  mu_ymd2jd_noleap(&jd, yy, mo, dd);
  mu_jdhms2js(js, jd, hh, mn, ss);
  *js -= MU_JDATE_MAX;
  return;
}

void mu_ymdhms2js(long long* js, 
                  int yy, int mo, int dd, int hh, int mn, int ss) {
  if (ignore_leap_year) {
    mu_ymdhms2js_noleap(js, yy,mo, dd, hh, mn, ss);
  } else {
    mu_ymdhms2js_leap(js, yy, mo, dd, hh, mn, ss);
  }
  return;
}

/*
 * Compute date time equivalent of number of julian seconds
 * Ignores leap seconds
 * if js < 0: consider ignore_leap_year=1 and use js1 = js + jsmax
 */
void mu_js2ymdhms_leap(long long js, 
                  int* yy, int* mo, int* dd, int* hh, int* mn, int* ss) {
  long jd;
  mu_js2jdhms(js, &jd, hh, mn, ss);
  mu_jd2ymd_leap(jd, yy, mo, dd);
  return;
}

void mu_js2ymdhms_noleap(long long js, 
                  int* yy, int* mo, int* dd, int* hh, int* mn, int* ss) {
  long jd;
  long long js0;
  //TODO: if (js >= 0) {ERROR}
  js0 = js + MU_JDATE_MAX;
  mu_js2jdhms(js0, &jd, hh, mn, ss);
  mu_jd2ymd_noleap(jd, yy, mo, dd);
  return;
}

void mu_js2ymdhms(long long js, 
                  int* yy, int* mo, int* dd, int* hh, int* mn, int* ss) {
  if (js >= 0) {
    mu_js2ymdhms_leap(js,  yy, mo, dd, hh, mn, ss);
  } else {
    mu_js2ymdhms_noleap(js, yy, mo, dd, hh, mn, ss);
  }
  return;
}
