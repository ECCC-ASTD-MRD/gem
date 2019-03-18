/* RMNLIB - Library of useful routines for C and FORTRAN programming
 * Copyright (C) 1975-2001  Division de Recherche en Prevision Numerique
 *                          Environnement Canada
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "ezscint.h"
#include "ez_funcdef.h"

void f77name(ez_calcdist)(float *distance, float *lat1, float *lon1, float *lat2, float *lon2);
void f77name(ez_calcarea)(float *area, float lats[], float lons[]);
void f77name(ez_calcarea2)(float *area, float lats[], float lons[]);
void c_ez_calcdist(float *distance, float lat1, float lon1, float lat2, float lon2);
void c_ez_calcdist2(double *distance, float lat1, float lon1, float lat2, float lon2);
void c_ez_calcarea_rect(float *area,     float lat1, float lon1, float lat2, float lon2);
void c_ez_calcarea(float *area, float lats[], float lons[]);
void c_ez_calcarea2(double *area, float lats[], float lons[]);

void f77name(ez_calcdist)(float *distance, float *lat1, float *lon1, float *lat2, float *lon2)
{
   c_ez_calcdist(distance, *lat1, *lon1, *lat2, *lon2);
}

void f77name(ez_calcdist2)(double *distance, float *lat1, float *lon1, float *lat2, float *lon2)
{
   c_ez_calcdist2(distance, *lat1, *lon1, *lat2, *lon2);
}

void f77name(ez_calcarea_rect)(float *area, float *lat1, float *lon1, float *lat2, float *lon2)
{
   c_ez_calcarea_rect(area, *lat1, *lon1, *lat2, *lon2);
}

void f77name(ez_calcarea)(float *area, float lats[], float lons[])
{
   c_ez_calcarea(area, lats, lons);
}

void f77name(ez_calcarea2)(float *area, float lats[], float lons[])
{
   double area2;
   c_ez_calcarea2(&area2, lats, lons);
   *area = (float)area2;
}
/*
c_ez_calcdist :
   This function computes the distance between 2 latlon points on the sphere.
   Source of the formula : http://mathworld.wolfram.com/GreatCircle.html
   Latitudes and longitudes are in degrees.
*/


void c_ez_calcdist(float *distance, float lat1, float lon1, float lat2, float lon2)
{
   double degre_a_radian;
   double radlat1, radlat2, radlon1, radlon2;
   double dist;
   double earth_radius = 6370997.;

   degre_a_radian = M_PI / 180.0;

   radlat1 = lat1 * degre_a_radian;
   radlat2 = lat2 * degre_a_radian;
   radlon1 = lon1 * degre_a_radian;
   radlon2 = lon2 * degre_a_radian;

   *distance = (float)(earth_radius*acos(cos(radlat1)*cos(radlat2)*cos(radlon1-radlon2)+sin(radlat1)*sin(radlat2)));
}

void c_ez_calcdist2(double *distance, float lat1, float lon1, float lat2, float lon2)
{
   double degre_a_radian;
   double radlat1, radlat2, radlon1, radlon2;
   double dist;
   double earth_radius = 6370997.;
   double a,c,dlat,dlon,sindlat, sindlon;

   degre_a_radian = M_PI / 180.0;

   radlat1 = lat1 * degre_a_radian;
   radlat2 = lat2 * degre_a_radian;
   radlon1 = lon1 * degre_a_radian;
   radlon2 = lon2 * degre_a_radian;

   dlon = radlon2 - radlon1;
   dlat = radlat2 - radlat1;
   sindlat = sin(dlat*0.5);
   sindlon = sin(dlon*0.5);
   a = sindlat*sindlat + cos(radlat1) * cos(radlat2) * sindlon*sindlon;
   c = 2.0 * atan2(sqrt(a),sqrt(1.0-a));
   c = 2.0 * asin(sqrt(a));
   *distance = (float)(c*earth_radius);

/*   *distance = (float)(earth_radius*acos(cos(radlat1)*cos(radlat2)*cos(radlon1-radlon2)+sin(radlat1)*sin(radlat2)));*/
}
/*
  c_ez_calcarea :
   This function computes the area of the solid rectangle formed by  2 latlon points on the sphere.
   Source of the formula : http://mathworld.wolfram.com/GreatCircle.html
                           http://mathworld.wolfram.com/SphericalTrigonometry.html
                           http://mathworld.wolfram.com/SphericalTriangle.html
   Latitudes and longitudes are in degrees.

   The computation is done by splitting the solid rectangle formed by the latlon points into 2 triangles,
   compute the area of each triangle and add the two areas

   If we have
                                    seg_d
                     *--------------------------------------* (lat2, lon2)
                     *                         -----------  *
             seg_e   *          ---------------             * seg_b
                     *----------     seg_c                  *
      (lat1, lon1)   *--------------------------------------*
                                     seg_a

  We have the corresponding angles for triangle 1

                                                          a * (lat2, lon2)
                                               -----------  *
                                ---------------             *
                     *----------                         c  *
      (lat1, lon1)   *-b -----------------------------------*

  and for triangle 2


                     *------------------------------------e-* (lat2, lon2)
                     * c                       -----------  *
                     *          ---------------
                     *-d--------
      (lat1, lon1)   *
*/



void c_ez_calcarea_rect(float *area, float lat1, float lon1, float lat2, float lon2)
{
   double degre_a_radian;
   double radlat1, radlat2, radlon1, radlon2;
   double dist;
   double earth_radius = 6370997.;
   double a,b,c,d,e;
   double seg_a, seg_b, seg_c, seg_d, seg_e;
   double area1, area2;

   degre_a_radian = M_PI / 180.0;

   radlat1 = lat1 * degre_a_radian;
   radlat2 = lat2 * degre_a_radian;
   radlon1 = lon1 * degre_a_radian;
   radlon2 = lon2 * degre_a_radian;

   /*
   seg_a = distance pt(1,1) - pt(2,1)
   seg_b = distance pt(2,1) - pt(2,2)
   seg_c = distance pt(1,1) - pt(2,2)
   seg_d = distance pt(1,2) - pt(2,2)
   seg_e = distance pt(1,1) - pt(1,2)
   */

   /* These computations have not been optimized */

   seg_a = acos(cos(radlat1)*cos(radlat1)*cos(radlon1-radlon2)+sin(radlat1)*sin(radlat1));
   seg_b = acos(cos(radlat1)*cos(radlat2)*cos(radlon2-radlon2)+sin(radlat1)*sin(radlat2));
   seg_c = acos(cos(radlat1)*cos(radlat2)*cos(radlon1-radlon2)+sin(radlat1)*sin(radlat2));
   seg_d = acos(cos(radlat2)*cos(radlat2)*cos(radlon1-radlon2)+sin(radlat2)*sin(radlat2));
   seg_e = acos(cos(radlat1)*cos(radlat2)*cos(radlon1-radlon1)+sin(radlat1)*sin(radlat2));
/*   printf("seg_a:%f seg_b:%f seg_c:%f seg_d:%f seg_e:%f\n", seg_a, seg_b, seg_c, seg_d, seg_e);*/

   c = acos((cos(seg_c)-cos(seg_b)*cos(seg_a))/(sin(seg_b)*sin(seg_a)));
   b = asin(sin(seg_b)*sin(c)/sin(seg_c));
   a = asin(sin(seg_a)*sin(c)/sin(seg_c));
   area1 = ((a + b + c) - M_PI) * earth_radius * earth_radius;
/*   printf("a:%f b:%f c:%f area:%f\n", a, b, c, area1);*/

/*   printf("seg_a:%f seg_b:%f seg_c:%f seg_d:%f seg_e:%f\n", seg_a, seg_b, seg_c, seg_d, seg_e);*/
   c = acos((cos(seg_c)-cos(seg_e)*cos(seg_d))/(sin(seg_e)*sin(seg_d)));
   e = asin(sin(seg_e)*sin(c)/sin(seg_c));
   d = asin(sin(seg_d)*sin(c)/sin(seg_c));
   area2 = ((d + e + c) - M_PI) * earth_radius * earth_radius;
/*   printf("d:%f e:%f c:%f area:%f\n", d, e, c, area2);*/
   *area = (float)(area1 + area2);

}



void c_ez_calcarea(float *area, float lats[], float lons[])
{
   double degre_a_radian;
   double radlat[4], radlon[4];
   double dist;
   double earth_radius = 6370997.;
   double a,b,c,d,e;
   double seg_a, seg_b, seg_c, seg_d, seg_e;
   double area1, area2;
   int i;

   degre_a_radian = M_PI / 180.0;

   for (i=0; i < 4; i++)
      {
      radlat[i] = lats[i] * degre_a_radian;
      radlon[i] = lons[i] * degre_a_radian;
      }


   /*
   seg_a = distance pt(1,1) - pt(2,1)
   seg_b = distance pt(2,1) - pt(2,2)
   seg_c = distance pt(1,1) - pt(2,2)
   seg_d = distance pt(1,2) - pt(2,2)
   seg_e = distance pt(1,1) - pt(1,2)
   */

   /* These computations have not been optimized */

   seg_a = acos(cos(radlat[0])*cos(radlat[1])*cos(radlon[0]-radlon[1])+sin(radlat[0])*sin(radlat[1]));
   seg_b = acos(cos(radlat[1])*cos(radlat[2])*cos(radlon[1]-radlon[2])+sin(radlat[1])*sin(radlat[2]));
   seg_c = acos(cos(radlat[0])*cos(radlat[2])*cos(radlon[0]-radlon[2])+sin(radlat[0])*sin(radlat[2]));
   seg_d = acos(cos(radlat[2])*cos(radlat[3])*cos(radlon[2]-radlon[3])+sin(radlat[2])*sin(radlat[3]));
   seg_e = acos(cos(radlat[0])*cos(radlat[3])*cos(radlon[0]-radlon[3])+sin(radlat[0])*sin(radlat[3]));
/*   printf("seg_a:%f seg_b:%f seg_c:%f seg_d:%f seg_e:%f\n", seg_a, seg_b, seg_c, seg_d, seg_e);*/

   c = acos((cos(seg_c)-cos(seg_b)*cos(seg_a))/(sin(seg_b)*sin(seg_a)));
   b = asin(sin(seg_b)*sin(c)/sin(seg_c));
   a = asin(sin(seg_a)*sin(c)/sin(seg_c));
   area1 = ((a + b + c) - M_PI) * earth_radius * earth_radius;

/*   printf("seg_a:%f seg_b:%f seg_c:%f seg_d:%f seg_e:%f\n", seg_a, seg_b, seg_c, seg_d, seg_e);*/
   c = acos((cos(seg_c)-cos(seg_e)*cos(seg_d))/(sin(seg_e)*sin(seg_d)));
   e = asin(sin(seg_e)*sin(c)/sin(seg_c));
   d = asin(sin(seg_d)*sin(c)/sin(seg_c));
   area2 = ((d + e + c) - M_PI) * earth_radius * earth_radius;
   *area = (float)(area1 + area2);
   if (*area < 0.0)
      {
      printf("a:%f b:%f c:%f area:%f\n", a, b, c, area1);
      printf("d:%f e:%f c:%f area:%f\n", d, e, c, area2);
      printf("area:%f\n", *area);
      }

}

void c_ez_calcarea2(double *area, float lats[], float lons[])
{
   double degre_a_radian;
   double radlat[4], radlon[4];
   double dist;
   double earth_radius = 6370997.;
   double a,b,c,d,e;
   double seg_a, seg_b, seg_c, seg_d, seg_e;
   double area1, area2;
   double s, excess;
   int i;

   degre_a_radian = M_PI / 180.0;

   for (i=0; i < 4; i++)
      {
      radlat[i] = lats[i] * degre_a_radian;
      radlon[i] = lons[i] * degre_a_radian;
      }


   /*
   seg_a = distance pt(1,1) - pt(2,1)
   seg_b = distance pt(2,1) - pt(2,2)
   seg_c = distance pt(1,1) - pt(2,2)
   seg_d = distance pt(1,2) - pt(2,2)
   seg_e = distance pt(1,1) - pt(1,2)
   */

   /* These computations have not been optimized */

   c_ez_calcdist2(&seg_a, lats[0], lons[0], lats[1], lons[1]);
   c_ez_calcdist2(&seg_b, lats[1], lons[1], lats[2], lons[2]);
   c_ez_calcdist2(&seg_c, lats[0], lons[0], lats[2], lons[2]);
   c_ez_calcdist2(&seg_d, lats[2], lons[2], lats[3], lons[3]);
   c_ez_calcdist2(&seg_e, lats[3], lons[3], lats[0], lons[0]);
/*
   printf("seg_a:%f seg_b:%f seg_c:%f seg_d:%f seg_e:%f\n", seg_a, seg_b, seg_c, seg_d, seg_e);
*/

   seg_a /= earth_radius;
   seg_b /= earth_radius;
   seg_c /= earth_radius;
   seg_d /= earth_radius;
   seg_e /= earth_radius;

   s = 0.5*(seg_a+seg_b+seg_c);
   excess = sqrt(tan(0.5*s)*tan(0.5*(s-seg_a))*tan(0.5*(s-seg_b))*tan(0.5*(s-seg_c)));
   excess = atan(excess);
   excess = 4.0 * excess;
   area1 = earth_radius * earth_radius * excess;
/*   printf("s: %f excess: %f area1: %f\n", s, excess, area1);*/

   s = 0.5*(seg_c+seg_d+seg_e);
   excess = sqrt(tan(0.5*s)*tan(0.5*(s-seg_c))*tan(0.5*(s-seg_d))*tan(0.5*(s-seg_e)));
   excess = atan(excess);
   excess *= 4.0;
   area2 = earth_radius * earth_radius * excess;
/*   printf("s: %f excess: %f area2: %f\n", s, excess, area2);*/

   *area = (float)(area1 + area2);
   if (*area < 0.0)
      {
      printf("area1:%f\n", area1);
      printf("area2:%f\n", area2);
      }

}
/*
main()
{
float distance, area, lat1, lon1, lat2, lon2;
int i;
lat1=0.0;
lat2=1.0;
lon1=0.0,
lon2=1.0;

c_ez_calcdist(&distance, lat1, lon1, lat2, lon2);
printf("%f\n", distance);
c_ez_calcarea(&area, lat1, lon1, lat2, lon2);
printf("%f\n", area/1.0e6);

for (i=-90; i < 89; i++)
   {
   lat1 = (float) i+0.01;
   lat2 = lat1 + 1.0;
   c_ez_calcarea(&area, lat2, lon2, lat1, lon1);
   printf("Aire: (%f %f) (%f %f) %f\n", lat1, lon1, lat2, lon2, area/1.0e6);
   }

}*/
