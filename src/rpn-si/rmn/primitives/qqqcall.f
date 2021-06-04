*/* RMNLIB - Library of useful routines for C and FORTRAN programming
* * Copyright (C) 1975-2001  Division de Recherche en Prevision Numerique
* *                          Environnement Canada
* *
* * This library is free software; you can redistribute it and/or
* * modify it under the terms of the GNU Lesser General Public
* * License as published by the Free Software Foundation,
* * version 2.1 of the License.
* *
* * This library is distributed in the hope that it will be useful,
* * but WITHOUT ANY WARRANTY; without even the implied warranty of
* * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* * Lesser General Public License for more details.
* *
* * You should have received a copy of the GNU Lesser General Public
* * License along with this library; if not, write to the
* * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
* * Boston, MA 02111-1307, USA.
* */
***fonction qqqcall
      integer function qqqcall(extrn,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,
     % p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,
     % p25,p26,p27,p28,p29,p30,p31,p32,p33,p34,p35,p36,p37,p38,
     %p39,p40,p41)
      external extrn
      integer  extrn
*notes
*     ceci est un morceau ( le reste dans hlp_c) du mecanisme d'appel
*     a une fonction fortran par le biais de pointeurs renvoyes par
*     locf.
*     qqqcall est appele par la fonction c rmtcall
**
      qqqcall=extrn(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,
     % p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,
     % p25,p26,p27,p28,p29,p30,p31,p32,p33,p34,p35,p36,p37,p38,
     %p39,p40,p41)
      return
      end
