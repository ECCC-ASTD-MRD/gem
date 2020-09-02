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


#define powerSpan 65
static double powerOf2s[powerSpan];
static int powerOf2sInitialized = 0;

#define isDouble 1
#define FLOAT_4_8 double
#define compact_IEEEblock_FLOAT_4_8 compact_IEEEblock_double
#include "compact_IEEEblock_h.c"


#define isDouble 0
#define FLOAT_4_8 float
#define compact_IEEEblock_FLOAT_4_8 compact_IEEEblock_float
#include "compact_IEEEblock_h.c"
