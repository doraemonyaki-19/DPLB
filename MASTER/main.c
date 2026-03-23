/*  Main routine for suspension code  */

/***********************************************************************
 * ASLB : Lattice-Boltzmann simulation code for deformable particle+polymer+
 * lattice Boltzmann fluid
 * Copyright (C) 2019 Yeng-Long Chen
 *
 * based on 
 * Susp3D: Lattice-Boltzmann simulation code for particle-fluid suspensions
 * Copyright (C) 2003 Tony Ladd
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * ***********************************************************************/

#include "header.h"
#include <time.h>

int  num_proc, n_proc;

int main (int argc, char **argv)
{
  double ntime;
  double start, finish;
  char  *work_dir;
  time_t t1, t2; 
 
  init_procs (&argc, &argv);                          /* Initialization */
  num_proc = 1; 
  n_proc   = 0;

  work_dir = argv[1];                                 /* Work directory */
 
  if(argc == 1)
    work_dir=".";
 
  if (n_proc == 0)
    {
      fprintf (stdout, "Running on %d processors\n", num_proc);
      fflush  (stdout);
    }
  
  t1=time(NULL);
  start = wclock();
  driver (work_dir);
  finish = wclock();
  t2=time(NULL);  

  fini_procs ();                                      /* Wrap up */
  
  fprintf (stdout, "Elapsed total time on proc %3d: %le (%le %le); real time = %d seconds\n", n_proc, finish-start, start, finish, (int)(t2-t1));
  fflush  (stdout);
  return (0);
}
