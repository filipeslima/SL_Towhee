#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2005-2014 Marcus G. Martin                       *
c     * see the file license.gpl for the full license information      *
c     *                                                                *
c     * This program is free software; you can redistribute it and/or  *
c     * modify it under the terms of the GNU General Public License    *
c     * as published by the Free Software Foundation; either version 2 *
c     * of the License, or (at your option) any later version.         *
c     *                                                                *
c     * This program is distributed in the hope that it will be useful,*
c     * but WITHOUT ANY WARRANTY; without even the implied warranty of *
c     * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the  *
c     * GNU General Public License for more details.                   *
c     *                                                                *
c     * You should have received a copy of the GNU General Public      *
c     * License along with this program; if not, write to the Free     *
c     * Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,*
c     * MA  02111-1307, USA.                                           *
c     ******************************************************************
      subroutine twh_writeintro()
c     ******************************************************************
c     * writes out the introductory information about Towhee           *
c     *                                                                *
c     * split from mainloop 05-13-2005 by M.G. Martin                  *
c     * last modified 07-25-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
      character*7 version
      character*20 date

      version = '7.1.0'
      call twh_version(GLB_SET,version)
      date = '(July 25 2014)'
c     --- Output some information about the code
      write(6,'(1x,a22,1x,a,1x,a)') 'MCCCS Towhee - Version',version
     &     ,date
c     --- 2014 copyright plus code contributors
      write(6,'(1x,a)') 'Copyright (C) <2014> Marcus G. Martin'
      write(6,'(3(1x,a))')
     &     'Code Contributors for 2014:'
     &     ,'Marcus G. Martin;'
c     --- 2013 copyright plus code contributors
      write(6,'(1x,a)') 'Copyright (C) <2013> Marcus G. Martin'
      write(6,'(3(1x,a))')
     &     'Code Contributors for 2013:'
     &     ,'Wesley Emeneker;'
     &     ,'Marcus G. Martin;'
c     --- 2012 copyright plus code contributors
      write(6,'(1x,a)') 'Copyright (C) <2012> Marcus G. Martin'
      write(6,'(2(1x,a))')
     &     'Code Contributors for 2012:'
     &     ,'Marcus G. Martin;'
c     --- 2011 copyright plus code contributors
      write(6,'(1x,a)') 'Copyright (C) <2011> Marcus G. Martin'
      write(6,'(4(1x,a))')
     &     'Code Contributors for 2011:'
     &     ,'Fred James;'
     &     ,'Marcus G. Martin;'
     &     ,'Suleiman Oloriegbe;'
c     --- 2010 copyright plus code contributors
      write(6,'(1x,a)') 'Copyright (C) <2010> Marcus G. Martin'
      write(6,'(5(1x,a))')
     &     'Code Contributors for 2010:'
     &     ,'Rene Haber;'
     &     ,'Andrej Lajovic;'
     &     ,'Marcus G. Martin;'
     &     ,'Loukas Peristeras;'
c     --- 2009 copyright plus code contributors
      write(6,'(1x,a)') 'Copyright (C) <2009> Marcus G. Martin'
      write(6,'(4(1x,a))')
     &     'Code Contributors for 2009:'
     &     ,'Rene Haber;'
     &     ,'Marcus G. Martin;'
     &     ,'Loukas Peristeras;'
c     --- 2008 copyright plus code contributors
      write(6,'(1x,a)') 'Copyright (C) <2008> Marcus G. Martin'
      write(6,'(6(1x,a))')
     &     'Code Contributors for 2008:'
     &     ,'Rene Haber;'
     &     ,'Iyad A. Hijazi;'
     &     ,'Marcus G. Martin;'
     &     ,'Loukas Peristeras;'
     &     ,'Craig Tenney;'
c     --- 2007 copyright plus code contributors
      write(6,'(1x,a)') 'Copyright (C) <2007> Marcus G. Martin'
      write(6,'(7(1x,a))')
     &     'Code Contributors for 2007:'
     &     ,'Bernhard Eckl;'
     &     ,'Frank Heilmann;'
     &     ,'Jianhui Li;'
     &     ,'Marcus G. Martin;'
     &     ,'Craig M. Tenney;'
     &     ,'Ozgur Yazaydin'
c     --- 2006 copyright plus code contributors
      write(6,'(1x,a)') 'Copyright (C) <2006> Marcus G. Martin'
      write(6,'(9(1x,a))')
     &     'Code Contributors for 2006:'
     &     ,'Alan A. Chen;'
     &     ,'Arben Jusufi;'
     &     ,'Christian D. Lorenz;'
     &     ,'Marcus G. Martin;'
     &     ,'Jeffrey J. Potoff;'
     &     ,'Vincent K. Shen;'
     &     ,'Matthew A. Wyczalkowski;'
     &     ,'Ozgur Yazaydin'
c     --- 2005 copyright plus code contributors
      write(6,'(1x,a)') 'Copyright (C) <2005> Marcus G. Martin'
      write(6,'(5(1x,a))')
     &     'Code Contributors for 2005:'
     &     ,'Marcus G. Martin;'
     &     ,'Alan A. Chen;'
     &     ,'Christian D. Lorenz;'
     &     ,'Matthew A. Wyczalkowski'
c     --- 2004 copyright plus code contributors
      write(6,'(1x,a)') 'Copyright (C) <2004> Marcus G. Martin'
      write(6,'(9(1x,a))') 
     &     'Code Contributors for 2004:'
     &     ,'Marcus G. Martin;'
     &     ,'Christian D. Lorenz;'
     &     ,'Peter A. Schultz;'
     &     ,'Aidan P. Thompson;'
     &     ,'Jeffrey J. Potoff;'
     &     ,'Steve R. Lustig;'
     &     ,'Paul S. Crozier;'
     &     ,'Nicholas du Preez'
c     --- 2003 copyright plus code contributors
      write(6,'(1x,a)') 'Copyright (C) <2003> Marcus G. Martin'
      write(6,'(7(1x,a))') 
     &     'Code Contributors for 2003:'
     &     ,'Marcus G. Martin;'
     &     ,'Christian D. Lorenz;'
     &     ,'David B. Ritch;'
     &     ,'Michael L. Greenfield;'
     &     ,'Yogesh Srivastava;'
     &     ,'Evangelos A Coutsias'
c     --- 2002
      write(6,'(1x,a)') 'Copyright (C) <2002> Marcus G. Martin'
      write(6,'(5(1x,a))') 
     &     'Code Contributors for 2002:'
     &     ,'Marcus G. Martin;'
     &     ,'Christian D. Lorenz;'
     &     ,'Aidan P. Thompson;'
     &     ,'Allen G. Sault'
c     --- 2000-2001
      write(6,'(1x,a)') 'Copyright (C) <2000-2001> Marcus G. Martin'
c     --- 1999
      write(6,'(1x,a,1x,a)') 
     &     'Copyright (C) <1999> J. Ilja Siepmann and' 
     &     ,'Marcus G. Martin'
      write(6,'(6(1x,a))')
     &     'Code Contributors for 1999:'
     &     ,'J. Ilja Siepmann;'
     &     ,'Marcus G. Martin;'
     &     ,'Bin Chen;'
     &     ,'Collin D. Wick;'
     &     ,'John Stubbs'
      write(6,*)
c     --- links to the web manual
      write(6,'(1x,a,1x,a)') 
     &     'Direct comments about this code to Marcus G. Martin'
     &     ,'(marcus_martin@users.sourceforge.net)'
      write(6,*) 
     &     'See the MCCCS Towhee web site for more information '
     &     ,'and the users manual.'
      write(6,*) 'http://towhee.sourceforge.net'
      write(6,*) 
c     --- thank folks who gave us money for code development
      write(6,'(1x,a,1x,a)')
     &     'We would like to thank the following agencies for'
     &     ,'providing funding used to develop MCCCS Towhee'
      write(6,*) '-Department of Energy Industrial Technologies '
     &     ,'Program (MGM)'
      write(6,*) '-National Science Foundation (JIS)'
      write(6,*) 
     &     '-Department of Energy Computational Science Graduate '
     &     ,'Fellowship (MGM,CDW)'
      write(6,*)
c     --- GPL disclaimers
      write(6,'(1x,a,a)') 
     &     'This program is free software; you can redistribute '
     &     ,'it and/or modify'
      write(6,'(1x,a,a)') 
     &     'it under the terms of the GNU General Public License '
     &     ,'as published by'
      write(6,'(1x,a,a)') 
     &     'the Free Software Foundation; either version 2 of the'
     &     ,' License, or'
      write(6,'(1x,a)') '(at your option) any later version.'
      write(6,'(1x,a,a)') 
     &     'This program is distributed in the hope that it will'
     &     ,' be useful,'
      write(6,'(1x,a,a)') 
     &     'but WITHOUT ANY WARRANTY; without even the implied '
     &     ,'warranty of'
      write(6,'(1x,a,a)') 
     &     'MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.'
     &     ,' See the'
      write(6,'(1x,a)') 'GNU General Public License for more details.'
      write(6,*)
      write(6,'(1x,a,a)') 
     &     'You should have received a copy of the GNU General '
     &     ,'Public License'
      write(6,'(1x,a,a)') 
     &     'along with this program; if not, write to the Free '
     &     ,'Software'
      write(6,'(1x,a,a)') 
     &     'Foundation, Inc., 59 Temple Place, Suite 330, Boston,'
     &     ,' MA  02111-1307  USA'
      write(6,*)

      return
      end
