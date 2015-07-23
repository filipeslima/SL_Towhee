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
      subroutine twh_rwcbmc( lfailure,rwflag, rwfile)
c     ******************************************************************
c     * reads in or writes out the variables related to the settings of*
c     * the configurational-bias moves.                                *
c     * rwflag is either READ_FLAG or WRITE_FLAG                       *
c     *                                                                *
c     * rewritten from writetowhee and readtowhee on 07-29-2005        *
c     *           by M.G. Martin                                       *
c     * last modified 07-18-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_CHECK_LABEL
#define FUNCTION_GET_VIBRANG
#define FUNCTION_READ_LABELED_FLOAT
#define FUNCTION_READ_LABELED_INTEGER
#define FUNCTION_READ_LABELED_LOGICAL
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer rwflag
      integer rwfile
      logical lfailure
c     --- local variables
c     --- logical scalars
      logical lfound,lcreatemap
c     --- character string scalars
      character*15 idname,rwformat
      character*30 ctemp,cbmc_setting_style
      character*30 cbmc_bend_generation,cbmc_bond_generation
      character*30 cbmc_fit_strategy
      character*30 cbmc_dihedral_generation,cbmc_analysis
      character*50 goodstring,twobbias_style,twobbias_init
      character*50 threebbias_init,threebbias_style
      character*60 cbmc_formulation
c     --- character string arrays
      character*30 ctemparray
      dimension ctemparray(NTMAX)
c     --- integer scalars
      integer imolty,ibox,count,ivalue,mapmolty,potentialstyle
      integer nmolty,cubex,cubey,cubez,numboxes,gencode
      integer cbmc_form_code
c     --- integer arrays
      integer cube
      dimension cube(3)
      integer iarray
      dimension iarray(NTMAX)
c     --- double precision scalars
      double precision onepi,dvalue,radtodeg
c     --- double precision arrays
      double precision dtemp
      dimension dtemp(2)

      idname = 'rwcbmc'

c     --- retrieve constants
      call twh_constant_radtodeg(GLB_GET,radtodeg)
      call twh_nmolty(GLB_GET,nmolty)
      call twh_potentialstyle(GLB_GET,potentialstyle)
      call twh_numboxes(GLB_GET,numboxes)

c     --- make sure we have a valid rwflag
      if ( rwflag .ne. READ_FLAG .and. rwflag .ne. WRITE_FLAG ) then
         write(6,*) 'RWCBMC: invalid rwflag'
         write(6,*) rwflag
         lfailure = .true.
         return
      endif

      goodstring = 'cbmc_analysis'
      rwformat = '(a13)'
      if ( rwflag .eq. READ_FLAG ) then
         if ( twh_check_label(lfailure,rwfile,goodstring) ) then
            call twh_read_labeled_string(lfailure,rwfile,goodstring
     &           ,.false.,idname,cbmc_analysis)
            if ( lfailure ) return
         else
c           --- if not specified then default
            cbmc_analysis = 'normal'
         endif
         write(6,'(a14,1x,a30)') 'cbmc_analysis:',cbmc_analysis
         call twh_cbmc_analysis(GLB_SET,cbmc_analysis)
c        --- error check
         if ( cbmc_analysis .ne. 'normal'
     &        .and. cbmc_analysis .ne. 'full' ) then
            write(6,*) 'RWCBMC: unknown cbmc_analysis'
            lfailure = .true.
            return
         endif
      else
         write(rwfile,rwformat) goodstring
         call twh_cbmc_analysis(GLB_GET,cbmc_analysis)
         write(rwfile,'(a1,a30,a1)') char(39),cbmc_analysis,char(39)
      endif

c     --- cbmc formulation
      goodstring = 'cbmc_formulation'
      rwformat = '(a22)'
      if ( rwflag .eq. READ_FLAG ) then
         if ( twh_check_label(lfailure,rwfile,goodstring) ) then
            call twh_read_labeled_string(lfailure,rwfile,goodstring
     &           , .false. ,idname,cbmc_formulation)
            if ( lfailure ) return
         else
c           --- if not specified then set default
            cbmc_formulation = 'Martin and Frischknecht 2006'
         endif

         write(6,'(a17,1x,a60)') 'cbmc_formulation:',cbmc_formulation
         if ( cbmc_formulation .eq. 
     &        'Martin and Siepmann 1999 + Martin and Thompson 2004'
     &        ) then
            write(6,*)'Original Coupled-Decoupled formulation from'
            write(6,*) 'M.G. Martin; J.I. Siepmann; J. Phys. Chem. B'
     &           ,'103 2977-2980 (1999)'
            write(6,*) 'with a decoupled bond length selection from'
            write(6,*) 'M.G. Martin; A.P. Thompson; Fluid Phase'
     &           ,' Equilib. 217 105-110 (2004)'
         elseif ( cbmc_formulation .eq. 'Martin and Frischknecht 2006'
     &           ) then
            write(6,*) 'Coupled to pre-nonbond formulation from'
            write(6,*) 'M.G. Martin; A.L. Frischknecht; Mol. Phys.'
     &           ,' 104 2439-2456 (2006)'
         else
            write(6,*) 'RWCBMC: unknown cbmc_formulation'
            write(6,*) cbmc_formulation
            lfailure = .true.
            return
         endif
c        --- if we passed the safety checks then set the value
         call twh_cbmc_formulation(GLB_SET,cbmc_formulation)
      else
         write(rwfile,rwformat) goodstring
         call twh_cbmc_formulation(GLB_GET
     &        ,cbmc_formulation)
         write(rwfile,'(a1,a60,a1)') char(39)
     &        ,cbmc_formulation,char(39)
      endif
c     --- retrieve the cbmc formulation code for use later
      call twh_cbmc_form_code(GLB_GET,cbmc_form_code)

c     --- cbmc setting style
      goodstring = 'cbmc_setting_style'
      rwformat = '(a18)'
      if ( rwflag .eq. READ_FLAG ) then
         if ( twh_check_label(lfailure,rwfile,goodstring) ) then
            call twh_read_labeled_string(lfailure,rwfile,goodstring
     &           , .false. ,idname,cbmc_setting_style)
            if ( lfailure ) return
         else
c           --- if not specified then default
            cbmc_setting_style = 'Martin and Frischknecht'
         endif
         write(6,'(a19,1x,a30)') 'cbmc_setting_style:'
     &        ,cbmc_setting_style
         call twh_cbmc_setting_style(GLB_SET,cbmc_setting_style)
      else
         write(rwfile,rwformat) goodstring
         call twh_cbmc_setting_style(GLB_GET,cbmc_setting_style)
         write(rwfile,'(a1,a30,a1)') char(39),cbmc_setting_style
     &        ,char(39)
      endif

      if ( cbmc_setting_style .eq. 'explicit' ) then
         goodstring = 'cbmc_fit_strategy'
         rwformat = '(a17)'
         if ( rwflag .eq. READ_FLAG ) then
            if ( twh_check_label(lfailure,rwfile,goodstring) ) then
               call twh_read_labeled_string(lfailure,rwfile,goodstring
     &              , .false. ,idname,cbmc_fit_strategy)
               if ( lfailure ) return
            else
c              --- if not specified then default
               cbmc_fit_strategy = 'equilibrium'
            endif
            write(6,'(a18,1x,a30)') 'cbmc_fit_strategy:'
     &           ,cbmc_fit_strategy
            call twh_cbmc_fit_strategy(GLB_SET,cbmc_fit_strategy)
         else
            write(rwfile,rwformat) goodstring
            call twh_cbmc_fit_strategy(GLB_GET,cbmc_fit_strategy)
            write(rwfile,'(a1,a30,a1)') char(39),cbmc_fit_strategy
     &           ,char(39)
         endif

         goodstring = 'cbmc_nb_one_generation'
         rwformat = '(a23)'
         count = nmolty
         if ( rwflag .eq. READ_FLAG ) then
            ibox = 1
            call twh_read_labeled_string_array(lfailure,rwfile
     &           ,goodstring,count,.true.,idname,ctemparray)
            if ( lfailure ) return
            do imolty = 1,nmolty
               call twh_cbmc_nb_one_generation(GLB_SET,ibox,imolty
     &              ,ctemparray(imolty))
            enddo
            do ibox = 2,numboxes
               call twh_read_string_array(lfailure,rwfile,goodstring
     &              ,count,.true.,ctemparray)
               if ( lfailure ) return
               do imolty = 1,nmolty
                  call twh_cbmc_nb_one_generation(GLB_SET,ibox,imolty
     &                 ,ctemparray(imolty))
               enddo
            enddo
c           --- set all of the gencodes
            do ibox = 1,numboxes
               do imolty = 1,nmolty
                  call twh_cbmc_nb_one_generation(GLB_GET,ibox,imolty
     &                 ,ctemp)
                  if ( ctemp .eq. 'uniform' ) then
                     gencode = DIST_UNIFORM
                  elseif ( ctemp .eq. 'energy bias' ) then
                     gencode = DIST_ENERGY_BIAS
                  else
                     write(6,*) 'RWCBMC: unknown cbmc_nb_one_generation'
                     write(6,*) 'ibox,imolty,gen:',ibox,imolty,ctemp
                     lfailure = .true.
                     return
                  endif
                  call twh_cbmc_nb_one_gencode(GLB_SET,ibox,imolty
     &                 ,gencode)
               enddo
            enddo
         else
            write(rwfile,rwformat) goodstring
            do ibox = 1,numboxes
               do imolty = 1,nmolty
                  call twh_cbmc_nb_one_generation(GLB_GET,ibox,imolty
     &                 ,ctemparray(imolty))
               enddo
               write(rwfile,91)
     &              (char(39),ctemparray(imolty),char(39)
     &              ,imolty=1,nmolty)
 91            format(NTMAX (a1,a30,a1,1x))
            enddo
         endif

         do ibox = 1,numboxes
            lfound = .false.
            do imolty = 1,nmolty
               call twh_cbmc_nb_one_generation(GLB_GET,ibox,imolty
     &              ,ctemp)
               if ( ctemp .eq. 'energy bias' ) lfound = .true.
            enddo
            if ( lfound ) then
c              --- get all the variables related with energy biasing  
               if ( rwflag .eq. READ_FLAG ) then
                  write(6,*) 'energy biasing requested'
c                 --- check potentialstyle
                  if ( potentialstyle .ne. POT_INTERNAL ) then
                     write(6,*)'RWCBMC: potential style must be'
     &                    ,' internal' 
                     write(6,*)'when energybiasing requested'
                     lfailure = .true.
                     return
                  endif
               endif
               if ( rwflag .eq. READ_FLAG ) then
c                 --- read in molecule type to be mapped
                  mapmolty = twh_read_labeled_integer(lfailure
     &                 ,rwfile,'mapmolty',.true.,idname)
                  if ( lfailure ) return
                  call twh_mapmolty(GLB_SET,mapmolty)
c                 --- error check mapmolty
                  if ( mapmolty .lt. 1 .or. mapmolty .gt. nmolty ) then
                     write(6,*)'RWCBMC: invalid mapmolty' 
                     write(6,*)'must be between 1 and nmolty'
                     lfailure = .true.
                     return
                  endif
               else
c                 --- output mapmolty
                  write(rwfile,'(a8)') 'mapmolty'
                  call twh_mapmolty(GLB_GET,mapmolty)
                  write(rwfile,*) mapmolty
               endif

               if ( rwflag .eq. READ_FLAG ) then
c                 --- read in lcreatemap logical
                  lcreatemap = twh_read_labeled_logical(lfailure,rwfile
     &                 ,'lcreatemap',.true.,idname)
                  if ( lfailure ) return
                  call twh_lcreatemap(GLB_SET,lcreatemap)
               else
c                 --- output createmap logical
                  write(rwfile,'(a10)') 'lcreatemap'
                  call twh_lcreatemap(GLB_GET,lcreatemap)
                  write(rwfile,*) lcreatemap
               endif

               if ( lcreatemap ) then
                  if ( rwflag .eq. READ_FLAG ) then
c                    --- read number of cubelets in each direction
                     call twh_read_label(lfailure,rwfile
     &                    ,'cubex cubey cubez',.false.,idname)
                     if ( lfailure ) return
                     call twh_read_integer_array(lfailure,rwfile
     &                    ,'cubex cubey cubez',3,.false., cube)
                     if ( lfailure ) return
                     write(6,'(1x,a18,3(i10,1x))') 
     &                    ,'cubex,cubey,cubez:'
     &                    ,cube(1),cube(2),cube(3)
                     cubex=cube(1)
                     cubey=cube(2)
                     cubez=cube(3)
                     call twh_cubelet(GLB_SET,cubex,cubey,cubez)
c                    --- error check
                     if ( cubex .lt. 0 ) then
                        write(6,*) 'RWCBMC: cubex must be non-negative'
                        lfailure = .true.
                        return
                     endif
                     if ( cubey .lt. 0 ) then
                        write(6,*) 'RWCBMC: cubey must be non-negative'
                        lfailure = .true.
                        return
                     endif
                     if ( cubez .lt. 0 ) then
                        write(6,*) 'RWCBMC: cubez must be non-negative'
                        lfailure = .true.
                        return
                     endif
                  else
                     write(rwfile,'(a17)')'cubex cubey cubez'
                     call twh_cubelet(GLB_GET,cubex,cubey,cubez)
                     write(rwfile,*)cubex,cubey,cubez
                  endif
               endif
            endif
         enddo

         goodstring = 'nch_nb_one'
         rwformat = '(a10)'
         if ( rwflag .eq. READ_FLAG ) then
            call twh_read_labeled_integer_array(lfailure,rwfile
     &           ,goodstring,nmolty,.true.,idname,iarray)
            if ( lfailure ) return
c           --- set and error check
            do imolty = 1,nmolty
               call twh_nch_nb_one(GLB_SET,imolty,iarray(imolty))
               if ( iarray(imolty) .le. 0 ) then
                  write(6,*) 'RWCBMC: all values of nch_nb_one must be'
     &                 ,' positive'
                  lfailure = .true.
                  return
               endif
            enddo
         else
            write(rwfile,rwformat) goodstring
            do imolty=1,nmolty
               call twh_nch_nb_one(GLB_GET,imolty,iarray(imolty))
            enddo
            write(rwfile,*) (iarray(imolty),imolty=1,nmolty)
         endif

         goodstring = 'nch_nb'
         rwformat = '(a6)'
         if ( rwflag .eq. READ_FLAG ) then
            call twh_read_labeled_integer_array(lfailure,rwfile
     &           ,goodstring,nmolty, .true., idname, iarray)
            if ( lfailure ) return
c           --- set and error check
            do imolty = 1,nmolty
               call twh_nch_nb(GLB_SET,imolty,iarray(imolty))
               if ( iarray(imolty) .le. 0 ) then
                  write(6,*) 'RWCBMC: all values of nch_nb must be'
     &                 ,' positive'
                  lfailure = .true.
                  return
               endif
            enddo
         else
            write(rwfile,rwformat) goodstring
            do imolty = 1,nmolty
               call twh_nch_nb(GLB_GET,imolty,iarray(imolty))
            enddo
            write(rwfile,*) (iarray(imolty),imolty=1,nmolty)
         endif
         
         if ( cbmc_form_code .eq. CBMC_FORM_MF2006 )then
            goodstring = 'nch_pre_nb'
            rwformat = '(a10)'
            if ( rwflag .eq. READ_FLAG ) then
               call twh_read_labeled_integer_array(lfailure,rwfile
     &              ,goodstring,nmolty, .true., idname, iarray)
               if ( lfailure ) return
c              --- set and error check
               do imolty = 1,nmolty
                  call twh_nch_pre_nb(GLB_SET,imolty,iarray(imolty))
                  if ( iarray(imolty) .le. 0 ) then
                     write(6,*) 'RWCBMC: all values of nch_pre_nb must'
     &                    ,' be positive'
                     lfailure = .true.
                     return
                  endif
               enddo
            else
               write(rwfile,rwformat) goodstring
               do imolty = 1,nmolty
                  call twh_nch_pre_nb(GLB_GET,imolty,iarray(imolty))
               enddo
               write(rwfile,*) (iarray(imolty),imolty=1,nmolty)
            endif
         endif

c        --- torsion cb prob distribution style
         goodstring = 'cbmc_dihedral_generation'
         rwformat = '(a24)'
         if ( rwflag .eq. READ_FLAG ) then
            call twh_read_labeled_string(lfailure,rwfile,goodstring
     &           , .false. ,idname,cbmc_dihedral_generation)
            if ( lfailure ) return
            call twh_cbmc_dihedral_generation(GLB_SET
     &           ,cbmc_dihedral_generation)
            write(6,'(a25,1x,a30)') 'cbmc_dihedral_generation:'
     &           ,cbmc_dihedral_generation
            if ( cbmc_dihedral_generation .eq. 'ideal' ) then
               call twh_cbmc_dihedgen_code(GLB_SET,CBMC_IDEAL)
               write(6,*) 'Use true distributions to generate torsions'
     &              ,' in config-bias'
            elseif ( cbmc_dihedral_generation .eq. 'global gaussian'
     &              ) then
               call twh_cbmc_dihedgen_code(GLB_SET,CBMC_GLOBAL_GAUSSIAN)
               write(6,*) 'Use a global gaussian distribution to'
     &              ,' generate torsions in config-bias'
               write(6,*)'and correct this in the rosenbluth weights'
            elseif ( cbmc_dihedral_generation .eq. 'autofit gaussian'
     &              ) then
               call twh_cbmc_dihedgen_code(GLB_SET
     &              ,CBMC_AUTOFIT_GAUSSIAN)
               write(6,*) 'Use individual gaussian distributions to'
     &              ,' generate torsions in config-bias'
               write(6,*)'and correct this in the rosenbluth weights'
            elseif ( cbmc_dihedral_generation .eq. 
     &              'ideal + autofit gaussian' ) then
               call twh_cbmc_dihedgen_code(GLB_SET
     &              ,CBMC_IDEAL_PLUS_AUTOFIT_GAUSSIAN)
               write(6,*) 'Use a linear combination of the ideal'
     &              ,' distribution plus'
               write(6,*) 'automatically fit individual gaussian'
     &              ,' distrubtions to generate torsions'
               write(6,*)'and correct this in the rosenbluth weights'
            else
               write(6,*) 'RWCBMC: invalid cbmc_dihedral_generation'
               lfailure = .true.
               return
            endif
         else
            write(rwfile,rwformat) goodstring
            call twh_cbmc_dihedral_generation(GLB_GET
     &           ,cbmc_dihedral_generation)
            write(rwfile,'(a1,a30,a1)') char(39)
     &           ,cbmc_dihedral_generation,char(39)
         endif

         if ( cbmc_dihedral_generation .eq. 'global gaussian' ) then
c          --- standard deviations need for this type
            goodstring = 'sdevtor'
            rwformat = '(a7)'
            if ( rwflag .eq. READ_FLAG ) then
               dvalue = twh_read_labeled_float(lfailure,rwfile
     &              ,goodstring, .false.,idname)
               if ( lfailure ) return
               call twh_sdevtor(GLB_SET,dvalue)
               write(6,'(a8,f10.5)') 'sdevtor:',dvalue
c              --- error check
               if ( dvalue .lt. 0.0d0 ) then
                  write(6,*) 'RWCBMC: sdevtor must be non-negative'
                  lfailure = .true.
                  return
               endif
               if ( dvalue .gt. 360.0d0 ) then
                  write(6,*) 'RWCBMC: sdvetor must be less than 360.0'
                  lfailure = .true.
                  return
               endif
            else
               write(rwfile,rwformat) goodstring
c              --- convert back to degrees
               call twh_sdevtor(GLB_GET,dvalue)
               dvalue = dvalue*radtodeg
               write(rwfile,*) dvalue
            endif
         endif

         if ( cbmc_dihedral_generation .eq. 'autofit gaussian' .or.
     &        cbmc_dihedral_generation .eq. 'ideal + autofit gaussian'
     &        ) then
c           --- dihedral peak weight style
            goodstring = 'dihedral_peak_weight_style'
            rwformat = '(a26)'
            if ( rwflag .eq. READ_FLAG ) then
               call twh_read_labeled_string(lfailure,rwfile,goodstring
     &              ,.false.,idname,ctemp)
               if ( lfailure ) return
               call twh_dihedral_peak_weight_style(GLB_SET,ctemp)
            else
               write(rwfile,rwformat) goodstring
               call twh_dihedral_peak_weight_style(GLB_GET,ctemp)
               write(rwfile,'(a1,a30,a1)') char(39),ctemp,char(39)
            endif

c           --- standard deviation multiplier
            goodstring = 'dihedral_sdev_multiplier'
            rwformat = '(a24)'
            if ( rwflag .eq. READ_FLAG ) then
               dvalue = twh_read_labeled_float(lfailure,rwfile
     &              ,goodstring, .false.,idname)
               if ( lfailure ) return
               call twh_dihedral_sdev_multiplier(GLB_SET,dvalue)
               write(6,'(a25,f10.5)') 'dihedral_sdev_multiplier:'
     &              ,dvalue
c              --- error check
               if ( dvalue .le. 0.0d0 ) then
                  write(6,*) 'RWCBMC: dihedral_sdev_multiplier must be'
     &                 ,' positive'
                  lfailure = .true.
                  return
               endif
            else
               write(rwfile,rwformat) goodstring
               call twh_dihedral_sdev_multiplier(GLB_GET,dvalue)
               write(rwfile,*) dvalue
            endif
         endif

         if ( cbmc_dihedral_generation .eq. 'ideal + autofit gaussian'
     &        ) then
c           --- ideal linear combination fraction
            goodstring = 'dihedral_ideal_fraction'
            rwformat = '(a23)'
            if ( rwflag .eq. READ_FLAG ) then
               dvalue = twh_read_labeled_float(lfailure,rwfile
     &              ,goodstring, .false.,idname)
               if ( lfailure ) return
               call twh_dihedral_ideal_fraction(GLB_SET,dvalue)
               write(6,'(a24,1x,f6.4)') 'dihedral_ideal_fraction:'
     &              ,dvalue
c              --- error check
               if ( dvalue .lt. 0.0d0 .or. dvalue .gt. 1.0d0 ) then
                  write(6,*) 'RWCBMC: dihedral_ideal_fraction must be'
     &                 ,' between 0.0 and 1.0 (inclusive)'
                  lfailure = .true.
                  return
               endif
            else
               write(rwfile,rwformat) goodstring
               call twh_dihedral_ideal_fraction(GLB_GET,dvalue)
               write(rwfile,*) dvalue
            endif
         endif

c        --- information for CBMC inner torsion growth
         goodstring = 'nch_tor'
         rwformat = '(a7)'
         if ( rwflag .eq. READ_FLAG ) then
            call twh_read_labeled_integer_array(lfailure,rwfile
     &           ,goodstring,nmolty, .true., idname, iarray)
            if ( lfailure ) return
c           --- set nch_tor
            do imolty = 1,nmolty
               ivalue = iarray(imolty)
               call twh_nch_tor(GLB_SET,imolty,ivalue)
            enddo
c           --- error check
            do imolty = 1,nmolty
               call twh_nch_tor(GLB_GET,imolty,ivalue)
               if ( ivalue .le. 0 ) then
                  write(6,*) 'RWCBMC: all values of nch_tor must be'
     &                 ,' positive'
                  lfailure = .true.
                  return
               endif
            enddo
         else
            write(rwfile,rwformat) goodstring
            do imolty = 1,nmolty
               call twh_nch_tor(GLB_GET,imolty,iarray(imolty))
            enddo
            write(rwfile,*) (iarray(imolty),imolty=1,nmolty)
         endif

c        --- information for CBMC torsion growth for fixed endpoint
         goodstring = 'nch_tor_connect'
         rwformat = '(a15)'
         if ( rwflag .eq. READ_FLAG ) then
            call twh_read_labeled_integer_array(lfailure,rwfile
     &           ,goodstring,nmolty, .true., idname, iarray)
            if ( lfailure ) return
c           --- set and error check
            do imolty = 1,nmolty
               call twh_nch_tor_connect(GLB_SET,imolty,iarray(imolty))
               if ( iarray(imolty) .le. 0 ) then
                  write(6,*) 'RWCBMC: nch_tor_connect must be positive'
                  lfailure = .true.
                  return
               endif
            enddo
         else
            write(rwfile,rwformat) goodstring
            do imolty = 1,nmolty
               call twh_nch_tor_connect(GLB_GET,imolty,iarray(imolty))
            enddo
            write(rwfile,*) (iarray(imolty),imolty=1,nmolty)
         endif

c        --- bend cb prob distribution style
         goodstring = 'cbmc_bend_generation'
         rwformat = '(a20)'
         if ( rwflag .eq. READ_FLAG ) then
            call twh_read_labeled_string(lfailure,rwfile,goodstring
     &           , .false. ,idname,cbmc_bend_generation)
            if ( lfailure ) return
            call twh_cbmc_bend_generation(GLB_SET,cbmc_bend_generation)
            write(6,'(a21,1x,a30)') 'cbmc_bend_generation:'
     &           ,cbmc_bend_generation
            if ( cbmc_bend_generation .eq. 'ideal' ) then
               call twh_cbmc_bendgen_code(GLB_SET,CBMC_IDEAL)
               write(6,*) 'Use ideal distributions to generate angles'
     &              ,' in configurational-bias'
            elseif ( cbmc_bend_generation .eq. 'global gaussian' ) then
               call twh_cbmc_bendgen_code(GLB_SET,CBMC_GLOBAL_GAUSSIAN)
               write(6,*) 'Use a globally defined gaussian distribution'
     &              ,' in configurational-bias'
               write(6,*) 'and correct this in the rosenbluth weights'
               write(6,*) 'Bend style A mean = equilibrium angle'
               write(6,*) 'Bend style A standard deviation ='
     &              ,' sdev_bend_a'
               write(6,*) 'Bend style B has multiple regions'
               write(6,*) 'Bend style B mean = function of central'
     &              ,' atom bondpatt'
               write(6,*) 'Bend style B standard deviation ='
     &              ,' sdev_bend_b'
            elseif ( cbmc_bend_generation .eq. 'autofit gaussian' ) then
               call twh_cbmc_bendgen_code(GLB_SET,CBMC_AUTOFIT_GAUSSIAN)
               write(6,*) 'Use an individually defined gaussian'
     &              ,' distribution in configurational-bias'
               write(6,*) 'and correct this in the rosenbluth weights'
               write(6,*) 'Bend style A mean and standard deviation fit'
     &              ,' to cos(theta) exp(-beta U_bend)'
               write(6,*) 'Bend style B has multiple regions'
               write(6,*) 'Bend style B means and standard deviations'
     &              ,' fit exp(-beta U_bend)'
            elseif ( cbmc_bend_generation 
     &              .eq. 'ideal + autofit gaussian' ) then
               call twh_cbmc_bendgen_code(GLB_SET
     &              ,CBMC_IDEAL_PLUS_AUTOFIT_GAUSSIAN)
               write(6,*) 'Use a linear combination of the ideal'
     &              ,' distributions' 
               write(6,*) 'and individually defined gaussian'
     &              ,' distribution in configurational-bias'
               write(6,*) 'and correct this in the rosenbluth weights'
               write(6,*) 'Bend style A mean and standard deviation fit'
     &              ,' to cos(theta) exp(-beta U_bend)'
               write(6,*) 'Bend style B has multiple regions'
               write(6,*) 'Bend style B means and standard deviations'
     &              ,' fit exp(-beta U_bend)'
            else
               write(6,*) 'RWCBMC: invalid cbmc_bend_generation'
               lfailure = .true.
               return
            endif
         else
            write(rwfile,rwformat) goodstring
            call twh_cbmc_bend_generation(GLB_GET,cbmc_bend_generation)
            write(rwfile,'(a1,a30,a1)') char(39),cbmc_bend_generation
     &           ,char(39)
         endif

         if ( cbmc_bend_generation .eq. 'global gaussian' ) then
c           --- standard deviation for use with this style
            goodstring = 'sdevbena'
            rwformat = '(a8)'
            if ( rwflag .eq. READ_FLAG ) then
               dvalue =  twh_read_labeled_float(lfailure,rwfile
     &              ,goodstring, .false.,idname)
               if ( lfailure ) return
               call twh_sdevbena(GLB_SET,dvalue)
               write(6,'(a9,f10.5)') 'sdevbena:',dvalue
c              --- error check
               if ( dvalue .lt. 0.0d0 ) then
                  write(6,*) 'RWCBMC: sdevbena must be non-negative'
                  lfailure = .true.
                  return
               endif
               if ( dvalue .gt. 180.0d0 ) then
                  write(6,*) 'RWCBMC: sdevbena must be less than 180.0'
                  lfailure = .true.
                  return
               endif
            else
               write(rwfile,rwformat) goodstring
c              --- convert back to degrees
               call twh_sdevbena(GLB_GET,dvalue)
               write(rwfile,*) dvalue*radtodeg
            endif
c           --- standard deviation for use with this style
            goodstring = 'sdevbenb'
            rwformat = '(a8)'
            if ( rwflag .eq. READ_FLAG ) then
               dvalue =  twh_read_labeled_float(lfailure,rwfile
     &              ,goodstring, .false.,idname)
               if ( lfailure ) return
               call twh_sdevbenb(GLB_SET,dvalue)
               write(6,'(a9,f10.5)') 'sdevbenb:',dvalue
c              --- error check
               if ( dvalue .lt. 0.0d0 ) then
                  write(6,*) 'RWCBMC: sdevbenb must be non-negative'
                  lfailure = .true.
                  return
               endif
               if ( dvalue .gt. 360.0 ) then
                  write(6,*) 'RWCBMC: sdevbenb must be less than 360.0'
                  lfailure = .true.
                  return
               endif
            else
               write(rwfile,rwformat) goodstring
c              --- convert back to degrees
               call twh_sdevbenb(GLB_GET,dvalue)
               write(rwfile,*) dvalue*radtodeg
            endif

         elseif ( cbmc_bend_generation .eq. 'autofit gaussian' .or.
     &           cbmc_bend_generation .eq. 'ideal + autofit gaussian'
     &           ) then
c           --- need standard deviation multiplier
            goodstring = 'bend_a_sdev_multiplier'
            rwformat = '(a22)'
            if ( rwflag .eq. READ_FLAG ) then
               dvalue = twh_read_labeled_float(lfailure,rwfile
     &              ,goodstring, .false.,idname)
               if ( lfailure ) return
               call twh_bend_a_sdev_multiplier(GLB_SET,dvalue)
               write(6,'(a23,1x,f10.5)') 'bend_a_sdev_multiplier:'
     &              ,dvalue
               if ( dvalue .le. 0.0d0 ) then
                  write(6,*) 'RWCBMC: bend_a_sdev_multiplier must be'
     &                 ,' positive'
                  lfailure = .true.
                  return
               endif
            else
               write(rwfile,rwformat) goodstring
               call twh_bend_a_sdev_multiplier(GLB_GET,dvalue)
               write(rwfile,*) dvalue
            endif

c           --- need standard deviation multiplier
            goodstring = 'bend_b_sdev_multiplier'
            rwformat = '(a22)'
            if ( rwflag .eq. READ_FLAG ) then
               dvalue = twh_read_labeled_float(lfailure,rwfile
     &              ,goodstring, .false.,idname)
               if ( lfailure ) return
               call twh_bend_b_sdev_multiplier(GLB_SET,dvalue)
               write(6,'(a23,1x,f10.5)') 'bend_b_sdev_multiplier:'
     &              ,dvalue
               if ( dvalue .le. 0.0d0 ) then
                  write(6,*) 'RWCBMC: bend_b_sdev_multiplier must be'
     &                 ,' positive'
                  lfailure = .true.
                  return
               endif
            else
               write(rwfile,rwformat) goodstring
               call twh_bend_b_sdev_multiplier(GLB_GET,dvalue)
               write(rwfile,*) dvalue
            endif
         endif

         if ( cbmc_bend_generation .eq. 'ideal + autofit gaussian') then
c           --- need ideal fractions for bend a and b
            goodstring = 'bend_a_ideal_fraction'
            rwformat = '(a21)'
            if ( rwflag .eq. READ_FLAG ) then
               dvalue = twh_read_labeled_float(lfailure,rwfile
     &              ,goodstring, .false.,idname)
               if ( lfailure ) return
               call twh_bend_a_ideal_fraction(GLB_SET,dvalue)
               write(6,'(a22,1x,f10.5)') 'bend_a_ideal_fraction:'
     &              ,dvalue
               if ( dvalue .lt. 0.0d0 .or. dvalue .gt. 1.0d0 ) then
                  write(6,*) 'RWCBMC: bend_a_ideal_fraction must be'
     &                 ,' in the range [0.0,1.0] inclusive'
                  lfailure = .true.
                  return
               endif
            else
               write(rwfile,rwformat) goodstring
               call twh_bend_a_ideal_fraction(GLB_GET,dvalue)
               write(rwfile,*) dvalue
            endif
            
            goodstring = 'bend_b_ideal_fraction'
            rwformat = '(a21)'
            if ( rwflag .eq. READ_FLAG ) then
               dvalue = twh_read_labeled_float(lfailure,rwfile
     &              ,goodstring, .false.,idname)
               if ( lfailure ) return
               call twh_bend_b_ideal_fraction(GLB_SET,dvalue)
               write(6,'(a22,1x,f10.5)') 'bend_b_ideal_fraction:'
     &              ,dvalue
               if ( dvalue .lt. 0.0d0 .or. dvalue .gt. 1.0d0 ) then
                  write(6,*) 'RWCBMC: bend_b_ideal_fraction must be'
     &                 ,' in the range [0.0,1.0] inclusive'
                  lfailure = .true.
                  return
               endif
            else
               write(rwfile,rwformat) goodstring
               call twh_bend_b_ideal_fraction(GLB_GET,dvalue)
               write(rwfile,*) dvalue
            endif
         endif

c        --- information for CBMC bond angle growth
         goodstring = 'nch_bend_a'
         rwformat = '(a10)'
         if ( rwflag .eq. READ_FLAG ) then
            call twh_read_labeled_integer_array(lfailure,rwfile
     &           ,goodstring, nmolty, .true., idname, iarray)
            if ( lfailure ) return
c           --- set and error check
            do imolty = 1,nmolty
               call twh_nch_bend_a(GLB_SET,imolty,iarray(imolty))
               if ( iarray(imolty) .le. 0 ) then
                  write(6,*) 'RWCBMC: nch_bend_a must be positive'
                  lfailure = .true.
                  return
               endif
            enddo
         else
            write(rwfile,rwformat) goodstring
            do imolty = 1,nmolty
               call twh_nch_bend_a(GLB_GET,imolty,iarray(imolty))
            enddo
            write(rwfile,*) (iarray(imolty),imolty=1,nmolty)
         endif

         goodstring = 'nch_bend_b'
         rwformat = '(a10)'
         if ( rwflag .eq. READ_FLAG ) then
            call twh_read_labeled_integer_array(lfailure,rwfile
     &           ,goodstring,nmolty, .true., idname, iarray)
            if ( lfailure ) return
c           --- set anderror check
            do imolty = 1,nmolty
               call twh_nch_bend_b(GLB_SET,imolty,iarray(imolty))
               if ( iarray(imolty) .le. 0 ) then
                  write(6,*) 'RWCBMC: nch_bend_b must be positive'
                  lfailure = .true.
                  return
               endif
            enddo
         else
            write(rwfile,rwformat) goodstring
            do imolty = 1,nmolty
               call twh_nch_bend_b(GLB_GET,imolty,iarray(imolty))
            enddo
            write(rwfile,*) (iarray(imolty),imolty=1,nmolty)
         endif

c        --- vibration cb prob distribution style
         goodstring = 'cbmc_bond_generation'
         rwformat = '(a20)'
         if ( rwflag .eq. READ_FLAG ) then
            call twh_read_labeled_string(lfailure,rwfile,goodstring
     &           , .false. ,idname,cbmc_bond_generation)
            if ( lfailure ) return
            call twh_cbmc_bond_generation(GLB_SET,cbmc_bond_generation)
            write(6,'(a21,1x,a30)') 'cbmc_bond_generation:'
     &           ,cbmc_bond_generation

            if ( cbmc_bond_generation .eq. 'r^2 with bounds' ) then
               call twh_cbmc_bondgen_code(GLB_SET,CBMC_R_SQ_WITH_BOUNDS)
               write(6,*) 'Use a bounded r^2 distribution to generate'
     &              ,' vibrations in config-bias'
            elseif ( cbmc_bond_generation .eq. 'global gaussian' ) then
               call twh_cbmc_bondgen_code(GLB_SET,CBMC_GLOBAL_GAUSSIAN)
               write(6,*) 'Use gaussian distribution to generate'
     &              ,' vibrations in config-bias'
               write(6,*)'   and correct this in the rosenbluth weights'
               write(6,*) 'mean = equilibrium bond length'
               write(6,*) 'standard deviation = sdevvib'
            elseif ( cbmc_bond_generation .eq. 'autofit gaussian' ) then
               call twh_cbmc_bondgen_code(GLB_SET,CBMC_AUTOFIT_GAUSSIAN)
               write(6,*) 'Use gaussian distribution to generate'
     &              ,' vibrations in config-bias'
               write(6,*)'   and correct this in the rosenbluth weights'
               write(6,*) 'mean and standard deviations from a gaussian'
     &              ,' to r^2 exp(-beta U_bond)'
            elseif ( cbmc_bond_generation .eq. 'ideal' ) then
               call twh_cbmc_bondgen_code(GLB_SET,CBMC_IDEAL)
               write(6,*) 'Use the ideal r^2 distribution to generate'
     &              ,' vibrations in config-bias'
            else
               write(6,*) 'RWCBMC: unknown cbmc_bond_generation'
               lfailure = .true.
               return
            endif
         else
            write(rwfile,rwformat) goodstring
            call twh_cbmc_bond_generation(GLB_GET,cbmc_bond_generation)
            write(rwfile,'(a1,a30,a1)') char(39),cbmc_bond_generation
     &           ,char(39)
         endif
         
         if ( cbmc_bond_generation .eq. 'r^2 with bounds' ) then
c           --- need vibrang for this type
            goodstring = 'vibrang'
            rwformat = '(a7)'
            if ( rwflag .eq. READ_FLAG ) then
               call twh_read_labeled_float_array(lfailure,rwfile
     &              ,goodstring,2, .false., idname, dtemp)
               if ( lfailure ) return
               call twh_vibrang(GLB_SET,1,dtemp(1))
               call twh_vibrang(GLB_SET,2,dtemp(2))
               write(6,'(a8,2(f10.5,1x))') 'vibrang:',twh_get_vibrang(1)
     &              ,twh_get_vibrang(2)
               if ( twh_get_vibrang(2) .le. twh_get_vibrang(1) ) then
                  write(6,*) 'RWCBMC: vibrang(2) must be > vibrang(1)'
                  lfailure = .true.
                  return
               endif
            else
               write(rwfile,rwformat) goodstring
c              --- output these values
               write(rwfile,*) twh_get_vibrang(1),twh_get_vibrang(2)
            endif
         elseif ( cbmc_bond_generation .eq. 'global gaussian' ) then
c           --- need sdevvib for this type
            goodstring = 'sdevvib'
            rwformat = '(a7)'
            if ( rwflag .eq. READ_FLAG ) then
               dvalue =  twh_read_labeled_float(lfailure,rwfile
     &              ,goodstring, .false.,idname)
               if ( lfailure ) return
               call twh_sdevvib(GLB_SET,dvalue)
               write(6,'(a8,f10.5)') 'sdevvib:',dvalue
c              --- error check
               if ( dvalue .lt. 0.0d0 ) then
                  write(6,*) 'RWCBMC: sdevvib must be non-negative'
                  lfailure = .true.
                  return
               endif
            else
               write(rwfile,rwformat) goodstring
               call twh_sdevvib(GLB_GET,dvalue)
               write(rwfile,*) dvalue
            endif
         elseif ( cbmc_bond_generation .eq. 'autofit gaussian' ) then
c           --- need bond_sdev_multiplier for this type
            goodstring = 'bond_sdev_multiplier'
            rwformat = '(a20)'
            if ( rwflag .eq. READ_FLAG ) then
               dvalue =  twh_read_labeled_float(lfailure,rwfile
     &              ,goodstring, .false.,idname)
               if ( lfailure ) return
               call twh_bond_sdev_multiplier(GLB_SET,dvalue)
               write(6,'(a21,f10.5)') 'bond_sdev_multiplier:',dvalue
c              --- error check
               if ( dvalue .le. 0.0d0 ) then
                  write(6,*) 'RWCBMC: bond_sdev_multiplier must be'
     &                 ,' positive'
                  lfailure = .true.
                  return
               endif
            else
               write(rwfile,rwformat) goodstring
               call twh_bond_sdev_multiplier(GLB_GET,dvalue)
               write(rwfile,*) dvalue
            endif
         endif

         goodstring = 'nch_vib'
         rwformat = '(a7)'
         if ( rwflag .eq. READ_FLAG ) then
            call twh_read_labeled_integer_array(lfailure,rwfile
     &           ,goodstring,nmolty, .true., idname, iarray)
            if ( lfailure ) return
c           --- set and Error checking
            do imolty = 1,nmolty
               call twh_nch_vib(GLB_SET,imolty,iarray(imolty))
               if ( iarray(imolty) .le. 0 ) then
                  write(6,*) 'RWCBMC: nch_vib must be positive'
                  lfailure = .true.
                  return
               endif
            enddo
         else
            write(rwfile,rwformat) goodstring
            do imolty = 1,nmolty
               call twh_nch_vib(GLB_GET,imolty,iarray(imolty))
            enddo
            write(rwfile,*) (iarray(imolty),imolty=1,nmolty)
         endif

c        --- fixed endpoint bias styles
         goodstring = 'two_bond_fixed_endpoint_bias_style'
         rwformat = '(a34)'
         if ( rwflag .eq. READ_FLAG ) then
            call twh_read_labeled_string(lfailure,rwfile,goodstring
     &           , .false. ,idname,twobbias_style)
            if ( lfailure ) return
            call twh_two_bond_fe_bias_style(GLB_SET,twobbias_style)
            write(6,'(a34,a1,1x,a50)') goodstring,':',twobbias_style
         else
            write(rwfile,rwformat) goodstring
            call twh_two_bond_fe_bias_style(GLB_GET,twobbias_style)
            write(rwfile,'(a1,a50,a1)') char(39),twobbias_style,char(39)
         endif

         if ( twobbias_style .eq. 'none' ) then
c           --- this style does not require further input or output
         elseif ( twobbias_style .eq.
     &           'analytic Boltzmann using angles' ) then
c           --- this style does not require further input or output
         elseif ( twobbias_style .eq. 
     &           'analytic Boltzmann dihedral energy sum' ) then
c           --- this style does not require further input or output
         elseif ( twobbias_style .eq.
     &           'autofit gaussian' ) then
            goodstring = 'two_bond_bias_sdev_multiplier'
            rwformat = '(a30)'
            if ( rwflag .eq. READ_FLAG ) then
               dvalue =  twh_read_labeled_float(lfailure,rwfile
     &              ,goodstring, .true.,idname)
               if ( lfailure ) return
               call twh_two_bond_bias_sdev_multiplier(GLB_SET,dvalue)
            else
               write(rwfile,rwformat) goodstring
               call twh_two_bond_bias_sdev_multiplier(GLB_GET,dvalue)
               write(rwfile,*) dvalue
            endif
            goodstring = 'two_bond_bias_vibrange'
            rwformat = '(a22)'
            if ( rwflag .eq. READ_FLAG ) then
               call twh_read_labeled_float_array(lfailure,rwfile
     &              ,goodstring,2, .true., idname, dtemp)
               if ( lfailure ) return
               call twh_two_bond_bias_vibrange(GLB_SET,1,dtemp(1))
               call twh_two_bond_bias_vibrange(GLB_SET,2,dtemp(2))
            else
               write(rwfile,rwformat) goodstring
               call twh_two_bond_bias_vibrange(GLB_GET,1,dtemp(1))
               call twh_two_bond_bias_vibrange(GLB_GET,2,dtemp(2))
               write(rwfile,*) dtemp(1),dtemp(2)
            endif
         elseif ( twobbias_style .eq.
     &           'self adapting gaussian using 1-3 distance'
     &           ) then
            goodstring = 'two_bond_bias_sdev_multiplier'
            rwformat = '(a30)'
            if ( rwflag .eq. READ_FLAG ) then
               dvalue = twh_read_labeled_float(lfailure,rwfile
     &              ,goodstring,.true.,idname)
               if ( lfailure ) return
               call twh_two_bond_bias_sdev_multiplier(GLB_SET,dvalue)
            else
               write(rwfile,rwformat) goodstring
               call twh_two_bond_bias_sdev_multiplier(GLB_GET,dvalue)
               write(rwfile,*) dvalue
            endif
c           --- need an update frequency, initial values, and mixing
c           --- ratio
c           --- initial values
            goodstring = 'two_bond_bias_initial_value'
            rwformat = '(a27)'
            if ( rwflag .eq. READ_FLAG ) then
               call twh_read_labeled_string(lfailure,rwfile,goodstring
     &              , .true. ,idname,twobbias_init)
               if ( lfailure ) return
               call twh_two_bond_bias_initial_value(GLB_SET
     &              ,twobbias_init)
            else
               write(rwfile,rwformat) goodstring
               call twh_two_bond_bias_initial_value(GLB_GET
     &              ,twobbias_init)
               write(rwfile,'(a1,a50,a1)') char(39),twobbias_init
     &              ,char(39)
            endif
c           --- compute frequency
            goodstring = 'two_bond_bias_compute_frequency'
            rwformat = '(a31)'
            if ( rwflag .eq. READ_FLAG ) then
               ivalue = twh_read_labeled_integer(lfailure,rwfile
     &              , goodstring, .true.,idname)
               if ( lfailure ) return
               call twh_two_bond_bias_compute_frequency(GLB_SET,ivalue)
            else
               write(rwfile,rwformat) goodstring
               call twh_two_bond_bias_compute_frequency(GLB_GET,ivalue)
               write(rwfile,*) ivalue
            endif
c           --- update frequency
            goodstring = 'two_bond_bias_update_frequency'
            rwformat = '(a30)'
            if ( rwflag .eq. READ_FLAG ) then
               ivalue = twh_read_labeled_integer(lfailure,rwfile
     &              , goodstring, .true.,idname)
               if ( lfailure ) return
               call twh_two_bond_bias_update_frequency(GLB_SET,ivalue)
            else
               write(rwfile,rwformat) goodstring
               call twh_two_bond_bias_update_frequency(GLB_GET,ivalue)
               write(rwfile,*) ivalue
            endif
c           --- old fraction
            goodstring = 'two_bond_bias_old_fraction'
            rwformat = '(a26)'
            if ( rwflag .eq. READ_FLAG ) then
               dvalue = twh_read_labeled_float(lfailure,rwfile
     &              , goodstring, .true.,idname)
               if ( lfailure ) return
               call twh_two_bond_bias_old_fraction(GLB_SET,dvalue)
            else
               write(rwfile,rwformat) goodstring
               call twh_two_bond_bias_old_fraction(GLB_GET,dvalue)
               write(rwfile,*) dvalue
            endif
         else
c           --- unknown style
            lfailure = .true.
            write(6,*) 'RWCBMC: unknown'
     &           ,' two_bond_fe_bias_style'
            write(6,*) twobbias_style
            return
         endif
         
         goodstring = 'three_bond_fixed_endpoint_bias_style'
         rwformat = '(a36)'
         if ( rwflag .eq. READ_FLAG ) then
            call twh_read_labeled_string(lfailure,rwfile,goodstring
     &           ,.false.,idname,threebbias_style)
            if ( lfailure ) return
            call twh_three_bond_fe_bias_style(GLB_SET,threebbias_style)
            write(6,'(a36,a1,1x,a50)') goodstring,':',threebbias_style
         else
            write(rwfile,rwformat) goodstring
            call twh_three_bond_fe_bias_style(GLB_GET,threebbias_style)
            write(rwfile,'(a1,a50,a1)') char(39),threebbias_style
     &           ,char(39)
         endif

         if ( threebbias_style .eq. 'none' ) then
c           --- no other variables required
         elseif ( threebbias_style
     &           .eq. 'analytic using max and min 2-4 distance' ) then
c           --- no other variables required
         elseif ( threebbias_style .eq.
     &           'autofit gaussian using max and min 2-4 distance' 
     &           ) then
            goodstring = 'three_bond_bias_sdev_multiplier'
            rwformat = '(a31)'
            if ( rwflag .eq. READ_FLAG ) then
               dvalue = twh_read_labeled_float(lfailure,rwfile
     &              ,goodstring,.true. ,idname)
               if ( lfailure ) return
               call twh_three_bond_bias_sdev_multiplier(GLB_SET,dvalue)
            else
               write(rwfile,rwformat) goodstring
               call twh_three_bond_bias_sdev_multiplier(GLB_GET,dvalue)
               write(rwfile,*) dvalue
            endif 

         elseif ( threebbias_style .eq. 
     &           'self adapting gaussian using 1-4 distance' ) then
            goodstring = 'three_bond_bias_sdev_multiplier'
            rwformat = '(a31)'
            if ( rwflag .eq. READ_FLAG ) then
               dvalue = twh_read_labeled_float(lfailure,rwfile
     &              ,goodstring,.true.,idname)
               if ( lfailure ) return
               call twh_three_bond_bias_sdev_multiplier(GLB_SET,dvalue)
            else
               write(rwfile,rwformat) goodstring
               call twh_three_bond_bias_sdev_multiplier(GLB_GET,dvalue)
               write(rwfile,*) dvalue
            endif 
c           --- initial value
            goodstring = 'three_bond_bias_initial_value'
            rwformat = '(a30)'
            if ( rwflag .eq. READ_FLAG ) then
               call twh_read_labeled_string(lfailure,rwfile, goodstring
     &              ,.true. ,idname,threebbias_init)
               if ( lfailure ) return
               call twh_three_bond_bias_initial_value(GLB_SET
     &              ,threebbias_init)
            else
               write(rwfile,rwformat) goodstring
               call twh_three_bond_bias_initial_value(GLB_GET
     &              ,threebbias_init)
               write(rwfile,'(a1,a50,a1)') char(39)
     &              ,threebbias_init,char(39)
            endif
c           --- compute frequency
            goodstring = 'three_bond_bias_compute_frequency'
            rwformat = '(a33)'
            if ( rwflag .eq. READ_FLAG ) then
               ivalue = twh_read_labeled_integer(lfailure,rwfile
     &              ,goodstring,.true.,idname)
               if ( lfailure ) return
               call twh_three_bond_bias_compute_freq(GLB_SET,ivalue)
            else
               write(rwfile,rwformat) goodstring
               call twh_three_bond_bias_compute_freq(GLB_GET,ivalue)
               write(rwfile,*) ivalue
            endif
c           --- update frequency
            goodstring = 'three_bond_bias_update_frequency'
            rwformat = '(a33)'
            if ( rwflag .eq. READ_FLAG ) then
               ivalue = twh_read_labeled_integer(lfailure,rwfile
     &              ,goodstring, .true.,idname)
               if ( lfailure ) return
               call twh_three_bond_bias_update_freq(GLB_SET,ivalue)
            else
               write(rwfile,rwformat) goodstring
               call twh_three_bond_bias_update_freq(GLB_GET,ivalue)
               write(rwfile,*) ivalue
            endif
c           --- old fraction
            goodstring = 'three_bond_bias_old_fraction'
            rwformat = '(a28)'
            if ( rwflag .eq. READ_FLAG ) then
               dvalue = twh_read_labeled_float(lfailure,rwfile
     &              , goodstring, .true.,idname)
               if ( lfailure ) return
               call twh_three_bond_bias_old_fraction(GLB_SET,dvalue)
            else
               write(rwfile,rwformat) goodstring
               call twh_three_bond_bias_old_fraction(GLB_GET,dvalue)
               write(rwfile,*) dvalue
            endif
         else
            write(6,*) 'RWCBMC: unknown ',goodstring
            write(6,*) threebbias_style
            lfailure = .true.
            return
         endif
         
      elseif ( cbmc_setting_style .eq. 'Martin and Thompson FPE 2004'
     &        ) then
         if ( rwflag .eq. READ_FLAG ) then
c           --- use the default setting and styles from that paper
            call twh_cbmc_fit_strategy(GLB_SET,'equilibrium')
            call twh_cbmc_nb_one_generation(GLB_INIT,numboxes,nmolty
     &           ,'uniform')
            call twh_cbmc_nb_one_gencode(GLB_INIT,numboxes,nmolty
     &           ,DIST_UNIFORM)
            call twh_cbmc_bond_generation(GLB_SET,'r^2 with bounds')
            call twh_cbmc_bondgen_code(GLB_SET,CBMC_R_SQ_WITH_BOUNDS)
            call twh_cbmc_bend_generation(GLB_SET,'ideal')
            call twh_cbmc_bendgen_code(GLB_SET,CBMC_IDEAL)
            call twh_cbmc_dihedral_generation(GLB_SET,'ideal')
            call twh_cbmc_dihedgen_code(GLB_SET,CBMC_IDEAL)
            call twh_two_bond_fe_bias_style(GLB_SET
     &           ,'analytic Boltzmann dihedral energy sum')
            call twh_three_bond_fe_bias_style(GLB_SET
     &           ,'analytic using max and min 2-4 distance')
            call twh_nch_nb_one(GLB_INIT,nmolty,10)
            call twh_nch_nb(GLB_INIT,nmolty,10)
            call twh_nch_tor(GLB_INIT,nmolty,360)
            call twh_nch_tor_connect(GLB_INIT,nmolty,360)
            call twh_nch_bend_a(GLB_INIT,nmolty,100)
            call twh_nch_bend_b(GLB_INIT,nmolty,100)
            call twh_nch_vib(GLB_INIT,nmolty,1000)
            call twh_nch_pre_nb(GLB_INIT,nmolty,1)
            call twh_vibrang(GLB_SET,1,0.85d0)
            call twh_vibrang(GLB_SET,2,1.15d0)
         endif
      elseif ( cbmc_setting_style .eq. 'default ideal'
     &        ) then
         if ( rwflag .eq. READ_FLAG ) then
c           --- use all of the default setting and styles that were
c           --- suggested by the manual at one point
            call twh_cbmc_fit_strategy(GLB_SET,'equilibrium')
            call twh_cbmc_nb_one_generation(GLB_INIT,numboxes,nmolty
     &           ,'uniform')
            call twh_cbmc_nb_one_gencode(GLB_INIT,numboxes,nmolty
     &           ,DIST_UNIFORM)
            call twh_cbmc_bond_generation(GLB_SET,'r^2 with bounds')
            call twh_cbmc_bondgen_code(GLB_SET,CBMC_R_SQ_WITH_BOUNDS)
            call twh_cbmc_bend_generation(GLB_SET,'ideal')
            call twh_cbmc_bendgen_code(GLB_SET,CBMC_IDEAL)
            call twh_cbmc_dihedral_generation(GLB_SET,'ideal')
            call twh_cbmc_dihedgen_code(GLB_SET,CBMC_IDEAL)
            call twh_two_bond_fe_bias_style(GLB_SET
     &           ,'analytic Boltzmann dihedral energy sum')
            call twh_three_bond_fe_bias_style(GLB_SET
     &           ,'analytic using max and min 2-4 distance')
            call twh_nch_nb_one(GLB_INIT,nmolty,10)
            call twh_nch_nb(GLB_INIT,nmolty,10)
            call twh_nch_tor(GLB_INIT,nmolty,360)
            call twh_nch_tor_connect(GLB_INIT,nmolty,360)
            call twh_nch_bend_a(GLB_INIT,nmolty,1000)
            call twh_nch_bend_b(GLB_INIT,nmolty,1000)
            call twh_nch_vib(GLB_INIT,nmolty,1000)
            call twh_nch_pre_nb(GLB_INIT,nmolty,1)
            call twh_vibrang(GLB_SET,1,0.85d0)
            call twh_vibrang(GLB_SET,2,1.15d0)
         endif
      elseif ( cbmc_setting_style .eq. 'widom ideal'
     &        ) then
         if ( rwflag .eq. READ_FLAG ) then
c           --- use the default ideal intramolecular generation, but
c           --- set all nb terms to 1 so that the widom insertion to
c           --- essentially remove the poor biasing seen in some 
c           ---- molecules related to the coulombic terms
            call twh_cbmc_fit_strategy(GLB_SET,'equilibrium')
            call twh_cbmc_nb_one_generation(GLB_INIT,numboxes,nmolty
     &           ,'uniform')
            call twh_cbmc_nb_one_gencode(GLB_INIT,numboxes,nmolty
     &           ,DIST_UNIFORM)
            call twh_cbmc_bond_generation(GLB_SET,'r^2 with bounds')
            call twh_cbmc_bondgen_code(GLB_SET,CBMC_R_SQ_WITH_BOUNDS)
            call twh_cbmc_bend_generation(GLB_SET,'ideal')
            call twh_cbmc_bendgen_code(GLB_SET,CBMC_IDEAL)
            call twh_cbmc_dihedral_generation(GLB_SET,'ideal')
            call twh_cbmc_dihedgen_code(GLB_SET,CBMC_IDEAL)
            call twh_two_bond_fe_bias_style(GLB_SET,'none')
            call twh_three_bond_fe_bias_style(GLB_SET,'none')
            call twh_nch_nb_one(GLB_INIT,nmolty,1)
            call twh_nch_nb(GLB_INIT,nmolty,1)
            call twh_nch_tor(GLB_INIT,nmolty,360)
            call twh_nch_tor_connect(GLB_INIT,nmolty,360)
            call twh_nch_bend_a(GLB_INIT,nmolty,1000)
            call twh_nch_bend_b(GLB_INIT,nmolty,1000)
            call twh_nch_vib(GLB_INIT,nmolty,1000)
            call twh_nch_pre_nb(GLB_INIT,nmolty,1)
            call twh_vibrang(GLB_SET,1,0.85d0)
            call twh_vibrang(GLB_SET,2,1.15d0)
         endif
      elseif ( cbmc_setting_style .eq. 'default autofit gaussian'
     &        ) then
         if ( rwflag .eq. READ_FLAG ) then
            call twh_cbmc_fit_strategy(GLB_SET,'equilibrium')
            call twh_cbmc_dihedral_generation(GLB_SET
     &           ,'autofit gaussian')
            call twh_cbmc_dihedgen_code(GLB_SET,CBMC_AUTOFIT_GAUSSIAN)
            call twh_dihedral_peak_weight_style(GLB_SET,'uniform')
            call twh_dihedral_sdev_multiplier(GLB_SET,1.0d0)
            call twh_cbmc_bend_generation(GLB_SET,'autofit gaussian')
            call twh_cbmc_bendgen_code(GLB_SET,CBMC_AUTOFIT_GAUSSIAN)
            call twh_bend_a_sdev_multiplier(GLB_SET,1.0d0)
            call twh_bend_b_sdev_multiplier(GLB_SET,1.0d0)
            call twh_cbmc_nb_one_generation(GLB_INIT,numboxes,nmolty
     &           ,'uniform')
            call twh_cbmc_nb_one_gencode(GLB_INIT,numboxes,nmolty
     &           ,DIST_UNIFORM)
            call twh_cbmc_bond_generation(GLB_SET,'autofit gaussian')
            call twh_cbmc_bondgen_code(GLB_SET,CBMC_AUTOFIT_GAUSSIAN)
            call twh_bond_sdev_multiplier(GLB_SET,1.0d0)
            call twh_two_bond_fe_bias_style(GLB_SET,'autofit gaussian')
            call twh_two_bond_bias_sdev_multiplier(GLB_SET,1.0d0)
            call twh_two_bond_bias_vibrange(GLB_SET,1,0.5d0)
            call twh_two_bond_bias_vibrange(GLB_SET,2,1.5d0)
            call twh_three_bond_fe_bias_style(GLB_SET 
     &           ,'autofit gaussian using max and min 2-4 distance')
            call twh_three_bond_bias_sdev_multiplier(GLB_SET,1.0d0)
            call twh_nch_nb_one(GLB_INIT,nmolty,10)
            call twh_nch_nb(GLB_INIT,nmolty,10)
            call twh_nch_bend_a(GLB_INIT,nmolty,1)
            call twh_nch_bend_b(GLB_INIT,nmolty,1)
            call twh_nch_vib(GLB_INIT,nmolty,1)
            if ( cbmc_form_code .eq. CBMC_FORM_MS1999_MT2004 ) then
               call twh_nch_tor(GLB_INIT,nmolty,10)
               call twh_nch_tor_connect(GLB_INIT,imolty,10)
               call twh_nch_pre_nb(GLB_INIT,nmolty,1)
            elseif ( cbmc_form_code .eq. CBMC_FORM_MF2006 ) then
               call twh_nch_tor(GLB_INIT,nmolty,1)
               call twh_nch_tor_connect(GLB_INIT,nmolty,1)
               call twh_nch_pre_nb(GLB_INIT,nmolty,10)
            else
               write(6,*) 'RWCBMC: unknown cbmc_form_code'
     &              ,' for this default setting'
               lfailure = .true.
               return
            endif
            call twh_vibrang(GLB_SET,1,0.85d0)
            call twh_vibrang(GLB_SET,2,1.15d0)
         endif
      elseif ( cbmc_setting_style .eq. 'Martin and Frischknecht'
     &       ) then
         if ( rwflag .eq. READ_FLAG ) then
            call twh_cbmc_fit_strategy(GLB_SET,'equilibrium')
            call twh_cbmc_dihedral_generation(GLB_SET
     &           ,'ideal + autofit gaussian')
            call twh_cbmc_dihedgen_code(GLB_SET
     &           ,CBMC_IDEAL_PLUS_AUTOFIT_GAUSSIAN)
            call twh_dihedral_peak_weight_style(GLB_SET,'uniform')
            call twh_dihedral_sdev_multiplier(GLB_SET,1.0d0)
            call twh_dihedral_ideal_fraction(GLB_SET,0.01d0)
            call twh_cbmc_bend_generation(GLB_SET
     &           ,'ideal + autofit gaussian')
            call twh_cbmc_bendgen_code(GLB_SET
     &           ,CBMC_IDEAL_PLUS_AUTOFIT_GAUSSIAN)
            call twh_bend_a_sdev_multiplier(GLB_SET,1.0d0)
            call twh_bend_b_sdev_multiplier(GLB_SET,1.0d0)
            call twh_bend_a_ideal_fraction(GLB_SET,0.01d0)
            call twh_bend_b_ideal_fraction(GLB_SET,0.01d0)
            call twh_cbmc_nb_one_generation(GLB_INIT,numboxes,nmolty
     &           ,'uniform')
            call twh_cbmc_nb_one_gencode(GLB_INIT,numboxes,nmolty
     &           ,DIST_UNIFORM)
            call twh_cbmc_bond_generation(GLB_SET,'autofit gaussian')
            call twh_cbmc_bondgen_code(GLB_SET,CBMC_AUTOFIT_GAUSSIAN)
            call twh_bond_sdev_multiplier(GLB_SET,1.0d0)
            call twh_two_bond_fe_bias_style(GLB_SET
     &           ,'analytic Boltzmann using angles')
            call twh_three_bond_fe_bias_style(GLB_SET
     &           ,'analytic using max and min 2-4 distance')
            call twh_nch_nb_one(GLB_INIT,nmolty,10)
            call twh_nch_nb(GLB_INIT,nmolty,10)
            call twh_nch_vib(GLB_INIT,nmolty,1)
            call twh_nch_bend_a(GLB_INIT,nmolty,1)
            if ( cbmc_form_code .eq. CBMC_FORM_MS1999_MT2004 ) then
               call twh_nch_pre_nb(GLB_INIT,nmolty,1)
               call twh_nch_tor(GLB_INIT,nmolty,100)
               call twh_nch_tor_connect(GLB_INIT,nmolty,100)
               call twh_nch_bend_b(GLB_INIT,nmolty,100)
            elseif ( cbmc_form_code .eq. CBMC_FORM_MF2006 ) then
               call twh_nch_pre_nb(GLB_INIT,nmolty,100)
               call twh_nch_tor(GLB_INIT,nmolty,1)
               call twh_nch_tor_connect(GLB_INIT,nmolty,1)
               call twh_nch_bend_b(GLB_INIT,nmolty,1)
            else
               write(6,*) 'RWCBMC: unknown cbmc_form_code'
     &              ,' for this default setting'
               lfailure = .true.
               return
            endif
         endif
      else
         write(6,*) 'RWCBMC: unknown cbmc_setting_style'
         write(6,*) cbmc_setting_style
         lfailure = .true.
         return
      endif

      if ( rwflag .eq. READ_FLAG ) then
         call twh_constant_pi(GLB_GET,onepi)
c        --- convert to radians
         call twh_sdevtor(GLB_GET,dvalue)
         dvalue = dvalue*onepi/180.0d0
         call twh_sdevtor(GLB_SET,dvalue)
c        --- convert to radians
         call twh_sdevbena(GLB_GET,dvalue)
         dvalue = dvalue*onepi/180.0d0
         call twh_sdevbena(GLB_SET,dvalue)
c        --- convert to radians 
         call twh_sdevbenb(GLB_GET,dvalue)
         dvalue = dvalue*onepi/180.0d0
         call twh_sdevbenb(GLB_SET,dvalue)
      endif

      return
      end
