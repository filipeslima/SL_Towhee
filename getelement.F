#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2002-2008 Marcus G. Martin                       *
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
      subroutine twh_getelement(lfailure,flag,testmass,testelement
     &,ielement )
c     ******************************************************************
c     * checks our list of elemental masses and returns that element   *
c     * if its mass is within the tolerance of testmass                *
c     * alternatively returns a mass given an element name             *
c     * or, returns a mass given an atom number                        *
c     *                                                                *
c     * originally written 03-26-2002 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- debugging logical
#define DEBUG_GETELEMENT 0     
c     --- variables passed to the function
      logical lfailure
      character*2 testelement
      character*(*) flag
      integer ielement
      double precision testmass
c     --- local variables
      integer totalmasses,imass
      parameter ( totalmasses = 110)
      double precision toler
      parameter (toler = 0.01d0)
      double precision masslist
      dimension masslist(totalmasses)
      character*2 elelist
      dimension elelist(totalmasses)

#if DEBUG_GETELEMENT
      write(6,*) 'GETELEMENT: flag,testmass,testelement,ielement'
     &,flag,testmass,testelement,ielement
#endif

      masslist(1) = 1.00794d0
      masslist(2) = 4.002602d0
      masslist(3) = 6.941d0
      masslist(4) = 9.012182d0
      masslist(5) = 10.811d0
      masslist(6) = 12.0107d0
      masslist(7) = 14.0067d0
      masslist(8) = 15.9994d0
      masslist(9) = 18.9984032d0
      masslist(10) = 20.1797d0
      masslist(11) = 22.989770d0
      masslist(12) = 24.3050d0
      masslist(13) = 26.981538d0
      masslist(14) = 28.0855d0
      masslist(15) = 30.973761d0
      masslist(16) = 32.065d0
      masslist(17) = 35.453d0
      masslist(18) = 39.948d0
      masslist(19) = 39.0983d0
      masslist(20) = 40.078d0
      masslist(21) = 44.955910d0
      masslist(22) = 47.867d0
      masslist(23) = 50.9415d0
      masslist(24) = 51.9961d0
      masslist(25) = 54.938049d0
      masslist(26) = 55.845d0
      masslist(27) = 58.933200d0
      masslist(28) = 58.6934d0
      masslist(29) = 63.546d0
      masslist(30) = 65.39d0
      masslist(31) = 69.723d0
      masslist(32) = 72.64d0
      masslist(33) = 74.92160d0
      masslist(34) = 78.96d0
      masslist(35) = 79.904d0
      masslist(36) = 83.80d0
      masslist(37) = 85.4678d0
      masslist(38) = 87.62d0
      masslist(39) = 88.90585d0
      masslist(40) = 91.224d0
      masslist(41) = 92.90638d0
      masslist(42) = 95.94d0
      masslist(43) = 98.0d0
      masslist(44) = 101.07d0
      masslist(45) = 102.90550d0
      masslist(46) = 106.42d0
      masslist(47) = 107.8682d0
      masslist(48) = 112.411d0
      masslist(49) = 114.818d0
      masslist(50) = 118.710d0
      masslist(51) = 121.760d0
      masslist(52) = 127.60d0
      masslist(53) = 126.90447d0
      masslist(54) = 131.293d0
      masslist(55) = 132.90545d0
      masslist(56) = 137.327d0
      masslist(57) = 138.9055d0
      masslist(58) = 140.116d0
      masslist(59) = 140.90765d0
      masslist(60) = 144.24d0
      masslist(61) = 145.0d0
      masslist(62) = 150.36d0
      masslist(63) = 151.964d0
      masslist(64) = 157.25d0
      masslist(65) = 158.92534d0
      masslist(66) = 162.50d0
      masslist(67) = 164.93032d0
      masslist(68) = 167.259d0
      masslist(69) = 168.93421d0
      masslist(70) = 173.04d0
      masslist(71) = 174.967d0
      masslist(72) = 178.49d0
      masslist(73) = 180.9479d0
      masslist(74) = 183.84d0
      masslist(75) = 186.207d0
      masslist(76) = 190.23d0
      masslist(77) = 192.217d0
      masslist(78) = 195.078d0
      masslist(79) = 196.96655d0
      masslist(80) = 200.59d0
      masslist(81) = 204.3833d0
      masslist(82) = 207.2d0
      masslist(83) = 208.98038d0
      masslist(84) = 209.d0
      masslist(85) = 210.d0
      masslist(86) = 222.d0
      masslist(87) = 223.d0
      masslist(88) = 226.d0
      masslist(89) = 227.d0
      masslist(90) = 232.0381d0
      masslist(91) = 231.03588d0
      masslist(92) = 238.02891d0
      masslist(93) = 237.d0
      masslist(94) = 244.d0
      masslist(95) = 243.d0
      masslist(96) = 247.d0
      masslist(97) = 247.d0
      masslist(98) = 251.d0
      masslist(99) = 252.d0
      masslist(100) = 257.d0
      masslist(101) = 258.d0
      masslist(102) = 259.d0
      masslist(103) = 262.d0
      masslist(104) = 261.d0
      masslist(105) = 262.d0
      masslist(106) = 266.d0
      masslist(107) = 264.d0
      masslist(108) = 277.d0
      masslist(109) = 268.d0
      masslist(110) = 1.00794d0


      elelist(1) = ' H'
      elelist(2) = 'He'
      elelist(3) = 'Li'
      elelist(4) = 'Be'
      elelist(5) = ' B'
      elelist(6) = ' C'
      elelist(7) = ' N'
      elelist(8) = ' O'
      elelist(9) = ' F'
      elelist(10) = 'Ne'
      elelist(11) = 'Na'
      elelist(12) = 'Mg'
      elelist(13) = 'Al'
      elelist(14) = 'Si'
      elelist(15) = ' P'
      elelist(16) = ' S'
      elelist(17) = 'Cl'
      elelist(18) = 'Ar'
      elelist(19) = ' K'
      elelist(20) = 'Ca'
      elelist(21) = 'Sc'
      elelist(22) = 'Ti'
      elelist(23) = ' V'
      elelist(24) = 'Cr'
      elelist(25) = 'Mn'
      elelist(26) = 'Fe'
      elelist(27) = 'Co'
      elelist(28) = 'Ni'
      elelist(29) = 'Cu'
      elelist(30) = 'Zn'
      elelist(31) = 'Ga'
      elelist(32) = 'Ge'
      elelist(33) = 'As'
      elelist(34) = 'Se'
      elelist(35) = 'Br'
      elelist(36) = 'Kr'
      elelist(37) = 'Rb'
      elelist(38) = 'Sr'
      elelist(39) = ' Y'
      elelist(40) = 'Zr'
      elelist(41) = 'Nb'
      elelist(42) = 'Mo'
      elelist(43) = 'Tc'
      elelist(44) = 'Ru'
      elelist(45) = 'Rh'
      elelist(46) = 'Pd'
      elelist(47) = 'Ag'
      elelist(48) = 'Cd'
      elelist(49) = 'In'
      elelist(50) = 'Sn'
      elelist(51) = 'Sb'
      elelist(52) = 'Te'
      elelist(53) = ' I'
      elelist(54) = 'Xe'
      elelist(55) = 'Cs'
      elelist(56) = 'Ba'
      elelist(57) = 'La'
      elelist(58) = 'Ce'
      elelist(59) = 'Pr'
      elelist(60) = 'Nd'
      elelist(61) = 'Pm'
      elelist(62) = 'Sm'
      elelist(63) = 'Eu'
      elelist(64) = 'Gd'
      elelist(65) = 'Tb'
      elelist(66) = 'Dy'
      elelist(67) = 'Ho'
      elelist(68) = 'Er'
      elelist(69) = 'Tm'
      elelist(70) = 'Yb'
      elelist(71) = 'Lu'
      elelist(72) = 'Hf'
      elelist(73) = 'Ta'
      elelist(74) = ' W'
      elelist(75) = 'Re'
      elelist(76) = 'Os'
      elelist(77) = 'Ir'
      elelist(78) = 'Pt'
      elelist(79) = 'Au'
      elelist(80) = 'Hg'
      elelist(81) = 'Tl'
      elelist(82) = 'Pb'
      elelist(83) = 'Bi'
      elelist(84) = 'Po'
      elelist(85) = 'At'
      elelist(86) = 'Rn'
      elelist(87) = 'Fr'
      elelist(88) = 'Ra'
      elelist(89) = 'Ac'
      elelist(90) = 'Th'
      elelist(91) = 'Pa'
      elelist(92) = ' U'
      elelist(93) = 'Np'
      elelist(94) = 'Pu'
      elelist(95) = 'Am'
      elelist(96) = 'Cm'
      elelist(97) = 'Bk'
      elelist(98) = 'Cf'
      elelist(99) = 'Es'
      elelist(100) = 'Fm'
      elelist(101) = 'Md'
      elelist(102) = 'No'
      elelist(103) = 'Lr'
      elelist(104) = 'Rf'
      elelist(105) = 'Db'
      elelist(106) = 'Sg'
      elelist(107) = 'Bh'
      elelist(108) = 'Hs'
      elelist(109) = 'Mt'
      elelist(110) = 'H '

      if ( flag .eq. 'getelement' ) then
c        --- see if we match any of these masses
         do imass = 1,totalmasses
            if ( (masslist(imass) - testmass)**2 .lt. toler**2 ) then
c              --- this is a match
               testelement = elelist(imass)
               return
            endif
         enddo
c        --- if we got here there was no match
         write(6,*) 'GETELELMENT: no match found for mass',testmass
         testelement = '--'
      elseif ( flag .eq. 'getmass' ) then
c        --- see if we match any of these element names
         do imass = 1,totalmasses
            if ( testelement .eq. elelist(imass) ) then
c              --- this is a match
               testmass = masslist(imass)
               return
            endif
         enddo
c        --- if we got here then there was no match
         write(6,*) 'GETELEMENT: no match found for element',testelement
         testmass = 1.0d0
      elseif ( flag .eq. 'get int' ) then
c        --- match with element names and return the element integer
         do imass = 1,totalmasses
            if ( testelement .eq. elelist(imass) ) then
c              --- this is a match
               ielement = imass
               return
            endif
         enddo
c        --- if we got here then there was no match
         write(6,*) 'GETELEMENT: no match found for element',testelement
         ielement = 0
      elseif ( flag .eq. 'get e and m' ) then
         testmass = masslist(ielement)
         testelement = elelist(ielement)
      else
         write(6,*) 'GETELEMENT: invalid flag'
         write(6,*) flag
         lfailure = .true.
         return
      endif

      return
      end
