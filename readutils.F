#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2001-2011 Marcus G. Martin                       *
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

c     ******************************************************************
c     * This file contains a collection of utilities to simplify and   *
c     * organize reading of parameters from files.                     *
c     * Generally, towhee input files consist of a label               *
c     * (eg. 'temperature')                                            *
c     * followed by one or more values.  Such values may be on one line*
c     * (for the case of a single number, or an array) or multiple     *
c     * lines.                                                         *
c     * While the input format for Towhee tends to be strict, with     *
c     * the proper labels required in the proper places, this set of   *
c     * functions allows for:                                          *
c     *   1) comments (defined as lines beginning with the symbol # )  *
c     *   2) optional labels, which may be checked for and read only if*
c     *      they are known to exist.                                  *
c     *                                                                *
c     * Function types:                                                *
c     * (TTT is one of integer, float, logical or string )             *
c     * twh_read_labeled_TTT -- top level convenience function, reads  *
c     *     label and following data.  Error if label incorrect.       *
c     * twh_read_label -- reads and confirms the label, error if label *
c     *     incorrect                                                  *
c     * twh_read_TTT -- reads given data type                          *
c     * twh_read_dir_string -- special-purpose version of              *
c     *    twh_read_string, for reading directory strings.             *
c     * twh_check_label -- reads and checks label, no error if label   *
c     *    incorrect                                                   *
c     * twh_peek_label -- reads and returns label name, then rewinds so*
c     *    it may be read again with twh_read_labeled_TTT functions    *
c     * twh_skip_comments -- reads lines repeatedly until non-comment  *
c     *    is found                                                    *
c     * twh_skip_lines -- skips data and/or label line(s).  comments   *
c     *    ignored.                                                    *
c     *                                                                *
c     * Arrays are read with the following subroutine:                 *
c     * twh_read_labeled_TTT_array -- top level convenience function,  *
c     *     reads label and array of specified size of following data  *
c     * twh_read_TTT_array - reads one line containing an array of data*
c     ******************************************************************
#define STDOUT 6
      function twh_read_labeled_integer(lfailure,io_unit,label,
     &    print_value,idname)
c     ******************************************************************
c     * twh_read_labeled_TTT(io_unit, label, print_value, idname)      *
c     *                                                                *
c     * Read label of given name and the associated value value.       *
c     * io_unit specifies the file number to read from                 *
c     * label is the label we expect                                   *
c     * Print_value specifies whether the label and value output to    *
c     * screen.                                                        *
c     * idname is the calling subroutine, printed only in case of error*
c     *                                                                *
c     * If label does not match the read label, an error is printed    *
c     * and execution halts.                                           *
c     * If the read value is not of the type expected (e.g. integer),  *
c     * an error is printed and execution halts.                       *
c     *                                                                *
c     * Lines beginning with the '#' symbol are ignored.               *
c     *                                                                *
c     * last modified 08-05-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_READ_INTEGER
#define FUNCTION_READ_LABELED_INTEGER
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure,print_value
      character*(*) label, idname
      integer io_unit
c     --- local vars
      integer value

c     --- Confirm correct label
      call twh_read_label(lfailure,io_unit,label,.false.,idname)
      if ( lfailure ) goto 66
      value = twh_read_integer(lfailure,io_unit,label,.false.)
      if ( lfailure ) goto 66
      if (print_value) then
         write(STDOUT,'(1x,a,a,1x,i12)') label,':',value
      endif

      twh_read_labeled_integer = value
      return
 66   write(STDOUT,*) 'failure in TWH_READ_LABELED_INTEGER'
      twh_read_labeled_integer = 0
      return
      end


      function twh_read_labeled_float(lfailure,io_unit, label
     &     ,print_value, idname)
c     ******************************************************************
c     * last modified 08-05-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_READ_FLOAT
#define FUNCTION_READ_LABELED_FLOAT
#include "functions.h"
c     --- variables passed to/from the function
      logical lfailure,print_value
      character*(*) label, idname
      integer io_unit
c     --- local vars
      double precision value

c     --- Confirm correct label
      call twh_read_label(lfailure,io_unit,label,.false.,idname)
      if ( lfailure ) goto 77
      value = twh_read_float(lfailure,io_unit,label,.false.)
      if ( lfailure ) goto 77

      if (print_value) then
         write(STDOUT,'(1x,a,a,1x,g12.6)') label,':',value
      endif

      twh_read_labeled_float = value
      return
 77   write(STDOUT,*) 'failure in TWH_READ_LABELED_FLOAT'
      twh_read_labeled_float = 0.0
      return
      end


      function twh_read_labeled_logical(lfailure,io_unit,label
     &     ,print_value,idname)
c     ******************************************************************
c     * last modified 08-05-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_READ_LABELED_LOGICAL
#define FUNCTION_READ_LABEL
#define FUNCTION_READ_LOGICAL
#include "functions.h"
c     --- variables passed to/from the function
      logical lfailure,print_value
      character*(*) label, idname
      integer io_unit
c     --- local vars
      logical value

c     --- Confirm correct label
      call twh_read_label(lfailure,io_unit,label,.false.,idname)
      if ( lfailure ) goto 78
      value = twh_read_logical(lfailure,io_unit,label,.false.)
      if ( lfailure ) goto 78

      if (print_value) then
         write(STDOUT,'(1x,a,a,1x,l1)') label,': ',value
      endif
      twh_read_labeled_logical = value
      return
 78   write(STDOUT,*) 'failure in TWH_READ_LABELED_LOGICAL'
      twh_read_labeled_logical = .false.
      return
      end

      subroutine twh_read_labeled_string(lfailure,io_unit,label
     &     ,print_value,idname,cvalue)
c     ******************************************************************
c     * Any leading or trailing quotes from the string are stripped.   *
c     * i.e.  'hello' becomes hello                                    *
c     *                                                                *
c     * last modified 08-05-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the function
      logical lfailure,print_value
      character*(*) label,idname,cvalue
      integer io_unit
c     --- local variables

c     --- Confirm correct label
      call twh_read_label(lfailure,io_unit,label,.false.,idname)
      if ( lfailure ) return
c     --- read in the string
      call twh_read_string(lfailure,io_unit,label,.false.,cvalue)
      if ( lfailure ) return

      if (print_value) then
         write(STDOUT,'(1x,a,a,1x,a)') label,':',cvalue
      endif
      return
      end

      function twh_read_integer(lfailure,io_unit,label,print_value)
c     ******************************************************************
c     * twh_read_TTT(io_unit, label, print_value)                      *
c     *                                                                *
c     * Read and return single value.                                  *
c     * io_unit specifies the file number to read from                 *
c     * label is the label we are reading, output only in case of error*
c     *    it may be blank ('')                                        *
c     * Print_value specifies whether the label and value output to    *
c     *   screen.                                                      *
c     *                                                                *
c     * If the read value is not of the type expected (e.g. integer),  *
c     *   an error is printed and execution halts.                     *
c     *                                                                *
c     * Lines beginning with the '#' symbol are ignored.               *
c     *                                                                *
c     * last modified 08-05-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_READ_INTEGER
#include "functions.h"
c     --- variables passed to/from the function
      integer io_unit
      logical lfailure
      logical print_value
      character*(*) label


      call twh_skip_comments(lfailure,io_unit)
      if ( lfailure ) goto 100
c     --- read the value
      read(io_unit,*,ERR=100) twh_read_integer

      if (print_value) then
         write(STDOUT,'(1x,i12)') twh_read_integer
      endif

      return
c     --- get here in case of read error
100   if ( lfailure ) then
         write(STDOUT,*) 'Error skipping comments'
      elseif (label .eq. '') then
         write(STDOUT, *) "Error reading integer value in unit ",io_unit
      else
         write(STDOUT, *) "Error reading integer value for label ",label
      endif
      lfailure = .true.
      return
      end

      function twh_read_float(lfailure,io_unit,label,print_value)
c     ******************************************************************
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_READ_FLOAT
#include "functions.h"
c     --- variables passed to/from the function
      integer io_unit
      logical lfailure,print_value
      character*(*) label

      call twh_skip_comments(lfailure,io_unit)
      if ( lfailure ) goto 101
c     --- read the value
      read(io_unit,*,ERR=101) twh_read_float

      if (print_value) then
         write(STDOUT,'(1x,g12.6)') twh_read_float
      endif

      return
c     --- get here in case of read error
101   if ( lfailure ) then
         write(6,*) 'Error skipping comments'
      elseif (label .eq. '') then
         write(STDOUT, *) "Error reading float value in unit ",io_unit
      else
         write(STDOUT, *) "Error reading float value for label ",label
      endif
      lfailure = .true.
      return
      end

      function twh_read_logical(lfailure,io_unit,label,print_value)
c     ******************************************************************
c     * last modified 08-05-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_READ_LOGICAL
#include "functions.h"
c     --- variables passed to/from the function
      integer io_unit
      logical lfailure,print_value
      character*(*) label

      call twh_skip_comments(lfailure,io_unit)
      if ( lfailure ) goto 102
c     --- read the value
      read(io_unit,*,ERR=102) twh_read_logical

      if (print_value) then
         write(STDOUT,'(1x,L1)') twh_read_logical
      endif

      return
c     --- get here in case of read error
102   if ( lfailure ) then
         write(STDOUT,*) 'Error skipping comments'
      elseif (label .eq. '') then
         write(STDOUT, *) "Error reading logical value in unit ",io_unit
      else
         write(STDOUT, *) "Error reading logical value for label ",label
      endif
      lfailure = .true.
      return
      end

      subroutine twh_read_string(lfailure,io_unit,label,print_value
     &     ,cvalue)
c     ******************************************************************
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the function
      integer io_unit
      logical lfailure,print_value
      character*(*) label,cvalue
c     --- local variables
      character*100 readvalue

      call twh_skip_comments(lfailure,io_unit)
      if ( lfailure ) return
c     --- read the value
      read(io_unit,"(a60)",ERR=103) readvalue
      call twh_strip_quotes(readvalue)
c     --- set cvalue to readvalue, required in the case of stripping
c     --- quotes from a maximum sized string
      cvalue = readvalue

      if (print_value) then
         write(STDOUT,'(1x,a)') cvalue
      endif

      return
c     --- get here in case of read error.
103   if (label .eq. '') then
         write(STDOUT, *) "Error reading string value in unit ",io_unit
      else
         write(STDOUT, *) "Error reading string value for label ",label
      endif
      lfailure = .true.
      return
      end

      function twh_read_dir_string(lfailure,io_unit,label,print_value)
c     ******************************************************************
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_READ_DIR_STRING
#include "functions.h"
c     --- variables passed to/from the function
      integer io_unit
      logical lfailure,print_value
      character*(*) label

c     --- local vars
      character*15 control_list

      call twh_skip_comments(lfailure,io_unit)
      if ( lfailure ) return
c     --- create control list
c     --- Note: this is sensitive to size of MAXDIRLENGTH -- be sure
c     --- control_list is what you expect if reads fail.
      write(control_list,"(a2,i3,a1)"),"(a",MAXDIRLENGTH,")"
c     --- read the value
      read(io_unit,control_list,ERR=104) twh_read_dir_string

      if (print_value) then
         write(STDOUT,'(a)') twh_read_dir_string
      endif

      return
c     --- get here in case of read error.
104   if (label .eq. '') then
         write(STDOUT, *) "Error reading string value in unit ",io_unit
      else
         write(STDOUT, *) "Error reading string value for label ",label
      endif
      lfailure = .true.
      return
      end

      subroutine twh_read_labeled_integer_array(lfailure,io_unit,label,
     &    array_length, print_values, idname, array_values)
c     ******************************************************************
c     * twh_read_labeled_TTT_array(io_unit, label, array_length,       *
c     *    print_values, idname, array_values)                         *
c     *                                                                *
c     * Read label of given name and the associated array value.       *
c     * io_unit specifies the file number to read from                 *
c     * label is the label we expect                                   *
c     * array_length is the length of the array we will be reading     *
c     * print_values specifies whether the label and values output to  *
c     *  screen.                                                       *
c     * idname is the calling subroutine, printed only in case of error*
c     * array_values is the actual array we will read into.  Its value *
c     *    is modified                                                 *
c     *                                                                *
c     * If label does not match the read label, an error is printed    *
c     * and execution halts.                                           *
c     * If the read value is not of the type expected (e.g. integer),  *
c     *  an error is printed and execution halts.                      *
c     * If the number of entries in the input line differs from        *
c     *  array_length, error is printed and execution halts.           *
c     *                                                                *
c     * Lines beginning with the '#' symbol are ignored.               *
c     *                                                                *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure,print_values
      character*(*) label, idname
      integer io_unit, array_length
      integer array_values(array_length) 

c     --- local vars
      integer i

c     --- Confirm correct label
      call twh_read_label(lfailure,io_unit,label,.false.,idname)
      if ( lfailure ) return
      call twh_read_integer_array(lfailure,io_unit,label,array_length
     &     ,.false.,array_values)
      if ( lfailure ) return

      if (print_values) then
         write(STDOUT,'(1x,a,a,10(1x,i12))') label,':',
     &      (array_values(i), i=1,array_length)
      endif

      return
      end

      subroutine twh_read_labeled_float_array(lfailure,io_unit,label,
     &    array_length, print_values, idname, array_values)
c     ******************************************************************
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure,print_values
      character*(*) label, idname
      integer io_unit, array_length
      double precision array_values(array_length) 

c     --- local vars
      integer i

c     --- Confirm correct label
      call twh_read_label(lfailure,io_unit,label,.false.,idname)
      if ( lfailure ) return
      call twh_read_float_array(lfailure,io_unit,label,array_length
     &     ,.false.,array_values)
      if ( lfailure ) return

      if (print_values) then
         write(STDOUT,'(1x,a,a,1x,10(g12.6))') label,':',
     &      (array_values(i), i=1,array_length)
      endif

      return
      end


      subroutine twh_read_labeled_logical_array(lfailure,io_unit,label
     &     ,array_length, print_values, idname, array_values)
c     ******************************************************************
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure,print_values
      character*(*) label, idname
      integer io_unit, array_length
      logical array_values(array_length) 

c     --- local vars
      integer i

c     --- Confirm correct label
      call twh_read_label(lfailure,io_unit,label,.false.,idname)
      if ( lfailure ) return
      call twh_read_logical_array(lfailure,io_unit,label,array_length
     &     ,.false.,array_values)
      if ( lfailure ) return

      if (print_values) then
         write(STDOUT,'(1x,a,a,10(1x,i12))') label,':',
     &      (array_values(i), i=1,array_length)
      endif

      return
      end



      subroutine twh_read_labeled_string_array(lfailure,io_unit,label,
     &    array_length, print_values, idname, array_values)
c     ******************************************************************
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure,print_values
      character*(*) label, idname
      integer io_unit, array_length
      character*(*) array_values(array_length) 

c     --- local vars
      integer i

c     --- Confirm correct label
      call twh_read_label(lfailure,io_unit,label,.false.,idname)
      if ( lfailure ) return
      call twh_read_string_array(lfailure,io_unit,label,array_length
     &     ,.false.,array_values)
      if ( lfailure ) return

      if (print_values) then
         write(STDOUT,'(1x,a,a,1x,10(a))') label,':',
     &      (array_values(i), i=1,array_length)
      endif

      return
      end


      subroutine twh_read_integer_array(lfailure,io_unit,label
     &     ,array_length,print_values, array_values)
c     ******************************************************************
c     * twh_read_TTT_array(io_unit, label, array_length, print_values, *
c     *   array_values)                                                *
c     *                                                                *
c     * Read and return an array of given type.                        *
c     * io_unit specifies the file number to read from                 *
c     * label is the label we are reading, output only in case of error*
c     *    it may be blank ('')                                        *
c     * array_length is the length of the array we will be reading     *
c     * print_values specifies whether the label and values output to  *
c     *    screen.                                                     *
c     * array_values is the actual array we will read into.  Its value *
c     *    is modified                                                 *
c     *                                                                *
c     * If the read value is not of the type expected (e.g. integer),  *
c     *   an error is printed and execution halts.                     *
c     * If the number of entries in the input line differs from        *
c     *   array_length, error is printed and execution halts.          *
c     *                                                                *
c     * Lines beginning with the '#' symbol are ignored.               *
c     *                                                                *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer io_unit, array_length
      character*(*) label
      logical print_values
      integer array_values(array_length) 
c     --- local variables
      integer i

      call twh_skip_comments(lfailure,io_unit)
      if ( lfailure ) return
c     --- read the value
      read(io_unit,*,ERR=100) (array_values(i), i=1,array_length)

      if (print_values) then
         write(STDOUT,'(10(1x,i12))')
     &        (array_values(i), i=1,array_length)
      endif

      return
c     --- get here in case of read error
100   if (label .eq. '') then
         write(STDOUT, *) "Error reading integer array in unit ",io_unit
      else
         write(STDOUT, *) "Error reading integer array for label ",label
      endif
      lfailure = .true.
      return
      end


      subroutine twh_read_float_array(lfailure,io_unit,label
     &     ,array_length,print_values, array_values)
c     ******************************************************************
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer io_unit, array_length
      character*(*) label
      logical lfailure,print_values
      double precision array_values(array_length) 

c     --- local variables
      integer i

      call twh_skip_comments(lfailure,io_unit)
      if ( lfailure ) return
c     --- read the value
      read(io_unit,*,ERR=100) (array_values(i), i=1,array_length)

      if (print_values) then
         write(STDOUT,'(10(g12.6))') (array_values(i), i=1,array_length)
      endif

      return
c     --- get here in case of read error
100   if (label .eq. '') then
         write(STDOUT, *) "Error reading float array in unit ",io_unit
      else
         write(STDOUT, *) "Error reading float array for label ",label
      endif
      lfailure = .true.
      return
      end


      subroutine twh_read_logical_array(lfailure,io_unit,label
     &     ,array_length,print_values, array_values)
c     ******************************************************************
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer io_unit, array_length
      character*(*) label
      logical lfailure,print_values
      logical array_values(array_length) 

c     --- local variables
      integer i

      call twh_skip_comments(lfailure,io_unit)
      if ( lfailure ) return
c     --- read the value
      read(io_unit,*,ERR=100) (array_values(i), i=1,array_length)

      if (print_values) then
         write(STDOUT,'(10(1x,L1))') (array_values(i), i=1,array_length)
      endif

      return
c     --- get here in case of read error
100   if (label .eq. '') then
         write(STDOUT, *) "Error reading logical array in unit ",io_unit
      else
         write(STDOUT, *) "Error reading logical array for label ",label
      endif
      lfailure = .true.
      return
      end

      subroutine twh_read_string_array(lfailure,io_unit,label
     &     ,array_length,print_values, array_values)
c     ******************************************************************
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer io_unit, array_length
      character*(*) label
      logical lfailure,print_values
      character*(*) array_values(array_length) 

c     --- local variables
      integer i

      call twh_skip_comments(lfailure,io_unit)
      if ( lfailure ) return
c     --- read the value
      read(io_unit,*,ERR=100) (array_values(i), i=1,array_length)

      if (print_values) then
         write(STDOUT,'(10(1x,a))') (array_values(i), i=1,array_length)
      endif

      return
c     --- get here in case of read error
100   if (label .eq. '') then
         write(STDOUT, *) "Error reading string array in unit ",io_unit
      else
         write(STDOUT, *) "Error reading string array for label ",label
      endif
      lfailure = .true.
      return
      end


      subroutine twh_read_label(lfailure,io_unit,label,print_value
     &     ,idname)
c     ******************************************************************
c     * Read label of given name.
c     * In case of a mismatched label this function halts execution.
c     * (Use twh_check_label to verify the expected label if necessary).
c     * Lines beginning with the '#' symbol are ignored.
c     * Leading blanks are ignored
c     *                                                                *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      character*(*) label, idname
      integer io_unit
      logical lfailure,print_value

      call twh_read_label_formatted(lfailure,io_unit,label,print_value,
     &      idname, '(a50)')
      if ( lfailure ) return

      return
      end


      subroutine twh_read_label_formatted(lfailure,io_unit,label
     &     ,print_value,idname, format_s)
c     ******************************************************************
c     * Read label of given name and specified label format (useful for*
c     * instances where part of label is variable).                    *
c     * In case of a mismatched label this function halts execution.   *
c     * (Use twh_check_label to verify the expected label if necessary)*
c     * Lines beginning with the '#' symbol are ignored.               *
c     * Leading blanks are ignored                                     *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none

c     --- variables passed to/from the subroutine
      character*(*) label, idname, format_s
      integer io_unit
      logical lfailure,print_value

c     -- local variables
      character*50 teststring

      call twh_skip_comments(lfailure,io_unit)
      if ( lfailure ) return
 
      read(io_unit,format_s,ERR=501), teststring

      call twh_strip_blank(teststring)

      if ( label .ne. teststring ) then
         write(STDOUT,*) 'Mismatch of strings in ',idname
         write(STDOUT,*) '  Expected String: ',label
         write(STDOUT,*) '  Found String: ', teststring
         write(STDOUT,*) 'Fatal error, quitting. '
         lfailure = .true.
         return
      endif
      if (print_value) then
         write(STDOUT,*) label,":"
      endif
      return

501   write(STDOUT,*) 'Error reading label in input file unit ',io_unit
      lfailure = .true.
      return
      end



      function twh_check_label(lfailure,io_unit,label)
c     ******************************************************************
c     * Check whether string in file matches label, and returns        *
c     * true or false accordingly.  When the funciton returns, the file*
c     * position remains unchanged -- meaning that this function can be*
c     * called repeatedly.                                             *
c     * Lines beginning with the '#' symbol are ignored.               *
c     *                                                                *
c     * last modified 08-05-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_CHECK_LABEL
#define FUNCTION_PEEK_LABEL
#include "functions.h"
c     --- variables passed to/from the subroutine
      character*(*) label
      logical lfailure
      integer io_unit
c     --- local variables
      character*128 s

      s = twh_peek_label(lfailure,io_unit)
      if ( lfailure ) goto 345
      if (label .eq. s) then
         twh_check_label = .true.
      else
         twh_check_label = .false.
      endif
      return
 345  write(6,*) 'failure in TWH_CHECK_LABEL'
      twh_check_label = .false.
      return
      end

      function twh_peek_label(lfailure,io_unit)
c     ******************************************************************
c     * Read and return label.  When the funciton returns, the file    *
c     * position remains unchanged -- meaning that normal label/value  *
c     * reading functions (e.g. twh_read_labeled_integer) can then be  *
c     * called.                                                        *
c     * Lines beginning with the '#' symbol are ignored.               *
c     *                                                                *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_PEEK_LABEL
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer io_unit

      call twh_skip_comments(lfailure,io_unit)
      if ( lfailure ) return
      read(io_unit,*,ERR=505) twh_peek_label
      backspace(io_unit)
      return

505   write(STDOUT,*) 'Error reading input file unit ',io_unit
      lfailure = .true.
      return
      end


      subroutine twh_skip_comments(lfailure,io_unit)
c     ******************************************************************
c     * Advances to the first line which is not a comment.             *
c     * last modified 08-04-2001 by M.G. Martin                        *
c     ******************************************************************
      implicit none
      logical lfailure,twh_check_comment

c     --- variables passed to/from the subroutine
      integer io_unit

      do while ( twh_check_comment(lfailure,io_unit) )
         if ( lfailure ) return
      enddo

c900   if (twh_check_comment(io_unit)) then
c         goto 900
c      endif

      return
      end

      subroutine twh_skip_lines(lfailure,io_unit,num_lines)
c     ******************************************************************
c     * Reads a number of lines (given by num_lines) which are label   *
c     * and/or data lines.  Comment lines are ignored and do not count *
c     * in the num_lines count.                                        *
c     *                                                                *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer io_unit, num_lines
c     -- local variables
      integer l
      character*50 stringread

      do l=1,num_lines
         call twh_skip_comments(lfailure,io_unit)
         call twh_read_string(lfailure,io_unit, '', .false.,stringread)
         if ( lfailure ) return
      enddo

      return
      end

      subroutine twh_strip_quotes(s)
c     ******************************************************************
c     * Strips leading and trailing quotes (') from the string s if    *
c     * they exist.  If quotes do not exist, or if they're internal    *
c     * (neither first nor last character) then they are left as is.   *
c     *                                                                *
c     * last modified 03-08-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_STRING_LENGTH
#include "functions.h"
c     --- variables passed to/from the subroutine
      character*(*) s
c     --- local variables
      integer i

c     --- find and strip leading quote
      if (s(1:1) .eq. '''') then
         s(1:len(s)-1) = s(2:len(s))
         s(len(s):len(s)) = ' '
      endif

c     --- find length of remaining string, ignoring padding
      i = twh_get_string_length(s)

c     --- i now contains the length of non-padded string
      if (i .gt. 0) then 
         if (s(i:i) .eq. '''') then
            s(i:i) = ' '
         endif
      endif

      return
      end

      function twh_get_string_length(s)
c     ******************************************************************
c     * Determines string length, ignoring blank padding               *
c     * integer returned is index of last character.                   *
c     *                                                                *
c     * last modified 03-08-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_STRING_LENGTH
#include "functions.h"
c     --- variable passed to/from this function
      character*(*) s
c     --- local variable
      integer i

c     --- Code courtesy of
c     --- http://www.star.le.ac.uk/~cgp/prof77.html#tth_sEc7.1i
      do 150, i = len(s), 1, -1 
         if(s(i:i) .ne. ' ' .and. s(i:i) .ne. '\0') go to 151 
150   continue 
151   continue
      twh_get_string_length = i
      return
      end


      subroutine twh_strip_blank(s)
c     ******************************************************************
c     * this routine first converts Tabs (\t) into a single space.     *
c     * It then strips all leading space from a string (including those*
c     * that were Tabs.                                                *
c     * Internal spaces are left intact                                *
c     * Data returned in variable s                                    *
c     * last modified 01-24-2009 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      character*(*) s
c     --- local variables
      integer i
      character*100 stmp

c     --- first, convert all tabs to spaces (these lead to weirdness
c     --- in string comparison). 
      do 155, i = 1, len(s)
         if(s(i:i) .eq. '\t' .or. s(i:i) .eq. char(9) ) then
            s(i:i) = ' '
         endif
155   continue 

c     --- count the number of leading spaces
      do 160, i = 1, len(s)
         if(s(i:i) .ne. ' ') go to 161 
160   continue 

c     --- i now contains the length of initial padding
161   if (i .gt. 1) then 
         stmp = s(i:len(s))
         s = stmp
      endif

      return
      end


      function twh_check_comment(lfailure,io_unit)
c     ******************************************************************
c     * Checks whether the current line begins with the '#' symbol, and 
c     * returns true if it does.
c     * The file position is rewound if a comment is not found; hence,
c     * repeated calls will skip through comments, but will not move
c     * past data.
c     *                                                                *
c     * last modified 08-05-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_CHECK_COMMENT
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer io_unit
c     --- local variables
      character*1 s

      read(io_unit,'(a1)',ERR=510) s

      if (s .eq. '#') then
         twh_check_comment = .true.
         return
      else
         backspace(io_unit)
         twh_check_comment = .false.
         return
      endif

510   write(STDOUT,*) 'Error reading input file unit ',io_unit
      twh_check_comment = .false.
      lfailure = .true.
      return
      end



