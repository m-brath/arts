/* Copyright (C) 2003 Oliver Lemke <olemke@uni-bremen.de>

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */


////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////
/*!
  \file   bifstream.cc
  \author Oliver Lemke <olemke@uni-bremen.de>
  \date   2003-01-23

  \brief This file contains the class implementation of bifstream.

*/

#include <fstream>
#include "bifstream.h"
using namespace std;

void bifstream::seek(long pos, Offset offs)
{
  if(!in) { err = NotOpen; return; }

  switch(offs) {
  case Set: this->seekg(pos, ios::beg); break;
  case Add: this->seekg(pos, ios::cur); break;
  case End: this->seekg(pos, ios::end); break;
  }
}

long bifstream::pos()
{
  if(!in) { err = NotOpen; return 0; }
  return (long)this->tellg();
}

bifstream::Byte bifstream::getByte()
{
  int read;

  if(this->good ()) {
        read = this->get ();
        if(read == EOF) err |= Eof;
        return (Byte)read;
  } else {
        err |= NotOpen;
        return 0;
  }
}


/* Overloaded input operators */
bifstream &operator>> (bifstream &bif, double &n)
{ n = (double)bif.readFloat (binio::Double); return (bif); }

bifstream &operator>> (bifstream &bif, float &n)
{ n = (float)bif.readFloat (binio::Double); return (bif); }

bifstream &operator>> (bifstream &bif, long &n)
{ n = (long)bif.readInt (8); return (bif); }

bifstream &operator>> (bifstream &bif, int &n)
{ n = (int)bif.readInt (8); return (bif); }

