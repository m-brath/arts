/* Copyright (C) 2000 Stefan Buehler <sbuehler@uni-bremen.de>

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

// Stuff related to the calculation of absorption coefficients.

#include "arts.h"
#include "vecmat.h"
#include "messages.h"
#include "file.h"
#include "absorption.h"
#include "wsv.h"
#include "md.h"
#include "math_funcs.h"
#include "atm_funcs.h"

void linesReadFromHitran(// WS Output:
                         ARRAYofLineRecord& lines,
                          // Control Parameters:
                         const string& filename,
                         const Numeric& fmin,
                         const Numeric& fmax)
{
  ifstream is;

  out2 << "  Reading file: " << filename << '\n';
  open_input_file(is, filename);

  bool go_on = true;
  while ( go_on )
    {
      LineRecord lr;
      if ( lr.ReadFromHitranStream(is) )
	{
	  // If we are here the read function has reached eof and has
	  // returned no data.
	  go_on = false;
	}
      else
	{
	  if ( fmin <= lr.F() )
	    {
	      if ( lr.F() <= fmax )
		lines.push_back(lr);
	      else
		go_on = false;
	    }
	}
    }
  out2 << "  Read " << lines.size() << " lines.\n";
}


void linesReadFromMytran2(// WS Output:
			  ARRAYofLineRecord& lines,
                          // Control Parameters:
			  const string& filename,
			  const Numeric& fmin,
			  const Numeric& fmax)
{
  ifstream is;

  out2 << "  Reading file: " << filename << '\n';
  open_input_file(is, filename);

  bool go_on = true;
  while ( go_on )
    {
      LineRecord lr;
      if ( lr.ReadFromMytran2Stream(is) )
	{
	  // If we are here the read function has reached eof and has
	  // returned no data.
	  go_on = false;
	}
      else
	{
	  // lines are not necessarily frequency sorted 
	  if ( fmin <= lr.F() )
	    if ( lr.F() <= fmax )
	      lines.push_back(lr);
	}
    }
  out2 << "  Read " << lines.size() << " lines.\n";
}

void linesReadFromJpl(// WS Output:
		      ARRAYofLineRecord& lines,
		      // Control Parameters:
		      const string& filename,
		      const Numeric& fmin,
		      const Numeric& fmax)
{
  ifstream is;

  out2 << "  Reading file: " << filename << '\n';
  open_input_file(is, filename);

  bool go_on = true;
  while ( go_on )
    {
      LineRecord lr;
      if ( lr.ReadFromJplStream(is) )
	{
	  // If we are here the read function has reached eof and has
	  // returned no data.
	  go_on = false;
	}
      else
	{
	  // we expect lines to be sorted
	  if ( fmin <= lr.F() )
	    {
	      if ( lr.F() <= fmax )
		lines.push_back(lr);
	      else
		go_on = false;
	    }
	}
    }
  out2 << "  Read " << lines.size() << " lines.\n";
}


void lines_per_tgCreateFromLines(// WS Output:
                                  ARRAYofARRAYofLineRecord& lines_per_tg,
                                  // WS Input:
                                  const ARRAYofLineRecord&   lines,
                                  const TagGroups&           tag_groups)
{
  // The species lookup data:
  extern const ARRAY<SpeciesRecord> species_data;

  // As a safety feature, we will watch out for the case that a
  // species is included in the calculation, but not all lines are
  // used. For this we need an array to flag the used species:
  
  // For some weird reason, ARRAYs of bool do not work, although all
  // other types seem to work fine. So in this single case, I'll use
  // the stl vector directly. The other place where this is done is in
  // the function executor in main.cc.
  // FIXME: Fix this when ARRAY<bool> works.
  std::vector<bool> species_used (species_data.size(),false);
      
  // Make lines_per_tg the right size:
  lines_per_tg = ARRAYofARRAYofLineRecord(tag_groups.size());

  // Unfortunately, MTL conatains a bug that leads to all elements of
  // the outer ARRAY of an ARRAY<ARRAY>> pointing to the same data
  // after creation. So we need to fix this explicitly:
  for ( size_t i=0; i<lines_per_tg.size(); ++i )
    lines_per_tg[i] = ARRAYofLineRecord();

  // Loop all lines in the input line list:
  for ( size_t i=0; i<lines.size(); ++i )
    {
      // Get a convenient reference to the current line:
      const LineRecord& this_line = lines[i];

      // We have to test each tag group in turn (in the order in which
      // they appear in the controlfile). The line is assigned to the
      // first tag group that fits.

      // The flag found is used to break the for loops when the right
      bool found = false;

      // We need to define j here, since we need the value outside the
      // for loop:
      size_t j;

      // Loop the tag groups:
      for ( j=0; j<tag_groups.size() && !found ; ++j ) 
	{
	  // A tag group can contain several tags:
	  for ( size_t k=0; k<tag_groups[j].size() && !found; ++k )
	    {
	      // Get a reference to the current tag (not really
	      // necessary, but makes for nicer notation below):
	      const OneTag& this_tag = tag_groups[j][k];

	      // Now we will test different attributes of the line
	      // against matching attributes of the tag. If any
	      // attribute does not match, we continue with the next tag
	      // in the tag group. (Exception: Species, see just below.)

	      // Test species. If this attribute does not match we don`t
	      // have to test the other tags in this group, since all
	      // tags must belong to the same species.
	      if ( this_tag.Species() != this_line.Species() ) break;

	      // Test isotope. The isotope can either match directly, or
	      // the Isotope of the tag can be one larger than the
	      // number of isotopes, which means `all'. Test the second
	      // condition first, since this will probably be more often
	      // used.
	      if ( this_tag.Isotope() != this_line.SpeciesData().Isotope().size() )
		if ( this_tag.Isotope() != this_line.Isotope() )
		  continue;

	    // Test frequncy range we take both the lower (Lf) and the
	    // upper (Uf) border frequency to include the `equal' case.
	    // Both Lf and Uf can also be negative, which means `no limit'

	    // Take the lower limit first:
	      if ( this_tag.Lf() >= 0 )
		if ( this_tag.Lf() > this_line.F() )
		  continue;

	    // Then the upper limit:
	      if ( this_tag.Uf() >= 0 )
		if ( this_tag.Uf() < this_line.F() )
		  continue;

	    // When we get here, this_tag has survived all tests. That
	    // means it matches the line perfectly!
	      found = true;
	    }
	}

      // If a matching tag was found, this line can be used in the
      // calculation. Add it to the line list for this tag group.
      if (found)
	{
	  // We have to use j-1 here, since j was still increased by
	  // one after the matching tag has been found.
	  lines_per_tg[j-1].push_back(this_line);

	  // Flag this species as used, if not already done:
	  if ( !species_used[this_line.Species()] )
	    species_used[this_line.Species()] = true;
	}
      else
	{
	  // Safety feature: Issue a warning messages if the lines for a
	  // species are only partly covered by tags.
	  if ( species_used[this_line.Species()] )
	    {
	      out0 << "Your tags include other lines of species "
		   << this_line.SpeciesData().Name()
		   << ",\n"
		   << "why do you not include line "
		   << i
		   << " (at "
		   << this_line.F()
		   << " Hz)?\n";
	    }
	}
   }

  // Write some information to the lowest priority output stream.
  for (size_t i=0; i<tag_groups.size(); ++i)
    {
	out3 << "  " << i << ":";

	for (size_t s=0; s<tag_groups[i].size(); ++s)
	  out3 << " " << tag_groups[i][s].Name();

	out3 << ": " << lines_per_tg[i].size() << " lines\n";
    }

}


void linesWriteToFile(// WS Input:
		      const ARRAYofLineRecord& lines,
		      // Control Parameters:
		      const string& f)
{
  string filename = f;
  
  // If filename is empty then set the default filename:
  if ( "" == filename )
    {
      extern const string basename;                       
      filename = basename+".lines.al";
    }

  ofstream os;

  out2 << "  Writing file: " << filename << '\n';
  open_output_file(os, filename);

  write_lines_to_stream(os,lines);
}


void lines_per_tgWriteToFile(// WS Input:
			      const ARRAYofARRAYofLineRecord& lines_per_tg,
			      // Control Parameters:
			      const string& f)
{
  string filename = f;
  
  // If filename is empty then set the default filename:
  if ( "" == filename )
    {
      extern const string basename;                       
      filename = basename+".lines_per_tg.al";
    }

  ofstream os;

  out2 << "  Writing file: " << filename << '\n';
  open_output_file(os, filename);

  os << lines_per_tg.size() << '\n';

  for ( size_t i=0; i<lines_per_tg.size(); ++i )
    {
      const ARRAYofLineRecord& lines = lines_per_tg[i];
      os << lines.size() << '\n';
      write_lines_to_stream( os, lines );
    }
}


void tag_groupsDefine(// WS Output:
                      TagGroups& tag_groups,
                      // Control Parameters:
                      const ARRAYofstring& tags)
{
  tag_groups = TagGroups(tags.size());

  // Each element of the array of strings tags defines one tag
  // group. Let's work through them one by one.
  for ( size_t i=0; i<tags.size(); ++i )
    {
      // There can be a comma separated list of tag definitions, so we
      // need to break the string apart at the commas.
      ARRAY<string> tag_def;
      bool go_on = true;
      string these_tags = tags[i];
      while (go_on)
	{
	  size_t n = these_tags.find(',');
	  if ( n >= these_tags.size() )
	    {
	      // There are no more commas.
	      tag_def.push_back(these_tags);
	      go_on = false;
	    }
	  else
	    {
	      tag_def.push_back(these_tags.substr(0,n));
	      these_tags.erase(0,n+1);
	    }
	}

      // tag_def now holds the different tag strings for this group.
      //      cout << "tag_def =\n" << tag_def << endl;


      // Unfortunately, MTL conatains a bug that leads to all elements of
      // the outer ARRAY of an ARRAY<ARRAY>> pointing to the same data
      // after creation. So we need to fix this explicitly:
      tag_groups[i] = ARRAY<OneTag>();

      for ( size_t s=0; s<tag_def.size(); ++s )
	{
	  // Remove leading whitespace, if there is any:
	  while ( ' '  == tag_def[s][0] ||
		  '\t' == tag_def[s][0]    )	tag_def[s].erase(0,1);

	  OneTag this_tag(tag_def[s]);

	  // Safety check: For s>0 check that the tags belong to the same species.
	  if (s>0)
	    if ( tag_groups[i][0].Species() != this_tag.Species() )
	      throw runtime_error("Tags in a tag group must belong to the same species.");

	  tag_groups[i].push_back(this_tag);
	}
    }

  // Print list of tag groups to the most verbose output stream:
  out3 << "  Defined tag groups:";
  for ( size_t i=0; i<tag_groups.size(); ++i )
    {
      out3 << "\n  " << i << ":";
      for ( size_t s=0; s<tag_groups[i].size(); ++s )
	{
	  out3 << " " << tag_groups[i][s].Name();
	}
    }
  out3 << '\n';

//  cout << endl << endl << tag_groups << endl;
}



void lineshapeDefine(// WS Output:
		     ARRAYofsizet&    lineshape,
		     ARRAYofsizet&    lineshape_norm,
		     // WS Input:
		     const TagGroups& tag_groups,
		     const string&    shape,
		     const string&    normalizationfactor)
{
  // Make lineshape and normalization factor data visible:
  extern const ARRAY<LineshapeRecord> lineshape_data;
  extern const ARRAY<LineshapeNormRecord> lineshape_norm_data;


  // generate the right number of elements
  size_t tag_sz = tag_groups.size();
  lineshape.resize(tag_sz);
  lineshape_norm.resize(tag_sz);

  // Is this lineshape available?
  int found0=-1;
  for ( size_t i=0; i<lineshape_data.size() && (found0 == -1) ; ++i )
    {
      const string& str = lineshape_data[i].Name();
      if (str == shape) 
	{
	  out2 << "  Selected lineshape: " << str << "\n";
	  found0=i;
	}
    }

  // Is this normalization to the lineshape available?
  int found1=-1;
  for ( size_t i=0; i<lineshape_norm_data.size() && (found1 == -1); ++i )
    {
      const string& str = lineshape_norm_data[i].Name();
      if (str == normalizationfactor) 
	{
	  out2 << "  Selected normalization factor: " << normalizationfactor << "\n";
	  found1=i;
	}
    }


  // did we find the lineshape and normalization factor?
  if (found0 == -1)
    throw runtime_error("Selected lineshape not available.");
  if (found1 == -1)
    throw runtime_error("Selected normalization to lineshape not available.");


  // now set the lineshape and lineshape_norm workspace variables 
  for (size_t i=0; i<tag_sz; i++)
    {
      lineshape[i]=(size_t) found0;
      lineshape_norm[i]=(size_t) found1;
    }	  
}

void lineshape_per_abs_tagDefine(// WS Output:
				 ARRAYofsizet&         lineshape,
				 ARRAYofsizet&         lineshape_norm,
				 // WS Input:
				 const TagGroups&      tag_groups,
				 const ARRAYofstring&  shape,
				 const ARRAYofstring&  normalizationfactor)
{
  // Make lineshape and normalization factor data visible:
  extern const ARRAY<LineshapeRecord> lineshape_data;
  extern const ARRAY<LineshapeNormRecord> lineshape_norm_data;

  // check that the number of elements are equal
  size_t tag_sz = tag_groups.size();
  if ( (tag_sz != shape.size()) ||
       (tag_sz != normalizationfactor.size()) )
    {
      ostringstream os;
      os << "lineshape_per_abs_tagDefine: number of elements does\n"
	 << "not match the number of tags defined.";
      throw runtime_error(os.str());
    }
      

  // generate the right number of elements
  lineshape.resize(tag_sz);
  lineshape_norm.resize(tag_sz);

  // Is this lineshape available?
  for (size_t k=0; k<tag_sz; ++k)
    {
      int found0=-1;
      for ( size_t i=0; i<lineshape_data.size() && (found0 == -1); ++i )
	{
	  const string& str = lineshape_data[i].Name();
	  if (str == shape[k]) 
	    {
	      out2 << "  Tags: [";
	      for (size_t s=0; s<tag_groups[k].size()-1; ++s)
		out2 << tag_groups[k][s].Name() << ", "; 
	      out2 << tag_groups[k][tag_groups[k].size()-1].Name() << "]\n";
	      out2 << "  Selected lineshape: " << str << "\n";
	      found0=i;
	    }
	}

      // Is this normalization to the lineshape available?
      int found1=-1;
      for ( size_t i=0; i<lineshape_norm_data.size() && (found1 == -1); ++i )
	{
	  const string& str = lineshape_norm_data[i].Name();
	  if (str == normalizationfactor[k]) 
	    {
	      out2 << "  Selected normalization factor: " << normalizationfactor[k] << "\n";
	      found1=i;
	    }
	}


      // did we find the lineshape and normalization factor?
      if (found0 == -1)
	{
	  ostringstream os;
	  os << "Selected lineshape not available: "<< shape[k] <<"\n";
	  throw runtime_error(os.str());
	}
      if (found1 == -1)
	{
	  ostringstream os;
	  os << "Selected normalization to lineshape not available: "<< 
	    normalizationfactor[k] <<"\n";
	  throw runtime_error(os.str());
	}

      // now set the lineshape and lineshape_norm workspace variables 
      lineshape[k]=(size_t) found0;
	  lineshape_norm[k]=(size_t) found1;
    }
}



void raw_vmrs_1dReadFromScenario(// WS Output:
                                 ARRAYofMATRIX&   raw_vmrs_1d,
                                 // WS Input:
                                 const TagGroups& tag_groups,
                                 // Control Parameters:
                                 const string&    basename)
{
  // The species lookup data:
  extern const ARRAY<SpeciesRecord> species_data;

  // We need to read one profile for each tag group.
  for ( size_t i=0; i<tag_groups.size(); ++i )
    {
      // Determine the name.
      string name =
	basename + "." +
	species_data[tag_groups[i][0].Species()].Name() + ".am";

      // Add an element for this tag group to the vmr profiles:
      raw_vmrs_1d.push_back(MATRIX());

      // Read the VMR:
      // (We use the workspace method MatrixReadAscii for this.)
      MatrixReadAscii(raw_vmrs_1d[i],"",name);
    }
}

void Atm2dFromRaw1D(// WS Output:
                    ARRAYofVECTOR& 	 t_abs_2d,
                    ARRAYofVECTOR& 	 z_abs_2d,
                    ARRAYofMATRIX& 	 vmrs_2d,
                    // WS Input:      
                    const VECTOR&  	 p_abs,
                    const MATRIX&  	 raw_ptz_1d,
                    const ARRAYofMATRIX& raw_vmrs_1d)
{
  // FIXME: This function is terrible! Make this better using MTL functionality.

  // This function uses a lot of copying. Rather inefficient. The
  // problem is that the raw matrices do not directly fit the
  // interpolation routines. If this turns out to be too slow, it
  // should be replace by a routine using element-wise
  // interpolation. This can be also efficient, if the search for the
  // right interpolation point is done efficiently (see function
  // `hunt' in numerical recipies).

  // Also, I'm sure the copying could be done more elegantly, but I
  // just wasted a lot of time trying to do this with matrix / vector
  // notation. 


  //---------------< 1. Interpolation of temperature and altitude >---------------
  {  
    // Safety check: Make sure that raw_ptz_1d really is a [x,3] matrix:
    if ( 3 != raw_ptz_1d.dim(2) )
      {
	ostringstream os;
	os << "The variable raw_ptz_1d does not have the right dimensions,\n"
	   << "dim(2) should be 3, but is actually "<< raw_ptz_1d.dim(2);
	throw runtime_error(os.str());
      }

    // Break up raw_ptz_1d in p_raw, tz_raw.
    // The reason why we take tz_raw as a matrix is that the
    // interpolation can then be done simultaneously, hence slightly
    // more efficient.

    // p_raw is column 1:
    VECTOR p_raw;
    col( p_raw, 1, raw_ptz_1d );

    // tz_raw is column 2-3:
    MATRIX tz_raw;
    col( tz_raw, 2, 3, raw_ptz_1d );

    // Now interpolate tz_raw to p_abs grid:
    MATRIX tz_intp;
    interp_lin_col( tz_intp,
		    p_raw, tz_raw, p_abs );

    // Extract t_abs_2d:
    t_abs_2d.resize(0);
    t_abs_2d.push_back(VECTOR());
    col( t_abs_2d[0], 1, tz_intp );

    // Extract z_abs_2d:
    z_abs_2d.resize(0);
    z_abs_2d.push_back(VECTOR());
    col( z_abs_2d[0], 2, tz_intp );
  }

  //---------------< 2. Interpolation of VMR profiles >---------------
  {
    // We will write everything to the first array element of vmrs_2d
    // (more array elements would only be used in a 2D calculation).

    // Get room for our results:
    vmrs_2d.resize(0);
    vmrs_2d.push_back(MATRIX());
  
    // Get a convenient reference:
    MATRIX& intp = vmrs_2d[0];

    // Set dimensions.
    // The first dimension is the number of profiles (= the number of
    // tag groups). The second dimension is the dimension of the new
    // pressure grid.
    intp.resize( raw_vmrs_1d.size() , p_abs.size() );
  
    // We need this for each profile, therefore we define it here:
    VECTOR target;

    // For sure, we need to loop through all VMR profiles:
    for ( size_t j=0; j<raw_vmrs_1d.size(); ++j )
      {
	// Get a reference to the profile we are concerned with
	const MATRIX& raw = raw_vmrs_1d[j];

	// Raw should be a matrix with dimension [x,2], the first column
	// is the raw pressure grid, the second column the VMR values.
      
	// Safety check to ensure this:
	if ( 2 != raw.dim(2) )
	  {
	    ostringstream os;
	    os << "The variable raw_vmrs_1d("
	       << j
	       << ") does not have the right dimensions,\n"
	       << "dim(2) should be 2, but is actually "<< raw.dim(2);
	    throw runtime_error(os.str());
	  }

	// Extract p_raw and vmr_raw:
	VECTOR p_raw;
	col( p_raw, 1, raw );

	VECTOR vmr_raw;
	col( vmr_raw, 2, raw );

	// Interpolate:
	interp_lin( target,
		    p_raw, vmr_raw, p_abs );

	// Put the result in the apropriate row of intp:
	for ( size_t i=0; i<p_abs.size(); ++i )
	  {
	    intp(j+1,i+1) = target[i];
	  }
      }
  }
}

void AtmFromRaw1D(// WS Output:
		  VECTOR& 	 t_abs,
		  VECTOR& 	 z_abs,
		  ARRAYofVECTOR& vmrs,
		  // WS Input:      
		  const VECTOR&  	 p_abs,
		  const MATRIX&  	 raw_ptz_1d,
		  const ARRAYofMATRIX& raw_vmrs_1d)
{
  
  //---------------< 1. Interpolation of temperature and altitude >---------------
  {  
    // Safety check: Make sure that raw_ptz_1d really is a [x,3] matrix:
    if ( 3 != raw_ptz_1d.dim(2) )
      {
	ostringstream os;
	os << "The variable raw_ptz_1d does not have the right dimensions,\n"
	   << "dim(2) should be 3, but is actually "<< raw_ptz_1d.dim(2);
	throw runtime_error(os.str());
      }

    // Break up raw_ptz_1d in p_raw, tz_raw.
    // The reason why we take tz_raw as a matrix is that the
    // interpolation can then be done simultaneously, hence slightly
    // more efficient.

    // p_raw is column 1:
    VECTOR p_raw;
    col( p_raw, 1, raw_ptz_1d );

    // tz_raw is column 2-3:
    MATRIX tz_raw;
    col( tz_raw, 2, 3, raw_ptz_1d );

    // Now interpolate tz_raw to p_abs grid:
    MATRIX tz_intp;
    interp_lin_col( tz_intp,
		    p_raw, tz_raw, p_abs );

    // Extract t_abs:
    col( t_abs, 1, tz_intp );

    // Extract z_abs:
    col( z_abs, 2, tz_intp );
  }

  //---------------< 2. Interpolation of VMR profiles >---------------
  {
    // Make vmrs the right size:
    vmrs.resize(raw_vmrs_1d.size());
    
    // For sure, we need to loop through all VMR profiles:
    for ( size_t j=0; j<raw_vmrs_1d.size(); ++j )
      {
	// Get a reference to the profile we are concerned with:
	const MATRIX& raw = raw_vmrs_1d[j];

	// Get a reference to the place where we want to put the
	// interpolated profile:
	VECTOR& this_vmr = vmrs[j];

	// Raw should be a matrix with dimension [x,2], the first column
	// is the raw pressure grid, the second column the VMR values.
      
	// Safety check to ensure this:
	if ( 2 != raw.dim(2) )
	  {
	    ostringstream os;
	    os << "The variable raw_vmrs_1d("
	       << j
	       << ") does not have the right dimensions,\n"
	       << "dim(2) should be 2, but is actually "<< raw.dim(2);
	    throw runtime_error(os.str());
	  }

	// Extract p_raw and vmr_raw:
	VECTOR p_raw;
	col( p_raw, 1, raw );

	VECTOR vmr_raw;
	col( vmr_raw, 2, raw );

	// Interpolate:
	interp_lin( this_vmr,
		    p_raw, vmr_raw, p_abs );
      }
  }
}



// Algorithm based on equation from my (PE) thesis (page 274) and the book
// Meteorology today for scientists and engineers by R.B. Stull (pages 9-10). 
//
void z_absHydrostatic(
          VECTOR&    z_abs,
    const VECTOR&    p_abs,
    const VECTOR&    t_abs,
    const VECTOR&    h2o_abs,
    const Numeric&   g0,
    const Numeric&   pref,
    const Numeric&   zref,
    const int&       niter )
{
  const size_t   np = p_abs.size();
        size_t   i;                     // altitude index
        Numeric  g;                     // gravitational acceleration
        Numeric  r;                     // water mixing ratio in gram/gram
        Numeric  tv;                    // virtual temperature
        Numeric  dz;                    // step geometrical altitude
        VECTOR   ztmp(np);              // temporary storage for z_abs
  extern const Numeric EARTH_RADIUS;

  if ( (z_abs.size()!=np) || (t_abs.size()!=np) || (h2o_abs.size()!=np) )
    throw runtime_error("The vectors p_abs, t_abs, z_abs and h2o_abs do not all have the same length.");

  if ( niter < 1 )
    throw runtime_error("The number of iterations must be > 0.");

  for ( int iter=1; iter<=niter; iter++ )
  {
    // Init ztmp
    ztmp(1) = z_abs(1);

    // Calculate new altitudes (relative z_abs(1)) and store in ztmp
    for ( i=1; i<np; i++ )
    {
      // Calculate g 
      g  = ( g_of_z(EARTH_RADIUS,g0,z_abs(i)) + 
             g_of_z(EARTH_RADIUS,g0,z_abs(i+1)) ) / 2.0;

      // Calculate weight mixing ratio for water assuming constant average
      // molecular weight of the air
      r  = 18/28.96 * (h2o_abs(i)+h2o_abs(i+1))/2;

      // The virtual temperature (no liquid water)
      tv = (1+0.61*r) * (t_abs(i)+t_abs(i+1))/2;

      // The change in vertical altitude from i to i+1 
      dz = 287.053*tv/g * log( p_abs(i)/p_abs(i+1) );
      ztmp(i+1) = ztmp(i) + dz;
    }

    // Match the altitude of the reference point
    dz = interpp( p_abs, ztmp, pref ) - zref;
    z_abs = ztmp - dz;
  }
}



void absCalc(// WS Output:
             MATRIX&        		     abs,
             ARRAYofMATRIX& 		     abs_per_tg,
             // WS Input:		  
             const VECTOR&  		     f_mono,
             const VECTOR&  		     p_abs,
             const VECTOR&  		     t_abs,           
             const ARRAYofVECTOR&            vmrs,
             const ARRAYofARRAYofLineRecord& lines_per_tg,
	     const ARRAYofsizet&             lineshape,
	     const ARRAYofsizet&             lineshape_norm)
{
  // Check that vmrs and lines_per_tg really have the
  // same array dimension:
  if ( vmrs.size() != lines_per_tg.size() )
    {
      ostringstream os;
      os << "Variable vmrs must have the same dimension as lines_per_tg.\n"
	 << "vmrs.size() = " << vmrs.size() << '\n'
	 << "lines_per_tg.size() = " << lines_per_tg.size();
      throw runtime_error(os.str());
    }
  
  // Initialize abs and abs_per_tg. The array dimension of abs_per_tg
  // is the same as that of lines_per_tag.
  abs.resize(f_mono.size(), p_abs.size());
  abs = 0;
  abs_per_tg.resize(0);
  abs_per_tg.resize(lines_per_tg.size());

  // Call abs_species for each tag group.
  for ( size_t i=0; i<lines_per_tg.size(); ++i )
    {
      out2 << "  Tag group " << i << '\n';
      
      // Make this element of abs_per_tg the right size:
      abs_per_tg[i].resize(f_mono.size(), p_abs.size());
      abs_per_tg[i] = 0;

      abs_species( abs_per_tg[i],
		   f_mono,
		   p_abs,
		   t_abs,
		   vmrs[i],
		   lines_per_tg[i],
		   lineshape[i],
		   lineshape_norm[i]);
      
      // Add up to the total absorption:
      abs = abs + abs_per_tg[i];
    }
}



//// refr_indexBoudourisDryAir ///////////////////////////////////////////////
/**
   Calculates the refractive index for dry air at microwave frequncies 
   following Boudouris 1963.

   The expression is also found in Chapter 5 of the Janssen book.

   The atmosphere is assumed to have no water vapour.

   \retval   refr        refractive index
   \param    p_abs       absorption pressure grid
   \param    t_abs       temperatures at p_abs

   \author Patrick Eriksson
   \date   2000-09-30
*/
void refr_indexBoudourisDryAir (
                    VECTOR&          refr_index,
              const VECTOR&          p_abs,
              const VECTOR&          t_abs )
{
  refr_index = 1.0 + 0.77593e-6*ediv(p_abs,t_abs);
}
