/**
 * @file Resource.h
 * @author Kishalay Kundu <kishalay.kundu@gmail.com>
 * @section LICENSE
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details at
 * http://www.gnu.org/copyleft/gpl.html
 *
 * @section DESCRIPTION
 * The virtual resource class for used by Plugin's to store and manipulate data.
 * Function pointers in Resource need to be implemented by derived classes.
 */

#pragma once

#include <string>
#include <boost/shared_ptr.hpp>

using namespace std;
using namespace boost;

namespace SF {

	class Resource {

	public:
		boost::shared_ptr < string > _name;
		boost::shared_ptr < string > _owner;

	public:
		Resource ();
		virtual ~Resource () { }

		Resource (const Resource&);
		Resource& operator = (const Resource&);

		// methods to be custom-designed by derived classes
		void (* draw) (Resource&);
		void (* touch) (Resource&);
		void (* transform) (Resource&);
		void (* reprogram) (Resource&);
	};
}
