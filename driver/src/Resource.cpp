/**
 * @file Resource.cpp
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

#include "Resource.h"

namespace SF {

	// default constructor
	Resource::Resource ()
	: draw (NULL), touch (NULL), transform (NULL), reprogram (NULL)
	{ }

	// copy constructor
	Resource::Resource (const Resource& r)
	: _name (r._name), _owner (r._owner),
	  draw (r.draw), touch (r.touch), transform (r.transform), reprogram (r.reprogram)
	{ }

	// assignment operator
	Resource&
	Resource::operator = (const Resource& r)
	{
		_name = r._name;
		_owner = r._owner;

		draw = r.draw;
		touch = r.touch;
		transform = r.transform;
		reprogram = r.reprogram;

		return *this;
	}
}
