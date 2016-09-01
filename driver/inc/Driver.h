/**
 * @file Driver.h
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
 * The driver for the Simulate framework.  This is the central authority
 * that loads and stores all the plugins. Additionally, it also stores all
 * the resources that used by the plugins. It is through the driver that
 * all the plugins exchange information and synchronize with  each other.
 * The Driver also contains the display class that is in charge of
 * controlling rendering related behavior.
 */

#pragma once

#include <vector>
#include <boost/shared_ptr.hpp>

#include "Plugin.h"

#include "Preprocess.h"

using namespace std;
using namespace boost;

namespace SF {

	class GL_Window;
	class Resource;

	class Driver {

	public:

		// display-feedback specific data
		boost::shared_ptr< GL_Window > _display;

		// resources
		vector< boost::shared_ptr< Resource > > _resources;

		// plugin-library-specific variables
	protected:
		vector< Plugin* > _plugins;
		vector< PluginDestructor* > _pluginDestructors;

	public:
		Driver (int &argc, char **argv);
		~Driver ();

	public:
		void run ();
		void cleanup ();
	};
}
