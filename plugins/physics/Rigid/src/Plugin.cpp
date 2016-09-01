/**
 * @file Plugin.cpp
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
 * The plugin class for the driver in the Simulate framework. The actual
 * implementation is left to the individual implementation libraries.
 * The Plugin in the Simulate framework acts as the manager of a set of
 * resources. The resources may be in the form of the library's choosing.
 * The synchronicity/ data-safety of the threads with respect to the
 * resources is assured by having the plugin threads access these resources
 * using pre-defined boost semaphores.
 */

#include <cassert>
#include <cstring>

extern "C" {
#include <libxml/tree.h>
#include <libxml/parser.h>
}

#include <vector>
#include <string>

#include <boost/thread.hpp>
#include <boost/shared_ptr.hpp>

#include "Preprocess.h"
#include "aabb.h"
#include "ThreadControl.h"
#include "Plugin.h"
#include "Driver.h"
#include "Display.h"

#include "Mesh.h"

using namespace std;
using namespace boost;

namespace SF {

	// function to parse configuration file
	static bool
	parse (const string &cfgFile, vector <string> &configs)
	{
		assert (!cfgFile.empty ());
		assert (configs.empty ());

		// get a pointer to xml file
		xmlDocPtr doc = xmlReadFile (cfgFile.c_str (), NULL, 0);
		if (!doc){
			PRINT ("error: could not read %s\n", cfgFile.c_str ());
			return false;
		}

		// get document root element
		xmlNodePtr node = xmlDocGetRootElement (doc);
		if (strcmp (reinterpret_cast< const char* > (node->name), "SFRigidConfig")){
			PRINT ("error: root element in %s in not of SFRigidConfig type", cfgFile.c_str ());
			xmlFreeDoc (doc);
			return false;
		}

		// get children nodes
		node = node->children;
		node = node->next;
		while (node){

			if (!strcmp (reinterpret_cast <const char *> (node->name), "configFile")){
				char *fname = reinterpret_cast <char *> (xmlGetProp (node, reinterpret_cast <const xmlChar *> ("name")));
				configs.push_back (string (fname));
				free (fname); fname = NULL;
			}
			node = node->next;
			node = node->next;
		}
		assert (!configs.empty ());

		// clean up and leave
		xmlFreeDoc (doc);
		xmlCleanupParser ();

		return true;
	}

  // plugin constructor
  EXPORT
  Plugin::Plugin (const string &config, Driver &driver)
  {
		// parse input configuration files
		vector <string> configFiles;

		if (!parse (config, configFiles)){
			PRINT ("error parsing %s....aborting\n", config.c_str ());
			exit (EXIT_FAILURE);
		}
		assert (!configFiles.empty ());

		_resources.reserve (configFiles.size ());
		for (size_t i = 0; i < configFiles.size (); ++i){

			Resource *mesh = new RM::Mesh (configFiles.at (i), driver);
			_resources.push_back (boost::shared_ptr <Resource> (mesh));

			driver._resources.push_back (_resources.at (i));
			driver._display.get ()->addDrawables (_resources.at (i));
			driver._display.get ()->addMoveables (_resources.at (i));
		}
  }

  // plugin destructor
  EXPORT
  Plugin::~Plugin () { }

	// function to bind constructor to shared library symbol list
	EXPORT
	Plugin *NewPlugin (const string &config, Driver &driver)
	{
		return new Plugin (config, driver);
	}

	// function to bind destructor to shared library symbol list
	EXPORT
	void DeletePlugin (Plugin *p)
	{
		delete p;
	}

	// method to synchronize own resources with those of other plugins
	void Plugin::synchronize (const string &config, const vector <boost::shared_ptr <Resource > > &resources) { }

	// run method
	void
	Plugin::run ()
	{
    _threads = new boost::thread [_resources.size ()];

    for (size_t i = 0; i < _resources.size (); ++i){
      boost::thread xm (&SF::RM::Mesh::run, dynamic_cast <SF::RM::Mesh *> (_resources.at (i).get ()));
      _threads [i] = boost::move (xm);
    }
    PRINT ("libRigid threads started\n");
	}

	// cleanup method
	void
	Plugin::cleanup ()
	{
	}
}
