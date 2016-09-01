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

#include "Common.h"
#include "Mesh.h"
#include "Scene.h"

using namespace std;
using namespace boost;

namespace SF {

  XFE::Scene _scene;

	// function to parse configuration file
	static int
	parse (const string &cfgFile, vector <string> &configs)
	{
		assert (!cfgFile.empty ());
		assert (configs.empty ());

		// get a pointer to xml file
		xmlDocPtr doc = xmlReadFile (cfgFile.c_str (), NULL, 0);
		if (!doc){
			PRINT ("error: could not read %s\n", cfgFile.c_str ());
			return 0;
		}

		// get document root element
		xmlNodePtr node = xmlDocGetRootElement (doc);
		if (strcmp (reinterpret_cast< const char* > (node->name), "SFXFEMConfig")){
			PRINT ("error: root element in %s in not of SFXFEMConfig type", cfgFile.c_str ());
			xmlFreeDoc (doc);
			return 0;
		}

		// get children nodes
		node = node->children;
		node = node->next;

		int result = 0;
		while (node){

			if (!strcmp (reinterpret_cast <const char *> (node->name), "configFile")){
				char *fname = reinterpret_cast <char *> (xmlGetProp (node, reinterpret_cast <const xmlChar *> ("name")));
				configs.push_back (string (fname));
				free (fname); fname = NULL;
			}
			else if (!strcmp (reinterpret_cast <const char *> (node->name), "threadpool")){
				char *sname = reinterpret_cast <char *> (xmlGetProp (node, reinterpret_cast <const xmlChar *> ("size")));
				for (unsigned int i = 0; i < strlen (sname); ++i){
					if (!isdigit (sname[i])){
						PRINT ("error: threadpool size \'%s\' is not a number", sname);
						free (sname);
						xmlFreeDoc (doc);
						xmlCleanupParser ();
						return 0;
					}
				}
				result = atoi (sname);
				free (sname); sname = NULL;
			}

			node = node->next;
			node = node->next;
		}
		assert (!configs.empty ());

		// clean up and leave
		xmlFreeDoc (doc);
		xmlCleanupParser ();

		return result;
	}

	// plugin constructor
	EXPORT
	Plugin::Plugin (const string &config, Driver &driver)
	{
		// parse input configuration files
		vector <string> configFiles;

    unsigned int numThreads = static_cast <unsigned int> (parse (config, configFiles));
		if (!numThreads){
			PRINT ("error parsing %s....aborting\n", config.c_str ());
			exit (EXIT_FAILURE);
		}
		assert (!configFiles.empty ());
		_scene.resizePool (numThreads);

		_resources.reserve (configFiles.size ());
		for (unsigned int i = 0; i < configFiles.size (); ++i){

			Resource *mesh = new XFE::Mesh (configFiles.at (i), driver);
			_resources.push_back (boost::shared_ptr <Resource> (mesh));
			_scene.addMesh (*mesh);

			driver._resources.push_back (_resources.at (i));
			driver._display.get ()->addDrawables (_resources.at (i));
		}
	}

	// plugin destructor
	EXPORT
	Plugin::~Plugin () { PRINT ("here\n"); }

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
	void
	Plugin::synchronize (const string &config, const vector <boost::shared_ptr <Resource > > &resources)
	{
	  vector <string> configFiles;
    parse (config, configFiles);

    Resource *r;
    string bladeOwner, bladeName;
    for (unsigned int i = 0; i < configFiles.size (); ++i){

      // register blade
      XFE::getConfigParameter (configFiles [i], "blade_name", bladeName);
      XFE::getConfigParameter (configFiles [i], "blade_owner", bladeOwner);

      if (!bladeName.empty () && !bladeOwner.empty ()){
        for (unsigned int j = 0; j < resources.size (); ++j){
          r = resources [j].get ();
          if (!r->_name.get ()->compare (bladeName) && !r->_owner.get ()->compare (bladeOwner)){
            _scene.addBlade (r);
            break;
          }
        }
      }
    }
	}

	// run method
	void
	Plugin::run ()
	{
    _threads = new boost::thread [_resources.size () + 1];

    // start scene thread
    boost::thread xs (&SF::XFE::Scene::run, dynamic_cast <SF::XFE::Scene *> (&_scene));
    _threads [0] = boost::move (xs);

    for (unsigned int i = 1; i <= _resources.size (); ++i){
      boost::thread xm (&SF::XFE::Mesh::run, dynamic_cast <SF::XFE::Mesh *> (_resources.at (i - 1).get ()));
      _threads [i] = boost::move (xm);
    }
    PRINT ("libCudaXfem threads started\n");
	}

	// cleanup method
	void
	Plugin::cleanup ()
	{
	  for (unsigned int i = 1; i <= _resources.size (); ++i){
      boost::thread xm (&SF::XFE::Mesh::cleanup, dynamic_cast <SF::XFE::Mesh *> (_resources.at (i - 1).get ()));
      _threads [i] = boost::move (xm);
	  }
	}

}
