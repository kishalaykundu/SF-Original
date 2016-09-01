/**
 * @file Driver.cpp
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
#include <cassert>
#include <cstring>

#include <vector>
#include <string>

extern "C" {
#include <dlfcn.h>
#include <libxml/tree.h>
#include <libxml/parser.h>
}

#include <boost/shared_ptr.hpp>

#include "Preprocess.h"
#include "aabb.h"

#include "Plugin.h"
#include "Display.h"
#include "Driver.h"

using namespace std;
using namespace boost;

namespace SF {

	// static function to parse and retrieve contents from a configuration file
	static bool
	parse (char* const cfgFile, vector <string>& types, vector <string>& properties, vector <string>& configs)
	{
		assert (cfgFile);

		// get pointer to xml file
		xmlDocPtr doc = xmlReadFile (cfgFile, NULL, 0);
		if (!doc){
			PRINT ("error: could not read %s\n", cfgFile);
			return false;
		}

		// get document root element
		xmlNodePtr node = xmlDocGetRootElement (doc);
		if (strcmp (reinterpret_cast <const char*> (node->name), "SFDriverConfig")){
			PRINT ("error: root element in %s in not of SFDriverConfig type", cfgFile);
			xmlFreeDoc (doc);
			return false;
		}

		// get children nodes
		node = node->children;
		node = node->next;

		while (node){

			// read in plugin-related properties
			if (!strcmp (reinterpret_cast <const char*> (node->name), "plugin")){

				char *pluginName = reinterpret_cast <char*> (xmlGetProp (node, reinterpret_cast <const xmlChar*> ("name")));
				assert (pluginName);
				char *configFile = reinterpret_cast <char*> (xmlGetProp (node, reinterpret_cast <const xmlChar*> ("config")));
				assert (configFile);

				types.push_back (string ("plugin"));
				properties.push_back (string (pluginName));
				configs.push_back (string (configFile));

				free (pluginName); pluginName = NULL;
				free (configFile); configFile = NULL;
			}
			// read in interface-related properties
			else if (!strcmp (reinterpret_cast <const char*> (node->name), "interface")){

				char* typeName = reinterpret_cast <char*> ( xmlGetProp (node, reinterpret_cast <const xmlChar*> ("type")));
				assert (typeName);
				char* configFile = reinterpret_cast <char*> ( xmlGetProp (node, reinterpret_cast <const xmlChar*> ("config")));
				assert (configFile);

				types.push_back (string ("interface"));
				properties.push_back (string (typeName));
				configs.push_back (string (configFile));

				free (typeName); typeName = NULL;
				free (configFile); configFile = NULL;
			}

			node = node->next; node = node->next;
		}

		// clean up and leave
		xmlFreeDoc (doc);
		xmlCleanupParser ();

		return true;
	}

	/**
	 * The only publicly available constructor for the Driver class.
	 * It reads a configuration file. Configuration file contains
	 * entries for each of the plugins that are to be loaded by the
	 * driver as well as interfaces (display or touch).
	 *
	 * Each line in configuration file contains three entries:
	 * 1. Module Type (plugin or interface),
	 * 2. Library/Interface-type name, and
	 * 3. Input configuration file.
	 *
	 * @param configFile Input XML styled configuration file
	 */
	Driver::Driver (int &argc, char **argv)
	{
		// sanity test
		char *cfgFile = argv [1];

		vector  <string> configFiles;
		vector  <string> pluginNames;

		// parse input configuration file and initialize display
		{
			vector  <string> moduleTypes;
			vector  <string> propertyNames;

			if (!parse (cfgFile, moduleTypes, propertyNames, configFiles)){
				PRINT ("error in parsing %s....Aborting\n", cfgFile);
				exit (EXIT_FAILURE);
			}
			assert (!moduleTypes.empty ());
			assert (!propertyNames.empty ());
			assert (!configFiles.empty ());

			// get display parameter and intialize display
			vector <string>::iterator miter (moduleTypes.begin ());
			vector <string>::iterator piter (propertyNames.begin ());
			vector <string>::iterator citer (configFiles.begin ());

/****************************************************************************************************/
/****************************************************************************************************/

			while (miter != moduleTypes.end ()){
				if (!miter->compare ("interface") && !piter->compare ("display")){
					_display = boost::shared_ptr <GL_Window> (new GL_Window (argc, argv, *citer));

					moduleTypes.erase (miter);
					propertyNames.erase (piter);
					configFiles.erase (citer);

					break;
				}
				else {
					++miter;
					++piter;
					++citer;
				}
			}
			assert (!moduleTypes.empty ());

			// copy plugin names and configuration files and erase the rest
			pluginNames = propertyNames;
		}

		// initialize plugins
		_plugins.reserve (pluginNames.size ());
		_pluginDestructors.reserve (pluginNames.size ());

		for (unsigned int i = 0; i  <pluginNames.size (); ++i){

			// open library
			void *handle = dlopen (pluginNames.at (i).c_str (), RTLD_NOW);
			if (!handle){
				PRINT ("error: could not open %s [%s]\n", pluginNames.at (i).c_str (), dlerror ());
				exit (EXIT_FAILURE);
			}

			// get plugin destructor function from plugin library
			PluginDestructor *destructor = reinterpret_cast <PluginDestructor*> (dlsym (handle, "DeletePlugin"));
			if (!destructor){
				PRINT ("error: could not find DeletePlugin in %s [%s]\n", pluginNames.at (i).c_str (), dlerror ());
				exit (EXIT_FAILURE);
			}

			// get plugin constructor function from plugin library
			PluginConstructor *constructor = reinterpret_cast <PluginConstructor*> (dlsym (handle, "NewPlugin"));
			if (!constructor){
				PRINT ("error: could not find NewPlugin in %s [%s]\n", pluginNames.at (i).c_str (), dlerror ());
				exit (EXIT_FAILURE);
			}

			// initialize plugin
			Plugin* plugin = constructor (configFiles.at (i), const_cast <Driver &> (* this));
			_plugins.push_back (plugin);
		}

		// allow plugins to synchronize between each other
		for (unsigned int i = 0; i  <_plugins.size (); ++i){
			_plugins.at (i)->synchronize (configFiles.at (i), _resources);
		}
	}

	// destructor
	Driver::~Driver ()
	{
		for (unsigned int i = 0; i  <_plugins.size(); ++i ){
			_pluginDestructors.at(i) (_plugins.at(i));
		}
	}

	// the run method: this starts the plugin threads and the display
	void
	Driver::run ()
	{
	  // register self to the display
	  _display.get ()->_parent = const_cast <Driver *> (this);

		for (unsigned int i = 0; i  < _plugins.size(); ++i){
			_plugins.at (i)->run ();
		}
		_display.get ()->run ();
	}

  // the cleanup method
	void
	Driver::cleanup ()
	{
	  for (unsigned int i = 0; i < _plugins.size (); ++i){
      _plugins.at (i)->cleanup ();
	  }
	}
}
