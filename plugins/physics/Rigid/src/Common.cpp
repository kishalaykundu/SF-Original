/**
 * @file Common.cpp
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
 * Common functions for the Rigid library.
 */
#include <cstring>
#include <string>

extern "C" {
#include <libxml/tree.h>
#include <libxml/parser.h>
}

#include "Preprocess.h"
#include "Common.h"

using namespace std;

namespace SF {
	namespace RM {

		bool getConfigParameter (const string &cfgFile, const char *param, string &result)
		{
			if (!result.empty ()){
				result.clear ();
			}

			// get a pointer to xml file
			xmlDocPtr doc = xmlReadFile (cfgFile.c_str (), NULL, 0);
			if (!doc){
				PRINT ("error: could not read %s\n", cfgFile.c_str ());
				return false;
			}

			// get document root element
			xmlNodePtr node = xmlDocGetRootElement (doc);
			if (strcmp (reinterpret_cast< const char* > (node->name), "SFRigidInfo")){
				PRINT ("error: root element in %s in not of SFXFEMInfo type", cfgFile.c_str ());
				xmlFreeDoc (doc);
				return false;
			}

			// get children nodes
			node = node->children;
			node = node->next;

			while (node){

				if (!strcmp (reinterpret_cast <const char *> (node->name), "dataInfo")){
					char *fname = reinterpret_cast <char *> (xmlGetProp (node, reinterpret_cast <const xmlChar *> (param)));
					if (fname){
            result = string (fname);
            free (fname); fname = NULL;
            break;
					}
				}

				node = node->next;
				node = node->next;
			}

			// clean up and leave
			xmlFreeDoc (doc);
			xmlCleanupParser ();

			return true;
		}
	}
}
