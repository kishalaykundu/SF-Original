/**
 * @file main.cpp
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
 * The driver of the Integrate Framework. This is where everything starts.
 */
#include <cstdlib>
#include <cstring>
#include <iostream>

#include "Preprocess.h"
#include "Driver.h"

using namespace std;
using namespace SF;

int
main (int argc, char* argv[])
{
	// sanity check
	if (argc < 2 || !strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")){
		cerr << "Usage: runsim <input configuration file>\n" << endl;
	}

	Driver _driver (argc, argv);

	_driver.run ();

	exit (EXIT_SUCCESS);
}
