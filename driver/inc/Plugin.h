/**
 * @file Plugin.h
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
#pragma once

#include <vector>
#include <string>
#include <boost/thread.hpp>
#include <boost/shared_ptr.hpp>

using namespace std;
using namespace boost;

namespace SF {

	class Driver;
	class Resource;

	class Plugin {

	private:
		boost::thread *_threads;

	public:
		std::vector <boost::shared_ptr <Resource> > _resources;

	private:
		Plugin ()
		: _threads (NULL)
		{ }

		Plugin (const Plugin &p)
		: _threads (p._threads), _resources (p._resources)
		{ }

		Plugin& operator = (const Plugin& p)
		{
			_threads = p._threads;
			_resources = p._resources;

			return *this;
		}

	public:
		Plugin (const string &, Driver &);
		virtual ~Plugin ();

		virtual void synchronize (const string &config, const vector <boost::shared_ptr <Resource> > &driverResources);
		virtual void run ();
		virtual void cleanup ();
	};

	// set of functions to get constructors and destructors from shared libs
	extern "C" Plugin* NewPlugin (const string &c, Driver &d);
	typedef Plugin* PluginConstructor (const string &, Driver &);

	extern "C" void DeletePlugin (Plugin *t);
	typedef void PluginDestructor (Plugin *);
}
