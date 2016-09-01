/**
 * @file ThreadControl.h
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
 * The thread control structure for plugins controlled by Driver.
 */

# pragma once

#include <cassert>
#include <vector>
#include <boost/interprocess/sync/interprocess_semaphore.hpp>

using namespace std;
using namespace boost::interprocess;

namespace SF {

	class ThreadControl {

	protected:
		vector< interprocess_semaphore* > _mutex;

	public:
		ThreadControl(){}
		~ThreadControl()
		{
			for(size_t i = 0; i < _mutex.size(); ++i ){
				delete _mutex.at(i);
			}
		}

		// method to check if empty
		inline bool empty(){ return _mutex.empty(); }

		// method to check size
		inline size_t size(){ return _mutex.size(); }

		// method to access individual members
		inline interprocess_semaphore&
		operator [] (unsigned int index)
		{
			assert (index < static_cast <unsigned int>(_mutex.size()));
			return *(_mutex[index]);
		}

		// method to push a semaphore
		inline void
		push_back(unsigned int number)
		{
			interprocess_semaphore* ism = new interprocess_semaphore (number);
			_mutex.push_back (ism);
		}
	};
}
