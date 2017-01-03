/*
 Boost Software License - Version 1.0 - August 17th, 2003

 Permission is hereby granted, free of charge, to any person or organization
 obtaining a copy of the software and accompanying documentation covered by
 this license (the "Software") to use, reproduce, display, distribute,
 execute, and transmit the Software, and to prepare derivative works of the
 Software, and to permit third-parties to whom the Software is furnished to
 do so, all subject to the following:

 The copyright notices in the Software and this entire statement, including
 the above license grant, this restriction and the following disclaimer,
 must be included in all copies of the Software, in whole or in part, and
 all derivative works of the Software, unless such copies or derivative
 works are solely in the form of machine-executable object code generated by
 a source language processor.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
 SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
 FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
 ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 DEALINGS IN THE SOFTWARE.
 */

#include "thread_pool.h"

// Just before listing 8.4
join_threads::join_threads(std::vector<std::thread>& threads_):
  threads(threads_)
{}

join_threads::~join_threads()
{
  for(unsigned long i=0;i<threads.size();++i)
  {
    if(threads[i].joinable())
      threads[i].join();
  }
}

// Listing 9.7:
void work_stealing_queue::push(data_type data)
{
  std::lock_guard<std::mutex> lock(the_mutex);
  the_queue.push_front(std::move(data));
}

bool work_stealing_queue::empty() const
{
  std::lock_guard<std::mutex> lock(the_mutex);
  return the_queue.empty();
}

bool work_stealing_queue::try_pop(data_type& res)
{
  std::lock_guard<std::mutex> lock(the_mutex);
  if(the_queue.empty())
  {
    return false;
  }

  res=std::move(the_queue.front());
  the_queue.pop_front();
  return true;
}

bool work_stealing_queue::try_steal(data_type& res)
{
  std::lock_guard<std::mutex> lock(the_mutex);
  if(the_queue.empty())
  {
    return false;
  }

  res=std::move(the_queue.back());
  the_queue.pop_back();
  return true;
}
